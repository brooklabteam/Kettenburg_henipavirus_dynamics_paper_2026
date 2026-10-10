rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(mclust)
library(ggh4x)
library(ggcorrplot)
library(factoextra)
library(paran)
library(hopkins)
library(cluster)
library(corrr) 
library(FactoMineR)
library(ggpmisc)
library(heatmaply)
library(seqinr)
library(scales)
library(ggdendroplot)

#use this script to calculate cutoffs for the luminex


#first record your wd for this repo (home) 
#cara adds hers here too - just comment out whichever you are not using
#homewd = "/Users/gwenddolenkettenburg/Desktop/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"


#now load in the merged luminex data from the "data" folder

#make sure this is latest catching file
dat <- read.csv(paste0(homewd, "/serology/data/lum_catch_merge.csv"), header = T, stringsAsFactors = F)
head(dat)


length(unique(dat$sampleid))
dat$sampleid[duplicated(dat$sampleid)]

test <- ddply(dat, .(bat_species, bat_sex), summarise, Ntot=length(unique(sampleid)))
test

#now for cutoffs, just slim to antigen data only and convert to long format
lumdat <- dplyr::select(dat, sampleid, bat_species, names(dat)[2:7])
head(lumdat)
lumdat.long <- melt(lumdat)
head(lumdat.long)
names(lumdat.long)[names(lumdat.long)=="variable"] <- "antigen"
names(lumdat.long)[names(lumdat.long)=="value"] <- "MFI"
unique(lumdat.long$antigen)
sort(unique(lumdat.long$MFI))
class(lumdat.long$MFI)
lumdat.long <- subset(lumdat.long, !is.na(MFI))

# # #convert any negative values to 1
# lumdat.long$MFI[lumdat.long$MFI<0] <- 1
# lumdat.long$log10MFI <- log10(lumdat.long$MFI)

lum.split <- dlply(lumdat.long,.(bat_species, antigen))

correct.negs <- function(df){
  abs.min <- abs(min(df$MFI))
  df$MFI <- df$MFI + abs.min + 1
  
  #then, also take the log10
  df$log10MFI <- log10(df$MFI)
  return(df)
}

lumdat.long <- data.table::rbindlist(lapply(lum.split, correct.negs))


#add family
lumdat.long$virus_family[lumdat.long$antigen== "NiV"| lumdat.long$antigen=="HeV"| lumdat.long$antigen=="CedV"| lumdat.long$antigen=="GhV"| lumdat.long$antigen=="AngV"| lumdat.long$antigen=="MojV"] <- "Paramyxovirus"
unique(lumdat.long$virus_family)

#and genus information
lumdat.long$virus_genus <- "Henipavirus"
unique(lumdat.long$virus_genus)

#factor antigens by relatedness
lumdat.long$antigen <- factor(lumdat.long$antigen, levels=c( "NiV", "HeV","CedV", "GhV", "AngV", "MojV"))

lumdat.long$bat_species <- factor(lumdat.long$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))



########################################################################################################
########################################################################################################
# Now, let's use dimensionality reduction to determine which antigens in the data we should retain:

# First, let's eliminate the Pteropus from consideration
# First, make a data subset of all bats with any MFI > 1000
dat.ind.split <- dlply(lumdat.long, .(sampleid, virus_family))
get.highMFI <- function(dat, cutoff){
  dat$MFI_high <- 0
  
  for(i in 1:nrow(dat)){
    if(dat$MFI[i]>cutoff){
      dat$MFI_high[i] <- 1
    }
  }
  if(sum(dat$MFI_high)>0){
    dat$MFI_high <- 1 
  }
  return(dat)
}

bat.df <- data.table::rbindlist(lapply(dat.ind.split, get.highMFI, cutoff=1000))
head(bat.df)
bat.high.df = subset(bat.df, MFI_high==1) #only 6 individual Pteropus with any MFI values >1000

#and plot just the Pteropus
p3 <- ggplot(subset(bat.high.df, bat_species=="Pteropus rufus")) + coord_cartesian(ylim = c(0,10000)) +
  facet_nested(bat_species~virus_family, scales = "free", space = "free") +
  geom_line(aes(x=antigen, y=MFI, color=sampleid, group=sampleid), size=.8) + theme_bw()+
  geom_point(aes(x=antigen, y=MFI, color=sampleid), size=2) + 
  #scale_color_manual(name="Individual Bat\nSample ID", values = c("tomato", "#00BFC4")) +
  geom_hline(aes(yintercept = 1000), color="red", linetype=2, size=.8) +
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        legend.position = c(.12,.84),
        axis.text.x = element_text(size=16), axis.title.y = element_text(size=16), 
        axis.text.y = element_text(size=13), strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=16, face="italic"), strip.text.x = element_text(size=16)) +
  guides(color=guide_legend(ncol=2)); p3 
# #so few bats, these are probably exposures from other bats
# 
# ggsave(file = paste0(homewd,"/serology/figures/pteropus-spillover-plot.png"),
#        plot = p3,
#        units="mm",  
#        width=80, 
#        height=40, 
#        scale=3.3, 
#        dpi=300)


# Now, focusing on the Eidolon and Rousettus let's decide which 
# antigens to focus on

# now try pca - using the untransformed data

# # # multiply MFI by 5 if 1:500 dilution is used
# lumdat.long$MFI = lumdat.long$MFI*5

# first, make subsets of the types of data you want to include in the PCA - 
# separate by virus genus and bat species
unique(lumdat.long$virus_genus) #"Henipavirus"  "Ebolavirus"   "Marburgvirus" "Dianlovirus"  "Cuevavirus"  
unique(lumdat.long$antigen[lumdat.long$virus_genus=="Henipavirus"]) # many!
unique(lumdat.long$MFI)
sum(is.na(lumdat.long$MFI))#0


#function to select antigens and evaluate
choose.antigens <- function(dat, species, genus, manual_cluster, do.save.plot, filename){
  #shrink dataset to area of focus
  df.shrink = subset(data.frame(dat), bat_species==species & virus_genus==genus)
  
  #select columns of relevance 
  df.slim <- dplyr::select(df.shrink, sampleid, antigen, MFI)
  
  #make wide
  df.wide = reshape(df.slim, direction="wide", idvar="sampleid",timevar = "antigen")
  
  #check for NAs in case some were missing a single antigen - replace with 1
  df.wide[is.na(df.wide)] <- 1
  
  #make a scaled and centered version 
  df.scaled <- scale(df.wide[,2:ncol(df.wide)], center = T, scale = T)
  
  #perform PCA - could do on scaled version but better to use original so you can plot by parts of it in the future
  pca.out <- prcomp(scale(df.wide[,2:ncol(df.wide)], center = T, scale = T), center = FALSE, scale. = FALSE) 
  pca.dat <-   pca.out$x
  #now, no scaling or centering needed in the actual PCA term
  sum.dat <- summary(pca.out) # how much cumulative variance per component
  #pca.out$rotation # loadings on each PCA
  p1 <- fviz_eig(pca.out, addlabels = TRUE) #scree plot
  #can see there is a strong pull from BdBV and a bit of overlap with EboV and the rest are just a grouping
  p2 <- fviz_cos2(pca.out, choice = "var", axes = 1) # contributions to Dim1
  p3 <- fviz_cos2(pca.out, choice = "var", axes = 2) # contributions to Dim2
  p4 <- fviz_cos2(pca.out, choice = "var", axes = 3) # contributions to Dim3
  #p5 <- fviz_cos2(pca.out, choice = "var", axes = 4) # contributions to Dim4
   
  eig.val<-get_eigenvalue(pca.out)
  #eig.val #retain number of dimensions in  where eigenvalue > 1
  
  
  dat.parallel <- paran(df.scaled, iterations=10000, quietly=FALSE,
        status=FALSE, all=TRUE, cfa=FALSE, graph=FALSE, seed=0, mat=NA, n=NA)
        #color=TRUE, col=c("black","red","blue"),
        #lty=c(1,2,3), lwd=1, legend=TRUE, file="",
        #width=640, height=640, grdevice="png", seed=0, mat=NA, n=NA)
  

  par.plot <- cbind.data.frame(components = rep(seq(1, length(c(dat.parallel$AdjEv)),1),3), eigenvalues = c(c(dat.parallel$AdjEv), c(dat.parallel$Ev), c(dat.parallel$RndEv)), status =rep(c("Adjusted Ev", "Unadjusted Ev", "Random Ev"), each = length(c(dat.parallel$AdjEv))))
  par.plot$status[par.plot$status=="Adjusted Ev"][(dat.parallel$Retained+1):length(c(dat.parallel$AdjEv))] <- "Adjusted Ev (unretained)"
  par.plot$status[par.plot$status=="Adjusted Ev"] <- "Adjusted Ev (retained)"
  
  #occasionally, if optimal of 1 component is selected, then you will have only 1 value
  #for the Adjusted Ev (retained). If that is the case, add one here
  
  if (length(par.plot$eigenvalues[par.plot$status=="Adjusted Ev (retained)"])==1){
    add <- subset(par.plot, status=="Adjusted Ev (retained)")
    add$components = add$components +1
    add$eigenvalues = par.plot$eigenvalues[par.plot$status=="Adjusted Ev (unretained)" & par.plot$components==unique(add$components)]
    par.plot <- rbind(par.plot, add)
  }
  
  colorz = c('Adjusted Ev (retained)' = "black", 'Adjusted Ev (unretained)' = "black", 'Unadjusted Ev' = "red", 'Random Ev' = "blue")
  fillz = c('Adjusted Ev (retained)' = "black", 'Adjusted Ev (unretained)' = "white", 'Unadjusted Ev' = "red", 'Random Ev' = "blue")
  
  p5 <- ggplot(data=par.plot) + theme_bw() + coord_cartesian(ylim=c(0,(max(par.plot$eigenvalues)+.1))) +
        geom_hline(aes(yintercept = 1), linetype=2, size=.8) +
        scale_color_manual(values=colorz) + scale_fill_manual(values=fillz) +
        geom_line(aes(x=components, y = eigenvalues, color=status), size=1) +
        geom_point(aes(x=components, y = eigenvalues, color=status, fill=status),   size=3, shape=21) +    
        geom_point(data = subset(par.plot, status=="Adjusted Ev (unretained)"), aes(x=components, y = eigenvalues, color=status, fill=status),   size=3, shape=21) +    
        theme(panel.grid = element_blank(), legend.title = element_blank(),
              legend.position = "inside", legend.position.inside = c(.78,.82))
  
  
  #hopkins stat for clusterability
  
  # Calculated values 0-0.3 indicate regularly-spaced data.
  # Values around 0.5 indicate random data.
  # Values 0.7-1 indicate clustered data.
  h.stat = hopkins(pca.dat, m=ceiling(nrow(pca.dat)/10))
  
  
  
  # Optimal cluster count by avg Silhouette width
  # kmeans
  p6a <- fviz_nbclust(pca.dat, FUNcluster=kmeans, k.max = (length(unique(df.shrink$antigen))+1)) + ggtitle("Optimal number of clusters by k-means")

  # Elbow method for kmeans - abandon
  # fviz_nbclust(pca.dat, kmeans, method = "wss") # not quantitative
  
  # clustering around medioids
  p6b <- fviz_nbclust(pca.dat, FUNcluster=cluster::pam, k.max = (length(unique(df.shrink$antigen))+1)) + ggtitle("Optimal number of clusters by k-medoids")
  
  
  if (is.na(manual_cluster)){
    pam1 <-eclust(pca.dat, "pam", k=as.numeric(as.character(p6b$data$clusters[p6b$data$y==max(p6b$data$y)])))   
  }else{
    pam1 <-eclust(pca.dat, "pam", k=as.numeric(manual_cluster))   
  }
  
  
  
  cluster <- pam1$cluster
  #and link to bat ID
  df.wide$cluster <- cluster
  
  # df.name.clust <- dplyr::select(df.wide, sampleid, cluster)
  # 
  # df.tmp <-  subset(dat, bat_species==species)
  # 
  # df.tmp <- merge(df.tmp, df.name.clust, by="sampleid")
  # df.tmp$cluster <- as.character(df.tmp$cluster)
  
  
  
  pca.out.new <- prcomp(scale(df.wide[,2:(ncol(df.wide)-1)], center = T, scale = T), center = FALSE, scale. = FALSE) 
  
  pca.out.new.dat <-  pca.out.new$x
  pca.out.new.plot <- cbind.data.frame(PC1=pca.out.new.dat[,1], PC2=pca.out.new.dat[,2])
  pca.out.new.plot$cluster <- cluster  
  pca.out.new.plot$cluster <- as.factor(pca.out.new.plot$cluster)
  df.wide$cluster <- as.factor(df.wide$cluster)
  
  
  p7 <- fviz_pca_biplot(pca.out.new, geom=c("point"), 
                  #habillage = df.wide$cluster,
                  fill.ind = df.wide$cluster, #select.ind = list(cos2=0.3),
                  col.var = c("cos2"),
                  #select.var = NULL, #c("MFI.AngV", "MFI.GhV", "MFI.NiV", "MFI.HeV", "MFI.MojV", "MFI.CedV"),
                  #col.var = c("cornflowerblue","mediumseagreen","cornflowerblue","cornflowerblue","mediumseagreen","tomato"),, 
                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                  #gradient.cols = c("black", "orange", "green"),
                  alpha.ind =0,
                  addEllipses=TRUE, ellipse.level=0.5,
                  repel = TRUE)
  
  
  p8 <- p7 + ggnewscale::new_scale_color() +  
        geom_point(data=  pca.out.new.plot, aes(x=PC1, y=PC2, color=cluster), size=2)
  
  
  #add to output
  pca.out$eigen <- eig.val
  pca.out$hopkins <-   h.stat
  
  
  
  #make plot
  out.plot.mid <- cowplot::plot_grid(p2,p3,p4, ncol=3,nrow=1, align = "hv", axis="tb", labels = c("C", "D", "E"), label_size = 22)
  out.plot.top <- cowplot::plot_grid(p1,p5,ncol=2,nrow = 1, align = "hv", axis="tb" , labels = c("A", "B"), label_size = 22)
  out.plot.top2 <- cowplot::plot_grid(p6a,p6b,ncol=2,nrow = 1, align = "hv", axis="tb", labels = c("F", "G"), label_size = 22)
  out<- cowplot::plot_grid(out.plot.top, out.plot.mid, out.plot.top2, p8, ncol=1, nrow = 4, labels = c("", "", "", "H"), label_size = 22) #+ theme(panel.background = element_rect(fill="white"), plot.background = element_rect(fill="white"))
  #print(out) 
  
  if(do.save.plot ==TRUE){

    ggsave(file = filename,
           plot = out,
           bg = "white",
           units="mm",  
           width=110, 
           height=140, 
           scale=3, 
           dpi=300)
    
  
  }
  #and return all the info
  return(pca.out)
  
  
  
}

PCA.Eid.hen <- choose.antigens(dat = lumdat.long, species = "Eidolon dupreanum", genus = "Henipavirus",  manual_cluster = NA,  do.save.plot = TRUE, filename=paste0(homewd,"/final_figures_tables/supplemental/FigSXX_PCA_plots_Eidolon.png"))
PCA.Rou.hen <- choose.antigens(dat = lumdat.long, species = "Rousettus madagascariensis", genus = "Henipavirus",  manual_cluster = NA, do.save.plot = TRUE, filename=paste0(homewd,"/final_figures_tables/supplemental/FigSXX_PCA_plots_Rousettus.png"))
PCA.Pter.hen <- choose.antigens(dat = lumdat.long, species = "Pteropus rufus", genus = "Henipavirus",  manual_cluster = NA, do.save.plot = TRUE,  filename=paste0(homewd,"/final_figures_tables/supplemental/FigSXX_PCA_plots_Pteropus.png"))

################################################################################
################################################################################

# Now, within each cluster, demonstrate that all MFIs correlate with the 
# most dominant antigen in the cluster


# load aa alignment for henipavirus
henipa.align <- read.alignment(file = paste0(homewd,"/serology/data/Henipavirus-Alignment.fasta"), format = "fasta")
# make distance matrix
dist.henipa <- dist.alignment(henipa.align)

hclust(dist.henipa)

plot.correlation <- function(dat, ant1, ant2, ant3, ant4, ant5, species, genus, dist.mat){
  
  
  #shrink dataset to area of focus
  df.shrink = subset(data.frame(dat), bat_species==species & virus_genus == genus)
  
  
  #select columns of relevance 
  df.slim <- dplyr::select(df.shrink, sampleid, antigen, MFI)
  
  #make wide
  df.wide = reshape(df.slim, direction="wide", idvar="sampleid",timevar = "antigen")

  
  #check for NAs in case some were missing a single antigen - replace with 1
  df.wide[is.na(df.wide)] <- 1
  names(df.wide)[2:ncol(df.wide)] <- gsub(pattern = "MFI.", replacement = "", x=names(df.wide)[2:ncol(df.wide)])
  
  if(genus=="Henipavirus"){
    df.wide <- dplyr::select(df.wide, "sampleid", "NiV", "HeV", "CedV", "GhV", "AngV", "MojV")  
  } else if(genus=="Orthoebolavirus"){
    df.wide <- dplyr::select(df.wide, "sampleid", "EboV", "BdbV", "BomV", "SudV", "RestV")  
  }
  
  
  
  p1 <-  ggcorrplot(cor(df.wide[,2:ncol(df.wide)], method = "spearman"), lab=T, lab_size = 3) #+ theme(plot.margin = unit(c(5.5,5.5,.1,.1), "points")) 
  
  
  hc.mat <- hclust(dist.mat)
  
  if(genus=="Orthoebolavirus"){
    hc.mat$order <- c(4,5,3,1,2)
    hc.mat$height = hc.mat$height-.5
    
    # p1add <- ggplot() + geom_dendro(hc.mat, size=1) + theme_bw() +
    #   scale_y_reverse() + scale_x_reverse()+ 
    #   theme(panel.grid = element_blank(), panel.border = element_blank(),
    #         axis.title = element_blank(), axis.text = element_blank(),
    #         plot.background = element_rect(fill="transparent"),
    #         panel.background = element_rect(fill="transparent"),
    #         axis.ticks = element_blank(), plot.margin = unit(c(0,3.3,0,1), "cm")) 
    
    p1add.left <- ggplot() + geom_dendro(hc.mat, size=1) + theme_bw() + coord_flip() +
      scale_y_reverse() + scale_x_reverse()+ 
      theme(panel.grid = element_blank(), panel.border = element_blank(),
            axis.title = element_blank(), axis.text = element_blank(),
            plot.background = element_rect(fill="transparent"),
            panel.background = element_rect(fill="transparent"),
            axis.ticks = element_blank(), plot.margin = unit(c(1.7,0, 2.8,.1), "cm")) 
    
    #p1.2 <- cowplot::plot_grid(p1,p1add, ncol = 1, nrow=2,  rel_heights = c(1,.2));  p1.2 
    p1.3 <- cowplot::plot_grid(p1add.left,p1, ncol = 2, nrow=1,  rel_widths = c(.2,1));  p1.3
  }else{
    
    #edit these a bit for ease of viewing
    hc.mat$height[2]<-hc.mat$height[2]-.3
    hc.mat$height[3]<-hc.mat$height[3]-.2
    hc.mat$height[4]<-hc.mat$height[4]-.1
    
    hc.mat$height =     hc.mat$height-.3
    
    # p1add <- ggplot() + geom_dendro(hc.mat, size=1) + theme_bw() +
    #   scale_y_reverse() + scale_x_reverse()+ 
    #   theme(panel.grid = element_blank(), panel.border = element_blank(),
    #         axis.title = element_blank(), axis.text = element_blank(),
    #         plot.background = element_rect(fill="transparent"),
    #         panel.background = element_rect(fill="transparent"),
    #         axis.ticks = element_blank(), plot.margin = unit(c(0,3.3,0,1), "cm")) 
    # 
    p1add.left <- ggplot() + geom_dendro(hc.mat, size=1) + theme_bw() + coord_flip() +
      scale_y_reverse() + scale_x_reverse()+ 
      theme(panel.grid = element_blank(), panel.border = element_blank(),
            axis.title = element_blank(), axis.text = element_blank(),
            plot.background = element_rect(fill="transparent"),
            panel.background = element_rect(fill="transparent"),
            axis.ticks = element_blank(), plot.margin = unit(c(1.7,0,2.7,.1), "cm")) 
    
    #p1.2 <- cowplot::plot_grid(p1,p1add, ncol = 1, nrow=2,  rel_heights = c(1,.2));  p1.2 
    p1.3 <- cowplot::plot_grid(p1add.left,p1, ncol = 2, nrow=1,  rel_widths = c(.2,1));  p1.3
    
  }
  
  
 
  
#   #plot
  p2 <- ggplot(df.wide,aes(get(ant1), get(ant2))) + theme_bw()  + theme(panel.grid = element_blank()) +
     stat_poly_line() + xlab(ant1) + ylab(ant2) +
     stat_poly_eq() + scale_y_log10() + scale_x_log10() +
     geom_point() +
     geom_smooth(method='lm')
#   
   p3 <- ggplot(df.wide,aes(get(ant1), get(ant3))) + theme_bw()  + theme(panel.grid = element_blank()) +
     stat_poly_line() + xlab(ant1) + ylab(ant3) +
     stat_poly_eq() + scale_y_log10() + scale_x_log10() +
     geom_point() +
     geom_smooth(method='lm')
#   
   p4 <- ggplot(df.wide,aes(get(ant1), get(ant4))) + theme_bw()  + theme(panel.grid = element_blank()) +
     stat_poly_line() + xlab(ant1) + ylab(ant4) +
     stat_poly_eq() + scale_y_log10() + scale_x_log10() +
     geom_point() +
     geom_smooth(method='lm')
   
   
   p5 <- ggplot(df.wide,aes(get(ant1), get(ant5))) + theme_bw()  + theme(panel.grid = element_blank()) +
     stat_poly_line() + xlab(ant1) + ylab(ant5) +
     stat_poly_eq() + scale_y_log10() + scale_x_log10() +
     geom_point() +
     geom_smooth(method='lm')
   #   
   
   p6 <- ggplot(df.wide,aes(get(ant2), get(ant3))) + theme_bw()  + theme(panel.grid = element_blank()) +
     stat_poly_line() + xlab(ant2) + ylab(ant3) +
     stat_poly_eq() + scale_y_log10() + scale_x_log10() +
     geom_point() +
     geom_smooth(method='lm')
   
   
   
   p.out <- cowplot::plot_grid(p1.3,p2,p3,p4,p5,p6, ncol = 3, nrow = 2, labels = c("A", "B", "C", "D", "E", "F"), label_size = 16)#, align = "hv", axis="tblr")

  family = unique(dat$virus_family[dat$antigen==ant1]) 
  ggsave(file = paste0(homewd,"/serology/figures/corr-plots-", gsub(" ", "-", species),"-",family, ".png"),
         plot = p.out,
         bg = "white",
         units="mm",  
         width=95, 
         height=55, 
         scale=3, 
         dpi=300)
  
  
  

}

plot.correlation(dat=lumdat.long, ant1 = "GhV", ant2 = "NiV", ant3 = "HeV", ant4 = "AngV", ant5 = "MojV", species = "Rousettus madagascariensis", genus = "Henipavirus", dist.mat = dist.henipa)
plot.correlation(dat=lumdat.long, ant1 = "GhV", ant2 = "NiV", ant3 = "HeV", ant4 = "AngV", ant5 = "MojV", species = "Eidolon dupreanum", genus = "Henipavirus", dist.mat = dist.henipa)
plot.correlation(dat=lumdat.long, ant1 = "GhV", ant2 = "NiV", ant3 = "HeV", ant4 = "AngV", ant5 = "MojV", species = "Pteropus rufus", genus = "Henipavirus", dist.mat = dist.henipa)


################################################################################
################################################################################
# Selection criteria suggests:

# Eidolon dupreanum : Henipavirus.
# Parallel analysis suggests maintaining 1 component, and PC1 has high loading 
# of both NiV and GhV. Eigenvalues are nonetheless close to 1 for PC2, PC3, even 
# PC4. PC2 is loaded almost entirely by AngV, while PC3 is dominated by HeV

# All 6 antigens have any MFI values > 1000 in at least one individual, 
# suggesting significant antigenic diversity for henipaviruses in this species.
# By Spearman, GhV is not highly correlated with any other antigen (r<.3), nor
# is AngV. As these two are important contributors to PC1 and PC2 (respectively),
# we will retain them in the dataset, especially since AngV is a known virus in 
# the system. NiV, HeV, CedV, and MojV are all correlated with one another,
# though HeV somewhat less so, and it does dominate PC3 independently. There
# are likely multiple other henipaviruses in this clade circulating. 


# Verdict: Seasonal analysis of seroprevalence to GhV, AngV

###############################################################
###############################################################
# Rousettus madagascariensis : Henipavirus.
# Parallel analysis suggests maintaining 1 component, and PC1 has high loading 
# of HeV and CedV and GhV. Eigenvalues are nonetheless close to 1 for PC2, which is 
# loaded almost entirely by AngV.  PC3 is dominated by NiV/MojV

# All 6 antigens have any MFI values > 1000 in multiple individuals, and all 
# are highly correlated (r> .4 and >.6 in most places). AngV and NiV are 
# somewhat less correlated with GhV but not hugely different -- we lack 
# justifiable cause though to act decisively beyond GhV. 2-3 clusters are 
# supported depending on the method.

# Verdict: Seasonal analysis of seroprevalence to GhV and AngV. 
# Clustering suggests the presence of a third HNV that is important
# in this system, though it does not clearly align with any of the antigens
# test. This could be a signature of AndV serology though we cannot say as this
# antigen was not successfully tested in this system.

# Verdict: just do seasonal analysis of seroprevalence to GhV and AngV in Rousettus.

###############################################################
###############################################################
# Pteropus rufus: Henipaviruses

# Not going to analyze as data were supportive only of a single, seronegative cluster
# and only 5 individuals had MFI values > 1000 for any HNV antigen. It is likely true 
# that there are circulating henipaviruses in P. rufus bats but we probably have not 
# identified and characterized these to sequence at this stage.

###############################################################
###############################################################
# Final decision:
# Eidolon: GhV, AngV, HeV, NiV (last two are most tenuous)
# Rousettus: GhV, RavV, BdbV 

# not going to explore Pteropus further.
