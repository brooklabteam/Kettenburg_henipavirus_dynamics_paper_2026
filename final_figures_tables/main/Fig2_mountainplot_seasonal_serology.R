#Script to make figure 2, which shows the mountain plots of henipavirus serology and the seasonal serology

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
library(ggimage)
library(ggtree)
library(ape)
library(ggnewscale)
library(gdata)
library(phylotools)
library(phylobase)
library(cowplot)

#use this script to calculate cutoffs for the luminex


#first record your wd for this repo (home) 
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/"

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
lumdat.long$virus_family <- "Paramyxovirus"
unique(lumdat.long$virus_family)

#and genus information
lumdat.long$virus_genus <- "Henipavirus"
lumdat.long$virus_genus[lumdat.long$antigen=="MojV"] <- "Parahenipavirus"
unique(lumdat.long$virus_genus)

#factor antigens by relatedness
lumdat.long$antigen <- factor(lumdat.long$antigen, levels=c( "NiV", "HeV","CedV", "GhV", "AngV", "MojV"))

lumdat.long$bat_species <- factor(lumdat.long$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))
 
lumdat.long.plot <- lumdat.long

#and heatmap MFI by species/antigen combination 
library(ggh4x)
p1 <- ggplot(lumdat.long.plot) + 
  geom_tile(aes(x=sampleid, y=antigen, fill=MFI)) + theme_bw()+
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(size=12), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        legend.key.width =  unit(c(1.5), "cm"),
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        legend.position = "bottom") +
  scale_fill_viridis_c(trans="log10") + facet_nested(virus_family+virus_genus~bat_species, scales = "free", space = "free");p1

# and also add a mountain plot for additional support
# here, with antigens arranged according to phylogenetic relatedness... 
# makes it a bit easier to see the relatedness between bats
# first pass is all antigens in order of true relatedness

# load aa alignment
henipa.align <- read.alignment(file = paste0(homewd,"serology/data/Henipavirus-Alignment.fasta"), format = "fasta")
# replace names
new.names <- strsplit(henipa.align$nam,"_")
length.names <- lapply(new.names, length)
new.names <- mapply(FUN='[', new.names, length.names)
henipa.align$nam <- new.names

# make distance matrix
dist.henipa <- dist.alignment(henipa.align)

colorz = c('AngV' = "#F8766D", 'HeV' ="#A58AFF", 
           'NiV' = "#DB72FB", 'GhV' = "#00C094", 
           'CedV'= "#53B400", "MojV" = "#619CFF")


#and put them together
image.df2 <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"), antigen=rep("NiV", 3), MFI=rep(21000,3), image=c(paste0(homewd,"/serology/images/eidolon_dupreanum_headshot.png"),paste0(homewd,"/serology/images/rousettus_madagascariensis_headshot.png"),paste0(homewd,"/serology/images/pteropus_rufus_headshot.png")))
image.df2$bat_species <- factor(image.df2$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))
image.df2$virus_family="Paramyxovirus"
image.df <- image.df2

lumdat.long$bat_species <- factor(lumdat.long$bat_species,  levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))

#some cutoffs were established by mixture modeling
cutoff.df <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 4), rep("Rousettus madagascariensis",3)), antigen=c("NiV", "HeV", "GhV", "AngV", "GhV", "BdbV", "RavV"))
cutoff.df$virus_family = "Paramyxovirus"
cutoff.df$MFI = c(936,2851,1141,1763,2642,5647,6617)
mean.cut <- ddply(cutoff.df,.(bat_species, virus_family), summarise, meanMFI = mean(MFI))

plot.cut <- cbind.data.frame(bat_species=rep(c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"),each=2),
                             virus_family=rep(c("Filovirus", "Paramyxovirus"), 3),
                             MFI= c(rep(1672,2), 2642, 6132, rep(1672,2))) #might come back and mixture model all of these


plot.cut$bat_species <- factor(plot.cut$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))

hc.henipa <- hclust(dist.henipa)

#edit these a bit for ease of viewing
hc.henipa$height[2]<-hc.henipa$height[2]-.3
hc.henipa$height[3]<-hc.henipa$height[3]-.2
hc.henipa$height[4]<-hc.henipa$height[4]-.1

p.h.1 <- ggplot() + geom_dendro(hc.henipa, size=1) + theme_bw() +
      scale_y_reverse() + scale_x_reverse()+ 
      theme(panel.grid = element_blank(), 
            panel.border = element_blank(),
            axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), plot.margin = unit(c(0,1,0,3.3), "cm")) 
p.h.1
p.h.2 <- ggplot() + geom_dendro(hc.henipa, size=1) + theme_bw() +
  scale_y_reverse() + scale_x_reverse()+
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), plot.margin = unit(c(0,1,0,1), "cm")) 
p.h.3 <- ggplot() + geom_dendro(hc.henipa, size=1) + theme_bw() +
  scale_y_reverse() + scale_x_reverse()+
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), plot.margin = unit(c(0,1,0,1), "cm")) 

p.h.row <- cowplot::plot_grid(p.h.1,p.h.2,p.h.3, ncol = 3, nrow = 1, rel_widths = c(1.22,1,1)) + theme(plot.margin = unit(c(0,0,0,0), "cm")) 

p2bsave <- ggplot(subset(lumdat.long)) + coord_cartesian(ylim = c(0,27000)) +
  facet_nested(.~bat_species, scales = "free", space = "free") +
  geom_line(aes(x=antigen, y=MFI, group=sampleid), show.legend=F) + theme_bw()+
  geom_point(aes(x=antigen, y=MFI, color=antigen), show.legend=F) + 
  geom_hline(data = subset(plot.cut, virus_family=="Paramyxovirus"), aes(yintercept = MFI), color="red", linetype=2, size=.8) +
  scale_color_manual(values=colorz) + ylab("mean fluorescence intensity (MFI)\n") +
  #geom_dendro(rowclus) +
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size=16), axis.title.y = element_text(size=16), 
        axis.text.y = element_text(size=16), 
        strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  geom_image(data=subset(image.df, virus_family=="Paramyxovirus"), aes(x=antigen, y=c(23000, 24000,23000), image = image), size=.4, nudge_x = .7); p2bsave
mountain <- cowplot::plot_grid(p2bsave, p.h.row, ncol = 1,nrow = 2, rel_heights = c(1,.2), rel_widths = c(1,1));mountain
mountain<-as.ggplot(mountain)

mountain <- mountain + theme(plot.margin = unit(c(0.1,0,0,2), "cm"))
mountain












###Put together final figure
Fig2<-plot_grid(mountain, NULL, tree1, labels=c("A","","B"),
                rel_widths=c(1,1,1), rel_heights = c(0.7,0.03,1),
                ncol=1,align="h", axis="l", label_size = 23) +
                theme(plot.background = element_rect(fill = "white"))
Fig2

Fig2<-plot_grid(tree1, NULL, mountain, labels=c("A","","B"),
                rel_widths=c(1,1,1), rel_heights = c(1,0.03,0.7),
                ncol=1,align="h", axis="l", label_size = 23) +
  theme(plot.background = element_rect(fill = "white"))
Fig2


#save trees
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/final_figures_tables/main/"
ggsave(file = paste0(homewd, "Fig2_mountainplot_seasonal_serology.png"),
       plot= Fig2,
       units="mm",  
       width=200, 
       height=160, 
       scale=2, 
       dpi=500)

