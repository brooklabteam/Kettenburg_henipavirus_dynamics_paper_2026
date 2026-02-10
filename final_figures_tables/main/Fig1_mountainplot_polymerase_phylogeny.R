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

#use this script to calculate cutoffs for the luminex


#first record your wd for this repo (home) 
#cara adds hers here too - just comment out whichever you are not using
homewd = "/Users/carabrook/Developer/luminex-serology-madabat"

#now load in the merged luminex data from the "data" folder

#make sure this is latest catching file
dat <- read.csv(paste0(homewd, "/data/lum_catch_merge.csv"), header = T, stringsAsFactors = F)
head(dat)


length(unique(dat$sampleid))
dat$sampleid[duplicated(dat$sampleid)]

test <- ddply(dat, .(bat_species, bat_sex), summarise, Ntot=length(unique(sampleid)))
test

#now for cutoffs, just slim to antigen data only and convert to long format
lumdat <- dplyr::select(dat, sampleid, bat_species, names(dat)[2:15])
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
lumdat.long$virus_family <- "Filovirus"
lumdat.long$virus_family[lumdat.long$antigen== "NiV"| lumdat.long$antigen=="HeV"| lumdat.long$antigen=="CedV"| lumdat.long$antigen=="GhV"| lumdat.long$antigen=="AngV"| lumdat.long$antigen=="MojV"] <- "Paramyxovirus"
unique(lumdat.long$virus_family)

#and genus information
lumdat.long$virus_genus <- "Henipavirus"
lumdat.long$virus_genus[lumdat.long$antigen=="EboV" |lumdat.long$antigen=="BdbV" | lumdat.long$antigen=="BomV" |lumdat.long$antigen=="SudV" | lumdat.long$antigen=="RestV"] <- "Orthoebolavirus"
lumdat.long$virus_genus[lumdat.long$antigen=="LloV"] <- "Cuevavirus"
lumdat.long$virus_genus[lumdat.long$antigen=="RavV"] <- "Orthomarburgvirus"
lumdat.long$virus_genus[lumdat.long$antigen=="MlaV"] <- "Dianlovirus"
unique(lumdat.long$virus_genus)

#factor antigens by relatedness
lumdat.long$antigen <- factor(lumdat.long$antigen, levels=c( "NiV", "HeV","CedV", "GhV", "AngV", "MojV",
                                                            "EboV", "BdbV", "BomV",  "SudV", "RestV",   "LloV",  "MlaV", "RavV"))

lumdat.long$bat_species <- factor(lumdat.long$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))
 
lumdat.long.plot <- lumdat.long

lumdat.long.plot$virus_genus[lumdat.long.plot$virus_genus=="Orthomarburgvirus"] <- "Ortho-\nmarburgvirus"
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

ggsave(file = paste0(homewd,"/figures/FigS1.png"),
       plot = p1,
       units="mm",  
       width=118, 
       height=130, 
       scale=3, 
       dpi=300)




# and also add a mountain plot for additional support
# here, with antigens arranged according to phylogenetic relatedness... 
# makes it a bit easier to see the relatedness between bats
# first pass is all antigens in order of true relatedness

# load aa alignment for henipavirus, filovirus, and orthoebolavirus
henipa.align <- read.alignment(file = paste0(homewd,"/data/phylo/Henipavirus-Alignment.fasta"), format = "fasta")
# replace names
new.names <- strsplit(henipa.align$nam,"_")
length.names <- lapply(new.names, length)
new.names <- mapply(FUN='[', new.names, length.names)
henipa.align$nam <- new.names
filo.align <- read.alignment(file = paste0(homewd,"/data/phylo/filovirus-Alignment.fasta"), format = "fasta")
new.names <- strsplit(filo.align$nam,"_")
length.names <- lapply(new.names, length)
new.names <- mapply(FUN='[', new.names, length.names)
filo.align$nam <- new.names
ebola.align <- read.alignment(file = paste0(homewd,"/data/phylo/Orthoebolavirus-Alignment.fasta"), format = "fasta")
new.names <- strsplit(ebola.align$nam,"_")
length.names <- lapply(new.names, length)
new.names <- mapply(FUN='[', new.names, length.names)
ebola.align$nam <- new.names

# make distance matrix
dist.henipa <- dist.alignment(henipa.align)
dist.ebola <-dist.alignment(ebola.align)
dist.filo <- dist.alignment(filo.align)


show_col(hue_pal()(length(unique(lumdat.long$antigen))))

colorz = c('AngV' = "#F8766D", 'HeV' ="#A58AFF", 
           'NiV' = "#DB72FB", 'GhV' = "#00C094", 
           'CedV'= "#53B400", "MojV" = "#619CFF",
           'EboV' = "#06A4FF", 'BdbV' = "#00BFC4",
           'BomV' = "#00C094", 'SudV' = "#00BC56", 
           'RestV' = "#53B400", 'LloV'=  "#A58AFF",
           'RavV' = "#E38900",'MlaV' = "#FF66A8")


#and put them together
image.df <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"), antigen=rep("EboV", 3), MFI=rep(21000,3), image=c(paste0(homewd,"/data/images/eidolon_dupreanum_headshot.png"),paste0(homewd,"/data/images/rousettus_madagascariensis_headshot.png"),paste0(homewd,"/data/images/pteropus_rufus_headshot.png")))
image.df$bat_species <- factor(image.df$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))
image.df$virus_family="Filovirus"

image.df2 <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"), antigen=rep("NiV", 3), MFI=rep(21000,3), image=c(paste0(homewd,"/data/images/eidolon_dupreanum_headshot.png"),paste0(homewd,"/data/images/rousettus_madagascariensis_headshot.png"),paste0(homewd,"/data/images/pteropus_rufus_headshot.png")))
image.df2$bat_species <- factor(image.df2$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))
image.df2$virus_family="Paramyxovirus"
image.df <- rbind(image.df, image.df2)

lumdat.long$bat_species <- factor(lumdat.long$bat_species,  levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))

#some cutoffs were established by mixture modeling
cutoff.df <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 4), rep("Rousettus madagascariensis",3)), antigen=c("NiV", "HeV", "GhV", "AngV", "GhV", "BdbV", "RavV"))
cutoff.df$virus_family = "Filovirus"
cutoff.df$virus_family[cutoff.df$antigen=="NiV"| cutoff.df$antigen== "HeV"| cutoff.df$antigen=="CedV"| cutoff.df$antigen== "GhV"| cutoff.df$antigen== "AngV"| cutoff.df$antigen== "MojV"] <- "Paramyxovirus"
cutoff.df$MFI = c(936,2851,1141,1763,2642,5647,6617)
mean.cut <- ddply(cutoff.df,.(bat_species, virus_family), summarise, meanMFI = mean(MFI))

plot.cut <- cbind.data.frame(bat_species=rep(c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"),each=2),
                             virus_family=rep(c("Filovirus", "Paramyxovirus"), 3),
                             MFI= c(rep(1672,2), 2642, 6132, rep(1672,2))) #might come back and mixture model all of these


plot.cut$bat_species <- factor(plot.cut$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))


p2 <- ggplot(lumdat.long) + coord_cartesian(ylim = c(0,27000)) +
      facet_nested(bat_species~virus_family, scales = "free", space = "free") +
      geom_line(aes(x=antigen, y=MFI, group=sampleid), show.legend=F) + theme_bw()+
      geom_point(aes(x=antigen, y=MFI, color=antigen), show.legend=F) + 
      scale_color_manual(values=colorz) +
      geom_hline(data = plot.cut, aes(yintercept = MFI), color="red", linetype=2, size=.8) +
      theme(axis.title.x = element_blank(), panel.grid = element_blank(),
            axis.text.x = element_text(size=16), axis.title.y = element_text(size=16), 
            axis.text.y = element_text(size=13), strip.background = element_rect(fill="white"),
            strip.text.y = element_text(size=16, face="italic"), strip.text.x = element_text(size=16)) +
      geom_image(data=image.df, aes(x=antigen, y=MFI, image = image), size=.5, nudge_x = .6); p2

hc.henipa <- hclust(dist.henipa)

#edit these a bit for ease of viewing
hc.henipa$height[2]<-hc.henipa$height[2]-.3
hc.henipa$height[3]<-hc.henipa$height[3]-.2
hc.henipa$height[4]<-hc.henipa$height[4]-.1



p.h.1 <- ggplot() + geom_dendro(hc.henipa, size=1) + theme_bw() +
      scale_y_reverse() + scale_x_reverse()+ 
      theme(panel.grid = element_blank(), panel.border = element_blank(),
            axis.title = element_blank(), axis.text = element_blank(),
            plot.background = element_rect(fill="transparent"),
            panel.background = element_rect(fill="transparent"),
            axis.ticks = element_blank(), plot.margin = unit(c(0,.2,.5,2.5), "cm")) 
p.h.2 <- ggplot() + geom_dendro(hc.henipa, size=1) + theme_bw() +
  scale_y_reverse() + scale_x_reverse()+
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        axis.ticks = element_blank(), plot.margin = unit(c(0,.6,.5,.8), "cm")) 

p.h.3 <- ggplot() + geom_dendro(hc.henipa, size=1) + theme_bw() +
  scale_y_reverse() + scale_x_reverse()+
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        axis.ticks = element_blank(), plot.margin = unit(c(0,1.2,.5,.3), "cm")) 

p.h.row <- cowplot::plot_grid(p.h.1,p.h.2,p.h.3, ncol = 3, nrow = 1, rel_widths = c(1.15,1,1)) + theme(plot.margin = unit(c(0,.5,.5,1), "cm")) 

p2b <- ggplot(subset(lumdat.long, virus_genus=="Henipavirus")) + coord_cartesian(ylim = c(0,27000)) +
  facet_nested(virus_genus~bat_species, scales = "free", space = "free") +
  geom_line(aes(x=antigen, y=MFI, group=sampleid), show.legend=F) + theme_bw()+
  geom_point(aes(x=antigen, y=MFI, color=antigen), show.legend=F) + 
  geom_hline(data = subset(plot.cut, virus_family=="Paramyxovirus"), aes(yintercept = MFI), color="red", linetype=2, size=.8) +
  scale_color_manual(values=colorz) + ylab("mean fluorescence intensity (MFI)\n") +
  #geom_dendro(rowclus) +
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        axis.text.x = element_text(size=16), axis.title.y = element_text(size=16), 
        axis.text.y = element_text(size=16), strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        plot.margin = unit(c(.5,.5,.1,.3), "cm")) +
  geom_image(data=subset(image.df, virus_family=="Paramyxovirus"), aes(x=antigen, y=c(23000, 24000,23000), image = image), size=.4, nudge_x = .7); p2b
p2c <- cowplot::plot_grid(p2b, p.h.row, ncol = 1,nrow = 2, rel_heights = c(1,.2));p2c


p2bsave <- ggplot(subset(lumdat.long, virus_genus=="Henipavirus")) + coord_cartesian(ylim = c(0,27000)) +
  facet_nested(virus_genus~bat_species, scales = "free", space = "free") +
  geom_line(aes(x=antigen, y=MFI, group=sampleid), show.legend=F) + theme_bw()+
  geom_point(aes(x=antigen, y=MFI, color=antigen), show.legend=F) + 
  geom_hline(data = subset(plot.cut, virus_family=="Paramyxovirus"), aes(yintercept = MFI), color="red", linetype=2, size=.8) +
  scale_color_manual(values=colorz) + ylab("mean fluorescence intensity (MFI)\n") +
  #geom_dendro(rowclus) +
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        axis.text.x = element_text(size=16), axis.title.y = element_text(size=16), 
        axis.text.y = element_text(size=16), strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        plot.margin = unit(c(.5,.5,.1,.3), "cm")) +
  geom_image(data=subset(image.df, virus_family=="Paramyxovirus"), aes(x=antigen, y=c(23000, 24000,23000), image = image), size=.3, nudge_x = .7); p2bsave
p2csave <- cowplot::plot_grid(p2bsave, p.h.row, ncol = 1,nrow = 2, rel_heights = c(1,.2));p2csave


ggsave(file = paste0(homewd,"/figures/archive/mountain-plot-henipa.png"),
       bg="transparent",
       plot = p2csave,
       units="mm",  
       width=110, 
       height=65, 
       scale=3, 
       dpi=300)


#edit 
p.h.row <- cowplot::plot_grid(p.h.1,p.h.2,p.h.3, ncol = 3, nrow = 1, rel_widths = c(1.15,1,1)) + theme(plot.margin = unit(c(0,.5,.1,1), "cm")) 
p2c <- cowplot::plot_grid(p2b, p.h.row, ncol = 1,nrow = 2, rel_heights = c(1,.16));p2c

hc.filo <- hclust(dist.filo)


# #edit these a bit for ease of viewing
hc.filo$height[1]<-hc.filo$height[1]-.3 
hc.filo$height[2]<-hc.filo$height[2]-.2
hc.filo$height[3]<-hc.filo$height[3]-.2
hc.filo$height[4]<-hc.filo$height[4]-.1
hc.filo$height[6]<-hc.filo$height[6]-.1

#change order
hc.filo$order[4:8] <- c(4,5,3,1,2)

  



p.f.1 <- ggplot() + geom_dendro(hc.filo, size=1) + theme_bw() +
  scale_y_reverse() + scale_x_reverse() + 
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        axis.ticks = element_blank(), plot.margin = unit(c(0,.2,.5,1.5), "cm")) 
p.f.2 <- ggplot() + geom_dendro(hc.filo, size=1) + theme_bw() +
  scale_y_reverse() + scale_x_reverse()+
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        axis.ticks = element_blank(), plot.margin = unit(c(0,.6,.5,.8), "cm")) 

p.f.3 <- ggplot() + geom_dendro(hc.filo, size=1) + theme_bw() +
  scale_y_reverse() + scale_x_reverse()+
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        axis.ticks = element_blank(), plot.margin = unit(c(0,1.2,.5,.3), "cm")) 

p.f.row <- cowplot::plot_grid(p.f.1,p.f.2,p.f.3, ncol = 3, nrow = 1, rel_widths = c(1.05,1,1) ) + theme(plot.margin = unit(c(0,.5,.5,1), "cm")) 



p2d <- ggplot(subset(lumdat.long, virus_family=="Filovirus")) + coord_cartesian(ylim = c(0,27000)) +
  facet_nested(virus_family~bat_species, scales = "free", space = "free") +
  geom_line(aes(x=antigen, y=MFI, group=sampleid), show.legend=F) + theme_bw()+
  geom_point(aes(x=antigen, y=MFI, color=antigen), show.legend=F) + 
  scale_color_manual(values=colorz) +
  geom_hline(data = subset(plot.cut, virus_family=="Filovirus"), aes(yintercept = MFI), color="red", linetype=2, size=.8) +
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        axis.text.x = element_text(size=16), axis.title.y = element_text(size=16), 
        axis.text.y = element_text(size=13), strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        plot.margin = unit(c(.5,.5,.1,.3), "cm")) +
  geom_image(data=subset(image.df, virus_family=="Filovirus"), aes(x=antigen, y=24000, image = image), size=.3, nudge_x = 6);p2d

#replot

p2e <- cowplot::plot_grid(p2d, p.f.row, ncol = 1,nrow = 2, rel_heights = c(1,.2));p2e


ggsave(file = paste0(homewd,"/figures/archive/mountain-plot-filo.png"),
       plot = p2e,
       bg="transparent",
       units="mm",  
       width=120, 
       height=65, 
       scale=3.3, 
       dpi=300)



p2f <- ggplot(subset(lumdat.long, virus_family=="Filovirus")) + coord_cartesian(ylim = c(0,27000)) +
  facet_nested(virus_family~bat_species, scales = "free", space = "free") +
  geom_line(aes(x=antigen, y=MFI, group=sampleid), show.legend=F) + theme_bw()+
  geom_point(aes(x=antigen, y=MFI, color=antigen), show.legend=F) + 
  scale_color_manual(values=colorz) + ylab("mean fluorescence intensity (MFI)\n") +
  geom_hline(data = subset(plot.cut, virus_family=="Filovirus"), aes(yintercept = MFI), color="red", linetype=2, size=.8) +
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        axis.text.x = element_text(size=16), axis.title.y = element_text(size=16), 
        axis.text.y = element_text(size=16), strip.background.x  = element_blank(),
        strip.background.y = element_rect(fill="white"), strip.text.y = element_text(size=16),
        strip.text.x = element_blank(),plot.margin = unit(c(.1,.5,.1,.3), "cm")); p2f #+
  #geom_image(data=subset(image.df, virus_family=="Filovirus"), aes(x=antigen, y=24000, image = image), size=.4, nudge_x = .6);p2f


p.f.1 <- ggplot() + geom_dendro(hc.filo, size=1) + theme_bw() +
  scale_y_reverse() + scale_x_reverse() + 
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        axis.ticks = element_blank(), plot.margin = unit(c(0,.15,.5,2.1), "cm")) 

p.f.row <- cowplot::plot_grid(p.f.1,p.f.2,p.f.3, ncol = 3, nrow = 1, rel_widths = c(1.05,1,1) ) + theme(plot.margin = unit(c(0,.5,.5,1), "cm")) 

p2g <- cowplot::plot_grid(p2f, p.f.row, ncol = 1,nrow = 2, rel_heights = c(1,.2));p2g


p2out <- cowplot::plot_grid(p2c, p2g, nrow = 2, ncol = 1, rel_heights = c(1.05,1),align = "hv", axis="lr",  labels = c("A", "B"), label_size = 28)

ggsave(file = paste0(homewd,"/figures/Fig1.png"),
       plot = p2out,
       bg="white",
       units="mm",  
       width=140, 
       height=110, 
       scale=3, 
       dpi=300)


########################################################################################################
########################################################################################################
# Now, let's use dimensionality reduction to determine which antigens in the data we should retain:
