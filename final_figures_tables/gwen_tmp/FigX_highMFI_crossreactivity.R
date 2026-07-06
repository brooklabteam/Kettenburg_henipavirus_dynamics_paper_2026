#exploratory high MFI cross reactivity

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
homewd = "/Users/gwenddolenkettenburg/Desktop/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
#homewd = "/Users/carabrook/Developer/luminex-serology-madabat"


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
bat.high.df = subset(bat.df, MFI_high==1) 

bat.high.df <- dplyr::select(bat.high.df, sampleid, bat_species, antigen, MFI)

bat.high.df.wide <- bat.high.df |> 
  pivot_wider(names_from = antigen, values_from = MFI)


numeric_cols <- bat.high.df.wide[, 3:8]
GhV_high <- bat.high.df.wide[max.col(numeric_cols, ties.method = "first") == 5, ]
GhV_high$highest_antigen="GhV"
GhV_high <- GhV_high %>% select(-GhV)
GhV_high.long <- melt(GhV_high)
head(GhV_high.long)
names(GhV_high.long)[names(GhV_high.long)=="variable"] <- "antigen"
names(GhV_high.long)[names(GhV_high.long)=="value"] <- "MFI"

HeV_high <- bat.high.df.wide[max.col(numeric_cols, ties.method = "first") == 1, ]
HeV_high$highest_antigen="HeV"
HeV_high <- HeV_high %>% select(-HeV)
HeV_high.long <- melt(HeV_high)
head(HeV_high.long)
names(HeV_high.long)[names(HeV_high.long)=="variable"] <- "antigen"
names(HeV_high.long)[names(HeV_high.long)=="value"] <- "MFI"

NiV_high <- bat.high.df.wide[max.col(numeric_cols, ties.method = "first") == 2, ]
NiV_high$highest_antigen="NiV"
NiV_high <- NiV_high %>% select(-NiV)
NiV_high.long <- melt(NiV_high)
head(NiV_high.long)
names(NiV_high.long)[names(NiV_high.long)=="variable"] <- "antigen"
names(NiV_high.long)[names(NiV_high.long)=="value"] <- "MFI"

CedV_high <- bat.high.df.wide[max.col(numeric_cols, ties.method = "first") == 3, ]
CedV_high$highest_antigen="CedV"
CedV_high <- CedV_high %>% select(-CedV)
CedV_high.long <- melt(CedV_high)
head(CedV_high.long)
names(CedV_high.long)[names(CedV_high.long)=="variable"] <- "antigen"
names(CedV_high.long)[names(CedV_high.long)=="value"] <- "MFI"

MojV_high <- bat.high.df.wide[max.col(numeric_cols, ties.method = "first") == 4, ]
MojV_high$highest_antigen="MojV"
MojV_high <- MojV_high %>% select(-MojV)
MojV_high.long <- melt(MojV_high)
head(MojV_high.long)
names(MojV_high.long)[names(MojV_high.long)=="variable"] <- "antigen"
names(MojV_high.long)[names(MojV_high.long)=="value"] <- "MFI"

AngV_high <- bat.high.df.wide[max.col(numeric_cols, ties.method = "first") == 6, ]
AngV_high$highest_antigen="AngV"
AngV_high <- AngV_high %>% select(-AngV)
AngV_high.long <- melt(AngV_high)
head(AngV_high.long)
names(AngV_high.long)[names(AngV_high.long)=="variable"] <- "antigen"
names(AngV_high.long)[names(AngV_high.long)=="value"] <- "MFI"

bat.high.df.merge <- rbind(GhV_high.long, HeV_high.long, NiV_high.long, CedV_high.long, 
                           MojV_high.long, AngV_high.long)
bat.high.df.merge$antigen<-factor(bat.high.df.merge$antigen, levels=c("HeV", "NiV", "GhV", "CedV", "MojV", "AngV"))

cutoff.df <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 6), rep("Rousettus madagascariensis",6), rep("Pteropus rufus",6)), antigen=c("HeV", "NiV", "GhV", "CedV", "MojV", "AngV"))
cutoff.df$MFI = c(2851,936,1141.25,777.5,676,1762.75,
                  964.5,2411.5, 2642, 757.75,1038.25, 871.75,
                  2033.5, 1687.5, 2231.5, 247, 1955, 198)
cutoff.df$antigen <- factor(cutoff.df$antigen, levels=c("HeV", "NiV", "GhV", "CedV", "MojV", "AngV"))
cutoff.df$highest_antigen <- factor(cutoff.df$antigen, levels=c("HeV", "NiV", "GhV", "CedV", "MojV", "AngV"))
cutoff.df$cutoff <- "mean"

cutoff.df2 <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 6), rep("Rousettus madagascariensis",6), rep("Pteropus rufus",6)), antigen=c("HeV", "NiV", "GhV", "CedV", "MojV", "AngV"))
cutoff.df2$MFI = c(2718,936,905.25,777.5,676,1667.25,
                   964.5,2383.5, 2642, 757.75,1038.25, 845,
                   2033.5, 1687.5, 2231.5, 247, 1955, 198)
cutoff.df2$antigen <- factor(cutoff.df2$antigen, levels=c("HeV", "NiV", "GhV", "CedV", "MojV", "AngV"))
cutoff.df2$highest_antigen <- factor(cutoff.df2$antigen, levels=c("HeV", "NiV", "GhV", "CedV", "MojV", "AngV"))
cutoff.df2$cutoff <- "lci"

cutoff.df3 <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 6), rep("Rousettus madagascariensis",6), rep("Pteropus rufus",6)), antigen=c("HeV", "NiV", "GhV", "CedV", "MojV", "AngV"))
cutoff.df3$MFI = c(3165.5,936,1619.25,777.5,676,2037.75,
                   1014,2484, 2642, 757.75,1038.25, 893,
                   2033.5, 1687.5, 2231.5, 247, 1955, 207)
cutoff.df3$antigen <- factor(cutoff.df3$antigen, levels=c("HeV", "NiV", "GhV", "CedV", "MojV", "AngV"))
cutoff.df3$highest_antigen <- factor(cutoff.df3$antigen, levels=c("HeV", "NiV", "GhV", "CedV", "MojV", "AngV"))
cutoff.df3$cutoff <- "uci"


cutoff.df <- rbind(cutoff.df, cutoff.df2, cutoff.df3)


p3 <- ggplot(subset(bat.high.df.merge)) +
  facet_nested(bat_species~highest_antigen, scales = "free", space = "free") +
  geom_line(aes(x=antigen, y=MFI, color=sampleid, group=sampleid), size=.8) + 
  #geom_jitter(aes(x=antigen, y=MFI, color=sampleid), width = .2,size=1) +
  #geom_violin(aes(x=antigen, y=MFI), quantiles = c(0,.25,.5,.75)) + 
  theme_bw()+
  geom_point(aes(x=antigen, y=MFI, color=sampleid), size=2) + 
  #scale_color_manual(name="Individual Bat\nSample ID", values = c("tomato", "#00BFC4")) +
  #geom_hline(aes(yintercept = 1000), color="red", linetype=2, size=.8) +
  geom_hline(data=subset(cutoff.df, cutoff=="mean"), aes(yintercept = MFI),  size=1, linetype=2)+
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size=16), axis.title.y = element_text(size=16), 
        axis.text.y = element_text(size=13), strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=16, face="italic"), strip.text.x = element_text(size=16)) +
  guides(color=guide_legend(ncol=2)); p3 

homewd = "/Users/gwenddolenkettenburg/Desktop/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
ggsave(file = paste0(homewd,"/final_figures_tables/gwen_tmp/FigX_highMFI_crossreactivity.png"),
       plot=p3,
       units="mm",  
       width=170, 
       height=80, 
       scale=3, 
       dpi=300)


