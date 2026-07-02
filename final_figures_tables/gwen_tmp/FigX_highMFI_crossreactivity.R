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
bat.high.df = subset(bat.df, MFI_high==1) #only 6 individual Pteropus with any MFI values >1000

p3 <- ggplot(subset(bat.high.df)) + #coord_cartesian(ylim = c(0,10000)) +
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



