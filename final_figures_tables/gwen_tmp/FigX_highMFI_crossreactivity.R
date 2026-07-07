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

#add seropos
dat2 <- read.csv(paste0(homewd, "/serology/data/mix_sub_data_75_90_96_all_wide.csv"), header = T, stringsAsFactors = F)
head(dat2)
dat2 <- dplyr::select(dat2, sampleid, seropos_NiV, seropos_HeV, seropos_CedV, seropos_GhV, seropos_AngV, seropos_MojV)

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


#make the seropos data long
dat2.long <- melt(dat2)
head(dat2.long)
names(dat2.long)[names(dat2.long)=="variable"] <- "antigen"
names(dat2.long)[names(dat2.long)=="value"] <- "seropos"
dat2.long$antigen <- gsub(pattern="seropos_", replacement = "", x=dat2.long$antigen, fixed = T)
dat2.long <- subset(dat2.long, !is.na(seropos))

#merge matching sampleids and antigens
bat.high.df.merge2 <- merge(bat.high.df.merge, dat2.long, by = c("sampleid", "antigen"))

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

# homewd = "/Users/gwenddolenkettenburg/Desktop/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
# ggsave(file = paste0(homewd,"/final_figures_tables/gwen_tmp/FigX_highMFI_crossreactivity.png"),
#        plot=p3,
#        units="mm",  
#        width=170, 
#        height=80, 
#        scale=3, 
#        dpi=300)

#now look at proportion of other antigens sample is seropositive to
colorz = c('AngV' = "#F8766D", 'HeV' ="#A58AFF", 
           'NiV' = "#DB72FB", 'GhV' = "#00C094", 
           'CedV'= "#53B400", "MojV" = "#619CFF")

#subset per bat species
ed.dat<-subset(bat.high.df.merge2, bat_species=="Eidolon dupreanum")
rm.dat<-subset(bat.high.df.merge2, bat_species=="Rousettus madagascariensis")
pr.dat<-subset(bat.high.df.merge2, bat_species=="Pteropus rufus")

ed.sampling<-ggplot(subset(ed.dat),aes(x = highest_antigen, y = seropos, fill = antigen, label=antigen)) +
  geom_col(
    position = "fill",
    stat = "identity") +
  #geom_text(size = 3, position = position_stack(vjust = 0.5))+
  facet_nested(bat_species~., scales="free", space="free",
               nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=colorz)+
  
  labs(
    x = "Highest antigen",
    y= "Proportion seropositive bats",
    fill="Antigen",
    title="")+
  theme_linedraw()+
  scale_y_continuous(expand=c(0,0))+
  guides(fill = guide_legend(ncol = 1))+
  #annotate("text", x="2018-02-15", y=5.2, label="*Lactation")+
  theme(plot.margin = margin(0, 1, 0, 30, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14),
        #axis.text.x = element_text(size=10),
        axis.text.x = element_text(angle = 90, size=14),
        legend.text = element_text(size=8, face="italic"),
        strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=14, face="italic", color="black"),
        strip.text.x = element_text(size=14, color="black"),
        legend.title = element_text(size=10),
        legend.position = "right");ed.sampling

rm.sampling<-ggplot(subset(rm.dat),aes(x = highest_antigen, y = seropos, fill = antigen, label=antigen)) +
  geom_col(
    position = "fill",
    stat = "identity") +
  #geom_text(size = 3, position = position_stack(vjust = 0.5))+
  facet_nested(bat_species~., scales="free", space="free",
               nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=colorz)+
  
  labs(
    x = "Highest antigen",
    y= "Proportion seropositive bats",
    fill="Antigen",
    title="")+
  theme_linedraw()+
  scale_y_continuous(expand=c(0,0))+
  guides(fill = guide_legend(ncol = 1))+
  #annotate("text", x="2018-02-15", y=5.2, label="*Lactation")+
  theme(plot.margin = margin(0, 1, 0, 30, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14),
        #axis.text.x = element_text(size=10),
        axis.text.x = element_text(angle = 90, size=14),
        legend.text = element_text(size=8, face="italic"),
        strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=14, face="italic", color="black"),
        strip.text.x = element_text(size=14, color="black"),
        legend.title = element_text(size=10),
        legend.position = "right");rm.sampling

pr.sampling<-ggplot(subset(pr.dat),aes(x = highest_antigen, y = seropos, fill = antigen, label=antigen)) +
  geom_col(
    position = "fill",
    stat = "identity") +
  #geom_text(size = 3, position = position_stack(vjust = 0.5))+
  facet_nested(bat_species~., scales="free", space="free",
               nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=colorz)+
  
  labs(
    x = "Highest antigen",
    y= "Proportion seropositive bats",
    fill="Antigen",
    title="")+
  theme_linedraw()+
  scale_y_continuous(expand=c(0,0))+
  guides(fill = guide_legend(ncol = 1))+
  #annotate("text", x="2018-02-15", y=5.2, label="*Lactation")+
  theme(plot.margin = margin(0, 1, 0, 30, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14),
        #axis.text.x = element_text(size=10),
        axis.text.x = element_text(angle = 90, size=14),
        legend.text = element_text(size=8, face="italic"),
        strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=14, face="italic", color="black"),
        strip.text.x = element_text(size=14, color="black"),
        legend.title = element_text(size=10),
        legend.position = "right");pr.sampling

plot4<-cowplot::plot_grid(ed.sampling, pr.sampling, rm.sampling, 
                          align="hv",axis='lr',ncol=1);plot4
plot<-as.ggplot(plot4)


#try making correlation plots
dat2 <- read.csv(paste0(homewd, "/serology/data/mix_sub_data_75_90_96_all_wide.csv"), header = T, stringsAsFactors = F)
head(dat2)
dat2 <- dplyr::select(dat2, sampleid, NiV,HeV,CedV,GhV,AngV,MojV, seropos_NiV, seropos_HeV, seropos_CedV, seropos_GhV, seropos_AngV, seropos_MojV)

# link this to the catching data
batwd = "/Users/gwenddolenkettenburg/Desktop/Developer/Madagascar-Bat-Data"
catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))
head(catch.dat)
# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)
#names(catch.dat)[names(catch.dat)=="net_date"] <- "processing_date"

merge.dat <- merge(dat2, catch.dat, by="sampleid")
head(merge.dat)

merge.dat$processing_date <- as.Date(merge.dat$processing_date)

#change roost site to just Moramanga region or Ankarana
unique(merge.dat$roost_site)
merge.dat$roost_site[merge.dat$roost_site=="Ambakoana_Pruf" | merge.dat$roost_site=="Marotsipohy_Pruf" | merge.dat$roost_site=="Analambotaka_Pruf" | merge.dat$roost_site=="Morarano_Pruf" | merge.dat$roost_site=="Marovitsika_Pruf" |  merge.dat$roost_site=="Mahialambo_Pruf"  | merge.dat$roost_site=="Lakato_Edup" | merge.dat$roost_site=="Angavobe_Edup" |merge.dat$roost_site=="Angavokely_Edup"| merge.dat$roost_site=="Maromizaha_Rmad"] <- "Moramanga"
merge.dat$roost_site[merge.dat$roost_site=="Mahabo_Edup" | merge.dat$roost_site=="Mahabo_Pruf"] <- "Mahabo"
merge.dat$roost_site[merge.dat$roost_site=="Makira_Rmad" | merge.dat$roost_site=="Makira_Pruf"] <- "Makira"
merge.dat$roost_site[merge.dat$roost_site!="Moramanga" & merge.dat$roost_site!="Makira" & merge.dat$roost_site!="Mahabo"] <- "Ankarana"
unique(merge.dat$roost_site)
merge.dat$processing_date <- as.Date(merge.dat$processing_date)
merge.dat$doy <- yday(merge.dat$processing_date)
head(merge.dat)

#first get midpoint date for sampling session
sampling.session <- ddply(merge.dat,.(sampling_session), summarise, min_date =min(processing_date), max_date=max(processing_date))
sampling.session$mid_date <- sampling.session$min_date + ((sampling.session$max_date-sampling.session$min_date)/2)
merge.dat <- merge(merge.dat, sampling.session, by="sampling_session")


#now get shedding prevalence
sero.prev <- ddply(merge.dat,.(bat_species, sampleid), summarise, 
                   NGhVpos = sum(seropos_GhV, na.rm=T), 
                   NAngVpos = sum(seropos_AngV, na.rm=T), 
                   NHeVpos = sum(seropos_HeV, na.rm=T), 
                   NNiVpos = sum(seropos_NiV, na.rm=T), 
                   NCedVpos = sum(seropos_CedV, na.rm=T), 
                   NMojVpos = sum(seropos_MojV, na.rm=T), 
                   Ntot=6, 
                   meanGhVMFI=mean(GhV, na.rm=T), 
                   meanAngVMFI=mean(AngV, na.rm=T),  
                   meanHeVMFI=mean(HeV, na.rm=T), 
                   meanNiVMFI=mean(NiV, na.rm=T),  
                   meanCedVMFI=mean(CedV, na.rm=T), 
                   meanMojVMFI=mean(MojV, na.rm=T),  
                   sdGhVMFI=sd(GhV, na.rm=T),
                   sdAngVMFI=sd(AngV, na.rm=T),
                   sdHeVMFI=sd(HeV, na.rm=T),
                   sdNiVMFI=sd(NiV, na.rm=T),
                   sdCedVMFI=sd(CedV, na.rm=T),
                   sdMojVMFI=sd(MojV, na.rm=T),
                   )
head(sero.prev)

# sero.prev.long <- melt(
#   sero.prev, 
#   id.vars = "sampleid", 
#   measure.vars = list(c("NGhVpos", "NAngVpos", "NHeVpos", "NNiVpos", "NCedVpos", "MojVpos"), 
#                       c("meanGhVMFI", "meanAngVMFI", "meanHeVMFI","meanNiVMFI","meanCedVMFI","meanMojVMFI"),
#                       c("sdGhVMFI", "sdAngVMFI", "sdHeVMFI","sdNiVMFI","sdCedVMFI","sdMojVMFI")),
#   value.name = c("Npos", "meanMFI", "sdMFI")
# )

sero.prev$GhVprev <- sero.prev$NGhVpos/sero.prev$Ntot
sero.prev$AngVprev <- sero.prev$NAngVpos/sero.prev$Ntot
sero.prev$HeVprev <- sero.prev$NHeVpos/sero.prev$Ntot
sero.prev$NiVprev <- sero.prev$NNiVpos/sero.prev$Ntot
sero.prev$CedVprev <- sero.prev$NCedVpos/sero.prev$Ntot
sero.prev$MojVprev <- sero.prev$NMojVpos/sero.prev$Ntot

sero.prev$NGhVneg <- sero.prev$Ntot - sero.prev$NGhVpos
sero.prev$NAngVneg <- sero.prev$Ntot - sero.prev$NAngVpos
sero.prev$NHeVneg <- sero.prev$Ntot - sero.prev$NHeVpos
sero.prev$NNiVneg <- sero.prev$Ntot - sero.prev$NNiVpos
sero.prev$NCedVneg <- sero.prev$Ntot - sero.prev$NCedVpos
sero.prev$NMojVneg <- sero.prev$Ntot - sero.prev$NMojVpos

sero.prev$GhVprevlci <- sero.prev$GhVprevuci <- sero.prev$AngVprevlci <- sero.prev$AngVprevuci <- 
  sero.prev$HeVprevlci <- sero.prev$HeVprevuci <- sero.prev$NiVprevlci <- sero.prev$NiVprevuci <- 
  sero.prev$CedVprevlci <- sero.prev$CedVprevuci <- sero.prev$MojVprevlci <- sero.prev$MojVprevuci <-NA

for(i in 1:nrow(sero.prev)){
  tmp <- binom.test(x=c(sero.prev$NGhVpos[i], sero.prev$NGhVneg[i]))$conf.int  
  tmp2 <- binom.test(x=c(sero.prev$NAngVpos[i], sero.prev$NAngVneg[i]))$conf.int  
  tmp3 <- binom.test(x=c(sero.prev$NHeVpos[i], sero.prev$NHeVneg[i]))$conf.int  
  tmp4 <- binom.test(x=c(sero.prev$NNiVpos[i], sero.prev$NNiVneg[i]))$conf.int  
  tmp5 <- binom.test(x=c(sero.prev$NCedVpos[i], sero.prev$NCedVneg[i]))$conf.int  
  tmp6 <- binom.test(x=c(sero.prev$NMojVpos[i], sero.prev$NMojVneg[i]))$conf.int  
  sero.prev$GhVprevlci[i] <- tmp[1]
  sero.prev$GhVprevuci[i] <- tmp[2]
  sero.prev$AngVprevlci[i] <- tmp2[1]
  sero.prev$AngVprevuci[i] <- tmp2[2]
  sero.prev$HeVprevlci[i] <- tmp3[1]
  sero.prev$HeVprevuci[i] <- tmp3[2]
  sero.prev$NiVprevlci[i] <- tmp4[1]
  sero.prev$NiVprevuci[i] <- tmp4[2]
  sero.prev$CedVprevlci[i] <- tmp5[1]
  sero.prev$CedVprevuci[i] <- tmp5[2]
  sero.prev$MojVprevlci[i] <- tmp6[1]
  sero.prev$MojVprevuci[i] <- tmp6[2]
}

head(sero.prev)
sero.prev$GhVMFI_uci <- sero.prev$meanGhVMFI+1.96*sero.prev$sdGhVMFI
sero.prev$GhVMFI_lci <- sero.prev$meanGhVMFI-1.96*sero.prev$sdGhVMFI
sero.prev$AngVMFI_uci <- sero.prev$meanAngVMFI+1.96*sero.prev$sdAngVMFI
sero.prev$AngVMFI_lci <- sero.prev$meanAngVMFI-1.96*sero.prev$sdAngVMFI
sero.prev$HeVMFI_uci <- sero.prev$meanHeVMFI+1.96*sero.prev$sdHeVMFI
sero.prev$HeVMFI_lci <- sero.prev$meanHeVMFI-1.96*sero.prev$sdHeVMFI
sero.prev$NiVMFI_uci <- sero.prev$meanNiVMFI+1.96*sero.prev$sdNiVMFI
sero.prev$NiVMFI_lci <- sero.prev$meanNiVMFI-1.96*sero.prev$sdNiVMFI
sero.prev$CedVMFI_uci <- sero.prev$meanCedVMFI+1.96*sero.prev$sdCedVMFI
sero.prev$CedVMFI_lci <- sero.prev$meanCedVMFI-1.96*sero.prev$sdCedVMFI
sero.prev$MojVMFI_uci <- sero.prev$meanMojVMFI+1.96*sero.prev$sdMojVMFI
sero.prev$MojVMFI_lci <- sero.prev$meanMojVMFI-1.96*sero.prev$sdMojVMFI



#subset and melt to make the sero.prev dataset long
part1 <- sero.prev[, c("sampleid","bat_species","GhVprev", "AngVprev","HeVprev","NiVprev","CedVprev","MojVprev")]
part1.long <- melt(part1)
names(part1.long)[names(part1.long)=="variable"] <- "antigen"
part1.long$antigen <- gsub(pattern="prev", replacement = "", x=part1.long$antigen, fixed = T)
names(part1.long)[names(part1.long)=="value"] <- "prev"

part2 <- sero.prev[, c("sampleid","bat_species","NGhVpos", "NAngVpos","NHeVpos","NNiVpos","NCedVpos","NMojVpos")]
part2.long <- melt(part2)
names(part2.long)[names(part2.long)=="variable"] <- "antigen"
part2.long$antigen <- gsub(pattern="pos", replacement = "", x=part2.long$antigen, fixed = T)
part2.long$antigen <- gsub(pattern="N", replacement = "", x=part2.long$antigen, fixed = T)
names(part2.long)[names(part2.long)=="value"] <- "pos"

part3 <- sero.prev[, c("sampleid","bat_species","meanGhVMFI", "meanAngVMFI","meanHeVMFI","meanNiVMFI","meanCedVMFI","meanMojVMFI")]
part3.long <- melt(part3)
names(part3.long)[names(part3.long)=="variable"] <- "antigen"
part3.long$antigen <- gsub(pattern="mean", replacement = "", x=part3.long$antigen, fixed = T)
part3.long$antigen <- gsub(pattern="MFI", replacement = "", x=part3.long$antigen, fixed = T)
names(part3.long)[names(part3.long)=="value"] <- "MFI"

part4 <- sero.prev[, c("sampleid","bat_species","GhVprevuci", "AngVprevuci","HeVprevuci","NiVprevuci","CedVprevuci","MojVprevuci")]
part4.long <- melt(part4)
names(part4.long)[names(part4.long)=="variable"] <- "antigen"
part4.long$antigen <- gsub(pattern="prevuci", replacement = "", x=part4.long$antigen, fixed = T)
names(part4.long)[names(part4.long)=="value"] <- "prevuci"

part5 <- sero.prev[, c("sampleid","bat_species","GhVprevlci", "AngVprevlci","HeVprevlci","NiVprevlci","CedVprevlci","MojVprevlci")]
part5.long <- melt(part5)
names(part5.long)[names(part5.long)=="variable"] <- "antigen"
part5.long$antigen <- gsub(pattern="prevlci", replacement = "", x=part5.long$antigen, fixed = T)
names(part5.long)[names(part5.long)=="value"] <- "prevlci"


sero.prev.long <- merge(part1.long, part2.long)
sero.prev.long <- merge(sero.prev.long, part3.long)
sero.prev.long <- merge(sero.prev.long, part4.long)
sero.prev.long <- merge(sero.prev.long, part5.long)


pMFIprev <- ggplot(data= subset(sero.prev.long, )) + 
  geom_point(aes(x=MFI, y=prev, color=antigen), show.legend = F) + theme_bw() +
  #geom_linerange(aes(x=MFI, ymin=prevlci, ymax=prevuci, color=antigen), linetype=3, alpha=.5, show.legend = F) +
  #geom_line(data=pred.df, aes(x=meanGhVMFI, y=prevalence)) +
  xlab("MFI") + ylab("Proportion of sample seropos to antigen") +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  scale_color_manual(name="Antigen", values = colorz) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=14), 
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.text = element_text(size=12),
        legend.position = "inside",
        legend.position.inside = c(.15,.73),
        legend.title = element_text(size=14));pMFIprev


