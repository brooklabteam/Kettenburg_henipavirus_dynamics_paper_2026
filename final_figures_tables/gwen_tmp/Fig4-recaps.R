rm(list=ls())
# these will be seasonal recapture data

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(mgcv)
library(lubridate)
library(lme4)
library(lmerTest)
library(ggh4x)
library(epitools)

#first record your wd for this repo (home) 
#homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
#batwd = "/Users/carabrook/Developer/Madagascar-Bat-Data"
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026"
batwd =  "/Users/gwenddolenkettenburg/Desktop/Developer/Madagascar-Bat-Data"

#load the recaptures
catch.dat <- read.csv(paste0(batwd,"/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)

head(catch.dat)
recap = subset(catch.dat, recap_new=="yes")
head(recap)
names(recap)
#slim to columns of interest
recap <- dplyr::select(recap, barcode_id, sampleid, roost_site, processing_date, bat_species, bat_age_class, bat_sex, mass_forearm_residual)
head(recap)

#and merge with shedding data too
shed.dat <- read.csv(file = paste0(homewd, "/virus_shedding/data/henipavirus_shedding_data.csv"),header = T, stringsAsFactors = F)
head(shed.dat)
names(shed.dat)

intersect(unique(recap$sampleid), unique(shed.dat$sampleid)) #164 shedding records for recaps

recap = subset(recap, processing_date < "2021-01-01")
length(recap$sampleid[recap$bat_species=="Eidolon dupreanum"]) #180
length(unique(recap$barcode_id[recap$bat_species=="Eidolon dupreanum"])) #88
length(recap$sampleid[recap$bat_species=="Rousettus madagascariensis"]) #115
length(unique(recap$barcode_id[recap$bat_species=="Rousettus madagascariensis"])) #57
length(recap$sampleid[recap$bat_species=="Pteropus rufus"]) #1
length(unique(recap$barcode_id[recap$bat_species=="Pteropus rufus"])) #1. only 1 sampling point for this bat, so cannot use it

#just take rousettus and eidolon
recap=subset(recap, bat_species!="Pteropus rufus")

#and merge with luminex data 

lum.dat <- read.csv(file = paste0(homewd, "/serology/data/mix_sub_data_75_90_96_long.csv"))
head(lum.dat)

lum.dat <- dplyr::select(lum.dat, sampleid, antigen, MFI, log10MFI, seropos, seropos_lci, seropos_uci)

recap <- merge(recap, lum.dat, by="sampleid")
head(recap)
recap$processing_date <- as.Date(recap$processing_date)



#let's just use the shedding status (not Ct) for now
shed.dat <- dplyr::select(shed.dat, sampleid, qPCR_AngV_multi_any_pos_ur, qPCR_GhV_multi_any_pos_ur)
head(shed.dat)

#and merge with recap
recap <- merge(recap, shed.dat, by="sampleid", all.x = T)
length(unique(recap$bat_species=="Eidolon dupreanum"))
recap$qPCR_AngV_multi_any_pos_ur[is.na(recap$qPCR_AngV_multi_any_pos_ur)] <- 2
recap$qPCR_GhV_multi_any_pos_ur[is.na(recap$qPCR_GhV_multi_any_pos_ur)] <- 2

length(unique(recap$sampleid[recap$qPCR_AngV_multi_any_pos_ur<2 |recap$qPCR_GhV_multi_any_pos_ur<2])) #101. not bad to start
head(recap)

#subset to the antigen/species of interest
recap$keep <- 0
recap$keep[recap$bat_species=="Eidolon dupreanum" & recap$antigen=="HeV"] <- 1
recap$keep[recap$bat_species=="Eidolon dupreanum" & recap$antigen=="GhV"] <- 1
recap$keep[recap$bat_species=="Eidolon dupreanum" & recap$antigen=="NiV"] <- 1
recap$keep[recap$bat_species=="Eidolon dupreanum" & recap$antigen=="AngV"] <- 1
recap$keep[recap$bat_species=="Rousettus madagascariensis" & recap$antigen=="GhV"] <- 1
recap$keep[recap$bat_species=="Rousettus madagascariensis" & recap$antigen=="AngV"] <- 1

recap = subset(recap, keep==1)
head(recap)

recap.sum <- ddply(recap, .(barcode_id, antigen), summarise, N=length(sampleid))
recap.sum <- ddply(recap.sum, .(barcode_id), summarise, Ncap=unique(N))

recap <- merge(recap, recap.sum, by="barcode_id")
head(recap)
recap$qPCR_GhV_multi_any_pos_ur[recap$qPCR_GhV_multi_any_pos_ur==1] <- "pos"
recap$qPCR_GhV_multi_any_pos_ur[recap$qPCR_GhV_multi_any_pos_ur==0] <- "neg"
recap$qPCR_GhV_multi_any_pos_ur[recap$qPCR_GhV_multi_any_pos_ur==2] <- "no data"

recap$qPCR_AngV_multi_any_pos_ur[recap$qPCR_AngV_multi_any_pos_ur==1] <- "pos"
recap$qPCR_AngV_multi_any_pos_ur[recap$qPCR_AngV_multi_any_pos_ur==0] <- "neg"
recap$qPCR_AngV_multi_any_pos_ur[recap$qPCR_AngV_multi_any_pos_ur==2] <- "no data"

recap$qPCR_GhV_multi_any_pos_ur <- factor(recap$qPCR_GhV_multi_any_pos_ur, levels=c("pos", "neg", "no data"))
recap$qPCR_AngV_multi_any_pos_ur <- factor(recap$qPCR_AngV_multi_any_pos_ur, levels=c("pos", "neg", "no data"))

Eid <- ggplot(data=subset(recap, bat_species=="Eidolon dupreanum" & Ncap>1),) + 
  facet_grid(antigen~.) + theme_bw() + 
  geom_point(aes(x=processing_date, y=MFI, color=barcode_id, shape = bat_sex), show.legend = F, size=3) +
  geom_line(aes(x=processing_date, y=MFI, color=barcode_id, group=barcode_id), show.legend = F, size=1); Eid

#or the shedding status
Eid_GhV_shed <- ggplot(data=subset(recap, bat_species=="Eidolon dupreanum" & Ncap>1 & antigen=="GhV"),) + 
  facet_grid(antigen~.) + theme_bw() + scale_color_manual(values=c("tomato", "cornflowerblue", "gray40"))+
  geom_point(aes(x=processing_date, y=barcode_id, color=qPCR_GhV_multi_any_pos_ur, shape = bat_sex),  size=3) +
  geom_line(aes(x=processing_date, y=barcode_id, color=qPCR_GhV_multi_any_pos_ur, group=barcode_id), show.legend = F, size=.8); Eid_GhV_shed


Eid_AngV_shed <- ggplot(data=subset(recap, bat_species=="Eidolon dupreanum" & Ncap>1 & antigen=="AngV"),) + 
  facet_grid(antigen~.) + theme_bw() + scale_color_manual(values=c("tomato", "cornflowerblue", "gray40"))+
  geom_point(aes(x=processing_date, y=barcode_id, color=qPCR_AngV_multi_any_pos_ur, shape = bat_sex),  size=3) +
  geom_line(aes(x=processing_date, y=barcode_id, color=qPCR_AngV_multi_any_pos_ur, group=barcode_id), show.legend = F, size=.8); Eid_AngV_shed


################################################################################
################################################################################
# and look up the old plots of recaps including prior status
rm(list=ls())

homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026"
batwd =  "/Users/gwenddolenkettenburg/Desktop/Developer/Madagascar-Bat-Data"
# homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
# batwd = "/Users/carabrook/Developer/Madagascar-Bat-Data"

setwd(homewd)
# re-merge the catching data

#load the recaptures
catch.dat <- read.csv(paste0(batwd,"/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)
head(catch.dat)
recap = subset(catch.dat, recap_new=="yes")
head(recap)
names(recap)
#slim to columns of interest
recap <- dplyr::select(recap, barcode_id, sampleid, roost_site, processing_date, bat_species, bat_age_class, bat_sex, mass_forearm_residual)
head(recap)

#just take rousettus and eidolon
recap=subset(recap, bat_species!="Pteropus rufus")
recap$processing_date <- as.Date(recap$processing_date)
min(recap$processing_date)#"2013-11-24"
max(recap$processing_date)#"2026-05-05"
recap <- subset(recap, processing_date < "2021-01-01")

#and merge with shedding data too
shed.dat <- read.csv(file = paste0(homewd, "/virus_shedding/data/henipavirus_shedding_data.csv"),header = T, stringsAsFactors = F)
head(shed.dat)
names(shed.dat)

intersect(unique(recap$sampleid), unique(shed.dat$sampleid)) #161 shedding records for recaps

#let's juts use the shedding status (not Ct) for now
shed.dat <- dplyr::select(shed.dat, sampleid, qPCR_AngV_multi_any_pos_ur, qPCR_GhV_multi_any_pos_ur)
head(shed.dat)

#and merge with recap
recap.shed <- merge(recap, shed.dat, by="sampleid")
recap.shed$qPCR_AngV_multi_any_pos_ur[is.na(recap.shed$qPCR_AngV_multi_any_pos_ur)] <- 2
recap.shed$qPCR_GhV_multi_any_pos_ur[is.na(recap.shed$qPCR_GhV_multi_any_pos_ur)] <- 2
length(unique(recap.shed$sampleid)) #133
length(unique(recap.shed$barcode_id)) #105

#and merge recaps with luminex data 
lum.dat <- read.csv(file = paste0(homewd, "/serology/data/mix_sub_data_75_90_96_long.csv"))
head(lum.dat)

lum.dat <- dplyr::select(lum.dat, sampleid, antigen, MFI, log10MFI, seropos, seropos_lci, seropos_uci)

recap.lum <- merge(recap, lum.dat, by="sampleid")
head(recap.lum)
recap.lum$processing_date <- as.Date(recap.lum$processing_date)
unique(recap.lum$barcode_id) #117

intersect(unique(recap.lum$barcode_id), unique(recap.shed$barcode_id)) #85 shedding bats records for recaps that also have serology
setdiff(unique(recap.lum$barcode_id), unique(recap.shed$barcode_id)) #34 bats with serology but no PCR
setdiff(unique(recap.shed$barcode_id), unique(recap.lum$barcode_id)) #20 bats with PCR but no serology 

#check dates
min(recap.shed$processing_date)#"2013-11-24"
max(recap.shed$processing_date)#2020-03-10"

head(recap.shed)

#we only want to compare with serology so just take after 2017
recap.shed <- subset(recap.shed, processing_date>"2017-01-01")
setdiff(unique(recap.shed$barcode_id), unique(recap.lum$barcode_id)) 
#now 6 bats with PCR but no serology 
setdiff(unique(recap.shed$sampleid), unique(recap.lum$sampleid)) #6 sampling events

head(recap.lum)

#slim recap.lum
#recap.lum <- dplyr::select(recap.lum, -(MFI), -(log10MFI))
head(recap.shed)
recap.shed$antigen <- "AngV"
recap.shed2 <- recap.shed
recap.shed2$antigen <- "GhV"
names(recap.shed)[9] <- "pos"
names(recap.shed2)[9] <- "pos"
recap.shed$pos #this is AngV pos
recap.shed2$pos <- recap.shed$qPCR_GhV_multi_any_pos_ur
recap.shed <- dplyr::select(recap.shed,-(qPCR_GhV_multi_any_pos_ur))
recap.shed2 <- dplyr::select(recap.shed2,-(qPCR_GhV_multi_any_pos_ur))
head(recap.shed)
head(recap.shed2)
recap.shed <- rbind(recap.shed, recap.shed2)
head(recap.shed)

#subset serology to the antigen/species of interest
recap.lum = subset(recap.lum, antigen=="GhV" | antigen=="AngV")
head(recap.lum)
recap.lum <- dplyr::select(recap.lum, -(seropos_lci), -(seropos_uci))
head(recap.lum)

#could bind them now but there will be repeated events
#let's give a different number system for PCR positives vs. seropositives
# you can be any of the following 8 categories:
# seropos + pcr pos
# seropos + pcr neg
# seropos + no pcr data
# seroneg + pcr pos
# seroneg + pcr neg
# seroneg + no pcr data
# no sero data + pcr pos
# no sero data + pcr neg

# merge but keep all sampleids
names(recap.lum)

sero.mol.merge <- merge(recap.lum, recap.shed, by=c("sampleid", "barcode_id", "roost_site", "processing_date", "bat_species", "bat_age_class", "bat_sex", "mass_forearm_residual", "antigen"), all = TRUE)
head(sero.mol.merge)

#and give the summary column
sero.mol.merge$status <- NA

sero.mol.merge$status[sero.mol.merge$seropos==0 & sero.mol.merge$pos==0] <- "seroneg + PCR neg"
sero.mol.merge$status[sero.mol.merge$seropos==0 & sero.mol.merge$pos==1] <- "seroneg + PCR pos"
sero.mol.merge$status[sero.mol.merge$seropos==0 & is.na(sero.mol.merge$pos)] <- "seroneg + no PCR data"
sero.mol.merge$status[sero.mol.merge$seropos==1 & sero.mol.merge$pos==1] <- "seropos + PCR pos"
sero.mol.merge$status[sero.mol.merge$seropos==1 & sero.mol.merge$pos==0] <- "seropos + PCR neg"
sero.mol.merge$status[sero.mol.merge$seropos==1 & is.na(sero.mol.merge$pos)] <- "seropos + no PCR data"
sero.mol.merge$status[is.na(sero.mol.merge$seropos) & sero.mol.merge$pos==1] <- "no sero data + PCR pos"
sero.mol.merge$status[is.na(sero.mol.merge$seropos) & sero.mol.merge$pos==0] <- "no sero data + PCR neg"

unique(sero.mol.merge$status)
sero.mol.merge$status <- factor(sero.mol.merge$status, levels = c("seroneg + PCR neg", "seroneg + PCR pos", "seroneg + no PCR data", "seropos + PCR neg", "seropos + PCR pos", "seropos + no PCR data", "no sero data + PCR neg", "no sero data + PCR pos"))

# split by antigen and barcode and make sure there are multiple timepoints. 
# color as bats that stay neg, seroconvert, serorevert or do both
id.recap <- function(df){
  if(nrow(df)>1){
    df.new <- arrange(df, processing_date)
    df.new$class <- NA
    if(sum(as.numeric(as.character(df$seropos)), na.rm = T)==0){
      df.new$class <- "always seronegative" # all seroneg
    }else if(sum(as.numeric(as.character(df$seropos)), na.rm = T)==length(df$seropos)){
      df.new$class <- "always seropositive" #all seropos
    }else{
      vec.sero <- unique(diff(as.numeric(as.character(df$seropos))))
      if(sum(unique(vec.sero), na.rm = T)==0){
        df.new$class <- "seroconversion-seroreversion" #both
      }else if((sum(unique(match(vec.sero,-1)),na.rm=T)==0) & (sum(unique(match(vec.sero,1)),na.rm=T)==1)){
        df.new$class <- "seroconversion" #both
      }else if((sum(unique(match(vec.sero,1)),na.rm=T)==0) & (sum(unique(match(vec.sero,-1)),na.rm=T)==1)){
        df.new$class <- "seroreversion" #both
      }
    }
    df.new$class2 <- "has positives"
    if (unique(((df.new$class=="always seronegative") & (sum(as.numeric(as.character(df.new$pos)), na.rm = T)==0)))){
      df.new$class2 <- "always seronegative & no PCR pos" 
    }
    return(df.new)
  }
  
}
recap.list <- dlply(sero.mol.merge, .(antigen, barcode_id))

recap.sum <- data.table::rbindlist(lapply(recap.list, id.recap))

recap.sum$seropos <- as.factor(recap.sum$seropos)
recap.sum$pos <- as.factor(recap.sum$pos)

#and plot
recap.sum$class2 <- as.factor(recap.sum$class2)

unique(recap.sum$class)
recap.sum$seropos <- as.character(recap.sum$seropos)
recap.sum$seropos[recap.sum$seropos==0] <- "seronegative"
recap.sum$seropos[recap.sum$seropos==1] <- "seropositive"
recap.sum$seropos[is.na(recap.sum$seropos)] <- "no sero data"
recap.sum$seropos <- factor(recap.sum$seropos, levels=c("seronegative", "seropositive", "no sero data"))

recap.sum$pos <- as.character(recap.sum$pos)
recap.sum$pos[recap.sum$pos==0] <- "PCR negative"
recap.sum$pos[recap.sum$pos==1] <- "PCR positive"
recap.sum$pos[is.na(recap.sum$pos)] <- "no PCR data"
recap.sum$pos <- factor(recap.sum$pos, levels=c("PCR negative", "PCR positive", "no PCR data"))

#filz <- c("seroneg + PCR neg" = "deepskyblue3", "seroneg + PCR pos" = "orange3", "seroneg + no PCR data" = "deepskyblue", "seropos + PCR neg" = "darkseagreen3", "seropos + PCR pos" = "tomato", "seropos + no PCR data" = "forestgreen" , "no sero data + PCR neg"= "cyan3", "no sero data + PCR pos" = "goldenrod2")
#filz = c('seronegative' = "navy",'seropositive'= "violetred")
filz <- c("PCR negative" = "cornflowerblue", "PCR positive" = "tomato", "no PCR data" = "gray65")

colorz1 = c('always seronegative' = "navy", 
            "always seropositive" = "violetred",
            "seroconversion-seroreversion"="goldenrod2", 
            "seroreversion" = "lightseagreen",
            "seroconversion" = "orangered")

shapez = c('seronegative' = 21, 'seropositive' = 24, 'no sero data'=22)
alphaz <- c('has positives' = 1, 'always seronegative & no PCR pos'=.1)

recap.sum$antigen <- factor(recap.sum$antigen, levels=c("GhV", "AngV"))
recap.sum$class <- factor(recap.sum$class, levels = c("always seronegative", "always seropositive", "seroconversion","seroreversion", "seroconversion-seroreversion"))

Eid.recap <- ggplot(subset(recap.sum, bat_species=="Eidolon dupreanum")) + theme_bw()+ 
  geom_line(aes(x=processing_date, y=barcode_id, group=barcode_id, color=class, alpha=class2), size=1.5) + 
  geom_point(aes(x=processing_date, y=barcode_id, group=barcode_id, shape=seropos, fill=pos, alpha=class2), size=4, color="gray65") + 
  scale_shape_manual(name="Serostatus", values = shapez) + 
  scale_fill_manual(name = "PCR status", values = filz) +
  scale_color_manual(name= "Serological\ntrajectory", values = colorz1) + ylab("Individual bat\n\n") +
  facet_nested(bat_sex~antigen, scales = "free", space = "free") + scale_alpha_manual(values = alphaz) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),  
        axis.title.y = element_text(size=16),
        legend.box = "vertical",
        legend.position = "right",
        #panel.spacing.y = unit(.1, "lines"),
        #panel.spacing.x = unit(1, "lines"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        plot.margin = unit(c(.3,.3,.3,2),"lines"),
        strip.text.x = element_text(size=16),
        strip.text.y = element_text(size=16),
        panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
  scale_x_continuous(breaks = c(as.Date("2014-07-01"),as.Date("2015-07-01"), as.Date("2016-07-01"), as.Date("2017-07-01"), as.Date("2018-07-01"), as.Date("2019-07-01")),labels = c(2014,2015,2016,2017,2018,2019)) +
  guides(alpha="none",  shape=guide_legend(order=1, override.aes = list(size=3, color="black")), fill=guide_legend(order=2, override.aes = list(shape=21,size=3)), color=guide_legend(order=3, ncol=1));Eid.recap

#leg <- get_legend(Eid.recap)

####################################################################################
####################################################################################
# expl
head(recap.sum)
#unique(recap.sum$roost_site[recap.sum$bat_species=="Eidolon dupreanum"])
recap.sum$antigen <- factor(recap.sum$antigen, levels=c("GhV", "AngV", "AndV", "GhVlike", "NiV","HeV"))
recap.sum$simple_class <- recap.sum$class
head(recap.sum)
levels(recap.sum$simple_class)
recap.sum$simple_class<-as.character(recap.sum$simple_class)
recap.sum$simple_class[recap.sum$simple_class!="always seronegative"] <- "dynamic"
recap.sum$simple_class <- as.factor(recap.sum$simple_class)
head(recap.sum)

alphaz2 = c("always seronegative" = .3, "dynamic" = 1 )
#recap.sum$antigen <- factor(recap.sum$antigen, levels=c( "GhV", "AngV"))

recap.sum$era <- "2014-2016"
recap.sum$era[recap.sum$processing_date>as.Date("2018-01-01")] <- "2018-2020"

# add <- recap.sum[1,]
# add$era <- "2014-2016"
# add$processing_date <- as.Date("2017-01-01")
# add$MFI <- NA
# add$antigen <- "NiV"
# add2 <- add
# add2$antigen <- "HeV"
# add <- rbind(add, add2)
# add2 <- add
# add2$era <- "2018-2020"
# add2$processing_date <- as.Date("2017-01-02")
# 
# recap.sum <- rbind(recap.sum, add, add2)

#shapez = c("female" = 21, "male" =24)

#and include plot of the cutoff line by species/antigen/era
# cutoff.df <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 6), rep("Rousettus madagascariensis",3)), antigen=c("NiV", "HeV", "NiV", "HeV", "GhV", "AngV", "GhV", "BdbV", "RavV"), era=c(rep("2014-2016",2), rep("2018-2020",7)))
# cutoff.df$MFI = c(757,215,936,2851,1141,1763,2642,5647,6617)
# cutoff.df$antigen <- factor(cutoff.df$antigen, levels=c("NiV", "HeV", "GhV", "AngV", "BdbV", "RavV"))

cutoff.df <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 2)), 
                              antigen=c("GhV", "AngV"), 
                              era=c(rep("2014-2016",2), rep("2018-2020",2)))
cutoff.df$MFI = c(1141,1763)
cutoff.df$antigen <- factor(cutoff.df$antigen, levels=c("GhV", "AngV"))

cutoff.df <- rbind(cutoff.df, c("Rousettus madagascariensis", "GhV", "2018-2020", 2642 ))
cutoff.df$MFI <- as.numeric(cutoff.df$MFI)

EidMFIall <- ggplot(data=subset(recap.sum, bat_species=="Eidolon dupreanum" )) + 
  facet_nested(bat_sex~antigen, scales = "free") + theme_bw() +
  scale_shape_manual(name="Serostatus", values = shapez) + 
  scale_alpha_manual(values = alphaz2) + ylab("mean fluorescence intensity (MFI)")+
  theme(axis.text.y = element_text(size=14),   
        #axis.ticks.y=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        legend.box = "vertical",
        legend.position = "left",
        #panel.spacing.y = unit(.1, "lines"),
        #panel.spacing.x = unit(1, "lines"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        plot.margin = unit(c(2,.3,0,2),"lines"),
        axis.text.x = element_text(size=14),  
        strip.text.x =element_blank(),
        strip.text.y = element_text(size=16),#, face="italic"),
        panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
  geom_hline(data = subset(cutoff.df, bat_species=="Eidolon dupreanum"), aes(yintercept = MFI), color="red", linetype=2) + 
  geom_line(aes(x=processing_date, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1.5)+
  geom_point(aes(x=processing_date, y=MFI, fill=barcode_id, shape = seropos, alpha = simple_class), show.legend = F, size=4) +
  scale_x_continuous(breaks = c(as.Date("2014-07-01"),as.Date("2015-07-01"), as.Date("2016-07-01"), as.Date("2017-07-01"), as.Date("2018-07-01"), as.Date("2019-07-01")),labels = c(2014,2015,2016,2017,2018,2019)) +
  guides(alpha="none",  shape=guide_legend(order=1));EidMFIall


#and the line up
Fig4 <- cowplot::plot_grid(Eid.recap, EidMFIall, 
                           labels = c("A", "B"),
                           ncol=1, nrow=2, rel_heights = c(1,1), align = "v", axis = "lr",label_size = 23)+
  theme(plot.margin = unit(c(0,2,0,0),"lines")); Fig4


Fig4.1 <- cowplot::plot_grid(Fig4, leg, ncol=2, nrow = 1, rel_widths = c(1,.2), align = "h", axis = "tb")+
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = unit(c(0,2,0,0),"lines")); Fig4.1



homewd = "/Users/gwenddolenkettenburg/Desktop/Developer/Kettenburg_henipavirus_dynamics_paper_2026/final_figures_tables/gwen_tmp"
ggsave(file = paste0(homewd,"/Fig4_recaps.png"),
       plot = Fig4,
       bg="white",
       units="mm",  
       width=120, 
       height=140, 
       scale=2, 
       dpi=500)


#and the Rousettus version for the supplement


Rou.recap <- ggplot(subset(recap.sum, bat_species=="Rousettus madagascariensis")) + theme_bw()+ 
  geom_line(aes(x=processing_date, y=barcode_id, group=barcode_id, color=class, alpha=class2), size=1.5) + 
  geom_point(aes(x=processing_date, y=barcode_id, group=barcode_id, shape=seropos, fill=pos, alpha=class2), size=4, color="gray65") + 
  scale_shape_manual(name="Serostatus", values = shapez) + 
  scale_fill_manual(name = "PCR status", values = filz) +
  scale_color_manual(name= "Serological\ntrajectory", values = colorz1) + ylab("Individual bat\n\n") +
  facet_nested(bat_sex~antigen, scales = "free", space = "free") + scale_alpha_manual(values = alphaz) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),  
        axis.title.y = element_text(size=16),
        legend.box = "vertical",
        legend.position = "right",
        #panel.spacing.y = unit(.1, "lines"),
        #panel.spacing.x = unit(1, "lines"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        plot.margin = unit(c(.3,.3,.3,2),"lines"),
        strip.text.x = element_text(size=16),
        strip.text.y = element_text(size=16),
        panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
  scale_x_continuous(breaks = c(as.Date("2014-07-01"),as.Date("2015-07-01"), as.Date("2016-07-01"), as.Date("2017-07-01"), as.Date("2018-07-01"), as.Date("2019-07-01")),labels = c(2014,2015,2016,2017,2018,2019)) +
  guides(alpha="none",  shape=guide_legend(order=1, override.aes = list(size=3, color="black")), fill=guide_legend(order=2, override.aes = list(shape=21,size=3)), color=guide_legend(order=3, ncol=1));Rou.recap


RouMFIall <- ggplot(data=subset(recap.sum, bat_species=="Rousettus madagascariensis" )) + 
  facet_nested(bat_sex~antigen, scales = "free") + theme_bw() +
  scale_shape_manual(name="Serostatus", values = shapez) + 
  scale_alpha_manual(values = alphaz2) + ylab("mean fluorescence intensity (MFI)")+
  theme(axis.text.y = element_text(size=14),   
        #axis.ticks.y=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        legend.box = "vertical",
        legend.position = "left",
        #panel.spacing.y = unit(.1, "lines"),
        #panel.spacing.x = unit(1, "lines"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        plot.margin = unit(c(2,.3,0,2),"lines"),
        axis.text.x = element_text(size=14),  
        strip.text.x =element_blank(),
        strip.text.y = element_text(size=16),#, face="italic"),
        panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
  geom_hline(data = subset(cutoff.df, bat_species=="Rousettus madagascariensis"), aes(yintercept = MFI), color="red", linetype=2) + 
  geom_line(aes(x=processing_date, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1.5)+
  geom_point(aes(x=processing_date, y=MFI, fill=barcode_id, shape = seropos, alpha = simple_class), show.legend = F, size=4) +
  scale_x_continuous(breaks = c(as.Date("2014-07-01"),as.Date("2015-07-01"), as.Date("2016-07-01"), as.Date("2017-07-01"), as.Date("2018-07-01"), as.Date("2019-07-01")),labels = c(2014,2015,2016,2017,2018,2019)) +
  guides(alpha="none",  shape=guide_legend(order=1));RouMFIall


#and the line up
Fig4supp <- cowplot::plot_grid(Rou.recap, RouMFIall, 
                           labels = c("A", "B"),
                           ncol=1, nrow=2, rel_heights = c(1,1), align = "v", axis = "lr",label_size = 23)+
  theme(plot.margin = unit(c(0,2,0,0),"lines")); Fig4supp

# 
# Fig4.1supp <- cowplot::plot_grid(Fig4supp, leg, ncol=2, nrow = 1, rel_widths = c(1,.2), align = "h", axis = "tb")+
#   theme(plot.background = element_rect(fill = "white", colour = NA),
#         plot.margin = unit(c(0,2,0,0),"lines")); Fig4.1supp

ggsave(file = paste0(homewd,"/final_figures_tables/gwen_tmp/Fig4_supp_recaps_Rousettus.png"),
       plot = Fig4supp,
       bg="white",
       units="mm",  
       width=80, 
       height=140, 
       scale=2, 
       dpi=500)
