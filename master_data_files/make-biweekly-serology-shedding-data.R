rm(list=ls())


library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(mgcv)
library(reshape2)
library(Hmisc)
library(epitools)


################################################################################
################################################################################

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/"

setwd(homewd)
# plot of seasonality of antigens for analysis
# also GAMs of predictors of seropositivity and/or seroprevalence
################################################################################
################################################################################
## Section 1. Prep the data

# First, load the luminex data after determining who is seropositive/
# seronegative in the other scripts
lum.dat <- read.csv(file = paste0(homewd, "master_data_files/mix_sub_data_75_90_96.csv"))
head(lum.dat)
lum.dat = subset(lum.dat, bat_species=="Eidolon dupreanum" & antigen=="GhV"  |bat_species=="Eidolon dupreanum" & antigen=="AngV" | bat_species=="Rousettus madagascariensis" & antigen=="GhV")

lum.dat <- dplyr::select(lum.dat, sampleid, antigen, MFI, log10MFI, seropos, seropos_lci, seropos_uci)
## merge the luminex data with the catching data


# link this to the catching data to plot by sampling session
catch.dat <- read.csv(file = paste0(homewd, "/master_data_files/catching_data.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))


head(catch.dat)
names(catch.dat)

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)

#and merge all three datasets
merge.dat <- merge(lum.dat, catch.dat, by="sampleid", all.x = T)
head(merge.dat)


#change roost site to just Moramanga region or Ankarana
unique(merge.dat$roost_site)
merge.dat$roost_site[merge.dat$roost_site=="Angavobe_Edup" |merge.dat$roost_site=="Angavokely_Edup"| merge.dat$roost_site=="Maromizaha_Rmad"] <- "Moramanga"
merge.dat$roost_site[merge.dat$roost_site!="Moramanga"] <- "Ankarana"

# #also add in quantitative age data where available:
# age.dat <- read.csv(file = paste0(batwd, "/Age-Data/aged_bats.csv"), header = T, stringsAsFactors = F)
#head(age.dat)

# #slim to just the columns you want to add
# age.dat <- dplyr::select(age.dat, sampleid, age)
# 
# #and merge
# merge.dat <- merge(merge.dat, age.dat, by="sampleid", all.x = T)
# head(merge.dat)

#add biweek based on birthday
merge.dat$processing_date <- as.Date(merge.dat$processing_date)

#there is one funky serological datapoint from 2014 that should not be in this dataset-- remove
merge.dat <- subset(merge.dat, processing_date>"2018-01-01")

#different sites using this function
age.from.birth <- function(df, day1){
  #df$yday_catch = yday(df$processing_date)
  df$day1 <- yday(as.Date(paste0(year(df$processing_date), "-", day1)))
  df$age_days = NA
  df$age_days[df$doy>df$day1] <- df$doy[df$doy>df$day1] - df$day1[df$doy>df$day1]
  
  df$diffday <- 365-df$day1
  df$age_days[is.na(df$age_days)] <- df$doy[is.na(df$age_days)]  + df$diffday[is.na(df$age_days)]
  
  #df$partial_age <- df$age_days/365
  age_days <-  df$age_days
  return(age_days)
}


merge.dat$doy <- yday(merge.dat$processing_date)


#and correct for birth pulse
merge.dat$doy_corrected <- merge.dat$doy

#Rousettus - assume same birth date in Ankarana and Moramanga because 
#we don't have strong evidence to do anything else (birthdate of "12-11" from Andrianiaina et al 2022)
merge.dat$doy_corrected[merge.dat$bat_species=="Rousettus madagascariensis" & !is.na(merge.dat$doy_corrected)] <- age.from.birth(merge.dat[merge.dat$bat_species=="Rousettus madagascariensis"& !is.na(merge.dat$doy),],"12-11")

# Eidolon  in Moramanga (birthdate of "11-15") . from Andrianiaina et al 2025
merge.dat$doy_corrected[merge.dat$bat_species=="Eidolon dupreanum" & !is.na(merge.dat$doy_corrected) & merge.dat$roost_site=="Moramanga"] <- age.from.birth(merge.dat[merge.dat$bat_species=="Eidolon dupreanum"& !is.na(merge.dat$doy) & merge.dat$roost_site=="Moramanga",], "11-15")

# Eidolon in Ankarana (birthdate of "09-26")
merge.dat$doy_corrected[merge.dat$bat_species=="Eidolon dupreanum" & !is.na(merge.dat$doy_corrected) & merge.dat$roost_site=="Ankarana"] <- age.from.birth(merge.dat[merge.dat$bat_species=="Eidolon dupreanum"& !is.na(merge.dat$doy) & merge.dat$roost_site=="Ankarana",], "09-26")

head(merge.dat)

#no need to correct the 4 pregnant ankarana bats because our model allows the births paced out for 2x extra biweeks in both directions of birth date

#now add a biweek to the doy_corrected data
biwk.df <- cbind.data.frame(doy=1:365, biwk=trunc(seq(from=1, to=26.99999, length.out=365)))
ddply(biwk.df,.(biwk), summarise, N=length(doy)) #14


#replace doy with the corrected one
merge.dat <- dplyr::select(merge.dat, -(doy))
names(merge.dat)[names(merge.dat)=="doy_corrected"] <- "doy"
head(merge.dat)

#now get biweek from this
merge.dat <- merge(merge.dat, biwk.df, by="doy")
head(merge.dat)

merge.dat$bat_age_class[merge.dat$young_of_year=="yes"] <- "J"
subset(merge.dat, is.na(bat_age_class))

#merge.dat$bat_age_class[merge.dat$pupmom_pair=="yes" & merge.dat$bat_age_class!="J" & merge.dat$bat_age_class!="L" & !is.na(merge.dat$sampleid)] <- "L"

#and drop the columns you don't need
merge.dat <- dplyr::select(merge.dat, sampleid, barcode_id, biwk, antigen, MFI, seropos, seropos_lci, seropos_uci, roost_site, bat_species, processing_date, bat_sex, bat_age_class)

ggplot(data = merge.dat) + geom_violin(aes(x=biwk, y=seropos,group=biwk)) + geom_jitter(aes(x=biwk, y=seropos,color=bat_age_class), height = 0.2, width = 0.2) + facet_grid(bat_species~antigen)


#this is the serological data by biweek of the year
#save
write.csv(merge.dat,  file = paste0(homewd, "/master_data_files/biweekly-serology-data.csv"), row.names = F)

# next, make virus shedding data, also by biweek
rm(list=ls())


homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/"

setwd(homewd)


# First, load the shedding data 
shed.dat <- read.csv(file = paste0(homewd, "/virus_shedding/data/henipavirus_shedding_data.csv"))
head(shed.dat)
names(shed.dat)

shed.dat <- dplyr::select(shed.dat, sampleid, qPCR_AngV_multi_N_ur, qPCR_AngV_multi_L_ur, qPCR_GhV_multi_M_ur, qPCR_GhV_multi_L_ur)

shed.dat$AngV_pos <- 0
shed.dat$GhV_pos <- 0
shed.dat$AngV_pos[shed.dat$qPCR_AngV_multi_L==1 &shed.dat$qPCR_AngV_multi_N==1] <- 1
shed.dat$GhV_pos[shed.dat$qPCR_GhV_multi_M==1 &shed.dat$qPCR_GhV_multi_L==1] <- 1

sum(shed.dat$AngV_pos)#31
sum(shed.dat$GhV_pos)#119

shed.dat <- dplyr::select(shed.dat, sampleid,AngV_pos, GhV_pos)


# link this to the catching data to plot by sampling session
catch.dat <- read.csv(file = paste0(homewd, "/master_data_files/catching_data.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))

head(catch.dat)
names(catch.dat)

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)

#and merge all three datasets
setdiff(shed.dat$sampleid, catch.dat$sampleid)

merge.dat <- merge(shed.dat, catch.dat, by="sampleid")
head(merge.dat)


merge.dat$processing_date <- as.Date(merge.dat$processing_date)

#change roost site to just Moramanga region or Ankarana
unique(merge.dat$roost_site)
merge.dat$roost_site[merge.dat$roost_site=="Ambakoana_Pruf" | merge.dat$roost_site=="Marotsipohy_Pruf" | merge.dat$roost_site=="Analambotaka_Pruf" | merge.dat$roost_site=="Morarano_Pruf" | merge.dat$roost_site=="Marovitsika_Pruf" |  merge.dat$roost_site=="Mahialambo_Pruf"  | merge.dat$roost_site=="Lakato_Edup" | merge.dat$roost_site=="Angavobe_Edup" |merge.dat$roost_site=="Angavokely_Edup"| merge.dat$roost_site=="Maromizaha_Rmad"] <- "Moramanga"
merge.dat$roost_site[merge.dat$roost_site=="Mahabo_Edup" | merge.dat$roost_site=="Mahabo_Pruf"] <- "Mahabo"
merge.dat$roost_site[merge.dat$roost_site=="Makira_Rmad" | merge.dat$roost_site=="Makira_Pruf"] <- "Makira"
merge.dat$roost_site[merge.dat$roost_site!="Moramanga" & merge.dat$roost_site!="Makira" & merge.dat$roost_site!="Mahabo"] <- "Ankarana"

merge.dat$processing_date <- as.Date(merge.dat$processing_date)
merge.dat$doy <- yday(merge.dat$processing_date)

head(merge.dat)


# #also add in quantitative age data where available:
# age.dat <- read.csv(file = paste0(batwd, "/Age-Data/aged_bats.csv"), header = T, stringsAsFactors = F)
# head(age.dat)
# 
# #slim to just the columns you want to add
# age.dat <- dplyr::select(age.dat, sampleid, age)
# 
# #and merge
# merge.dat <- merge(merge.dat, age.dat, by="sampleid", all.x = T)
# head(merge.dat)

#add biweek based on birthday

#different sites using this function
age.from.birth <- function(df, day1){
  #df$yday_catch = yday(df$processing_date)
  df$day1 <- yday(as.Date(paste0(year(df$processing_date), "-", day1)))
  df$age_days = NA
  df$age_days[df$doy>df$day1] <- df$doy[df$doy>df$day1] - df$day1[df$doy>df$day1]
  
  df$diffday <- 365-df$day1
  df$age_days[is.na(df$age_days)] <- df$doy[is.na(df$age_days)]  + df$diffday[is.na(df$age_days)]
  
  #df$partial_age <- df$age_days/365
  age_days <-  df$age_days
  return(age_days)
}


#and correct doy for birth pulse
merge.dat$doy_corrected <- merge.dat$doy

#Rousettus - assume same birth date in Ankarana and Moramanga because 
#we don't have strong evidence to do anything else ("12-11")
merge.dat$doy_corrected[merge.dat$bat_species=="Rousettus madagascariensis" & !is.na(merge.dat$doy_corrected)] <- age.from.birth(merge.dat[merge.dat$bat_species=="Rousettus madagascariensis"& !is.na(merge.dat$doy),],"12-11")

# Eidolon  in Moramanga (birthdate of "11-15")
merge.dat$doy_corrected[merge.dat$bat_species=="Eidolon dupreanum" & !is.na(merge.dat$doy_corrected) & merge.dat$roost_site=="Moramanga"] <- age.from.birth(merge.dat[merge.dat$bat_species=="Eidolon dupreanum"& !is.na(merge.dat$doy) & merge.dat$roost_site=="Moramanga",], "11-15")

# Eidolon in Ankarana (birthdate of "09-26")
merge.dat$doy_corrected[merge.dat$bat_species=="Eidolon dupreanum" & !is.na(merge.dat$doy_corrected) & merge.dat$roost_site=="Ankarana"] <- age.from.birth(merge.dat[merge.dat$bat_species=="Eidolon dupreanum"& !is.na(merge.dat$doy) & merge.dat$roost_site=="Ankarana",], "09-26")

head(merge.dat)

#now add a biweek to the doy_corrected data
biwk.df <- cbind.data.frame(doy=1:365, biwk=trunc(seq(from=1, to=26.99999, length.out=365)))
ddply(biwk.df,.(biwk), summarise, N=length(doy)) #14


#replace doy with the corrected one
merge.dat <- dplyr::select(merge.dat, -(doy))
names(merge.dat)[names(merge.dat)=="doy_corrected"] <- "doy"
head(merge.dat)

#now get biweek from this
merge.dat <- merge(merge.dat, biwk.df, by="doy")
head(merge.dat)

merge.dat$bat_age_class[merge.dat$young_of_year=="yes"] <- "J"
subset(merge.dat, is.na(bat_age_class))

#just keep GhV for now
#and drop the columns you don't need

merge.dat.GhV <- dplyr::select(merge.dat, sampleid, barcode_id, biwk,  GhV_pos, roost_site, bat_species, processing_date, bat_sex, bat_age_class)
merge.dat.GhV$antigen <- "GhV"
names(merge.dat.GhV)[names(merge.dat.GhV)=="GhV_pos"] <- "shed_pos"

merge.dat.AngV <- dplyr::select(merge.dat, sampleid, barcode_id, biwk,  AngV_pos, roost_site, bat_species, processing_date, bat_sex, bat_age_class)
merge.dat.AngV$antigen <- "AngV"
names(merge.dat.AngV)[names(merge.dat.AngV)=="AngV_pos"] <- "shed_pos"

merge.dat <- rbind(merge.dat.GhV, merge.dat.AngV)

merge.dat = subset(merge.dat, bat_species !="Pteropus rufus")


ggplot(data = merge.dat) + geom_violin(aes(x=biwk, y=shed_pos,group=biwk)) + geom_jitter(aes(x=biwk, y=shed_pos,color=bat_age_class), height = 0.2, width = 0.2) + facet_grid(bat_species~antigen)

#this is the serological data by biweek of the year
#save
write.csv(merge.dat,  file = paste0(homewd, "/master_data_files/biweekly-shedding-data.csv"), row.names = F)
