rm(list=ls())


library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(mgcv)
library(reshape2)
library(Hmisc)


################################################################################
################################################################################

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
batwd = "/Users/carabrook/Developer/Madagascar-Bat-Data"
################################################################################
################################################################################
## Section 1. Prep the data

# First, load the luminex data after determining who is seropositive/
# seronegative in the other scripts
lum.dat <- read.csv(file = paste0(homewd, "/master_data_files/mix_sub_data_75_90_96.csv"))
head(lum.dat)

lum.dat <- dplyr::select(lum.dat, sampleid, antigen, MFI, log10MFI, seropos, seropos_lci, seropos_uci)
## merge the luminex data with the catching data

# link this to the catching data to plot by sampling session
catch.dat <- read.csv(file = paste0(homewd, "/master_data_files/catching_data.csv"), header = T, stringsAsFactors = F)


head(catch.dat)
names(catch.dat)

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)

#and merge all three datasets
merge.dat <- merge(lum.dat, catch.dat, by="sampleid", all.x = T)
head(merge.dat)

#get the midpoint processing date per sampling session
merge.dat$processing_date <- as.Date(merge.dat$processing_date)

#there is one funky serological datapoint from 2014 that should not be in this dataset-- remove
merge.dat <- subset(merge.dat, processing_date>"2018-01-01")


#change roost site to just Moramanga region or Ankarana
unique(merge.dat$roost_site)
merge.dat$roost_site[merge.dat$roost_site=="Angavobe_Edup" |merge.dat$roost_site=="Angavokely_Edup"| merge.dat$roost_site=="Maromizaha_Rmad"] <- "Moramanga"
merge.dat$roost_site[merge.dat$roost_site!="Moramanga"] <- "Ankarana"

#now, correct doy year to represent the birth pulse for the 
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
#we don't have strong evidence to do anything else ("12-11")
merge.dat$doy_corrected[merge.dat$bat_species=="Rousettus madagascariensis" & !is.na(merge.dat$doy_corrected)] <- age.from.birth(merge.dat[merge.dat$bat_species=="Rousettus madagascariensis"& !is.na(merge.dat$doy),],"12-11")

#Pteropus - no pteropus bats in merged dataset (birthdate of "09-28")
#merge.dat$doy_corrected[merge.dat$bat_species=="Pteropus rufus" & !is.na(merge.dat$doy_corrected)] <- age.from.birth(merge.dat[merge.dat$bat_species=="Pteropus rufus"& !is.na(merge.dat$doy),], "09-28")

# Eidolon  in Moramanga (birthdate of "11-15")
merge.dat$doy_corrected[merge.dat$bat_species=="Eidolon dupreanum" & !is.na(merge.dat$doy_corrected) & merge.dat$roost_site=="Moramanga"] <- age.from.birth(merge.dat[merge.dat$bat_species=="Eidolon dupreanum"& !is.na(merge.dat$doy) & merge.dat$roost_site=="Moramanga",], "11-15")
# Eidolon in Ankarana (birthdate of "09-26")
merge.dat$doy_corrected[merge.dat$bat_species=="Eidolon dupreanum" & !is.na(merge.dat$doy_corrected) & merge.dat$roost_site=="Ankarana"] <- age.from.birth(merge.dat[merge.dat$bat_species=="Eidolon dupreanum"& !is.na(merge.dat$doy) & merge.dat$roost_site=="Ankarana",], "09-26")

head(merge.dat)


#also add in quantitative age data where available:
age.dat <- read.csv(file = paste0(batwd, "/Age-Data/aged_bats.csv"), header = T, stringsAsFactors = F)
head(age.dat)
unique(age.dat$notes)
#slim to just the columns you want to add
age.dat <- dplyr::select(age.dat, sampleid, age)

#and merge
merge.dat <- merge(merge.dat, age.dat, by="sampleid", all.x = T)
head(merge.dat)

#now slim to just moramanga region and just bats with age info
merge.dat = subset(merge.dat, !is.na(age))

#and get biweek of the year and attach to the data
biwk.df <- cbind.data.frame(doy=1:365, biwk=trunc(seq(from=1, to=26.99999, length.out=365)))
ddply(biwk.df,.(biwk), summarise, N=length(doy)) #14



#replace doy with the corrected one
merge.dat <- dplyr::select(merge.dat, -(doy))
names(merge.dat)[names(merge.dat)=="doy_corrected"] <- "doy"
#now get biweek from this

merge.dat <- merge(merge.dat, biwk.df, by="doy")
head(merge.dat)

#and slim to all you need 
merge.dat <- dplyr::select(merge.dat, roost_site, bat_species, processing_date, roost_site, biwk, sampleid, age, antigen, seropos, seropos_lci, seropos_uci, MFI, mass_forearm_residual)
head(merge.dat)
names(merge.dat)[names(merge.dat)=="seropos"] <- "prev"
names(merge.dat)[names(merge.dat)=="seropos_lci"] <- "prev_lci"
names(merge.dat)[names(merge.dat)=="seropos_uci"] <- "prev_uci"
names(merge.dat)[names(merge.dat)=="bat_species"] <- "species"
#and save
write.csv(merge.dat, file = paste0(homewd,"/master_data_files/age-seroprevalence-data.csv"), row.names = F)




