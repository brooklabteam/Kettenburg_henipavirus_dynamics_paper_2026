rm(list=ls())


library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(mgcv)
library(reshape2)
library(Hmisc)
library(exactci)

# Does shedding one HNV make you more likely to shed multiple???
# Answer both cumulatively and across species appears to be no, not really
# Some trends that testing + for AngV means that you are likely to also be positive for GhV 
# but none are significant

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/"
batwd = "/Users/carabrook/Developer/Madagascar-Bat-Data"

# First, load the shedding data 
shed.dat <- read.csv(file = paste0(homewd, "/virus_shedding/data/henipavirus_shedding_data.csv"))
head(shed.dat)
names(shed.dat)


shed.dat <- dplyr::select(shed.dat, sampleid, GhV_shedding_ur, AngV_shedding_ur, AndV_shedding_ur, GhV_shedding_ct_ur, AngV_shedding_ct_ur, AndV_shedding_ct_ur)

names(shed.dat)[names(shed.dat) == "GhV_shedding_ur"] <- "GhV_pos"
names(shed.dat)[names(shed.dat) == "AngV_shedding_ur"] <- "AngV_pos"
names(shed.dat)[names(shed.dat) == "AndV_shedding_ur"] <- "AndV_pos"

names(shed.dat)[names(shed.dat) == "GhV_shedding_ct_ur"] <- "GhV_ct"
names(shed.dat)[names(shed.dat) == "AngV_shedding_ct_ur"] <- "AngV_ct"
names(shed.dat)[names(shed.dat) == "AndV_shedding_ct_ur"] <- "AndV_ct"

head(shed.dat)
shed.dat$num_HNV_pos <- rowSums(cbind(shed.dat$GhV_pos, shed.dat$AngV_pos, shed.dat$AndV_pos))

# link this to the catching data to plot by sampling session
catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))

head(catch.dat)
#names(catch.dat)
#catch.dat$net_date[!is.na(catch.dat$net_date)] <- as.Date(catch.dat$net_date[!is.na(catch.dat$net_date)], format = "%m/%d/%Y")
#catch.dat$net_date[is.na(catch.dat$net_date)] <- catch.dat$processing_date[is.na(catch.dat$net_date)]

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)
#names(catch.dat)[names(catch.dat)=="net_date"] <- "processing_date"

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
unique(merge.dat$roost_site)
merge.dat$processing_date <- as.Date(merge.dat$processing_date)
merge.dat$doy <- yday(merge.dat$processing_date)

head(merge.dat)

#how many pos for each virus?
sum(merge.dat$GhV_pos)#130
sum(merge.dat$AngV_pos)#45
sum(merge.dat$AndV_pos)#23

#and prevalence of each by species - roughly equivalent across species for GhV (~10%) and AngV (~4%). AndV is low (7.5% in P. ruf but only 1.1% in E. dup and 0.6% in R. mad)
sum(merge.dat$GhV_pos[merge.dat$bat_species=="Pteropus rufus"])/length(merge.dat$GhV_pos[merge.dat$bat_species=="Pteropus rufus"])#0.1173709
sum(merge.dat$GhV_pos[merge.dat$bat_species=="Eidolon dupreanum"])/length(merge.dat$GhV_pos[merge.dat$bat_species=="Eidolon dupreanum"])#0.1588367
sum(merge.dat$GhV_pos[merge.dat$bat_species=="Rousettus madagascariensis"])/length(merge.dat$GhV_pos[merge.dat$bat_species=="Rousettus madagascariensis"])#0.1005917


sum(merge.dat$AngV_pos[merge.dat$bat_species=="Pteropus rufus"])/length(merge.dat$AngV_pos[merge.dat$bat_species=="Pteropus rufus"])#0.04694836
sum(merge.dat$AngV_pos[merge.dat$bat_species=="Eidolon dupreanum"])/length(merge.dat$AngV_pos[merge.dat$bat_species=="Eidolon dupreanum"])#0.049217
sum(merge.dat$AngV_pos[merge.dat$bat_species=="Rousettus madagascariensis"])/length(merge.dat$AngV_pos[merge.dat$bat_species=="Rousettus madagascariensis"])#0.03846154


sum(merge.dat$AndV_pos[merge.dat$bat_species=="Pteropus rufus"])/length(merge.dat$AndV_pos[merge.dat$bat_species=="Pteropus rufus"])#0.07511737
sum(merge.dat$AndV_pos[merge.dat$bat_species=="Eidolon dupreanum"])/length(merge.dat$AndV_pos[merge.dat$bat_species=="Eidolon dupreanum"])#0.01118568
sum(merge.dat$AndV_pos[merge.dat$bat_species=="Rousettus madagascariensis"])/length(merge.dat$AndV_pos[merge.dat$bat_species=="Rousettus madagascariensis"])#0.00591716


coinfections = subset(merge.dat, num_HNV_pos >1)
coinfections 
#11 bats out of 998 were coinfected with 2 HNVs (no bats were coinfected with all three):
# 2 Pteropus, 3 Rousettus, 6 Eidolon
# The two coinfected Pteropus were coinfected with GhV and AndV
# One Eidolon was coinfected with AngV and AndV
# All 3 Rousettus and 5 Eidolon were coinfected with GhV and AngV

# 2 of the Rousettus were from Ankarana and one from Maromizaha;
# 1 of the Pteropus was from Mahabo and one from Ambakoana;
# 3 of the Eidolon were from Angavobe/Angavokely and 3 from Ankarana

# 1 of the Pteropus was lactating, and two of the Eidolon were juveniles
# the rest were just adult bats

# 9/11 bats were from 2014
# most samples (7) were from Oct-Jan but we have 1 Pteropus from July
# and 2 Eidolon from March and 1 Rousettus from August

merge.dat$combo_class <- paste0(merge.dat$bat_sex, "_", merge.dat$bat_age_class)
merge.dat$combo_class[merge.dat$combo_class=="male_J" | merge.dat$combo_class=="female_J"] <- "juvenile"

merge.dat$combo_class <- factor(merge.dat$combo_class, levels = c("juvenile", "male_A", "female_NL", "female_P", "female_L"))
merge.dat$bat_species <- as.factor(merge.dat$bat_species)

#try the chi-squared approach from the Peel paper
tableGhVAngV <- table(merge.dat$GhV_pos, merge.dat$AngV_pos)
chisq.test(tableGhVAngV)
#data:  tableGhVAngV
#X-squared = 0.55131, df = 1, p-value = 0.4578

#observed vs. expected
n <- nrow(merge.dat)

pA <- mean(merge.dat$GhV_pos)
pB <- mean(merge.dat$AngV_pos)

expected_coinfection <- pA * pB * n
observed_coinfection <- sum(merge.dat$GhV_pos == 1 & merge.dat$AngV_pos == 1)

expected_coinfection #5.861723
observed_coinfection #8 

#Prob GhV given AngV: 0.06153846
P_GhV_given_AngV <- sum(merge.dat$AngV_pos == 1 & merge.dat$GhV_pos == 1) / sum(merge.dat$GhV_pos == 1)

#Prob AngV given GhV: 0.1777778
P_AngV_given_GhV <- sum(merge.dat$GhV_pos == 1 & merge.dat$AngV_pos == 1) / sum(merge.dat$AngV_pos == 1 )


# and do the same with AndV
tableGhVAndV <- table(merge.dat$GhV_pos, merge.dat$AndV_pos)
chisq.test(tableGhVAndV)
#data:  tableGhVAndV
#X-squared = 0.096637, df = 1, p-value = 0.7559


#Prob GhV given AndV: 0.01538462
P_GhV_given_AndV <- sum(merge.dat$AndV_pos == 1 & merge.dat$GhV_pos == 1) / sum(merge.dat$GhV_pos == 1)

#Prob AndV given GhV: 0.08695652
P_AndV_given_GhV <- sum(merge.dat$GhV_pos == 1 & merge.dat$AndV_pos == 1) / sum(merge.dat$AndV_pos == 1 )


# and do the same with AngV/AndV
tableAngVAndV <- table(merge.dat$AngV_pos, merge.dat$AndV_pos)
chisq.test(tableAngVAndV)
#data:  tableAngVAndV
#X-squared = 1.2799e-29, df = 1, p-value = 1


#Prob AngV given AndV: 0.02222222
P_AngV_given_AndV <- sum(merge.dat$AndV_pos == 1 & merge.dat$AngV_pos == 1) / sum(merge.dat$AngV_pos == 1)

#Prob AndV given AngV: 0.04347826
P_AndV_given_AngV <- sum(merge.dat$AngV_pos == 1 & merge.dat$AndV_pos == 1) / sum(merge.dat$AndV_pos == 1 )


### Now redo it all stratified by species
### Pteropus will look at GhV-AndV
### Eidolon and Rousettus will look at GhV-AngV
### To be thorough, Eidolon will also look at AngV-AndV


### Pteropus - GhV-AndV - not sig
# and do the same with AndV
tableGhVAndV <- table(merge.dat$GhV_pos[merge.dat$bat_species=="Pteropus rufus"], merge.dat$AndV_pos[merge.dat$bat_species=="Pteropus rufus"])
chisq.test(tableGhVAndV)
#data:  tableGhVAndV
#X-squared = 5.6811e-29, df = 1, p-value = 1


#Prob GhV given AndV: 0.08
P_GhV_given_AndV <- sum(merge.dat$AndV_pos[merge.dat$bat_species=="Pteropus rufus"] == 1 & merge.dat$GhV_pos[merge.dat$bat_species=="Pteropus rufus"] == 1) / sum(merge.dat$GhV_pos[merge.dat$bat_species=="Pteropus rufus"] == 1)

#Prob AndV given GhV: 0.125
P_AndV_given_GhV <- sum(merge.dat$GhV_pos[merge.dat$bat_species=="Pteropus rufus"] == 1 & merge.dat$AndV_pos[merge.dat$bat_species=="Pteropus rufus"] == 1) / sum(merge.dat$AndV_pos[merge.dat$bat_species=="Pteropus rufus"] == 1 )

### Rousettus - GhV-AngV
tableGhVAngV <- table(merge.dat$GhV_pos[merge.dat$bat_species=="Rousettus madagascariensis"], merge.dat$AngV_pos[merge.dat$bat_species=="Rousettus madagascariensis"])
chisq.test(tableGhVAngV)
#data:  tableGhVAngV
#X-squared = 1.257, df = 1, p-value = 0.2622


#Prob GhV given AngV:0.08823529
P_GhV_given_AngV <- sum(merge.dat$AngV_pos[merge.dat$bat_species=="Rousettus madagascariensis"] == 1 & merge.dat$GhV_pos[merge.dat$bat_species=="Rousettus madagascariensis"] == 1) / sum(merge.dat$GhV_pos[merge.dat$bat_species=="Rousettus madagascariensis"] == 1)

#Prob AngV given GhV: 0.2307692
P_AngV_given_GhV <- sum(merge.dat$GhV_pos[merge.dat$bat_species=="Rousettus madagascariensis"] == 1 & merge.dat$AngV_pos[merge.dat$bat_species=="Rousettus madagascariensis"] == 1) / sum(merge.dat$AngV_pos[merge.dat$bat_species=="Rousettus madagascariensis"] == 1 )


### Eidolon - GhV-AngV
tableGhVAngV <- table(merge.dat$GhV_pos[merge.dat$bat_species=="Eidolon dupreanum"], merge.dat$AngV_pos[merge.dat$bat_species=="Eidolon dupreanum"])
chisq.test(tableGhVAngV)
#data:  tableGhVAngV
#X-squared = 0.36183, df = 1, p-value = 0.5475


#Prob GhV given AngV:0.07042254
P_GhV_given_AngV <- sum(merge.dat$AngV_pos[merge.dat$bat_species=="Eidolon dupreanum"] == 1 & merge.dat$GhV_pos[merge.dat$bat_species=="Eidolon dupreanum"] == 1) / sum(merge.dat$GhV_pos[merge.dat$bat_species=="Eidolon dupreanum"] == 1)

#Prob AngV given GhV: 0.2272727
P_AngV_given_GhV <- sum(merge.dat$GhV_pos[merge.dat$bat_species=="Eidolon dupreanum"] == 1 & merge.dat$AngV_pos[merge.dat$bat_species=="Eidolon dupreanum"] == 1) / sum(merge.dat$AngV_pos[merge.dat$bat_species=="Eidolon dupreanum"] == 1 )


#and Eidolon AngV-AndV
tableAndVAngV <- table(merge.dat$AndV_pos[merge.dat$bat_species=="Eidolon dupreanum"], merge.dat$AngV_pos[merge.dat$bat_species=="Eidolon dupreanum"])
chisq.test(tableAndVAngV)
#data:  tableAndVAngV
#X-squared = 0.27867, df = 1, p-value = 0.5976


#Prob AndV given AngV:0.2
P_AndV_given_AngV <- sum(merge.dat$AngV_pos[merge.dat$bat_species=="Eidolon dupreanum"] == 1 & merge.dat$AndV_pos[merge.dat$bat_species=="Eidolon dupreanum"] == 1) / sum(merge.dat$AndV_pos[merge.dat$bat_species=="Eidolon dupreanum"] == 1)

#Prob AngV given AndV: 0.04545455
P_AngV_given_AndV <- sum(merge.dat$AndV_pos[merge.dat$bat_species=="Eidolon dupreanum"] == 1 & merge.dat$AngV_pos[merge.dat$bat_species=="Eidolon dupreanum"] == 1) / sum(merge.dat$AngV_pos[merge.dat$bat_species=="Eidolon dupreanum"] == 1 )

