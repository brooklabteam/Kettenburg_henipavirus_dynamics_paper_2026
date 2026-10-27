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
library(cowplot)
library(sjPlot)
library(lme4)
library(lmerTest)

################################################################################
################################################################################
 homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/"
 

# 
#homewd = "/Users/gwenddolenkettenburg/Desktop/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
#batwd = "/Users/gwenddolenkettenburg/Desktop/Developer/Madagascar-Bat-Data"

# plot of seasonality of antigens for analysis
# also GAMs of predictors of seropositivity and/or seroprevalence
################################################################################
################################################################################

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


sum(shed.dat$AngV_pos)#31
sum(shed.dat$GhV_pos)#119

# link this to the catching data to plot by sampling session
catch.dat <- read.csv(file = paste0(homewd, "/master_data_files/catching_data.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))


unique(catch.dat$mass_forearm_residual[catch.dat$young_of_year=="yes"]) #NA as it should be

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, pup_or_mom, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)

#and merge all three datasets
setdiff(shed.dat$sampleid, catch.dat$sampleid)


merge.dat <- merge(shed.dat, catch.dat, by="sampleid")
head(merge.dat)

length(unique(merge.dat$sampleid))
merge.dat$processing_date <- as.Date(merge.dat$processing_date)


#check prevalence overall
#GhV
sum(merge.dat$GhV_pos[merge.dat$bat_species=="Pteropus rufus"])/length(merge.dat$GhV_pos[merge.dat$bat_species=="Pteropus rufus"]) # 25/213 = 11.7%
sum(merge.dat$GhV_pos[merge.dat$bat_species=="Eidolon dupreanum"])/length(merge.dat$GhV_pos[merge.dat$bat_species=="Eidolon dupreanum"]) #83/447= 18.1%
sum(merge.dat$GhV_pos[merge.dat$bat_species=="Rousettus madagascariensis"])/length(merge.dat$GhV_pos[merge.dat$bat_species=="Rousettus madagascariensis"]) # 34/338 = 10.1%

#AndV
sum(merge.dat$AndV_pos[merge.dat$bat_species=="Pteropus rufus"])/length(merge.dat$AndV_pos[merge.dat$bat_species=="Pteropus rufus"]) # 16/213 = 7.5%
sum(merge.dat$AndV_pos[merge.dat$bat_species=="Eidolon dupreanum"])/length(merge.dat$AndV_pos[merge.dat$bat_species=="Eidolon dupreanum"]) #5/447= 1.1%
sum(merge.dat$AndV_pos[merge.dat$bat_species=="Rousettus madagascariensis"])/length(merge.dat$AndV_pos[merge.dat$bat_species=="Rousettus madagascariensis"]) # 2/338 = 0.6%

#AngV
sum(merge.dat$AngV_pos[merge.dat$bat_species=="Pteropus rufus"])/length(merge.dat$AngV_pos[merge.dat$bat_species=="Pteropus rufus"]) # 10/213 = 4.6%
sum(merge.dat$AngV_pos[merge.dat$bat_species=="Eidolon dupreanum"])/length(merge.dat$AngV_pos[merge.dat$bat_species=="Eidolon dupreanum"]) #22/447= 4.9%
sum(merge.dat$AngV_pos[merge.dat$bat_species=="Rousettus madagascariensis"])/length(merge.dat$AngV_pos[merge.dat$bat_species=="Rousettus madagascariensis"]) # 13/338 = 3.2%


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

#and sampling session
#summarise more broadly
sampling.date <- ddply(merge.dat, .(sampling_session), summarise, min_date = min(processing_date), max_date = max(processing_date))
sampling.date$mid_date <- (sampling.date$max_date-sampling.date$min_date)/2 + sampling.date$min_date

# #and slim to merge
sampling.merge <- dplyr::select(sampling.date, sampling_session, mid_date)
# 
merge.dat <- merge(merge.dat, sampling.merge, by = c( "sampling_session"))
merge.dat$year <- year(merge.dat$processing_date)
merge.dat$month <- as.character(month(merge.dat$processing_date))
merge.dat$month[merge.dat$month!="10"& merge.dat$month!="11" &merge.dat$month!="12"]<-paste0("0", merge.dat$month[merge.dat$month!="10"& merge.dat$month!="11" &merge.dat$month!="12"])
#merge.dat$mid_date <- as.Date(paste0(merge.dat$year, "-", merge.dat$month, "-15"))
head(merge.dat)

sampling.date <- ddply(merge.dat, .(bat_species, roost_site, mid_date), summarise)
sampling.date# examine

#
# #write over those that are super close
merge.dat$mid_date <- as.character(merge.dat$mid_date)
merge.dat$mid_date[merge.dat$bat_species=="Eidolon dupreanum" & merge.dat$roost_site=="Moramanga" & merge.dat$mid_date=="2018-02-15"] <- "2018-02-14"
merge.dat$mid_date[merge.dat$bat_species=="Rousettus madagascariensis" & merge.dat$roost_site=="Moramanga" & merge.dat$mid_date=="2019-01-09"] <- "2019-01-05"
merge.dat$mid_date <- as.Date(merge.dat$mid_date)

merge.dat$doy <- yday(merge.dat$processing_date)



tot.merge= merge.dat

head(tot.merge)
tot.merge$processing_date <- as.Date(tot.merge$processing_date)
tot.merge$month <- month(tot.merge$processing_date)
tot.merge$year <- year(tot.merge$processing_date)
#tot.merge = subset(tot.merge, bat_species=="Eidolon dupreanum" | bat_species =="Rousettus madagascariensis" )
# #tot.merge = subset(tot.merge, bat_species=="Eidolon dupreanum" | bat_species =="Rousettus madagascariensis" | bat_species =="Pteropus rufus")
# 

tot.merge$year <- as.factor(tot.merge$year)
tot.merge$roost_site <- as.factor(tot.merge$roost_site)
tot.merge$bat_sex <- as.factor(tot.merge$bat_sex)

#########
#########

# FIRST, we ask, what are the correlates of positive infection status with all three viruses in all three bat species?
### the ambiguous age class bats won't be included in the first models because they have NA values
length(tot.merge$bat_age_class[is.na(tot.merge$bat_age_class)])#46
#but they can be retained for seasonality

#verify juvenile status
unique(tot.merge$young_of_year)
unique(tot.merge$bat_age_class[tot.merge$young_of_year=="yes"]) 
tot.merge$bat_age_class[tot.merge$young_of_year=="yes"]<- "J" #all yoy are juveniles

# 
unique(tot.merge$bat_age_class)

# make a few extra variables
tot.merge$combo_class <- paste0(tot.merge$bat_sex, "-", tot.merge$bat_age_class)
#combine juveniles 
tot.merge$combo_class[tot.merge$bat_age_class=="J"] <- "juvenile"
tot.merge$combo_class[tot.merge$combo_class=="male-NA" | tot.merge$combo_class=="female-NA" | tot.merge$combo_class=="NA-NA"] <- NA
unique(tot.merge$combo_class) # don't use NAs in this analysis
tot.merge$MFR_class <- tot.merge$combo_class
tot.merge$MFR_class[tot.merge$MFR_class=="female-L" | tot.merge$MFR_class=="male-A" | tot.merge$MFR_class=="female-NL"] <- "adult"
tot.merge$MFR_class <- as.factor(tot.merge$MFR_class)
tot.merge$combo_class <- as.factor(tot.merge$combo_class)
unique(tot.merge$MFR_class) # don't use NAs in this analysis
tot.merge$mass_forearm_residual[tot.merge$MFR_class=="juvenile"] <- 0

#Eidolon first
gamGhV <- gam(GhV_pos ~ combo_class + mass_forearm_residual:MFR_class +
                     s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                     s(year, bs="re"),
                   family = "binomial",
                   method = "REML", 
                   gamma = 1,
                   select = TRUE,
                   data= subset(tot.merge, bat_species =="Eidolon dupreanum" ))
summary(gamGhV) #year and doy ankarana only  -- this leads us to explore the seasonal data more
#plot(gamGhV)


gamAngV <- gam(AngV_pos ~ combo_class + mass_forearm_residual:MFR_class +
                s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                s(year, bs="re"),
              family = "binomial",
              method = "REML", 
              gamma = 1,
              select = TRUE,
              data= subset(tot.merge, bat_species =="Eidolon dupreanum" ))
summary(gamAngV) #year and MFR adult (adult MFR slope is weakly positive - could just be coincidence). and ankarana weakly pos
#plot(gamAngV)

gamAndV <- gam(AndV_pos ~  combo_class + mass_forearm_residual:MFR_class +
                 s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                 s(year, bs="re"),
               family = "binomial",
               method = "REML", 
               gamma = 1,
               select = TRUE,
               data= subset(tot.merge, bat_species =="Eidolon dupreanum" ))
summary(gamAndV) #year only

### then Pteropus

gamPGhV <- gam(GhV_pos ~  combo_class + mass_forearm_residual:MFR_class +
                s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                s(year, bs="re"),
              family = "binomial",
              method = "REML", 
              gamma = 1,
              select = TRUE,
              data= subset(tot.merge, bat_species =="Pteropus rufus" ))
summary(gamPGhV) #year,  doy Moramanga, NL females borderline positively associated, MFR adult weakly positive
#plot(gamPGhV)


gamPAngV <- gam(AngV_pos ~  combo_class + mass_forearm_residual:MFR_class +
                 s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                 s(year, bs="re"),
               family = "binomial",
               method = "REML", 
               gamma = 1,
               select = TRUE,
               data= subset(tot.merge, bat_species =="Pteropus rufus" ))
summary(gamPAngV) #nothing 
#plot(gamPAngV)

gamPAndV <- gam(AndV_pos ~  combo_class + mass_forearm_residual:MFR_class +
                 s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                 s(year, bs="re"),
               family = "binomial",
               method = "REML", 
               gamma = 1,
               select = TRUE,
               data= subset(tot.merge, bat_species =="Pteropus rufus" ))
summary(gamPAndV) #year.

### and Rousettus


gamRGhV <- gam(GhV_pos ~  combo_class + mass_forearm_residual:MFR_class +
                 s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                 s(year, bs="re"),
               family = "binomial",
               method = "REML", 
               gamma = 1,
               select = TRUE,
               data= subset(tot.merge, bat_species =="Rousettus madagascariensis" ))
summary(gamRGhV) #year only
#plot(gamRGhV)


gamRAngV <- gam(AngV_pos ~  combo_class + mass_forearm_residual:MFR_class +
                  s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                  s(year, bs="re"),
                family = "binomial",
                method = "REML", 
                gamma = 1,
                select = TRUE,
                data= subset(tot.merge, bat_species =="Rousettus madagascariensis" ))
summary(gamRAngV) #moramanga weakly sig
#plot(gamRAngV)

gamRAndV <- gam(AndV_pos ~  combo_class + mass_forearm_residual:MFR_class +
                  s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                  s(year, bs="re"),
                family = "binomial",
                method = "REML", 
                gamma = 1,
                select = TRUE,
                data= subset(tot.merge, bat_species =="Rousettus madagascariensis" ))
summary(gamRAndV) #nothing
###

#make a csv file of all the model outputs
library(broom)
model_list<-list(gamGhV=gamGhV, gamAngV= gamAngV, gamAndV= gamAndV, 
                 gamPGhV=gamPGhV, gamPAngV= gamPAngV, gamPAndV= gamPAndV,
                 gamRGhV=gamRGhV, gamRAngV= gamRAngV, gamRAndV= gamRAndV)
all_summaries <- bind_rows(lapply(model_list, tidy), .id = "Model_Name")
write.csv(all_summaries, "intraannual_shed_gam_model_summaries.csv", row.names = FALSE)

# based on above, the only consistently significant results were the
# seasonal smoothers by site - so let's look at those. to make more 
# apparent, we chose to aggregate data by sampling session next


#####
session.sum = ddply(tot.merge, .(bat_species, sampling_session, mid_date, roost_site), summarise, GhVpos=sum(GhV_pos), AngVpos=sum(AngV_pos), AndVpos=sum(AndV_pos), Ntot=length(sampleid))
session.sum$doy <- yday(session.sum$mid_date)
session.sum$year <- as.factor(year(session.sum$mid_date))
session.sum$GhVneg <- session.sum$Ntot - session.sum$GhVpos
session.sum$AngVneg <- session.sum$Ntot - session.sum$AngVpos
session.sum$AndVneg <- session.sum$Ntot - session.sum$AndVpos

#Eidolon
gamSessionGhV <- gam(cbind(GhVpos, GhVneg) ~ 
                     s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                     s(year, bs="re"),
                   family = "binomial",
                   method = "REML", 
                   gamma = 1,
                   select = TRUE,
                   data= subset(session.sum, bat_species =="Eidolon dupreanum" ))
summary(gamSessionGhV) #year and ankarana 
#plot(gamSessionGhV)



gamSessionAngV <- gam(cbind(AngVpos, AngVneg) ~ 
                       s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                       s(year, bs="re"),
                     family = "binomial",
                     method = "REML", 
                     gamma = 1,
                     select = TRUE,
                     data= subset(session.sum, bat_species =="Eidolon dupreanum" ))
summary(gamSessionAngV) #year and ankarana
#plot(gamSessionAngV)

gamSessionAndV <- gam(cbind(AndVpos, AndVneg) ~ 
                        s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                        s(year, bs="re"),
                      family = "binomial",
                      method = "REML", 
                      gamma = 1,
                      select = TRUE,
                      data= subset(session.sum, bat_species =="Eidolon dupreanum" ))
summary(gamSessionAndV) #year 
#plot(gamSessionAndV)

#Pteropus

gamSessionPGhV <- gam(cbind(GhVpos, GhVneg) ~ 
                       s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                       s(year, bs="re"),
                     family = "binomial",
                     method = "REML", 
                     gamma = 1,
                     select = TRUE,
                     data= subset(session.sum, bat_species =="Pteropus rufus" ))
summary(gamSessionPGhV) #year and moramanga
#plot(gamSessionPGhV)



gamSessionPAngV <- gam(cbind(AngVpos, AngVneg) ~ 
                        s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                        s(year, bs="re"),
                      family = "binomial",
                      method = "REML", 
                      gamma = 1,
                      select = TRUE,
                      data= subset(session.sum, bat_species =="Pteropus rufus" ))
summary(gamSessionPAngV) #nothing
#plot(gamSessionPAngV)

gamSessionPAndV <- gam(cbind(AndVpos, AndVneg) ~ 
                        s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                        s(year, bs="re"),
                      family = "binomial",
                      method = "REML", 
                      gamma = 1,
                      select = TRUE,
                      data= subset(session.sum, bat_species =="Pteropus rufus" ))
summary(gamSessionPAndV) #year 
#plot(gamSessionPAndV)

# Rousettus
gamSessionRGhV <- gam(cbind(GhVpos, GhVneg) ~ 
                        s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                        s(year, bs="re"),
                      family = "binomial",
                      method = "REML", 
                      gamma = 1,
                      select = TRUE,
                      data= subset(session.sum, bat_species =="Rousettus madagascariensis" ))
summary(gamSessionRGhV) #year 
#plot(gamSessionRGhV)

gamSessionRAngV <- gam(cbind(AngVpos, AngVneg) ~ 
                         s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                         s(year, bs="re"),
                       family = "binomial",
                       method = "REML", 
                       gamma = 1,
                       select = TRUE,
                       data= subset(session.sum, bat_species =="Rousettus madagascariensis" ))
summary(gamSessionRAngV) #ankarana
#plot(gamSessionRAngV)

gamSessionRAndV <- gam(cbind(AndVpos, AndVneg) ~ 
                         s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                         s(year, bs="re"),
                       family = "binomial",
                       method = "REML", 
                       gamma = 1,
                       select = TRUE,
                       data= subset(session.sum, bat_species =="Rousettus madagascariensis" ))
summary(gamSessionRAndV) #nothing
#plot(gamSessionRAndV)

#again, save a csv of all the model outputs: 
library(broom)
model_list<-list(gamSessionGhV=gamSessionGhV, gamSessionAngV= gamSessionAngV, gamSessionAndV= gamSessionAndV, 
                 gamSessionPGhV=gamSessionPGhV, gamSessionPAngV= gamSessionPAngV, gamSessionPAndV= gamSessionPAndV,
                 gamSessionRGhV=gamSessionRGhV, gamSessionRAngV= gamSessionRAngV, gamSessionRAndV= gamSessionRAndV)
all_summaries <- bind_rows(lapply(model_list, tidy), .id = "Model_Name")
write.csv(all_summaries, "intraannual_session_shed_gam_model_summaries.csv", row.names = FALSE)


#####
# Great, let's plot the seasonal data

pred.df <- cbind.data.frame(mid_date = seq(from= as.Date("2019-01-01"), to = as.Date("2019-12-31"), by=1), doy=1:365)
pred.df$bat_species <- "Eidolon dupreanum"
pred.df$antigen <- "GhV"
pred.df2 <- pred.df
pred.df2$antigen <- "AngV"
pred.df3 <- pred.df
pred.df3$antigen <- "AndV"
pred.df <- rbind(pred.df, pred.df2, pred.df3)
pred.df$roost_site <- "Moramanga"
pred.df$year <- 2019

pred.df$prevalence <- NA
pred.df$prev_lci <- NA
pred.df$prev_uci <- NA


pred.df.P <- pred.df
pred.df.P$bat_species = "Pteropus rufus"

pred.df.R <- pred.df
pred.df.R$bat_species = "Rousettus madagascariensis"


pred.df <- rbind(pred.df, pred.df.P, pred.df.R)

pred.df$prevalence[pred.df$antigen=="GhV" & pred.df$bat_species=="Eidolon dupreanum"] <- predict(gamSessionGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="GhV" & pred.df$bat_species=="Eidolon dupreanum"]<- predict(gamSessionGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c("s(year)")) - 1.96*predict(gamSessionGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c("s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="GhV" & pred.df$bat_species=="Eidolon dupreanum"] <- predict(gamSessionGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c("s(year)")) + 1.96*predict(gamSessionGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c( "s(year)"), se.fit = TRUE)$se

pred.df$prevalence[pred.df$antigen=="AngV" & pred.df$bat_species=="Eidolon dupreanum"] <- predict(gamSessionAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="AngV" & pred.df$bat_species=="Eidolon dupreanum"]<- predict(gamSessionAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c("s(year)")) - 1.96*predict(gamSessionAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c("s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="AngV" & pred.df$bat_species=="Eidolon dupreanum"] <- predict(gamSessionAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c("s(year)")) + 1.96*predict(gamSessionAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c( "s(year)"), se.fit = TRUE)$se

pred.df$prevalence[pred.df$antigen=="AndV" & pred.df$bat_species=="Eidolon dupreanum"] <- predict(gamSessionAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="AndV" & pred.df$bat_species=="Eidolon dupreanum"]<- predict(gamSessionAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c("s(year)")) - 1.96*predict(gamSessionAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c("s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="AndV" & pred.df$bat_species=="Eidolon dupreanum"] <- predict(gamSessionAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c("s(year)")) + 1.96*predict(gamSessionAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Eidolon dupreanum",], exclude = c( "s(year)"), se.fit = TRUE)$se

pred.df$prevalence[pred.df$antigen=="GhV" & pred.df$bat_species=="Pteropus rufus"] <- predict(gamSessionPGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Pteropus rufus",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="GhV" & pred.df$bat_species=="Pteropus rufus"]<- predict(gamSessionPGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Pteropus rufus",], exclude = c("s(year)")) - 1.96*predict(gamSessionPGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Pteropus rufus",], exclude = c("s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="GhV" & pred.df$bat_species=="Pteropus rufus"] <- predict(gamSessionPGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Pteropus rufus",], exclude = c("s(year)")) + 1.96*predict(gamSessionPGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Pteropus rufus",], exclude = c( "s(year)"), se.fit = TRUE)$se

pred.df$prevalence[pred.df$antigen=="AngV" & pred.df$bat_species=="Pteropus rufus"] <- predict(gamSessionPAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Pteropus rufus",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="AngV" & pred.df$bat_species=="Pteropus rufus"]<- predict(gamSessionPAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Pteropus rufus",], exclude = c("s(year)")) - 1.96*predict(gamSessionPAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Pteropus rufus",], exclude = c("s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="AngV" & pred.df$bat_species=="Pteropus rufus"] <- predict(gamSessionPAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Pteropus rufus",], exclude = c("s(year)")) + 1.96*predict(gamSessionPAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Pteropus rufus",], exclude = c( "s(year)"), se.fit = TRUE)$se

pred.df$prevalence[pred.df$antigen=="AndV" & pred.df$bat_species=="Pteropus rufus"] <- predict(gamSessionPAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Pteropus rufus",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="AndV" & pred.df$bat_species=="Pteropus rufus"]<- predict(gamSessionPAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Pteropus rufus",], exclude = c("s(year)")) - 1.96*predict(gamSessionPAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Pteropus rufus",], exclude = c("s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="AndV" & pred.df$bat_species=="Pteropus rufus"] <- predict(gamSessionPAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Pteropus rufus",], exclude = c("s(year)")) + 1.96*predict(gamSessionPAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Pteropus rufus",], exclude = c( "s(year)"), se.fit = TRUE)$se


pred.df$prevalence[pred.df$antigen=="GhV" & pred.df$bat_species=="Rousettus madagascariensis"] <- predict(gamSessionRGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="GhV" & pred.df$bat_species=="Rousettus madagascariensis"]<- predict(gamSessionRGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c("s(year)")) - 1.96*predict(gamSessionRGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c("s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="GhV" & pred.df$bat_species=="Rousettus madagascariensis"] <- predict(gamSessionRGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c("s(year)")) + 1.96*predict(gamSessionRGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c( "s(year)"), se.fit = TRUE)$se

pred.df$prevalence[pred.df$antigen=="AngV" & pred.df$bat_species=="Rousettus madagascariensis"] <- predict(gamSessionRAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="AngV" & pred.df$bat_species=="Rousettus madagascariensis"]<- predict(gamSessionRAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c("s(year)")) - 1.96*predict(gamSessionRAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c("s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="AngV" & pred.df$bat_species=="Rousettus madagascariensis"] <- predict(gamSessionRAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c("s(year)")) + 1.96*predict(gamSessionRAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c( "s(year)"), se.fit = TRUE)$se

pred.df$prevalence[pred.df$antigen=="AndV" & pred.df$bat_species=="Rousettus madagascariensis"] <- predict(gamSessionRAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="AndV" & pred.df$bat_species=="Rousettus madagascariensis"]<- predict(gamSessionRAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c("s(year)")) - 1.96*predict(gamSessionRAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c("s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="AndV" & pred.df$bat_species=="Rousettus madagascariensis"] <- predict(gamSessionRAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c("s(year)")) + 1.96*predict(gamSessionRAndV, type = "response", newdata = pred.df[pred.df$antigen=="AndV" & pred.df$bat_species=="Rousettus madagascariensis",], exclude = c( "s(year)"), se.fit = TRUE)$se


pred.df$prev_lci[pred.df$prev_lci<0] <- 0
pred.df$prev_uci[pred.df$prev_uci>1] <- 1

session.sum$prevalence_GhV <- session.sum$GhVpos/session.sum$Ntot
session.sum$prevalence_AngV <- session.sum$AngVpos/session.sum$Ntot
session.sum$prevalence_AndV <- session.sum$AndVpos/session.sum$Ntot


session.sum$prevGhV_lci <- session.sum$prevGhV_uci <- session.sum$prevAngV_lci <- session.sum$prevAngV_uci <-  session.sum$prevAndV_lci <- session.sum$prevAndV_uci <- NA
for (i in 1:length(session.sum$prevalence_GhV)){
  session.sum$prevGhV_lci[i] <- binom.test(session.sum$GhVpos[i], session.sum$Ntot[i])$conf.int[1]
  session.sum$prevGhV_uci[i] <- binom.test(session.sum$GhVpos[i], session.sum$Ntot[i])$conf.int[2]
  session.sum$prevAngV_lci[i] <- binom.test(session.sum$AngVpos[i], session.sum$Ntot[i])$conf.int[1]
  session.sum$prevAngV_uci[i] <- binom.test(session.sum$AngVpos[i], session.sum$Ntot[i])$conf.int[2]
  session.sum$prevAndV_lci[i] <- binom.test(session.sum$AndVpos[i], session.sum$Ntot[i])$conf.int[1]
  session.sum$prevAndV_uci[i] <- binom.test(session.sum$AndVpos[i], session.sum$Ntot[i])$conf.int[2]
  
}


#and make the gestation database
gest.dat <- cbind.data.frame(bat_species=c("Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis"), start = c(as.Date("2019-07-07"), as.Date("2019-08-03"), as.Date("2019-09-11")), stop=c(as.Date("2019-09-28"), as.Date("2019-11-15"),  as.Date("2019-12-11")))
gest.dat$season = "gestation"
lact.dat <- cbind.data.frame(bat_species=c("Pteropus rufus","Pteropus rufus", "Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis",  "Rousettus madagascariensis"), start = c(as.Date("2019-09-29"), as.Date("2019-01-01"), as.Date("2019-11-16"), as.Date("2019-01-01"), as.Date("2019-12-12"),  as.Date("2019-01-01")), stop=c(as.Date("2019-12-31"), as.Date("2019-01-21"), as.Date("2019-12-31"), as.Date("2019-02-02"), as.Date("2019-12-31"), as.Date("2015-02-19")))
lact.dat$season = "lactation"

#and winter
cal.dat <-  cbind.data.frame(bat_species=c("Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis"), start = c(as.Date("2019-05-15"), as.Date("2019-05-15"),as.Date("2019-05-15")), stop=c(as.Date("2019-10-10"), as.Date("2019-10-10"),as.Date("2019-10-10")))
cal.dat$season <- "winter"
cal.dat <- rbind(cal.dat, gest.dat, lact.dat)
cal.dat$year <- year(cal.dat$stop)

cal.dat$doy_start <- yday(cal.dat$start)
cal.dat$doy_stop <- yday(cal.dat$stop)

sitez = c('Ankarana' ="#A58AFF",  'Moramanga' = "#00C094", 'Mahabo' = "#619CFF")


colorz = c('AngV' = "#F8766D", 'HeV' ="#A58AFF", 
           'NiV' = "#DB72FB", 'GhV' = "#00C094", 
           'CedV'= "#53B400", "MojV" = "#619CFF")


filz <- c('gestation' = "yellow", 'lactation' = "lightpink", 'winter' = "skyblue")

yearz <- c('2013'='3', '2014' = '7', '2015'='8', '2016'='18', '2018' = '16', '2019' = '17', '2020' = '15')
yearz <- c('2013'=3, '2014' = 7, '2015'=8, '2016'=18, '2018' = 16, '2019' = 17, '2020' = 15, '2021' = 4, '2023' = 6)

session.sum$bat_species <- factor(session.sum$bat_species, levels = c("Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis"))
pred.df$bat_species <- factor(pred.df$bat_species, levels = c("Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis"))
cal.dat$bat_species <- factor(cal.dat$bat_species, levels = c("Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis"))

head(session.sum)

session.wide <- dplyr::select(session.sum, bat_species, roost_site, doy,year, Ntot, prevalence_GhV, prevalence_AngV, prevGhV_lci, prevGhV_uci, prevAngV_lci, prevAngV_uci)

library(tidyr)
names(session.wide) <- c(names(session.wide)[1:5], "prevalence_GhV", "prevalence_AngV", "prev_lci_GhV", "prev_uci_GhV", "prev_lci_AngV", "prev_uci_AngV")
session.long <-session.wide %>%
  pivot_longer(
    cols = prevalence_GhV:prev_uci_AngV,
    names_to = c(".value", "antigen"),
    names_pattern = "(prevalence|prev_lci|prev_uci)_(GhV|AngV)"
  )
 head(session.long)

pIntrannualShedAll <- ggplot(data=session.long) +
  geom_rect(data=cal.dat, aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=doy, y=prevalence, size=Ntot, shape=year, color = roost_site)) + 
  geom_linerange(aes(x=doy, ymin=prev_lci, ymax=prev_uci, color=roost_site),linetype=2, alpha=.3, show.legend = F) +
  scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  scale_fill_manual(values=filz) +
  scale_size_continuous(name="N bats", breaks = c(10,20,30)) +
  scale_color_manual(values=sitez, name="roost site") +
  ggnewscale::new_scale_fill() + 
  ggnewscale::new_scale_color() + 
  geom_ribbon(data=subset(pred.df, antigen!="AndV"), aes(x=doy, ymin= prev_lci, ymax=prev_uci, fill=antigen),  alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df, antigen!="AndV"), aes(x=doy, y= prevalence, color=antigen), show.legend = F) +
  facet_grid(bat_species~antigen) + theme_bw() +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),, 
        strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16, face="italic"), 
        strip.text.x = element_text(size=16),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  scale_shape_manual(values=yearz) +
  coord_cartesian(xlim = c(0, 365), ylim = c(0, 1),expand = FALSE) + ylab("prevalence") + xlab("month of year") +
  scale_fill_manual(values=colorz) + scale_color_manual(values=colorz) + 
  guides(size=guide_legend(override.aes = list(color="black")), shape=guide_legend(override.aes = list(color="black")))
pIntrannualShedAll #this goes in supplement

dat_text <- data.frame(
  antigen = c("AngV","AngV", "AngV", "GhV", "GhV", "GhV"),
  bat_species = c("Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis",
         "Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis"),
 
  label = c(" p = NS Mahabo\n p = NS Makira\n p = NS Moramanga", 
            " p = 0.070 Ankarana\n p = NS Mahabo\n p = NS Moramanga", 
            " p = NS Ankarana\n p = NS Makira\n p = 0.088 Moramanga", 
            " p = NS Mahabo\n p = NS Makira\n p = 0.014 Moramanga", 
            " p = 0.008 Ankarana\n p = NS Mahabo\n p = NS Moramanga", 
            " p = NS Ankarana\n p = NS Makira\n p = NS Moramanga")
)

pIntrannualShedAll<-pIntrannualShedAll + geom_text(
  data    = dat_text,
  mapping = aes(x = 130, y = 0.5, label = label),
  hjust   = -0.1,
  vjust   = -1
);pIntrannualShedAll


ggsave(file = paste0(homewd,"/final_figures_tables/extended-data/ExtendedDataFig3.png"),
       plot=pIntrannualShedAll,
       units="mm",  
       width=80, 
       height=90, 
       scale=3, 
       dpi=300)



#and for the main text, plot just Eidolon data

pIntrannualShedEid <- ggplot(data=subset(session.sum, bat_species=="Eidolon dupreanum")) +
  geom_rect(data=subset(cal.dat, bat_species=="Eidolon dupreanum"), aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=doy, y=prevalence_GhV, size=Ntot, shape=year, color = roost_site)) + 
  geom_linerange(aes(x=doy, ymin=prevGhV_lci, ymax=prevGhV_uci, color=roost_site),linetype=2, alpha=.3, show.legend = F) +
  scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  scale_fill_manual(values=filz) +
  scale_color_manual(values=sitez, name="roost site") +
  #ggnewscale::new_scale_fill() + 
  geom_ribbon(data=subset(pred.df,antigen=="GhV"&bat_species=="Eidolon dupreanum"), aes(x=doy, ymin= prev_lci, ymax=prev_uci), fill = "#00C094", alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df,antigen=="GhV"&bat_species=="Eidolon dupreanum"), aes(x=doy, y= prevalence), color = "#00C094") +
  facet_grid(bat_species~.) + theme_bw() +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),, 
        strip.background = element_blank(), 
        strip.text.y = element_blank(), 
        #strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  scale_shape_manual(values=yearz) +
  coord_cartesian(xlim = c(0, 365), ylim = c(0, 1),expand = FALSE) + ylab("prevalence (GhV)") + xlab("month of year") +
  #scale_fill_manual(values=colorz) + scale_color_manual(values=colorz) + 
  guides(size=guide_legend(override.aes = list(color="black")), shape=guide_legend(override.aes = list(color="black")))
pIntrannualShedEid



Eid_Leg <- cowplot::get_legend(pIntrannualShedEid)

#annotate
pIntrannualShedEid_noLeg <- ggplot(data=subset(session.sum, bat_species=="Eidolon dupreanum")) +
  geom_rect(data=subset(cal.dat, bat_species=="Eidolon dupreanum"), aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3,show.legend = F) +
  geom_point(aes(x=doy, y=prevalence_GhV, size=Ntot, shape=year, color = roost_site),show.legend = F) + 
  geom_linerange(aes(x=doy, ymin=prevGhV_lci, ymax=prevGhV_uci, color=roost_site),linetype=2, alpha=.3, show.legend = F) +
  scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  scale_fill_manual(values=filz) +
  scale_color_manual(values=sitez, name="roost site") +
  #ggnewscale::new_scale_fill() + 
  geom_ribbon(data=subset(pred.df,antigen=="GhV"&bat_species=="Eidolon dupreanum"), aes(x=doy, ymin= prev_lci, ymax=prev_uci), fill = "#00C094", alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df,antigen=="GhV"&bat_species=="Eidolon dupreanum"), aes(x=doy, y= prevalence), color = "#00C094",show.legend = F) +
  facet_grid(bat_species~.) + theme_bw() +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),, 
        strip.background = element_blank(), 
        strip.text.y = element_blank(), 
        #strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  scale_shape_manual(values=yearz) +
  coord_cartesian(xlim = c(0, 365), ylim = c(0, 1),expand = FALSE) + ylab("prevalence (GhV)") + xlab("month of year") +
  #scale_fill_manual(values=colorz) + scale_color_manual(values=colorz) + 
  guides(size=guide_legend(override.aes = list(color="black")), shape=guide_legend(override.aes = list(color="black")))
pIntrannualShedEid_noLeg

dat_text <- data.frame(
  label = c(" p = 0.005 Ankarana\n p = NS Mahabo\n p = NS Moramanga"),
  cyl   = c("GhV")
)

pIntrannualShedEid_noLeg<-pIntrannualShedEid_noLeg + geom_text(
  data    = dat_text,
  mapping = aes(x = 210, y = 0.75, label = label),
  hjust   = -0.1,
  vjust   = -1
);pIntrannualShedEid_noLeg



#####

#Now, look at the effects of MFR on prevalence
tot.merge$MFR_cat <- NA
tot.merge$MFR_cat[tot.merge$mass_forearm_residual<0] <- "neg"
tot.merge$MFR_cat[tot.merge$mass_forearm_residual>0] <- "pos"
tot.merge$MFR_cat[tot.merge$MFR_class!="adult"] <- NA

#Eidolon 
t.test(x= tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$GhV_pos==1 & tot.merge$bat_species=="Eidolon dupreanum"], y=tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$GhV_pos==0 & tot.merge$bat_species=="Eidolon dupreanum"])
#GhV - no diff.  p-value = 0.8491

t.test(x= tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$AngV_pos==1 & tot.merge$bat_species=="Eidolon dupreanum"], y=tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$AngV_pos==0 & tot.merge$bat_species=="Eidolon dupreanum"])
#AngV - no diff. p-value = 0.1978

t.test(x= tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$AndV_pos==1 & tot.merge$bat_species=="Eidolon dupreanum"], y=tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$AndV_pos==0 & tot.merge$bat_species=="Eidolon dupreanum"])
#AndV - no diff. p-value = 0.581


#Pteropus
t.test(x= tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$GhV_pos==1 & tot.merge$bat_species=="Pteropus rufus"], y=tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$GhV_pos==0 & tot.merge$bat_species=="Pteropus rufus"])
#GhV - no diff. p-value = 0.05132

t.test(x= tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$AngV_pos==1 & tot.merge$bat_species=="Pteropus rufus"], y=tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$AngV_pos==0 & tot.merge$bat_species=="Pteropus rufus"])
#AngV - no diff. p-value = 0.9132

t.test(x= tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$AndV_pos==1 & tot.merge$bat_species=="Pteropus rufus"], y=tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$AndV_pos==0 & tot.merge$bat_species=="Pteropus rufus"])
#AndV - no diff. p-value = 0.8816


#Rousettus
t.test(x= tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$GhV_pos==1 & tot.merge$bat_species=="Rousettus madagascariensis"], y=tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$GhV_pos==0 & tot.merge$bat_species=="Rousettus madagascariensis"])
#GhV -  no diff. p-value = 0.1012

t.test(x= tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$AngV_pos==1 & tot.merge$bat_species=="Rousettus madagascariensis"], y=tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$AngV_pos==0 & tot.merge$bat_species=="Rousettus madagascariensis"])
#AngV - no diff. p-value = 0.2436

t.test(x= tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$AndV_pos==1 & tot.merge$bat_species=="Rousettus madagascariensis"], y=tot.merge$mass_forearm_residual[!is.na(tot.merge$MFR_cat) & tot.merge$AndV_pos==0 & tot.merge$bat_species=="Rousettus madagascariensis"])
#AndV - no diff.p-value = 0.4961

tot.merge$GhV_plot <- tot.merge$GhV_pos
tot.merge$GhV_plot[tot.merge$GhV_plot==0] <- "GhV-neg"
tot.merge$GhV_plot[tot.merge$GhV_plot==1] <- "GhV-pos"


tot.merge$AngV_plot <- tot.merge$AngV_pos
tot.merge$AngV_plot[tot.merge$AngV_plot==0] <- "AngV-neg"
tot.merge$AngV_plot[tot.merge$AngV_plot==1] <- "AngV-pos"

tot.merge$AndV_plot <- tot.merge$AndV_pos
tot.merge$AndV_plot[tot.merge$AndV_plot==0] <- "AndV-neg"
tot.merge$AndV_plot[tot.merge$AndV_plot==1] <- "AndV-pos"


#probably don't need the plots since nothing significant, but here they are
pMFRprevGhV <- ggplot(subset(tot.merge, !is.na(MFR_cat))) + theme_bw()+
            theme(panel.grid = element_blank(),
            strip.background = element_rect(fill="white"),
            strip.text.x = element_text(face="italic",size=16),
            strip.text.y = element_text(size=16),
            axis.text = element_text(size=14), axis.title = element_text(size=16),
            axis.title.x = element_blank()) +
            ylab("mass: forearm residual") +
            geom_violin(aes(x=GhV_plot, y=mass_forearm_residual, group=GhV_pos), draw_quantiles = c(0,.25,.5,.75)) + 
            geom_jitter(aes(x=GhV_plot, y=mass_forearm_residual, color=bat_sex, shape=roost_site), width = .2) +
            facet_grid(bat_sex~bat_species); pMFRprevGhV



pMFRprevAngV <- ggplot(subset(tot.merge, !is.na(MFR_cat))) + theme_bw()+
            theme(panel.grid = element_blank(),
            strip.background = element_rect(fill="white"),
            strip.text.x = element_text(face="italic",size=16),
            strip.text.y = element_text(size=16),
            axis.text = element_text(size=14), axis.title = element_text(size=16),
            axis.title.x = element_blank()) +
            ylab("mass: forearm residual") +
            geom_violin(aes(x=AngV_plot, y=mass_forearm_residual, group=AngV_pos), draw_quantiles = c(0,.25,.5,.75)) + 
            geom_jitter(aes(x=AngV_plot, y=mass_forearm_residual, color=bat_sex, shape=roost_site), width = .2) +
            facet_grid(bat_sex~bat_species); pMFRprevAngV



pMFRprevAndV <- ggplot(subset(tot.merge, !is.na(MFR_cat))) + theme_bw()+
              theme(panel.grid = element_blank(),
              strip.background = element_rect(fill="white"),
              strip.text.x = element_text(face="italic",size=16),
              strip.text.y = element_text(size=16),
              axis.text = element_text(size=14), axis.title = element_text(size=16),
              axis.title.x = element_blank()) +
              ylab("mass: forearm residual") +
              geom_violin(aes(x=AndV_plot, y=mass_forearm_residual, group=AndV_pos), draw_quantiles = c(0,.25,.5,.75)) + 
              geom_jitter(aes(x=AndV_plot, y=mass_forearm_residual, color=bat_sex, shape=roost_site), width = .2) +
              facet_grid(bat_sex~bat_species); pMFRprevAndV




####################
####################
#now let's look at predictors of Ct - low Ct is higher viral load
tot.merge$mass_cat <- NA
tot.merge$mass_cat[tot.merge$mass_forearm_residual<0 & !is.na(tot.merge$mass_forearm_residual)] <- "poor"
tot.merge$mass_cat[tot.merge$mass_forearm_residual>0 & !is.na(tot.merge$mass_forearm_residual)] <- "good"
tot.merge$mass_cat[tot.merge$combo_class=="juvenile"] <- "juvenile"
#tot.merge$mass_cat[tot.merge$combo_class=="subadult"] <- "subadult"
tot.merge$mass_cat[tot.merge$combo_class=="female-P"] <- "pregnant"
tot.merge$mass_cat <- as.factor(tot.merge$mass_cat)

#Eidolon first
#GhV
lmGhV <- lmer(GhV_ct ~ 
                mass_cat + (1| year),
              REML=T,
              data= subset(tot.merge, bat_species =="Eidolon dupreanum" &  GhV_ct>0)) #singular fit, no need for random effect

lmGhV <- lm(GhV_ct ~ 
                mass_cat,
              data= subset(tot.merge, bat_species =="Eidolon dupreanum" &  GhV_ct>0)) #singular fit, no need for random effect
summary(lmGhV) # same idea, juveniles, poor body condition adults all have low Ct

plot_model(lmGhV , type = "est")

#AngV
lmAngV <- lmer(AngV_ct ~ 
                mass_cat + (1| year),
              REML=T,
              data= subset(tot.merge, bat_species =="Eidolon dupreanum" &  AngV_ct>0)) #singular fit, no need for random effect

lmAngV <- lm(AngV_ct ~ 
              mass_cat,
            data= subset(tot.merge, bat_species =="Eidolon dupreanum" &  AngV_ct>0)) 
summary(lmAngV) # nothing

#AndV
lmAndV <- lmer(AndV_ct ~ 
                mass_cat + (1| year),
              REML=T,
              data= subset(tot.merge, bat_species =="Eidolon dupreanum" &  AndV_ct>0)) 

summary(lmAndV) # nothing. 

# Rousettus

#GhV
lmRGhV <- lmer(GhV_ct ~ 
                mass_cat + (1| year),
              REML=T,
              data= subset(tot.merge, bat_species =="Rousettus madagascariensis" &  GhV_ct>0)) #singular fit, no need for random effect

summary(lmRGhV) # pregnant females with low Ct

plot_model(lmRGhV , type = "est")

#AngV
lmRAngV <- lmer(AngV_ct ~ 
                 mass_cat + (1| year),
               REML=T,
               data= subset(tot.merge, bat_species =="Rousettus madagascariensis" &  AngV_ct>0)) #singular fit, no need for random effect

summary(lmRAngV) # nothing

#AndV
# lmRAndV <- lmer(AndV_ct ~ 
#                  mass_cat + (1| year),
#                REML=T,
#                data= subset(tot.merge, bat_species =="Rousettus madagascariensis" &  AndV_ct>0)) 

lmRAndV <- lm(AndV_ct ~ 
                  mass_cat,
                data= subset(tot.merge, bat_species =="Rousettus madagascariensis" &  AndV_ct>0)) 

summary(lmRAndV) # nothing. 

# Pteropus

#GhV
lmPGhV <- lmer(GhV_ct ~ 
                 mass_cat + (1| year),
               REML=T,
               data= subset(tot.merge, bat_species =="Pteropus rufus" &  GhV_ct>0)) 

summary(lmPGhV) # nothing



#AngV
# lmPAngV <- lmer(AngV_ct ~ 
#                   mass_cat + (1| year),
#                 REML=T,
#                 data= subset(tot.merge, bat_species =="Pteropus rufus" &  AngV_ct>0))

#summary(lmPAngV) 

#AndV
lmPAndV <- lmer(AndV_ct ~ 
                  mass_cat + (1| year),
                REML=T,
                data= subset(tot.merge, bat_species =="Pteropus rufus" &  AndV_ct>0)) #singular fit

lmPAndV <- lm(AndV_ct ~ mass_cat,
                data= subset(tot.merge, bat_species =="Pteropus rufus" &  AndV_ct>0)) 

summary(lmPAndV) # nothing. 

#save lm output to csv file
#install.packages("broom.mixed", type="source")
library(broom.mixed)
model_list<-list(lmGhV=lmGhV, lmAngV= lmAngV, lmAndV= lmAndV,
                 lmPGhV=lmPGhV,  lmPAndV= lmPAndV,
                 lmRGhV=lmRGhV, lmRAngV= lmRAngV, lmRAndV= lmRAndV)
all_summaries <- bind_rows(lapply(model_list, tidy), .id = "Model_Name")
write.csv(all_summaries, "bodycondition_ct_lmer_model_summaries.csv", row.names = FALSE)


tot.merge$cat_plot <- as.character(tot.merge$mass_cat)
tot.merge$cat_plot[tot.merge$cat_plot=="good"] <- "adult,\nMFR>0"
tot.merge$cat_plot[tot.merge$cat_plot=="poor"] <- "adult,\nMFR<0"
tot.merge$cat_plot[tot.merge$cat_plot=="pregnant"] <- "pregnant\nfemale"
tot.merge$cat_plot <- factor(tot.merge$cat_plot, levels = c("adult,\nMFR>0", "adult,\nMFR<0", "pregnant\nfemale", "juvenile"))

#and plot Ct by condition/juvenile status - leave out pregnant females since there are so few
pMFR_Ct_GhV <- ggplot(data= subset(tot.merge, GhV_ct>0 & !is.na(cat_plot) & cat_plot!="pregnant\nfemale")) +
  geom_violin(aes(x=cat_plot, y=GhV_ct), draw_quantiles = c(0,.25,.5,.75)) + 
  geom_jitter(aes(x=cat_plot, y=GhV_ct, color=roost_site, shape=as.factor(year)), width = .2, show.legend = F,size=3) +
  theme_bw() + ylab("Ct (GhV)") + facet_grid(bat_species~.)+
  scale_y_reverse() +
  scale_shape_manual(name="year", values=yearz) +
  scale_color_manual(values = sitez) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=14), 
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 16, face="italic"),
        axis.title.x = element_blank(), 
        #axis.title.x = element_text(size=16, colour = "transparent"),
        axis.title.y = element_text(size=16)); pMFR_Ct_GhV

#and for the main text, just Eidolon
pMFR_Ct_GhV_Eid <- ggplot(data= subset(tot.merge, GhV_ct>0 & !is.na(cat_plot) & cat_plot!="pregnant\nfemale" & bat_species=="Eidolon dupreanum")) +
  geom_violin(aes(x=cat_plot, y=GhV_ct), draw_quantiles = c(0,.25,.5,.75)) + 
  geom_jitter(aes(x=cat_plot, y=GhV_ct, color=roost_site, shape=as.factor(year)), width = .2, show.legend = F,size=3) +
  theme_bw() + ylab("Ct (GhV)") + facet_grid(bat_species~.)+
  scale_y_reverse() +
  scale_shape_manual(name="year", values=yearz) +
  scale_color_manual(values = sitez) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=14), 
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_blank(), 
        #axis.title.x = element_text(size=16, colour = "transparent"),
        axis.title.y = element_text(size=16)); pMFR_Ct_GhV_Eid

dat_text <- data.frame(
  label = c(" p = 0.02 adult, MFR<0\n  p = 0.001 juvenile")
)

pMFR_Ct_GhV_Eid<-pMFR_Ct_GhV_Eid + geom_text(
  data    = dat_text,
    mapping = aes(x = 2.01, y = 27, label = label),
  hjust   = -0.1,
  vjust   = -1
);pMFR_Ct_GhV_Eid



#these are pretty sparse
pMFR_Ct_AngV <- ggplot(data= subset(tot.merge, AngV_ct>0 & !is.na(cat_plot) & cat_plot!="pregnant\nfemale")) +
  geom_violin(aes(x=cat_plot, y=AngV_ct), draw_quantiles = c(0,.25,.5,.75)) + 
  geom_jitter(aes(x=cat_plot, y=AngV_ct, color=roost_site, shape=as.factor(year)), width = .2, show.legend = F,size=3) +
  theme_bw() + ylab("Ct (AngV)") + facet_grid(bat_species~.)+
  scale_y_reverse() +
  scale_shape_manual(name="year", values=yearz) +
  scale_color_manual(values = sitez) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=14), 
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 16, face="italic"),
        axis.title.x = element_blank(), 
        #axis.title.x = element_text(size=16, colour = "transparent"),
        axis.title.y = element_text(size=16)); pMFR_Ct_AngV


#too few datapoints here to say much
pMFR_Ct_AndV <- ggplot(data= subset(tot.merge, AndV_ct>0 & !is.na(cat_plot) & cat_plot!="pregnant\nfemale")) +
  geom_violin(aes(x=cat_plot, y=AndV_ct), draw_quantiles = c(0,.25,.5,.75)) + 
  geom_jitter(aes(x=cat_plot, y=AndV_ct, color=roost_site, shape=as.factor(year)), width = .2, show.legend = F,size=3) +
  theme_bw() + ylab("Ct (AndV)") + facet_grid(bat_species~.)+
  scale_y_reverse() +
  scale_shape_manual(name="year", values=yearz) +
  scale_color_manual(values = sitez) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=14), 
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 16, face="italic"),
        axis.title.x = element_blank(), 
        #axis.title.x = element_text(size=16, colour = "transparent"),
        axis.title.y = element_text(size=16)); pMFR_Ct_AndV


#supplementary plot - not including
# fig2_supp_ct <- cowplot::plot_grid(pMFR_Ct_GhV, pMFR_Ct_AngV, pMFR_Ct_AndV, ncol = 3, nrow=1, labels = c("A", "B", "C"), label_size = 22, label_x = -.02)
# 
# 
# ggsave(file = paste0(homewd,"/final_figures_tables/tmp/Fig2X-Ct-supplement-draft.png"),
#        plot=fig2_supp_ct,
#        units="mm",  
#        width=120, 
#        height=90, 
#        scale=3, 
#        dpi=300)


#poor body condition and juvenile status are associated with lower Ct


#finally, what about coinfection?
tot.merge$sum_infections <- rowSums(cbind(tot.merge$GhV_pos, tot.merge$AngV_pos, tot.merge$AndV_pos))

length(tot.merge$sum_infections[tot.merge$sum_infections>1 & tot.merge$bat_species=="Eidolon dupreanum"]) #6
length(tot.merge$sum_infections[tot.merge$sum_infections>1 & tot.merge$bat_species=="Pteropus rufus"]) #2
length(tot.merge$sum_infections[tot.merge$sum_infections>1 & tot.merge$bat_species=="Rousettus madagascariensis"]) #3

length(tot.merge$sum_infections[tot.merge$sum_infections>2 & tot.merge$bat_species=="Eidolon dupreanum"]) #0
length(tot.merge$sum_infections[tot.merge$sum_infections>2 & tot.merge$bat_species=="Pteropus rufus"]) #0
length(tot.merge$sum_infections[tot.merge$sum_infections>2 & tot.merge$bat_species=="Rousettus madagascariensis"]) #0
#Eidolon
glmCoInf <- glm(sum_infections ~ combo_class,
          family = "poisson",
          data= subset(tot.merge,  bat_species=="Eidolon dupreanum"))
summary(glmCoInf) #nothing

glmCoInfCat <- glm(sum_infections ~ mass_cat,
                family = "poisson",
                data= subset(tot.merge,  bat_species=="Eidolon dupreanum"))
summary(glmCoInfCat) #juveniles marginally more associated with more HNV infections for E.dup

#Pteropus 
glmCoInfP <- glm(sum_infections ~ combo_class,
                family = "poisson",
                data= subset(tot.merge,  bat_species=="Pteropus rufus"))
summary(glmCoInfP) #nothing

glmCoInfCatP <- glm(sum_infections ~ mass_cat,
                   family = "poisson",
                   data= subset(tot.merge,  bat_species=="Pteropus rufus"))
summary(glmCoInfCatP) #nothing

#Rousettus
glmCoInfR <- glm(sum_infections ~ combo_class,
                 family = "poisson",
                 data= subset(tot.merge,  bat_species=="Rousettus madagascariensis"))
summary(glmCoInfR) #nothing

glmCoInfCatR <- glm(sum_infections ~ mass_cat,
                    family = "poisson",
                    data= subset(tot.merge,  bat_species=="Rousettus madagascariensis"))
summary(glmCoInfCatR) #poor body condition marginally associated with more HNVs

#save a csv file with the glm output: 
model_list<-list(glmCoInf=glmCoInf, glmCoInfCat= glmCoInfCat,
                 glmCoInfP=glmCoInfP, glmCoInfCatP= glmCoInfCatP,
                 glmCoInfR=glmCoInfR, glmCoInfCatR= glmCoInfCatR)
all_summaries <- bind_rows(lapply(model_list, tidy), .id = "Model_Name")
#write.csv(all_summaries, "coinfection_shed_glm_model_summaries.csv", row.names = FALSE)


#and the Ct-prevalence plot

# First, load the shedding data again
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


# link this to the catching data to plot by sampling session
catch.dat <- read.csv(file = paste0(homewd, "/master_data_files/catching_data.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))

head(catch.dat)
#names(catch.dat)
#catch.dat$net_date[!is.na(catch.dat$net_date)] <- as.Date(catch.dat$net_date[!is.na(catch.dat$net_date)], format = "%m/%d/%Y")
#catch.dat$net_date[is.na(catch.dat$net_date)] <- catch.dat$processing_date[is.na(catch.dat$net_date)]

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)
#names(catch.dat)[names(catch.dat)=="net_date"] <- "processing_date"

#and merge all three datasets
setdiff(shed.dat$sampleid, catch.dat$sampleid)
shed.dat= subset(shed.dat, sampleid!="MIZ305" & sampleid!="MIZ546")

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

#subset to eidolon
#merge.eid <- subset(merge.dat, bat_species=="Eidolon dupreanum")

#now get prevalence for every sampling session
#first get midpoint date for sampling session

sampling.session <- ddply(merge.dat,.(sampling_session), summarise, min_date =min(processing_date), max_date=max(processing_date))
sampling.session$mid_date <- sampling.session$min_date + ((sampling.session$max_date-sampling.session$min_date)/2)

merge.eid <- merge(merge.dat, sampling.session, by="sampling_session")
merge.eid$GhV_ct[merge.eid$GhV_ct==0] <- NA
merge.eid$AngV_ct[merge.eid$AngV_ct==0] <- NA
merge.eid$AndV_ct[merge.eid$AndV_ct==0] <- NA

shed.prev <- ddply(merge.eid,.(sampling_session, bat_species, roost_site, mid_date), summarise, NGhVpos = sum(GhV_pos, na.rm=T), NAngVpos = sum(AngV_pos, na.rm=T), NAndVpos = sum(AndV_pos, na.rm=T), Ntot=length(GhV_pos), meanGhVCt=mean(GhV_ct, na.rm=T), meanAngVCt=mean(AngV_ct, na.rm=T),  meanAndVCt=mean(AndV_ct, na.rm=T), sdGhVCt=sd(GhV_ct, na.rm=T), sdAngVCt=sd(AngV_ct, na.rm=T), sdAndVCt=sd(AndV_ct, na.rm=T))
head(shed.prev)
shed.prev$GhVprev <- shed.prev$NGhVpos/shed.prev$Ntot
shed.prev$AngVprev <- shed.prev$NAngVpos/shed.prev$Ntot
shed.prev$AndVprev <- shed.prev$NAndVpos/shed.prev$Ntot

shed.prev$NGhVneg <- shed.prev$Ntot - shed.prev$NGhVpos
shed.prev$NAngVneg <- shed.prev$Ntot - shed.prev$NAngVpos
shed.prev$NAndVneg <- shed.prev$Ntot - shed.prev$NAndVpos

shed.prev$GhVprevlci <- shed.prev$GhVprevuci <- shed.prev$AngVprevlci <- shed.prev$AngVprevuci <- shed.prev$AndVprevlci <- shed.prev$AndVprevuci <-  NA

for(i in 1:nrow(shed.prev)){
  tmp <- binom.test(x=c(shed.prev$NGhVpos[i], shed.prev$NGhVneg[i]))$conf.int  
  tmp2 <- binom.test(x=c(shed.prev$NAngVpos[i], shed.prev$NAngVneg[i]))$conf.int  
  tmp3 <- binom.test(x=c(shed.prev$NAndVpos[i], shed.prev$NAndVneg[i]))$conf.int  
  shed.prev$GhVprevlci[i] <- tmp[1]
  shed.prev$GhVprevuci[i] <- tmp[2]
  shed.prev$AngVprevlci[i] <- tmp2[1]
  shed.prev$AngVprevuci[i] <- tmp2[2]
  shed.prev$AndVprevlci[i] <- tmp3[1]
  shed.prev$AndVprevuci[i] <- tmp3[2]
}

head(shed.prev)

##GhV
#Eidolon
Eid.shed.prev <- glm(cbind(NGhVpos, NGhVneg) ~ meanGhVCt:roost_site,
                family = "binomial",
                data = subset(shed.prev, bat_species=="Eidolon dupreanum"))

Eid.shed.prev.nosite <- glm(cbind(NGhVpos, NGhVneg) ~ meanGhVCt,
                     family = "binomial",
                     data = subset(shed.prev, bat_species=="Eidolon dupreanum"))

summary(Eid.shed.prev) #all three sites sig negatively associated with prevalence and Ct, Moramanga strongest
summary(Eid.shed.prev.nosite) #with no site, sign is even stronger

Pter.shed.prev <- glm(cbind(NGhVpos, NGhVneg) ~ meanGhVCt:roost_site,
                     family = "binomial",
                     data = subset(shed.prev, bat_species=="Pteropus rufus"))

summary(Pter.shed.prev) #moramanga sig but opposite pattern

Rou.shed.prev <- glm(cbind(NGhVpos, NGhVneg) ~ meanGhVCt:roost_site,
                      family = "binomial",
                      data = subset(shed.prev, bat_species=="Rousettus madagascariensis"))

summary(Rou.shed.prev) #nothing sig


#and as a sensitivity analysis, try versions that weight those 
#sessions with more individuals contributing to Ct score

Eid.shed.prev.a <- glm(cbind(NGhVpos, NGhVneg) ~ meanGhVCt:roost_site,
           family = binomial,
           weights = NGhVpos,
           data = subset(shed.prev, bat_species=="Eidolon dupreanum"))

summary(Eid.shed.prev.a) ##all three sites are significant and neg as expected


Pter.shed.prev.a <- glm(cbind(NGhVpos, NGhVneg) ~ meanGhVCt:roost_site,
                       family = binomial,
                       weights = NGhVpos,
                       data = subset(shed.prev, bat_species=="Pteropus rufus"))

summary(Pter.shed.prev.a) ##both sites sig but positive which is weird
#but only 10 degrees of freedom



Rou.shed.prev.a <- glm(cbind(NGhVpos, NGhVneg) ~ meanGhVCt:roost_site,
                        family = binomial,
                        weights = NGhVpos,
                        data = subset(shed.prev, bat_species=="Rousettus madagascariensis"))

summary(Rou.shed.prev.a) ##moramanga sig but pos
#but only 10 degrees of freedom


# AngV
#Eidolon
Eid.shed.prev.AngV <- glm(cbind(NAngVpos, NAngVneg) ~ meanAngVCt:roost_site,
                     family = "binomial",
                     data = subset(shed.prev, bat_species=="Eidolon dupreanum"))

summary(Eid.shed.prev.AngV) #nothing

Pter.shed.prev.AngV <- glm(cbind(NAngVpos, NAngVneg) ~ meanAngVCt,
                      family = "binomial",
                      data = subset(shed.prev, bat_species=="Pteropus rufus"))

summary(Pter.shed.prev.AngV) #nothing

Rou.shed.prev.AngV <- glm(cbind(NAngVpos, NAngVneg) ~ meanAngVCt:roost_site,
                     family = "binomial",
                     data = subset(shed.prev, bat_species=="Rousettus madagascariensis"))

summary(Rou.shed.prev.AngV) #nothing sig


#and as a sensitivity analysis, try versions that weight those 
#sessions with more individuals contributing to Ct score

Eid.shed.prev.a.AngV <- glm(cbind(NAngVpos, NAngVneg) ~ meanAngVCt:roost_site,
                       family = binomial,
                       weights = NAngVpos,
                       data = subset(shed.prev, bat_species=="Eidolon dupreanum"))

summary(Eid.shed.prev.a.AngV) #nothing


Pter.shed.prev.a.AngV <- glm(cbind(NAngVpos, NAngVneg) ~ meanAngVCt,
                        family = binomial,
                        weights = NAngVpos,
                        data = subset(shed.prev, bat_species=="Pteropus rufus"))

summary(Pter.shed.prev.a.AngV) ##nothing




Rou.shed.prev.a.AngV <- glm(cbind(NAngVpos, NAngVneg) ~ meanAngVCt:roost_site,
                       family = binomial,
                       weights = NAngVpos,
                       data = subset(shed.prev, bat_species=="Rousettus madagascariensis"))

summary(Rou.shed.prev.a.AngV) ##nothing
#but only 10 degrees of freedom

# AndV
#Eidolon
Eid.shed.prev.AndV <- glm(cbind(NAndVpos, NAndVneg) ~ meanAndVCt,
                          family = "binomial",
                          data = subset(shed.prev, bat_species=="Eidolon dupreanum"))

summary(Eid.shed.prev.AndV) #nothing

Pter.shed.prev.AndV <- glm(cbind(NAndVpos, NAndVneg) ~ meanAndVCt,
                           family = "binomial",
                           data = subset(shed.prev, bat_species=="Pteropus rufus"))

summary(Pter.shed.prev.AndV) #nothing

Rou.shed.prev.AndV <- glm(cbind(NAndVpos, NAndVneg) ~ meanAndVCt:roost_site,
                          family = "binomial",
                          data = subset(shed.prev, bat_species=="Rousettus madagascariensis"))

summary(Rou.shed.prev.AndV) #nothing sig


#and as a sensitivity analysis, try versions that weight those 
#sessions with more individuals contributing to Ct score

Eid.shed.prev.a.AndV <- glm(cbind(NAndVpos, NAndVneg) ~ meanAndVCt,
                            family = binomial,
                            weights = NAndVpos,
                            data = subset(shed.prev, bat_species=="Eidolon dupreanum"))

summary(Eid.shed.prev.a.AndV) #nothing


Pter.shed.prev.a.AndV <- glm(cbind(NAndVpos, NAndVneg) ~ meanAndVCt,
                             family = binomial,
                             weights = NAndVpos,
                             data = subset(shed.prev, bat_species=="Pteropus rufus"))

summary(Pter.shed.prev.a.AndV) ##nothing




Rou.shed.prev.a.AndV <- glm(cbind(NAndVpos, NAndVneg) ~ meanAndVCt:roost_site,
                            family = binomial,
                            weights = NAndVpos,
                            data = subset(shed.prev, bat_species=="Rousettus madagascariensis"))

summary(Rou.shed.prev.a.AndV) ##nothing


#save csv file of glm statistics
model_list<-list(Eid.shed.prev=Eid.shed.prev, Eid.shed.prev.a= Eid.shed.prev.a, Rou.shed.prev=Rou.shed.prev, Rou.shed.prev.a= Rou.shed.prev.a, Pter.shed.prev=Pter.shed.prev, Pter.shed.prev.a= Pter.shed.prev.a,
                 Eid.shed.prev.AngV=Eid.shed.prev.AngV, Eid.shed.prev.a.AngV= Eid.shed.prev.a.AngV, Rou.shed.prev.AngV=Rou.shed.prev.AngV, Rou.shed.prev.a.AngV= Rou.shed.prev.a.AngV, Pter.shed.prev.AngV=Pter.shed.prev.AngV, Pter.shed.prev.a.AngV= Pter.shed.prev.a.AngV,
                 Eid.shed.prev.AndV=Eid.shed.prev.AndV, Eid.shed.prev.a.AndV= Eid.shed.prev.a.AndV, Rou.shed.prev.AndV=Rou.shed.prev.AndV, Rou.shed.prev.a.AndV= Rou.shed.prev.a.AndV, Pter.shed.prev.AndV=Pter.shed.prev.AndV, Pter.shed.prev.a.AndV= Pter.shed.prev.a.AndV)
all_summaries <- bind_rows(lapply(model_list, tidy), .id = "Model_Name")
write.csv(all_summaries, "weighted_unweighted_shed_session_glm_model_summaries.csv", row.names = FALSE)


#now try a permutation test - shuffle Cts across sampling sessions within a site
#first drop mahabo and makira because they had only one visit


weighted_cor <- function(x, y, w) {
  #rx <- rank(x, ties.method = "average")
  #ry <- rank(y, ties.method = "average")
  cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]
} #pearson
perm_cor_weighted <- function(df, virus, n_perm = 10000, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  if(virus=="GhV"){
    
    
    df=subset(df, !is.na(meanGhVCt))
    w <- df$NGhVpos  # or n_bats
    obs <- weighted_cor(df$meanGhVCt, df$GhVprev, w)
    
    perm_stats <- replicate(n_perm, {
      
      df_perm <- df
      
      for (s in unique(df$roost_site)) {
        idx <- which(df$roost_site == s)
        df_perm$meanGhVCt[idx] <- sample(df$meanGhVCt[idx])
      }
      
      weighted_cor(df_perm$meanGhVCt, df_perm$GhVprev, w)
    })
    
  }else if (virus=="AngV"){
    
    df=subset(df, !is.na(meanAngVCt))
    w <- df$NAngVpos  # or n_bats
    obs <- weighted_cor(df$meanAngVCt, df$AngVprev, w)
    
    perm_stats <- replicate(n_perm, {
      
      df_perm <- df
      
      for (s in unique(df$roost_site)) {
        idx <- which(df$roost_site == s)
        df_perm$meanAngVCt[idx] <- sample(df$meanAngVCt[idx])
      }
      
      weighted_cor(df_perm$meanAngVCt, df_perm$AngVprev, w)
    })
    
  }else if (virus=="AndV"){
    
    df=subset(df, !is.na(meanAndVCt))
    w <- df$NAndVpos  # or n_bats
    obs <- weighted_cor(df$meanAndVCt, df$AndVprev, w)
    
    perm_stats <- replicate(n_perm, {
      
      df_perm <- df
      
      for (s in unique(df$roost_site)) {
        idx <- which(df$roost_site == s)
        df_perm$meanAndVCt[idx] <- sample(df$meanAndVCt[idx])
      }
      
      weighted_cor(df_perm$meanAndVCt, df_perm$AndVprev, w)
    })
    
  }
  
  perm_stats <- perm_stats[!is.na(perm_stats)]
  p_val <- mean(abs(perm_stats) >= abs(obs))
  
  list(observed = obs, p_value = p_val)
}
perm_cor_unweighted <- function(df, virus, n_perm = 10000, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  if(virus=="GhV"){
    
    
    df=subset(df, !is.na(meanGhVCt))
    #w <- df$NGhVpos  # or n_bats
    obs <- cor(df$meanGhVCt, df$GhVprev, method = "pearson")
    
    perm_stats <- replicate(n_perm, {
      
      df_perm <- df
      
      for (s in unique(df$roost_site)) {
        idx <- which(df$roost_site == s)
        df_perm$meanGhVCt[idx] <- sample(df$meanGhVCt[idx])
      }
      
      cor(df_perm$meanGhVCt, df_perm$GhVprev,method = "pearson")
    })
    
  }else if (virus=="AngV"){
    
    df=subset(df, !is.na(meanAngVCt))
    #w <- df$NAngVpos  # or n_bats
    obs <- cor(df$meanAngVCt, df$AngVprev, method = "pearson")
    
    perm_stats <- replicate(n_perm, {
      
      df_perm <- df
      
      for (s in unique(df$roost_site)) {
        idx <- which(df$roost_site == s)
        df_perm$meanAngVCt[idx] <- sample(df$meanAngVCt[idx])
      }
      
      cor(df_perm$meanAngVCt, df_perm$AngVprev, method = "pearson")
    })
    
  }else if (virus=="AndV"){
    
    df=subset(df, !is.na(meanAndVCt))
    #w <- df$NAngVpos  # or n_bats
    obs <- cor(df$meanAndVCt, df$AndVprev, method = "pearson")
    
    perm_stats <- replicate(n_perm, {
      
      df_perm <- df
      
      for (s in unique(df$roost_site)) {
        idx <- which(df$roost_site == s)
        df_perm$meanAndVCt[idx] <- sample(df$meanAndVCt[idx])
      }
      
      cor(df_perm$meanAndVCt, df_perm$AndVprev, method = "pearson")
    })
    
  }
  
  p_val <- mean(abs(perm_stats) >= abs(obs))
  
  list(observed = obs, p_value = p_val)
}

#Eidolon
perm_cor_weighted(df = subset(shed.prev, roost_site!="Mahabo" & bat_species=="Eidolon dupreanum"), virus="GhV",n_perm = 10000)
#$observed: -0.4468665; $p_value=  0.0402  significant
perm_cor_weighted(df = subset(shed.prev, roost_site!="Mahabo" & bat_species=="Eidolon dupreanum"), virus="AngV",n_perm = 10000)
#$observed: -0.3431113; $p_value= 0.6783. not significant
perm_cor_weighted(df = subset(shed.prev, roost_site!="Mahabo" & bat_species=="Eidolon dupreanum"), virus="AndV",n_perm = 10000)
#$observed: -0.9470351 $p_value = 0.0417, sig!

#Pteropus
perm_cor_weighted(df = subset(shed.prev, roost_site!="Mahabo" & roost_site!="Makira" & bat_species=="Pteropus rufus"), virus="GhV",n_perm = 10000)
#$observed: 0.408861; $p_value=  0.2261.  not significant 
perm_cor_weighted(df = subset(shed.prev, roost_site!="Mahabo" & roost_site!="Makira" & bat_species=="Pteropus rufus"), virus="AngV",n_perm = 10000)
#$observed: -0.9959159; $p_value= 0.3369. almost significant
perm_cor_weighted(df = subset(shed.prev, roost_site!="Mahabo" & roost_site!="Makira" & bat_species=="Pteropus rufus"), virus="AndV",n_perm = 10000)
#$observed: 0.08487054 $p_value = 0.8459, not sig


#Rousettus
perm_cor_weighted(df = subset(shed.prev, roost_site!="Mahabo" & roost_site!="Makira" & bat_species=="Rousettus madagascariensis"), virus="GhV",n_perm = 10000)
#$observed: 0.5521216; $p_value= 0.0292.  significant but positively correlated again
perm_cor_weighted(df = subset(shed.prev, roost_site!="Mahabo" & roost_site!="Makira" & bat_species=="Rousettus madagascariensis"), virus="AngV",n_perm = 10000)
#$observed: -0.8819563; $p_value= 0.1691. almost significant
perm_cor_weighted(df = subset(shed.prev, roost_site!="Mahabo" & roost_site!="Makira" & bat_species=="Rousettus madagascariensis"), virus="AndV",n_perm = 10000)
#NA



#and unweighted
#Eidolon
perm_cor_unweighted(df = subset(shed.prev, roost_site!="Mahabo" & bat_species=="Eidolon dupreanum"), virus="GhV",n_perm = 10000)
#$observed: -0.2951848; $p_value= 0.1552. not sig
perm_cor_unweighted(df = subset(shed.prev, roost_site!="Mahabo" & bat_species=="Eidolon dupreanum"), virus="AngV",n_perm = 10000)
#$observed: -0.3294097; $p_value=0.7137. not significant
perm_cor_unweighted(df = subset(shed.prev, roost_site!="Mahabo" & bat_species=="Eidolon dupreanum"), virus="AndV",n_perm = 10000)
#$observed: -0.94703511 $p_value = 0.044, sig!

#Pteropus
perm_cor_unweighted(df = subset(shed.prev, roost_site!="Mahabo" & roost_site!="Makira" & bat_species=="Pteropus rufus"), virus="GhV",n_perm = 10000)
#$observed: 0.4132313; $p_value= 0.186.  not significant 
perm_cor_unweighted(df = subset(shed.prev, roost_site!="Mahabo" & roost_site!="Makira" & bat_species=="Pteropus rufus"), virus="AngV",n_perm = 10000)
#$observed: -0.5273378; $p_value= 0.2253. almost significant
perm_cor_unweighted(df = subset(shed.prev, roost_site!="Mahabo" & roost_site!="Makira" & bat_species=="Pteropus rufus"), virus="AndV",n_perm = 10000)
#$observed: 0.1182298 $p_value = 0.7483, not sig


#Rousettus
perm_cor_unweighted(df = subset(shed.prev, roost_site!="Mahabo" & roost_site!="Makira" & bat_species=="Rousettus madagascariensis"), virus="GhV",n_perm = 10000)
#$observed: 0.4500231; $p_value= 0.1133.  significant but positively correlated again
perm_cor_unweighted(df = subset(shed.prev, roost_site!="Mahabo" & roost_site!="Makira" & bat_species=="Rousettus madagascariensis"), virus="AngV",n_perm = 10000)
#$observed: -0.6673865; $p_value= 0.1152. almost significant
perm_cor_unweighted(df = subset(shed.prev, roost_site!="Mahabo" & roost_site!="Makira" & bat_species=="Rousettus madagascariensis"), virus="AndV",n_perm = 10000)
#NA

# now make a plot for the actual paper (either supp or one component of main):
shed.prev$GhVCt_uci <- shed.prev$meanGhVCt+1.96*shed.prev$sdGhVCt
shed.prev$GhVCt_lci <- shed.prev$meanGhVCt-1.96*shed.prev$sdGhVCt

shed.prev$AngVCt_uci <- shed.prev$meanAngVCt+1.96*shed.prev$sdAngVCt
shed.prev$AngVCt_lci <- shed.prev$meanAngVCt-1.96*shed.prev$sdAngVCt

#these don't work because no sd because only 1 positive per capture
shed.prev$AndVCt_uci <- shed.prev$meanAndVCt+1.96*shed.prev$sdAndVCt
shed.prev$AndVCt_lci <- shed.prev$meanAndVCt-1.96*shed.prev$sdAndVCt

head(shed.prev)

pred.df <- cbind.data.frame(meanGhVCt = rep(25:39, 2), roost_site=rep(c("Moramanga", "Ankarana"), each=15))
pred.df$bat_species <- "Eidolon dupreanum"
pred.df1 <- pred.df
pred.df1$bat_species <- "Pteropus rufus"
pred.df1 = subset(pred.df1, roost_site=="Moramanga")
pred.df2 <- pred.df
pred.df2$bat_species <- "Rousettus madagascariensis"
pred.df <- rbind(pred.df, pred.df1, pred.df2)
pred.df$prevalence <- pred.df$prevalence_lci <- pred.df$prevalence_uci <- NA

pred.df$prevalence[pred.df$bat_species=="Eidolon dupreanum"] <- predict(Eid.shed.prev, newdata = pred.df[pred.df$bat_species=="Eidolon dupreanum",], type = "response")
pred.df$prevalence_lci[pred.df$bat_species=="Eidolon dupreanum"] <- predict(Eid.shed.prev, newdata = pred.df[pred.df$bat_species=="Eidolon dupreanum",], type = "response") -1.96*predict(Eid.shed.prev, newdata = pred.df[pred.df$bat_species=="Eidolon dupreanum",], type = "response", se.fit = T)$se
pred.df$prevalence_uci[pred.df$bat_species=="Eidolon dupreanum"] <- predict(Eid.shed.prev, newdata = pred.df[pred.df$bat_species=="Eidolon dupreanum",], type = "response") +1.96*predict(Eid.shed.prev, newdata = pred.df[pred.df$bat_species=="Eidolon dupreanum",], type = "response", se.fit = T)$se


pred.df$prevalence[pred.df$bat_species=="Pteropus rufus"] <- predict(Pter.shed.prev, newdata = pred.df[pred.df$bat_species=="Pteropus rufus",], type = "response")
pred.df$prevalence_lci[pred.df$bat_species=="Pteropus rufus"] <- predict(Pter.shed.prev, newdata = pred.df[pred.df$bat_species=="Pteropus rufus",], type = "response") -1.96*predict(Pter.shed.prev, newdata = pred.df[pred.df$bat_species=="Pteropus rufus",], type = "response", se.fit = T)$se
pred.df$prevalence_uci[pred.df$bat_species=="Pteropus rufus"] <- predict(Pter.shed.prev, newdata = pred.df[pred.df$bat_species=="Pteropus rufus",], type = "response") +1.96*predict(Pter.shed.prev, newdata = pred.df[pred.df$bat_species=="Pteropus rufus",], type = "response", se.fit = T)$se



pred.df$prevalence[pred.df$bat_species=="Rousettus madagascariensis"] <- predict(Rou.shed.prev, newdata = pred.df[pred.df$bat_species=="Rousettus madagascariensis",], type = "response")
pred.df$prevalence_lci[pred.df$bat_species=="Rousettus madagascariensis"] <- predict(Rou.shed.prev, newdata = pred.df[pred.df$bat_species=="Rousettus madagascariensis",], type = "response") -1.96*predict(Rou.shed.prev, newdata = pred.df[pred.df$bat_species=="Rousettus madagascariensis",], type = "response", se.fit = T)$se
pred.df$prevalence_uci[pred.df$bat_species=="Rousettus madagascariensis"] <- predict(Rou.shed.prev, newdata = pred.df[pred.df$bat_species=="Rousettus madagascariensis",], type = "response") +1.96*predict(Rou.shed.prev, newdata = pred.df[pred.df$bat_species=="Rousettus madagascariensis",], type = "response", se.fit = T)$se



pred.df$prevalence_lci[pred.df$prevalence_lci <0] <- 0
pred.df$prevalence_uci[pred.df$prevalence_uci >1] <- 1

shed.prev$year <- year(shed.prev$mid_date)


dat_text <- data.frame(
  label = c("p = 0 Ankarana\np = 0 Moramanga"),
  x     = c(Inf),
  y     = c(Inf))

pCtprev <- ggplot(data= subset(shed.prev, !is.na(meanGhVCt))) + 
  geom_point(aes(x=meanGhVCt, y=GhVprev, color=roost_site, size=Ntot, shape = as.factor(year)), show.legend = F) + theme_bw() +
  geom_linerange(aes(x=meanGhVCt, ymin=GhVprevlci, ymax=GhVprevuci, color=roost_site), linetype=3, alpha=.5, show.legend = F) +
  geom_linerange(aes(xmin=GhVCt_lci, xmax=GhVCt_uci, y=GhVprev, color=roost_site), linetype=3, alpha=.5, show.legend = F) +
  facet_grid(bat_species~.) +
  xlab("Ct (GhV)") + ylab("prevalence (GhV)") +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  scale_color_manual(name="roost site", values = sitez) +
  scale_shape_manual(values=yearz, name="year") +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=16, face="italic"),
        axis.text = element_text(size=14), 
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.text = element_text(size=12),
        legend.position = "inside",
        legend.position.inside = c(.15,.73),
        legend.title = element_text(size=14));pCtprev
#   
# #guides(color=guide_legend(ncol=1), shape="none", size="none")
# pCtprev<-pCtprev+geom_text(data= dat_text,
#                            mapping = aes(x = x, y = y, label = label, hjust=1.1, vjust=1.5, fontface="italic"));pCtprev

pCtprev_line <- ggplot(data= subset(shed.prev, !is.na(meanGhVCt))) + 
  geom_point(aes(x=meanGhVCt, y=GhVprev, color=roost_site, size=Ntot, shape = as.factor(year)), show.legend = F) + theme_bw() +
  geom_linerange(aes(x=meanGhVCt, ymin=GhVprevlci, ymax=GhVprevuci, color=roost_site), linetype=3, alpha=.5, show.legend = F) +
  geom_linerange(aes(xmin=GhVCt_lci, xmax=GhVCt_uci, y=GhVprev, color=roost_site), linetype=3, alpha=.5, show.legend = F) +
  geom_line(data=pred.df, aes(x=meanGhVCt, y=prevalence, color = roost_site)) +
  geom_ribbon(data=pred.df, aes(x=meanGhVCt, ymin= prevalence_lci, ymax=prevalence_uci, fill=roost_site), alpha=.3) +
  facet_grid(bat_species~.) +
  xlab("Ct (GhV)") + ylab("prevalence (GhV)") +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  scale_color_manual(name="roost site", values = sitez) +
  scale_fill_manual(name="roost site", values = sitez) +
  scale_shape_manual(values=yearz, name="year") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=14), 
        strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=16, face="italic"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.text = element_text(size=12),
        legend.position = "none",
        legend.position.inside = c(.15,.73),
        legend.title = element_text(size=14));pCtprev_line


#now, for Fig 1, focus on just the Eidolon data
pCtprev_line_Eid <- ggplot(data= subset(shed.prev, !is.na(meanGhVCt) & bat_species=="Eidolon dupreanum")) + 
  geom_point(aes(x=meanGhVCt, y=GhVprev, color=roost_site, size=Ntot, shape = as.factor(year)), show.legend = F) + theme_bw() +
  geom_linerange(aes(x=meanGhVCt, ymin=GhVprevlci, ymax=GhVprevuci, color=roost_site), linetype=3, alpha=.5, show.legend = F) +
  geom_linerange(aes(xmin=GhVCt_lci, xmax=GhVCt_uci, y=GhVprev, color=roost_site), linetype=3, alpha=.5, show.legend = F) +
  geom_line(data=subset(pred.df, bat_species=="Eidolon dupreanum"), aes(x=meanGhVCt, y=prevalence, color = roost_site)) +
  geom_ribbon(data=subset(pred.df, bat_species=="Eidolon dupreanum"), aes(x=meanGhVCt, ymin= prevalence_lci, ymax=prevalence_uci, fill=roost_site), alpha=.3) +
  facet_grid(bat_species~.) +
  xlab("Ct (GhV)") + ylab("prevalence (GhV)") +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  scale_color_manual(name="roost site", values = sitez) +
  scale_fill_manual(name="roost site", values = sitez) +
  scale_shape_manual(values=yearz, name="year") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=14), 
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.text = element_text(size=12),
        legend.position = "none",
        legend.position.inside = c(.15,.73),
        legend.title = element_text(size=14));pCtprev_line_Eid

dat_text <- data.frame(
  label = c(" p = 0.09 Ankarana\n  p = 0.07 Mahabo\n p = 0.04 Moramanga")
)

pCtprev_line_Eid<-pCtprev_line_Eid + geom_text(
  data    = dat_text,
  mapping = aes(x = 32, y = 0.75, label = label),
  hjust   = -0.1,
  vjust   = -1
);pCtprev_line_Eid



# 
# #guides(color=guide_legend(ncol=1), shape="none", size="none")
# pCtprev_line<-pCtprev_line+geom_text(data= dat_text,
#                            mapping = aes(x = x, y = y, label = label, hjust=1.1, vjust=1.5, fontface="italic"));pCtprev_line


#and all together for Fig 2
Fig2_left <- cowplot::plot_grid(pIntrannualShedEid_noLeg, pMFR_Ct_GhV_Eid, pCtprev_line_Eid,  ncol=3, nrow=1, labels = c("A", "B", "C", ""), rel_widths = c(1,.9,.9), label_size = 22, label_x = -.01); Fig2_left
Fig2 <- cowplot::plot_grid(Fig2_left, Eid_Leg, ncol=2, nrow=1, rel_widths = c(1,.1), label_size = 22); Fig2

Fig2 <- ggdraw(Fig2) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file = paste0(homewd,"/final_figures_tables/main/Fig2.png"),
       plot=Fig2,
       units="mm",  
       width=140, 
       height=50, 
       scale=3, 
       dpi=300)
