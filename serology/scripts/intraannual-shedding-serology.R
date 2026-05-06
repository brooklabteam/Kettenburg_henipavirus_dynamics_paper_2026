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
batwd = "/Users/carabrook/Developer/Madagascar-Bat-Data"

# 
# homewd = "/Users/gwenddolenkettenburg/Desktop/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
# batwd = "/Users/gwenddolenkettenburg/Desktop/Developer/Madagascar-Bat-Data"

# plot of seasonality of antigens for analysis
# also GAMs of predictors of seropositivity and/or seroprevalence
################################################################################
################################################################################
## Section 1. Prep the data

# First, load the luminex data after determining who is seropositive/
# seronegative in the other scripts
lum.dat <- read.csv(file = paste0(homewd, "/serology/data/mix_sub_data_75_90_96_long.csv"))
head(lum.dat)

lum.dat <- dplyr::select(lum.dat, sampleid, antigen, MFI, log10MFI, seropos, seropos_lci, seropos_uci)
## merge the luminex data with the catching data


# link this to the catching data to plot by sampling session
catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))
#catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/Historic-Catching-Data-Through2020/catching_data_historic.csv"))

head(catch.dat)
names(catch.dat)

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)

#and merge all three datasets
merge.dat <- merge(lum.dat, catch.dat, by="sampleid", all.x = T)
head(merge.dat)


#get the midpoint processing date per sampling session
merge.dat$processing_date <- as.Date(merge.dat$processing_date)
sampling.date <- ddply(merge.dat, .(bat_species, sampling_session), summarise, min_date = min(processing_date), max_date = max(processing_date))
sampling.date$mid_date <- (sampling.date$max_date-sampling.date$min_date)/2 + sampling.date$min_date

#and slim to merge
sampling.merge <- dplyr::select(sampling.date, bat_species, sampling_session, mid_date)

merge.dat <- merge(merge.dat, sampling.merge, by = c("bat_species", "sampling_session"))
merge.dat$year <- year(merge.dat$processing_date)
merge.dat$month <- as.character(month(merge.dat$processing_date))
merge.dat$month[merge.dat$month!="10"& merge.dat$month!="11" &merge.dat$month!="12"]<-paste0("0", merge.dat$month[merge.dat$month!="10"& merge.dat$month!="11" &merge.dat$month!="12"])
#merge.dat$mid_date <- as.Date(paste0(merge.dat$year, "-", merge.dat$month, "-15"))
head(merge.dat)

#just do the recent date
merge.dat <- subset(merge.dat, mid_date>"2018-01-01")


#change roost site to just Moramanga region or Ankarana
unique(merge.dat$roost_site)
merge.dat$roost_site[merge.dat$roost_site=="Angavobe_Edup" |merge.dat$roost_site=="Angavokely_Edup"| merge.dat$roost_site=="Maromizaha_Rmad"] <- "Moramanga"
merge.dat$roost_site[merge.dat$roost_site!="Moramanga"] <- "Ankarana"


#also add in quantitative age data where available:
age.dat <- read.csv(file = paste0(batwd, "/Age-Data/aged_bats.csv"), header = T, stringsAsFactors = F)
head(age.dat)
#slim to just the columns you want to add
age.dat <- dplyr::select(age.dat, sampleid, age)

#and merge
merge.dat <- merge(merge.dat, age.dat, by="sampleid", all.x = T)
head(merge.dat)

# #write over those that are super close
merge.dat$mid_date <- as.character(merge.dat$mid_date)
merge.dat$mid_date[merge.dat$bat_species=="Eidolon dupreanum" & merge.dat$roost_site=="Moramanga" & merge.dat$mid_date=="2018-02-15"] <- "2018-02-14"
merge.dat$mid_date[merge.dat$bat_species=="Rousettus madagascariensis" & merge.dat$roost_site=="Moramanga" & merge.dat$mid_date=="2019-01-09"] <- "2019-01-05"
merge.dat$mid_date <- as.Date(merge.dat$mid_date)

merge.dat$doy <- yday(merge.dat$processing_date)

####################################################################################################
####################################################################################################
head(merge.dat)

## add columns for day of year and year -- this time we use the processing date for the 
## specific bat rather than the mid_date of the sampling session


merge.dat$year <- year(merge.dat$processing_date)
merge.dat$year <- as.factor(merge.dat$year)
merge.dat$bat_sex <- as.factor(merge.dat$bat_sex)
merge.dat$roost_site <- as.factor(merge.dat$roost_site)

# remove two confusing bats from the dataset (duplicate sample id: MIZ540, MIZ 598)
merge.dat = subset(merge.dat, !is.na(bat_age_class))

# and remove that one juvenile for which sex is unknown
merge.dat = subset(merge.dat, bat_sex!="unknown")

# make a variable for repro class 

merge.dat$repro_class <- merge.dat$bat_age_class
unique(merge.dat$repro_class)


merge.dat$antigen <- as.factor(merge.dat$antigen)

#now also add an age class
merge.dat$age_class <- merge.dat$repro_class
merge.dat$age_class[merge.dat$age_class!="J"] <- "A"

#make into factors
merge.dat$repro_class <- as.factor(merge.dat$repro_class)
merge.dat$age_class <- as.factor(merge.dat$age_class)

#and the old status 
merge.dat$bat_age_class <- as.factor(merge.dat$bat_age_class)

#make combo class
merge.dat$combo_class <- paste0(merge.dat$bat_sex, "_", merge.dat$repro_class)

#juveniles will be lumped together
merge.dat$combo_class[merge.dat$combo_class=="female_J" | merge.dat$combo_class=="male_J"] <- "juvenile"
unique(merge.dat$combo_class)

merge.dat$combo_class <- as.factor(merge.dat$combo_class)


### there are around 1300 juveniles and YOY in the dataset that have no 
### mass/forearm residual because it is unfair to categorize them one way
### or the other so, we give them a score of '0' for this value
merge.dat$mass_forearm_residual[is.na(merge.dat$mass_forearm_residual)] <- 0


### droplevels on sex
merge.dat$bat_sex <- droplevels(merge.dat$bat_sex)

merge.dat$repro_class <- as.character(merge.dat$repro_class)
merge.dat$repro_class[!is.na(merge.dat$age) & merge.dat$age<1] <- "J"
merge.dat$repro_class <- as.factor(merge.dat$repro_class )

# now, that the data are ready, we can test predictors of seropositive
# status in the different bat species for each antigen

# Eidolon first 
merge.dat$roost_site
#GhV
# no seasonality by sex
m1GhVnoSeX <- gam(seropos ~ 
                    s(doy, k=10, bs="cc", by=roost_site) +
                    #s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex) +
                    #s(combo_class, bs="re") +
                    #s(roost_site, bs="re") +
                    s(year, bs = "re"), 
                  dat = subset(merge.dat, bat_species=="Eidolon dupreanum"  & antigen=="GhV"),
                  select = TRUE,
                  method = "REML",
                  gamma = 1,
                  family = "binomial") 
summary(m1GhVnoSeX)
#significant seasonality by doy for moramanga
plot(m1GhVnoSeX, pages = 1)

#with sex in seasonality
m1GhVSeX <- gam(seropos ~
                  s(doy, k=10, bs="cc", by=interaction(roost_site, bat_sex)) +
                  #s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                  #s(combo_class, bs="re") +
                  #s(age, bs="cr") +
                  #s(roost_site, bs="re") +
                  s(year, bs = "re"), 
                dat = subset(merge.dat, bat_species=="Eidolon dupreanum"  & antigen=="GhV"),
                select = TRUE,
                method = "REML",
                gamma = 1,
                family = "binomial") 
summary(m1GhVSeX)
#significant seasonality by doy in Moramanga males 
plot(m1GhVSeX, pages = 1)


AIC(m1GhVnoSeX, m1GhVSeX)
#simpler no sex model is favored though they are close

#AngV
m1AngVnoSeX <- gam(seropos ~
                     s(doy,k=10, bs="cc", by =roost_site)+
                     #s(doy,k=-1, bs="cc")+
                     #s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                     #s(combo_class, bs="re") +
                     #s(age, bs="cr") +
                     #s(roost_site, bs="re") +
                     s(year, bs = "re"), 
                   dat = subset(merge.dat, bat_species=="Eidolon dupreanum"  & antigen=="AngV"),
                   select = TRUE,
                   method = "REML",
                   gamma = 1,
                   family = "binomial") 
summary(m1AngVnoSeX) #no sig
plot(m1AngVnoSeX, pages = 1)

m1AngVSeX <- gam(seropos ~
                   s(doy,k=10, bs="cc", by=interaction(bat_sex, roost_site))+
                   #s(doy,k=-1, bs="cc")+
                   #s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                   #s(combo_class, bs="re") +
                   #s(age, bs="cr") +
                   #s(roost_site, bs="re") +
                   s(year, bs = "re"), 
                 dat = subset(merge.dat, bat_species=="Eidolon dupreanum"  & antigen=="AngV"),
                 select = TRUE,
                 method = "REML",
                 gamma = 1,
                 family = "binomial") 
summary(m1AngVSeX)# no sig
plot(m1AngVSeX, pages = 1)


AIC(m1AngVnoSeX, m1AngVSeX)
# sex improves model fit very slightly but not enough to justify


# and Rousettus models

#GhV
# no seasonality by sex
m2GhVnoSeX <- gam(seropos ~ 
                    s(doy,k=10, bs="cc", by=roost_site)+
                    #s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                    #s(combo_class, bs="re") +
                    #s(roost_site, bs="re") +
                    s(year, bs = "re"), 
                  dat = subset(merge.dat, bat_species=="Rousettus madagascariensis"  & antigen=="GhV"),
                  select = TRUE,
                  method = "REML",
                  gamma = 1,
                  family = "binomial") 
summary(m2GhVnoSeX)
#Nothing significant
plot(m2GhVnoSeX, pages = 1)

#with sex in seasonality
m2GhVSeX <- gam(seropos ~
                  s(doy,k=10, bs="cc", by=interaction(bat_sex, roost_site))+
                  #s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                  #s(combo_class, bs="re") +
                  #s(age, bs="cr") +
                  #s(roost_site, bs="re") +
                  s(year, bs = "re"), 
                dat = subset(merge.dat, bat_species=="Rousettus madagascariensis"  & antigen=="GhV"),
                select = TRUE,
                method = "REML",
                gamma = 1,
                family = "binomial") 
summary(m2GhVSeX)
#nothing significant
plot(m2GhVSeX, pages = 1)

AIC(m2GhVnoSeX, m2GhVSeX) #same but can't justify adding sex


# seasonal plotting

pred.df.Eid <- cbind.data.frame(doy=rep(1:365,3*4), year=rep(c(2018,2019,2020), each=(365*4)), antigen = rep(rep(c("GhV", "HeV", "AngV", "NiV"), each = 365), 3) )
#pred.df.Eid <- cbind.data.frame(doy=rep(1:365,3*4), year=rep(c(2018,2019,2020), each=(365*4)), antigen = rep(rep(c("GhV", "HeV", "AngV", "NiV"), each = 365), 3) )
pred.df.Eid$bat_sex = "female"
pred.df.tmp <- pred.df.Eid
pred.df.tmp$bat_sex <- "male"
pred.df.Eid <- rbind(pred.df.Eid, pred.df.tmp)
pred.df.Eid$roost_site="Moramanga"
pred.df.Eid$mass_forearm_residual = 0
pred.df.Eid$combo_class <- "male_A"
pred.df.Eid$combo_class[pred.df.Eid$bat_sex=="female"] <- "female_NL"

#predict without effect of year for plotting
#females
pred.df.Eid$predicted_seropos <- NA
pred.df.Eid$predicted_seropos_lci <- NA
pred.df.Eid$predicted_seropos_uci <- NA

#GhV
pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="GhV"] <- predict.gam(m1GhVnoSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="GhV",], type="response", exclude = "s(year)")
pred.df.Eid$predicted_seropos_lci[pred.df.Eid$antigen=="GhV"] <- pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="GhV"] - 1.96*predict.gam(m1GhVnoSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="GhV",], type="response", se.fit = T, exclude = "s(year)")$se
pred.df.Eid$predicted_seropos_uci[pred.df.Eid$antigen=="GhV"] <- pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="GhV"] + 1.96*predict.gam(m1GhVnoSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="GhV",], type="response", se.fit = T, exclude = "s(year)")$se

#AngV
pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="AngV"] <- predict.gam(m1AngVnoSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="AngV",], type="response", exclude = "s(year)")
pred.df.Eid$predicted_seropos_lci[pred.df.Eid$antigen=="AngV"] <- pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="AngV"] - 1.96*predict.gam(m1AngVnoSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="AngV",], type="response", se.fit = T, exclude = "s(year)")$se
pred.df.Eid$predicted_seropos_uci[pred.df.Eid$antigen=="AngV"] <- pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="AngV"] + 1.96*predict.gam(m1AngVnoSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="AngV",], type="response", se.fit = T, exclude = "s(year)")$se


pred.df.Eid$predicted_seropos_lci[pred.df.Eid$predicted_seropos_lci<0] <- 0
pred.df.Eid$predicted_seropos_uci[pred.df.Eid$predicted_seropos_uci>1] <- 1
pred.df.Eid$bat_species="Eidolon dupreanum"



#rousettus seasonal plotting
pred.df.Rou <- cbind.data.frame(doy=rep(1:365,3*1), year=rep(c(2018,2019,2020), each=(365*1)), antigen = rep(rep(c("GhV"), each = 365), 3))

pred.df.Rou$bat_sex = "female"
pred.df.tmp <- pred.df.Rou
pred.df.tmp$bat_sex <- "male"
pred.df.Rou <- rbind(pred.df.Rou, pred.df.tmp)
pred.df.Rou$roost_site="Moramanga"
pred.df.Rou$mass_forearm_residual = 0
pred.df.Rou$combo_class <- "male_A"
pred.df.Rou$combo_class[pred.df.Rou$bat_sex=="female"] <- "female_NL"


#predict without effect of year for plotting
#females
pred.df.Rou$predicted_seropos <- NA
pred.df.Rou$predicted_seropos_lci <- NA
pred.df.Rou$predicted_seropos_uci <- NA

pred.df.Rou$predicted_seropos[pred.df.Rou$antigen=="GhV"] <- predict.gam(m2GhVnoSeX, newdata = pred.df.Rou[pred.df.Rou$antigen=="GhV",], type="response", exclude = "s(year)")
pred.df.Rou$predicted_seropos_lci[pred.df.Rou$antigen=="GhV"] <- pred.df.Rou$predicted_seropos[pred.df.Rou$antigen=="GhV"] - 1.96*predict.gam(m2GhVnoSeX, newdata = pred.df.Rou[pred.df.Rou$antigen=="GhV",], type="response", se.fit = T, exclude = "s(year)")$se
pred.df.Rou$predicted_seropos_uci[pred.df.Rou$antigen=="GhV"] <- pred.df.Rou$predicted_seropos[pred.df.Rou$antigen=="GhV"] + 1.96*predict.gam(m2GhVnoSeX, newdata = pred.df.Rou[pred.df.Rou$antigen=="GhV",], type="response", se.fit = T, exclude = "s(year)")$se

pred.df.Rou$bat_species="Rousettus madagascariensis"


pred.df.supp <- rbind(pred.df.Eid, pred.df.Rou)# add in the timing of the reproductive season in Moramanga for both species
pred.df.supp$sig = 0
pred.df.supp$sig[pred.df.supp$bat_species=="Rousettus madagascariensis" &  pred.df.supp$antigen=="GhV"] <- 2
pred.df.supp$sig[pred.df.supp$bat_species=="Eidolon dupreanum" & pred.df.supp$antigen=="GhV"] <- 1
pred.df.supp$sig[pred.df.supp$bat_species=="Eidolon dupreanum" & pred.df.supp$antigen=="AngV"] <- 2



# add in the timing of the reproductive season in Moramanga for both species
# from angelo and santino's paper
# eidolon, gestation is from: Aug-3 - Nov 15 (104 days); lactation Nov 16 - Feb 2 (78 days)
# rousettus, gestation is from: Sept 11 - Dec 11 (91 days); lactation Dec 12 - Feb 19 (69 days)
# dry season = May 15 - October 10



pred.df <- cbind.data.frame(mid_date = seq(from= as.Date("2019-01-01"), to = as.Date("2019-12-31"), by=1), doy=1:365)
pred.df$bat_species <- "Eidolon dupreanum"
pred.df$antigen <- "GhV"
pred.df2 <- pred.df
pred.df2$antigen <- "AngV"
pred.df <- rbind(pred.df, pred.df2)
pred.df$roost_site <- "Moramanga"
pred.df$year <- 2019

pred.df$prevalence <- NA
pred.df$prev_lci <- NA
pred.df$prev_uci <- NA



pred.df$prevalence[pred.df$antigen=="GhV"] <- predict(m1GhVnoSeX, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="GhV"] <- predict(m1GhVnoSeX, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)")) - 1.96*predict(m1GhVnoSeX, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="GhV"] <- predict(m1GhVnoSeX, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)")) + 1.96*predict(m1GhVnoSeX, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)"), se.fit = TRUE)$se

pred.df$prevalence[pred.df$antigen=="AngV"] <- predict(m1AngVnoSeX, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="AngV"] <- predict(m1AngVnoSeX, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(year)")) - 1.96*predict(m1AngVnoSeX, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="AngV"] <- predict(m1AngVnoSeX, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(year)")) + 1.96*predict(m1AngVnoSeX, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(year)"), se.fit = TRUE)$se

pred.df$prev_lci[pred.df$prev_lci<0] <- 0
pred.df$prev_uci[pred.df$prev_uci>1] <- 1


#and make the gestation database
gest.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis"), start = c(as.Date("2019-08-03"), as.Date("2019-09-11")), stop=c(as.Date("2019-11-15"),  as.Date("2019-12-11")))
gest.dat$season = "gestation"
lact.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis",  "Rousettus madagascariensis"), start = c(as.Date("2019-11-16"), as.Date("2019-01-01"), as.Date("2019-12-12"),  as.Date("2019-01-01")), stop=c(as.Date("2019-12-31"), as.Date("2019-02-02"), as.Date("2019-12-31"), as.Date("2015-02-19")))
lact.dat$season = "lactation"

#and winter
cal.dat <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis"), start = c( as.Date("2019-05-15"),as.Date("2019-05-15")), stop=c(as.Date("2019-10-10"),as.Date("2019-10-10")))
cal.dat$season <- "winter"
cal.dat <- rbind(cal.dat, gest.dat, lact.dat)
cal.dat$year <- year(cal.dat$stop)

cal.dat$doy_start <- yday(cal.dat$start)
cal.dat$doy_stop <- yday(cal.dat$stop)




#collect data on prevalence just for plotting
head(merge.dat)
merge.dat$seropos <- as.numeric(as.character(merge.dat$seropos))
merge.dat$day_of_year_session <- yday(merge.dat$mid_date)
merge.plot <- ddply(merge.dat,.(bat_species, year, day_of_year_session, antigen), summarise, N=length(seropos), seropos=sum(seropos))
merge.plot = subset(merge.plot, antigen=="GhV" | antigen=="AngV")
merge.plot$seroprevalence <-   merge.plot$seropos/merge.plot$N
merge.plot$seroprevalence_lci <- merge.plot$seroprevalence_uci <- NA

for (i in 1:length(merge.plot$bat_species)){
  merge.plot$seroprevalence_lci[i] <-  binom.test(merge.plot$seropos[i],merge.plot$N[i])$conf.int[1]
  merge.plot$seroprevalence_uci[i] <-  binom.test(merge.plot$seropos[i],merge.plot$N[i])$conf.int[2]
  
}


merge.plot$sig = 1
# merge.plot$sig[merge.plot$bat_species=="Rousettus madagascariensis" & merge.plot$antigen=="GhV"] <- 1
# merge.plot$sig[merge.plot$bat_species=="Eidolon dupreanum" & merge.plot$antigen=="GhV"] <- 1

alphaz = c('0' = .3, '1'=.8)
alphaz.new <- c('0'=.4, '1'=.8)
pred.df.supp$sig <- 1
pred.df.supp$sig <- as.factor(pred.df.supp$sig)
merge.plot$sig <- as.factor(merge.plot$sig)

pred.df.supp$predicted_seropos_lci[pred.df.supp$predicted_seropos_lci<0] <- 0


colorz = c('AngV' = "#F8766D", 'HeV' ="#A58AFF", 
           'NiV' = "#DB72FB", 'GhV' = "#00C094", 
           'CedV'= "#53B400", "MojV" = "#619CFF")


colorz = c('AngV' = "#F8766D", 'GhV' = "#00C094")

filz <- c(colorz, 'gestation' = "yellow", 'lactation' = "lightpink", 'winter' = "skyblue")

filz1 <- filz[1:(length(filz)-3)]
filz2 <- filz[(length(filz)-2):length(filz)]

yearz <- c('2013'=3, '2014' = 7, '2015'=8, '2016'=18, '2018' = 16, '2019' = 17, '2020' = 15, '2021' = 4, '2023' = 6)


#slim to just 2019 for plotting purposes
pIntra <- ggplot(subset(pred.df.supp, sig==1 )) + 
  geom_point(data=merge.plot, aes(x=day_of_year_session, y=seroprevalence, size=N, color=antigen, shape=year, alpha = sig)) +
  geom_linerange(data=merge.plot, aes(x=day_of_year_session, ymin=seroprevalence_lci, ymax=seroprevalence_uci, color=antigen, alpha = sig, group=year), linetype=2, show.legend = F) +
  geom_rect(data=cal.dat, aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3) +
  scale_fill_manual(values=filz2) +
  scale_shape_manual(values=yearz) +
  ggnewscale::new_scale_fill() +
  scale_alpha_manual(values=alphaz.new) +
  geom_ribbon(data = subset(pred.df.supp, sig==2), aes(x=doy, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.15, show.legend = F) +
  geom_line(data = subset(pred.df.supp, sig==2),aes(x=doy, y=predicted_seropos, color=antigen), alpha=.3) +
  geom_ribbon(aes(x=doy, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.3, show.legend = F) +
  geom_line(aes(x=doy, y=predicted_seropos, color=antigen)) + 
  facet_grid(~bat_species) + theme_bw() +
  scale_color_manual(values=colorz) + scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  scale_fill_manual(values=filz1) +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  theme(axis.text = element_text(size=14), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16),
        strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16),panel.spacing = unit(1.5, "lines"),
        strip.text.x = element_text(size=16, face = "italic"), panel.grid = element_blank()) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + xlab("month of year")  + guides(alpha="none");pIntra


pIntraEid <- ggplot(subset(pred.df.supp, sig==1 & bat_species=="Eidolon dupreanum")) + 
  geom_point(data=subset(merge.plot, bat_species=="Eidolon dupreanum"), aes(x=day_of_year_session, y=seroprevalence, size=N, color=antigen, shape=year, alpha = sig)) +
  geom_linerange(data=subset(merge.plot, bat_species=="Eidolon dupreanum"), aes(x=day_of_year_session, ymin=seroprevalence_lci, ymax=seroprevalence_uci, color=antigen, alpha = sig, group=year), linetype=2, show.legend = F) +
  geom_rect(data=subset(cal.dat, bat_species=="Eidolon dupreanum"), aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3) +
  scale_fill_manual(values=filz2) +
  scale_size_continuous(name = "N bats", breaks = c(10,20,30)) +
  scale_shape_manual(values=yearz) +
  ggnewscale::new_scale_fill() +
  scale_alpha_manual(values=alphaz.new) +
  geom_ribbon(data = subset(pred.df.supp, sig==2 & bat_species=="Eidolon dupreanum" ), aes(x=doy, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.15, show.legend = F) +
  geom_line(data = subset(pred.df.supp, sig==2 & bat_species=="Eidolon dupreanum"),aes(x=doy, y=predicted_seropos, color=antigen), alpha=.3) +
  geom_ribbon(aes(x=doy, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.3, show.legend = F) +
  geom_line(aes(x=doy, y=predicted_seropos, color=antigen)) + 
  facet_grid(~bat_species) + theme_bw() +
  scale_color_manual(values=colorz) + 
  scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  scale_fill_manual(values=filz1) +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  theme(axis.text = element_text(size=14), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16), 
        strip.background = element_blank(),# element_rect(fill="white"), 
        #strip.text.y = element_text(size=16),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_blank(),
        #strip.text.x = element_text(size=16, face = "italic"), 
        panel.grid = element_blank()) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + xlab("month of year")  + guides(alpha="none");pIntraEid

#and no legend
pIntraEidnoLeg <- ggplot(subset(pred.df.supp, sig==1 & bat_species=="Eidolon dupreanum")) + 
  geom_point(data=subset(merge.plot, bat_species=="Eidolon dupreanum"), aes(x=day_of_year_session, y=seroprevalence, size=N, color=antigen, shape=year, alpha = sig), show.legend = F) +
  geom_linerange(data=subset(merge.plot, bat_species=="Eidolon dupreanum"), aes(x=day_of_year_session, ymin=seroprevalence_lci, ymax=seroprevalence_uci, color=antigen, alpha = sig, group=year), linetype=2, show.legend = F) +
  geom_rect(data=subset(cal.dat, bat_species=="Eidolon dupreanum"), aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3,show.legend = F) +
  scale_fill_manual(values=filz2) +
  scale_size_continuous(name = "N bats", breaks = c(10,20,30)) +
  scale_shape_manual(values=yearz) +
  ggnewscale::new_scale_fill() +
  scale_alpha_manual(values=alphaz.new) +
  geom_ribbon(data = subset(pred.df.supp, sig==2 & bat_species=="Eidolon dupreanum" ), aes(x=doy, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.15, show.legend = F) +
  geom_line(data = subset(pred.df.supp, sig==2 & bat_species=="Eidolon dupreanum"),aes(x=doy, y=predicted_seropos, color=antigen), alpha=.3, show.legend = F) +
  geom_ribbon(aes(x=doy, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.3, show.legend = F) +
  geom_line(aes(x=doy, y=predicted_seropos, color=antigen), show.legend = F) + 
  facet_grid(~bat_species) + theme_bw() +
  scale_color_manual(values=colorz) + 
  scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  scale_fill_manual(values=filz1) +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  theme(axis.text = element_text(size=14), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16), 
        strip.background = element_blank(),# element_rect(fill="white"), 
        #strip.text.y = element_text(size=16),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_blank(),
        #strip.text.x = element_text(size=16, face = "italic"), 
        panel.grid = element_blank()) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + xlab("month of year")  + guides(alpha="none");pIntraEidnoLeg



pIntraRou <- ggplot(subset(pred.df.supp, sig==1 & bat_species=="Rousettus madagascariensis")) + 
  geom_point(data=subset(merge.plot, bat_species=="Rousettus madagascariensis"), aes(x=day_of_year_session, y=seroprevalence, size=N, color=antigen, shape=year, alpha = sig)) +
  geom_linerange(data=subset(merge.plot, bat_species=="Rousettus madagascariensis"), aes(x=day_of_year_session, ymin=seroprevalence_lci, ymax=seroprevalence_uci, color=antigen, alpha = sig, group=year), linetype=2, show.legend = F) +
  geom_rect(data=subset(cal.dat, bat_species=="Rousettus madagascariensis"), aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3) +
  scale_fill_manual(values=filz2) +
  scale_size_continuous(name = "N bats", breaks = c(15,10,15)) +
  scale_shape_manual(values=yearz) +
  ggnewscale::new_scale_fill() +
  scale_alpha_manual(values=alphaz.new) +
  geom_ribbon(data = subset(pred.df.supp, sig==2 & bat_species=="Rousettus madagascariensis" ), aes(x=doy, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.15, show.legend = F) +
  geom_line(data = subset(pred.df.supp, sig==2 & bat_species=="Rousettus madagascariensis"),aes(x=doy, y=predicted_seropos, color=antigen), alpha=.3) +
  geom_ribbon(aes(x=doy, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.3, show.legend = F) +
  geom_line(aes(x=doy, y=predicted_seropos, color=antigen)) + 
  facet_grid(~bat_species) + theme_bw() +
  scale_color_manual(values=colorz) + 
  scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  scale_fill_manual(values=filz1) +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  theme(axis.text = element_text(size=14), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16), 
        strip.background = element_blank(),# element_rect(fill="white"), 
        #strip.text.y = element_text(size=16),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_blank(),
        #strip.text.x = element_text(size=16, face = "italic"), 
        panel.grid = element_blank()) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + xlab("month of year")  + guides(alpha="none");pIntraRou

#and no legend
pIntraRounoLeg <- ggplot(subset(pred.df.supp, sig==1 & bat_species=="Rousettus madagascariensis")) + 
  geom_point(data=subset(merge.plot, bat_species=="Rousettus madagascariensis"), aes(x=day_of_year_session, y=seroprevalence, size=N, color=antigen, shape=year, alpha = sig), show.legend = F) +
  geom_linerange(data=subset(merge.plot, bat_species=="Rousettus madagascariensis"), aes(x=day_of_year_session, ymin=seroprevalence_lci, ymax=seroprevalence_uci, color=antigen, alpha = sig, group=year), linetype=2, show.legend = F) +
  geom_rect(data=subset(cal.dat, bat_species=="Rousettus madagascariensis"), aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3,show.legend = F) +
  scale_fill_manual(values=filz2) +
  scale_size_continuous(name = "N bats", breaks = c(15,10,15)) +
  scale_shape_manual(values=yearz) +
  ggnewscale::new_scale_fill() +
  scale_alpha_manual(values=alphaz.new) +
  geom_ribbon(data = subset(pred.df.supp, sig==2 & bat_species=="Rousettus madagascariensis" ), aes(x=doy, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.15, show.legend = F) +
  geom_line(data = subset(pred.df.supp, sig==2 & bat_species=="Rousettus madagascariensis"),aes(x=doy, y=predicted_seropos, color=antigen), alpha=.3, show.legend = F) +
  geom_ribbon(aes(x=doy, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.3, show.legend = F) +
  geom_line(aes(x=doy, y=predicted_seropos, color=antigen), show.legend = F) + 
  facet_grid(~bat_species) + theme_bw() +
  scale_color_manual(values=colorz) + 
  scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  scale_fill_manual(values=filz1) +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  theme(axis.text = element_text(size=14), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16), 
        strip.background = element_blank(),# element_rect(fill="white"), 
        #strip.text.y = element_text(size=16),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_blank(),
        #strip.text.x = element_text(size=16, face = "italic"), 
        panel.grid = element_blank()) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + xlab("month of year")  + guides(alpha="none");pIntraRounoLeg


#and pair with intraannual shedding data


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


sum(shed.dat$AngV_pos)#45
sum(shed.dat$GhV_pos)#127

# link this to the catching data to plot by sampling session

#save sero dat
merge.sero.dat <- merge.dat

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


#tot.merge = subset(merge.dat, year<2022 & year>2016)
tot.merge= merge.dat

head(tot.merge)
tot.merge$processing_date <- as.Date(tot.merge$processing_date)
tot.merge$month <- month(tot.merge$processing_date)
tot.merge$year <- year(tot.merge$processing_date)
tot.merge = subset(tot.merge, bat_species=="Eidolon dupreanum" | bat_species =="Rousettus madagascariensis" )
# #tot.merge = subset(tot.merge, bat_species=="Eidolon dupreanum" | bat_species =="Rousettus madagascariensis" | bat_species =="Pteropus rufus")
# 

tot.merge$year <- as.factor(tot.merge$year)
tot.merge$roost_site <- as.factor(tot.merge$roost_site)
tot.merge$bat_sex <- as.factor(tot.merge$bat_sex)

gamGhVnoSex <- gam(GhV_pos ~ 
                     s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                     s(year, bs="re"),
                   #s(as.numeric(mid_date), bs="tp", k=6, by=interaction(roost_site, era)),
                   family = "binomial",
                   method = "REML", 
                   gamma = 1,
                   select = TRUE,
                   data= subset(tot.merge, bat_species =="Eidolon dupreanum" ))
summary(gamGhVnoSex) #year and ankarana
#plot(gamGhVnoSex)

gamGhVSex <- gam(GhV_pos ~ 
                   s(as.numeric(doy), bs="cc", k=10, by = interaction(roost_site, bat_sex)) +
                   s(year, bs="re"),
                 #s(as.numeric(mid_date), bs="tp", k=6, by = interaction(roost_site, bat_sex, era)),
                 family = "binomial",
                 method = "REML", 
                 gamma = 1,
                 select = TRUE,
                 data= subset(tot.merge, bat_species =="Eidolon dupreanum" ))
summary(gamGhVSex) # ankarana female, and year

AIC(gamGhVnoSex, gamGhVSex) #better with no sex

# #and AngV


gamAngVnoSex <- gam(AngV_pos ~ 
                      s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                      s(year, bs="re"),# +
                    #s(as.numeric(mid_date), bs="tp", k=6, by=interaction(roost_site, era)),
                    family = "binomial",
                    method = "REML", 
                    gamma = 1,
                    select = TRUE,
                    data= subset(tot.merge, bat_species =="Eidolon dupreanum" ))
summary(gamAngVnoSex) # year
#plot(gamAngVnoSex)

gamAngVSex <- gam(AngV_pos ~  
                    s(as.numeric(doy), bs="cc", k=10, by = interaction(roost_site, bat_sex)) +
                    #s(as.numeric(mid_date), bs="tp", k=6, by = interaction(roost_site, bat_sex, era)),
                    s(year, bs="re"),# +
                  family = "binomial",
                  method = "REML", 
                  gamma = 1,
                  select = TRUE,
                  data= subset(tot.merge, bat_species =="Eidolon dupreanum" ))
summary(gamAngVSex) #Anakarana female and year
AIC(gamAngVnoSex, gamAngVSex) #better with no sex

pred.df <- cbind.data.frame(mid_date = seq(from= as.Date("2019-01-01"), to = as.Date("2019-12-31"), by=1), doy=1:365)
pred.df$bat_species <- "Eidolon dupreanum"
pred.df$antigen <- "GhV"
pred.df2 <- pred.df
pred.df2$antigen <- "AngV"
pred.df <- rbind(pred.df, pred.df2)
pred.df$roost_site <- "Moramanga"
pred.df$year <- 2019

pred.df$prevalence <- NA
pred.df$prev_lci <- NA
pred.df$prev_uci <- NA



pred.df$prevalence[pred.df$antigen=="GhV"] <- predict(gamGhVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="GhV"] <- predict(gamGhVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)")) - 1.96*predict(gamGhVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="GhV"] <- predict(gamGhVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)")) + 1.96*predict(gamGhVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c( "s(year)"), se.fit = TRUE)$se

pred.df$prevalence[pred.df$antigen=="AngV"] <- predict(gamAngVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="AngV"] <- predict(gamAngVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c( "s(year)")) - 1.96*predict(gamAngVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c( "s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="AngV"] <- predict(gamAngVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c( "s(year)")) + 1.96*predict(gamAngVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(year)"), se.fit = TRUE)$se



pred.df$prev_lci[pred.df$prev_lci<0] <- 0
pred.df$prev_uci[pred.df$prev_uci>1] <- 1

#already made cal.dat above
#need to plot prevalence
prev.dat <- ddply(subset(tot.merge, !is.na(GhV_pos)),.(bat_species, roost_site, mid_date), summarise, GhVpos = sum(GhV_pos, na.rm = T), AngVpos=sum(AngV_pos, na.rm=T), Ntot = length(sampleid))
prev.dat.sex <- ddply(subset(tot.merge, !is.na(GhV_pos)),.(bat_species, bat_sex, roost_site, mid_date), summarise, GhVpos = sum(GhV_pos, na.rm = T), AngVpos=sum(AngV_pos, na.rm=T), Ntot = length(sampleid))
prev.dat$prevGhV <- prev.dat$GhVpos/prev.dat$Ntot
prev.dat$prevAngV <- prev.dat$AngVpos/prev.dat$Ntot
# 
prev.dat.sex$prevGhV <- prev.dat.sex$GhVpos/prev.dat.sex$Ntot
prev.dat.sex$prevAngV <- prev.dat.sex$AngVpos/prev.dat.sex$Ntot

# 
# #include era block to correct for the gaps in the time series
 prev.dat$era <- "II"
 prev.dat$era[prev.dat$mid_date<"2017-01-01"] <- "I"
 prev.dat$era[prev.dat$mid_date>"2021-01-01"] <- "III"
 prev.dat$era <- as.factor(prev.dat$era)
# 
# #and by month
 prev.dat$month <- month(prev.dat$mid_date)
 prev.dat$year <- year(prev.dat$mid_date)
 prev.dat <- prev.dat[complete.cases(prev.dat),]
 prev.dat$doy <- yday(prev.dat$mid_date)
 prev.dat$year <- as.factor(prev.dat$year)
# 
# #and by month
 prev.dat.sex$month <- month(prev.dat.sex$mid_date)
 prev.dat.sex$year <- year(prev.dat.sex$mid_date)
 prev.dat.sex <- prev.dat.sex[complete.cases(prev.dat.sex),]
 prev.dat.sex$doy <- yday(prev.dat.sex$mid_date)
 prev.dat.sex$year <- as.factor(prev.dat.sex$year)
# 
 prev.dat.sex$era <- "II"
 prev.dat.sex$era[prev.dat.sex$mid_date<"2017-01-01"] <- "I"
 prev.dat.sex$era[prev.dat.sex$mid_date>"2021-01-01"] <- "III"
 prev.dat.sex$era <- as.factor(prev.dat.sex$era)
# 
# 
# p2 <- ggplot(prev.dat) + facet_grid(~bat_species) +
#   ylab("prevalence GhV") + theme_bw() + #scale_shape_manual(values = c(21,24)) +
#   geom_point(aes(x=month, y=prevGhV, size=Ntot, shape=roost_site, color=as.character(year)));p2
# 
# # and gam to identify seasonal signatures
prev.dat$GhVneg <- prev.dat$Ntot-prev.dat$GhVpos
prev.dat.sex$GhVneg <- prev.dat.sex$Ntot-prev.dat.sex$GhVpos
prev.dat$roost_site <- as.factor(prev.dat$roost_site)
prev.dat$AngVneg <- prev.dat$Ntot-prev.dat$AngVpos
prev.dat.sex$AngVneg <- prev.dat.sex$Ntot-prev.dat.sex$AngVpos
prev.dat.sex$roost_site <- as.factor(prev.dat.sex$roost_site)

prev.dat$prevGhV_lci <- prev.dat$prevGhV_uci <- prev.dat$prevAngV_lci <- prev.dat$prevAngV_uci <- NA
for (i in 1:length(prev.dat$prevGhV)){
  prev.dat$prevGhV_lci[i] <- binom.test(prev.dat$GhVpos[i], prev.dat$Ntot[i])$conf.int[1]
  prev.dat$prevGhV_uci[i] <- binom.test(prev.dat$GhVpos[i], prev.dat$Ntot[i])$conf.int[2]
  prev.dat$prevAngV_lci[i] <- binom.test(prev.dat$AngVpos[i], prev.dat$Ntot[i])$conf.int[1]
  prev.dat$prevAngV_uci[i] <- binom.test(prev.dat$AngVpos[i], prev.dat$Ntot[i])$conf.int[2]
  
}


###rousettus

gamGhVnoSexRou <- gam(GhV_pos ~ 
                        s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                        #s(as.numeric(mid_date), bs="tp", k=6, by=roost_site),
                        s(year, bs="re"),
                      family = "binomial",
                      method = "REML", 
                      gamma = 1,
                      select = TRUE,
                      data= subset(tot.merge, bat_species =="Rousettus madagascariensis" ))
summary(gamGhVnoSexRou) #just year
#plot(gamGhVnoSex)

gamGhVSexRou <- gam(GhV_pos ~ 
                      s(as.numeric(doy), bs="cc", k=10, by = interaction(roost_site, bat_sex)) +
                      #s(as.numeric(mid_date), bs="tp", k=6, by = interaction(roost_site, bat_sex)),
                      s(year, bs="re"),
                    family = "binomial",
                    method = "REML", 
                    gamma = 1,
                    select = TRUE,
                    data= subset(tot.merge, bat_species =="Rousettus madagascariensis"))
summary(gamGhVSexRou) #just year

AIC(gamGhVnoSexRou, gamGhVSexRou) #better with sex but not by enough to matter


#and AngV

gamAngVnoSexRou <- gam(AngV_pos ~ 
                         s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                         #s(as.numeric(mid_date), bs="tp", k=6, by=roost_site),
                         s(year, bs="re"),
                       family = "binomial",
                       method = "REML", 
                       gamma = 1,
                       select = TRUE,
                       data= subset(tot.merge, bat_species =="Rousettus madagascariensis"))
summary(gamAngVnoSexRou) #year
#plot(gamAngVnoSex)

gamAngVSexRou <- gam(AngV_pos ~ 
                       s(as.numeric(doy), bs="cc", k=10, by = interaction(roost_site, bat_sex)) +
                       #s(as.numeric(mid_date), bs="tp", k=6, by = interaction(roost_site, bat_sex)),
                       s(year, bs="re"),
                     family = "binomial",
                     method = "REML", 
                     gamma = 1,
                     select = TRUE,
                     data= subset(tot.merge, bat_species =="Rousettus madagascariensis"))
summary(gamAngVSexRou) #year, moramanga male, ankarana female

AIC(gamAngVnoSexRou, gamAngVSexRou) #better with no sex

pred.df.Rou <- cbind.data.frame(mid_date = seq(from= as.Date("2019-01-01"), to = as.Date("2019-12-31"), by=1), doy=1:365)
pred.df.Rou$bat_species <- "Rousettus madagascariensis"
pred.df.Rou$antigen <- "GhV"
pred.df.Rou2 <- pred.df.Rou
pred.df.Rou2$antigen <- "AngV"
pred.df.Rou <- rbind(pred.df.Rou, pred.df.Rou2)
pred.df.Rou$roost_site <- "Moramanga"
pred.df.Rou$year <- 2019
# pred.df.Rou$era <- "I"
# pred.df.Rou$era[pred.df.Rou$mid_date < "2017-07-01"] <- "II"
# pred.df.Rou$era <- as.factor(pred.df.Rou$era)

pred.df.Rou$prevalence <- NA
pred.df.Rou$prev_lci <- NA
pred.df.Rou$prev_uci <- NA

pred.df.Rou$prevalence[pred.df.Rou$antigen=="GhV"] <- predict(gamGhVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="GhV",], exclude = c("s(year)"))
pred.df.Rou$prev_lci[pred.df.Rou$antigen=="GhV"] <- predict(gamGhVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="GhV",], exclude = c("s(year)")) - 1.96*predict(gamGhVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="GhV",], exclude = c( "s(year)"), se.fit = TRUE)$se
pred.df.Rou$prev_uci[pred.df.Rou$antigen=="GhV"] <- predict(gamGhVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="GhV",], exclude = c( "s(year)")) + 1.96*predict(gamGhVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="GhV",], exclude = c("s(year)"), se.fit = TRUE)$se

pred.df.Rou$prevalence[pred.df.Rou$antigen=="AngV"] <- predict(gamAngVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="AngV",], exclude = c( "s(year)"))
pred.df.Rou$prev_lci[pred.df.Rou$antigen=="AngV"] <- predict(gamAngVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="AngV",], exclude = c("s(year)")) - 1.96*predict(gamAngVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="AngV",], exclude = c( "s(year)"), se.fit = TRUE)$se
pred.df.Rou$prev_uci[pred.df.Rou$antigen=="AngV"] <- predict(gamAngVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="AngV",],exclude = c( "s(year)")) + 1.96*predict(gamAngVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="AngV",], exclude = c("s(year)"), se.fit = TRUE)$se


colorz = c('AngV' = "#F8766D", 'HeV' ="#A58AFF", 
           'NiV' = "#DB72FB", 'GhV' = "#00C094", 
           'CedV'= "#53B400", "MojV" = "#619CFF")


filz <- c('gestation' = "yellow", 'lactation' = "lightpink", 'winter' = "skyblue")

yearz <- c('2013'='3', '2014' = '7', '2015'='8', '2016'='18', '2018' = '16', '2019' = '17', '2020' = '15')
yearz <- c('2013'=3, '2014' = 7, '2015'=8, '2016'=18, '2018' = 16, '2019' = 17, '2020' = 15, '2021' = 4, '2023' = 6)

pIntrannualShed <- ggplot(data=subset(prev.dat)) +
  geom_rect(data=cal.dat, aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=doy, y=prevGhV, size=Ntot, shape=year), color="#00C094") + 
  geom_point(aes(x=doy, y=prevAngV, size=Ntot, shape=year), color="#F8766D") + 
  scale_fill_manual(values=filz) +
  ggnewscale::new_scale_fill() + scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  geom_linerange(aes(x=doy, ymin=prevGhV_lci, ymax=prevGhV_uci),color="#00C094", linetype=2, alpha=.3, show.legend = F) +
  geom_linerange(aes(x=doy, ymin=prevAngV_lci, ymax=prevAngV_uci),color="#F8766D", linetype=2, alpha=.3, show.legend = F) +
  geom_ribbon(data=pred.df, aes(x=doy, ymin= prev_lci, ymax=prev_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=pred.df, aes(x=doy, y= prevalence, color=antigen)) +
  geom_ribbon(data=pred.df.Rou, aes(x=doy, ymin= prev_lci, ymax=prev_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=pred.df.Rou, aes(x=doy, y= prevalence, color=antigen)) +
  facet_grid(~bat_species) + theme_bw() +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),, 
        strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  scale_shape_manual(values=yearz) +
  coord_cartesian(xlim = c(0, 365), ylim = c(0, 1),expand = FALSE) + ylab("prevalence") + xlab("month of year") +
  scale_fill_manual(values=colorz) + scale_color_manual(values=colorz) + 
  guides(size=guide_legend(override.aes = list(color="black")), shape=guide_legend(override.aes = list(color="black")))
pIntrannualShed


pIntrannualShedEid <- ggplot(data=subset(prev.dat, bat_species=="Eidolon dupreanum")) +
  geom_rect(data=subset(cal.dat, bat_species=="Eidolon dupreanum"), aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=doy, y=prevGhV, size=Ntot, shape=year), color="#00C094") + 
  geom_point(aes(x=doy, y=prevAngV, size=Ntot, shape=year), color="#F8766D") + 
  scale_fill_manual(values=filz) +
  scale_size_continuous(name = "N bats", breaks = c(10,20,30)) +
  ggnewscale::new_scale_fill() + 
  scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  geom_linerange(aes(x=doy, ymin=prevGhV_lci, ymax=prevGhV_uci),color="#00C094", linetype=2, alpha=.3, show.legend = F) +
  geom_linerange(aes(x=doy, ymin=prevAngV_lci, ymax=prevAngV_uci),color="#F8766D", linetype=2, alpha=.3, show.legend = F) +
  geom_ribbon(data=pred.df, aes(x=doy, ymin= prev_lci, ymax=prev_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=pred.df, aes(x=doy, y= prevalence, color=antigen)) +
  #geom_ribbon(data=pred.df.Rou, aes(x=doy, ymin= prev_lci, ymax=prev_uci, fill=antigen), alpha=.3, show.legend = F) +
  #geom_line(data=pred.df.Rou, aes(x=doy, y= prevalence, color=antigen)) +
  facet_grid(~bat_species) + theme_bw() +
  theme(axis.text = element_text(size=14), axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        strip.background = element_blank(),# element_rect(fill="white"), 
        strip.text = element_blank(),# .y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  scale_shape_manual(values=yearz) +
  coord_cartesian(xlim = c(0, 365), ylim = c(0, 1),expand = FALSE) + ylab("prevalence") + xlab("month of year") +
  scale_fill_manual(values=colorz) + scale_color_manual(values=colorz) + 
  guides(size=guide_legend(override.aes = list(color="black")), shape=guide_legend(override.aes = list(color="black")))
pIntrannualShedEid 

legEid <- get_legend(pIntrannualShedEid)

pIntrannualShedEidnoLeg <- ggplot(data=subset(prev.dat, bat_species=="Eidolon dupreanum")) +
  geom_rect(data=subset(cal.dat, bat_species=="Eidolon dupreanum"), aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3, show.legend = F) +
  geom_point(aes(x=doy, y=prevGhV, size=Ntot, shape=year), color="#00C094", show.legend = F) + 
  geom_point(aes(x=doy, y=prevAngV, size=Ntot, shape=year), color="#F8766D", show.legend = F) + 
  scale_fill_manual(values=filz) +
  scale_size_continuous(name = "N bats", breaks = c(10,20,30)) +
  ggnewscale::new_scale_fill() + 
  scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  geom_linerange(aes(x=doy, ymin=prevGhV_lci, ymax=prevGhV_uci),color="#00C094", linetype=2, alpha=.3, show.legend = F) +
  geom_linerange(aes(x=doy, ymin=prevAngV_lci, ymax=prevAngV_uci),color="#F8766D", linetype=2, alpha=.3, show.legend = F) +
  geom_ribbon(data=pred.df, aes(x=doy, ymin= prev_lci, ymax=prev_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=pred.df, aes(x=doy, y= prevalence, color=antigen), show.legend = F) +
  #geom_ribbon(data=pred.df.Rou, aes(x=doy, ymin= prev_lci, ymax=prev_uci, fill=antigen), alpha=.3, show.legend = F) +
  #geom_line(data=pred.df.Rou, aes(x=doy, y= prevalence, color=antigen)) +
  facet_grid(~bat_species) + theme_bw() +
  theme(axis.text = element_text(size=14), axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        strip.background = element_blank(),# element_rect(fill="white"), 
        strip.text = element_blank(),# .y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  scale_shape_manual(values=yearz) +
  coord_cartesian(xlim = c(0, 365), ylim = c(0, 1),expand = FALSE) + ylab("prevalence") + xlab("month of year") +
  scale_fill_manual(values=colorz) + scale_color_manual(values=colorz) + 
  guides(size=guide_legend(override.aes = list(color="black")), shape=guide_legend(override.aes = list(color="black")))


pIntrannualShedRou <- ggplot(data=subset(prev.dat, bat_species=="Rousettus madagascariensis")) +
  geom_rect(data=subset(cal.dat, bat_species=="Rousettus madagascariensis"), aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=doy, y=prevGhV, size=Ntot, shape=year), color="#00C094") + 
  geom_point(aes(x=doy, y=prevAngV, size=Ntot, shape=year), color="#F8766D") + 
  scale_fill_manual(values=filz) +
  scale_size_continuous(name = "N bats", breaks = c(15,10,15)) +
  ggnewscale::new_scale_fill() + scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  geom_linerange(aes(x=doy, ymin=prevGhV_lci, ymax=prevGhV_uci),color="#00C094", linetype=2, alpha=.3, show.legend = F) +
  geom_linerange(aes(x=doy, ymin=prevAngV_lci, ymax=prevAngV_uci),color="#F8766D", linetype=2, alpha=.3, show.legend = F) +
  #geom_ribbon(data=pred.df, aes(x=doy, ymin= prev_lci, ymax=prev_uci, fill=antigen), alpha=.3, show.legend = F) +
  #geom_line(data=pred.df, aes(x=doy, y= prevalence, color=antigen)) +
  geom_ribbon(data=pred.df.Rou, aes(x=doy, ymin= prev_lci, ymax=prev_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=pred.df.Rou, aes(x=doy, y= prevalence, color=antigen)) +
  facet_grid(~bat_species) + theme_bw() +
  theme(axis.text = element_text(size=14), axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        strip.background = element_blank(),# element_rect(fill="white"), 
        strip.text = element_blank(),# .y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  scale_shape_manual(values=yearz) +
  coord_cartesian(xlim = c(0, 365), ylim = c(0, 1),expand = FALSE) + ylab("prevalence") + xlab("month of year") +
  scale_fill_manual(values=colorz) + scale_color_manual(values=colorz) + 
  guides(size=guide_legend(override.aes = list(color="black")), shape=guide_legend(override.aes = list(color="black")))
pIntrannualShedRou

legRou <- get_legend(pIntrannualShedRou)

pIntrannualShedRounoLeg <- ggplot(data=subset(prev.dat, bat_species=="Rousettus madagascariensis")) +
  geom_rect(data=subset(cal.dat, bat_species=="Rousettus madagascariensis"), aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3, show.legend = F) +
  geom_point(aes(x=doy, y=prevGhV, size=Ntot, shape=year), color="#00C094", show.legend = F) + 
  geom_point(aes(x=doy, y=prevAngV, size=Ntot, shape=year), color="#F8766D", show.legend = F) + 
  scale_fill_manual(values=filz) +
  scale_size_continuous(name = "N bats", breaks = c(15,10,15)) +
  ggnewscale::new_scale_fill() + scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  geom_linerange(aes(x=doy, ymin=prevGhV_lci, ymax=prevGhV_uci),color="#00C094", linetype=2, alpha=.3, show.legend = F) +
  geom_linerange(aes(x=doy, ymin=prevAngV_lci, ymax=prevAngV_uci),color="#F8766D", linetype=2, alpha=.3, show.legend = F) +
  #geom_ribbon(data=pred.df, aes(x=doy, ymin= prev_lci, ymax=prev_uci, fill=antigen), alpha=.3, show.legend = F) +
  #geom_line(data=pred.df, aes(x=doy, y= prevalence, color=antigen)) +
  geom_ribbon(data=pred.df.Rou, aes(x=doy, ymin= prev_lci, ymax=prev_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=pred.df.Rou, aes(x=doy, y= prevalence, color=antigen), show.legend = F) +
  facet_grid(~bat_species) + theme_bw() +
  theme(axis.text = element_text(size=14), axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        strip.background = element_blank(),# element_rect(fill="white"), 
        strip.text = element_blank(),# .y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  scale_shape_manual(values=yearz) +
  coord_cartesian(xlim = c(0, 365), ylim = c(0, 1),expand = FALSE) + ylab("prevalence") + xlab("month of year") +
  scale_fill_manual(values=colorz) + scale_color_manual(values=colorz) + 
  guides(size=guide_legend(override.aes = list(color="black")), shape=guide_legend(override.aes = list(color="black")))
pIntrannualShedRounoLeg


#and put the two together
seroshed <- cowplot::plot_grid( pIntrannualShed, pIntra, ncol=1, nrow = 2, align = "hv")
seroshedEid <- cowplot::plot_grid(pIntrannualShedEidnoLeg, pIntraEidnoLeg, ncol=1, nrow = 2, align = "hv")
seroshedRou <- cowplot::plot_grid(pIntrannualShedRounoLeg, pIntraRounoLeg, ncol=1, nrow = 2, align = "hv")

seroshedEid1 <- cowplot::plot_grid(seroshedEid, legEid, ncol=2, nrow=1, rel_widths = c(1,.15))
seroshedRou1 <- cowplot::plot_grid(seroshedRou, legRou, ncol=2, nrow=1, rel_widths = c(1,.15))


ggsave(file = paste0(homewd,"/final_figures_tables/tmp/Intra-annual-sero-shed-Eid.png"),
       plot=seroshedEid1,
       units="mm",  
       width=60, 
       height=80, 
       scale=3, 
       dpi=300)


ggsave(file = paste0(homewd,"/final_figures_tables/tmp/Intra-annual-sero-shed-Rou.png"),
       plot=seroshedRou1,
       units="mm",  
       width=60, 
       height=80, 
       scale=3, 
       dpi=300)

#and combine with Ct

#keep Eid for main paper and Rousettus for the supplement

tot.merge$combo_class <- paste0(tot.merge$bat_sex, "-", tot.merge$bat_age_class)
tot.merge$combo_class[tot.merge$bat_age_class=="J"] <- "juvenile"
tot.merge = subset(tot.merge, combo_class!="NA-NA")


merge.sero.dat$combo_class <- paste0(merge.sero.dat$bat_sex, "-", merge.sero.dat$bat_age_class)
merge.sero.dat$combo_class[merge.sero.dat$bat_age_class=="J"] <- "juvenile"
merge.sero.dat = subset(merge.sero.dat, combo_class!="NA-NA")

tot.merge$carrying_pup <- 0
tot.merge$carrying_pup[tot.merge$pupmom_pair=="yes" & tot.merge$bat_age_class!="J"] <- 1
tot.merge$mass_forearm_residual[is.na(tot.merge$mass_forearm_residual)] <- 0
tot.merge$combo_class <- factor(tot.merge$combo_class, levels = c("female-NL", "female-P", "female-L", "male-A", "juvenile"))
merge.sero.dat$combo_class <- factor(merge.sero.dat$combo_class, levels = c("female-NL", "female-P", "female-L", "male-A", "juvenile"))
tot.merge$bat_sex
tot.merge$bat_age_class <- as.factor(tot.merge$bat_age_class)
tot.merge$roost_site <- factor(tot.merge$roost_site, levels=c(  "Mahabo", "Ankarana", "Moramanga"))
merge.sero.dat$roost_site <- factor(merge.sero.dat$roost_site, levels=c(  "Ankarana", "Moramanga"))
tot.merge$mass_cat <- NA
tot.merge$mass_cat[tot.merge$mass_forearm_residual<0] <- "poor"
tot.merge$mass_cat[tot.merge$mass_forearm_residual>0] <- "good"
tot.merge$mass_cat[tot.merge$mass_forearm_residua==0] <- "neutral"
tot.merge$mass_cat  <- factor(tot.merge$mass_cat, levels = c("neutral", "poor", "good"))


merge.sero.dat$mass_cat <- NA
merge.sero.dat$mass_cat[merge.sero.dat$mass_forearm_residual<0] <- "poor"
merge.sero.dat$mass_cat[merge.sero.dat$mass_forearm_residual>0] <- "good"
merge.sero.dat$mass_cat[merge.sero.dat$mass_forearm_residua==0] <- "neutral"
merge.sero.dat$mass_cat  <- factor(merge.sero.dat$mass_cat, levels = c("neutral", "poor", "good"))


  
m1 <- glm(GhV_ct~ combo_class,# +
            #roost_site,
            #family = "binomial",
            data= subset(tot.merge,  bat_species=="Eidolon dupreanum" & GhV_ct>0))
summary(m1)
plot_model(m1, type = "est")

#juveniles are associated with lower Ct values, so higher viral load

 
m2 <- glm(GhV_pos~ combo_class, #+
            #roost_site,
          family = "binomial",
          data= subset(tot.merge,  bat_species=="Eidolon dupreanum"))
summary(m2)
plot_model(m2)
#nothing is associated with pos status - some trending pos or neg but not definitive

m3 <- glm(MFI~ combo_class, #+
            #roost_site,
          #family = "binomial",
          data= subset(merge.sero.dat,  bat_species=="Eidolon dupreanum" & antigen=="GhV"))
summary(m3)
plot_model(m3, type = "est") 
# juveniles and lactating females are associated with higher MFI


m4 <- glm(seropos ~ combo_class,# +
           # roost_site,
          family = "binomial",
          data= subset(merge.sero.dat,  bat_species=="Eidolon dupreanum" & antigen=="GhV"))
summary(m4)
plot_model(m4)
#nothing associated with serostatus

#what about coinfection?
tot.merge$sum_infections <- rowSums(cbind(tot.merge$GhV_pos, tot.merge$AngV_pos, tot.merge$AndV_pos))

m5 <- glm(sum_infections ~ combo_class,
          family = "poisson",
          data= subset(tot.merge,  bat_species=="Eidolon dupreanum"))
summary(m5)
plot_model(m5, type = "est")
#juveniles are more likely to be infected with more HNVs



tot.sum <- ddply(tot.merge,.(bat_species, roost_site, combo_class), summarise, GhVpos = sum(GhV_pos), Ntot = length(GhV_pos))
tot.sum$prevalence <- tot.sum$GhVpos/tot.sum$Ntot
ggplot(tot.sum) + geom_point(aes(x=combo_class, y=prevalence)) + facet_grid(bat_species~roost_site)

ggplot(data= subset(merge.sero.dat, bat_species=="Eidolon dupreanum")) + 
  geom_boxplot(aes(x=roost_site, y=mass_forearm_residual)) #+ facet_grid(~bat_sex)


#and look at correlation MFI and MFR
ggplot(data= subset(merge.sero.dat, antigen=="GhV"  & bat_species=="Eidolon dupreanum")) + 
  geom_boxplot(aes(x=as.factor(seropos), y=mass_forearm_residual)) + facet_grid(~bat_sex)



ggplot(data= subset(merge.sero.dat, antigen=="GhV" & seropos==1 & bat_species=="Eidolon dupreanum" )) + 
            geom_point(aes(x=mass_forearm_residual, y=MFI, color=combo_class, shape=bat_sex), size=3)+ facet_grid(~roost_site)

ggplot(data=subset(merge.sero.dat, antigen=="GhV" & seropos==1 & bat_species=="Eidolon dupreanum")) + 
  geom_boxplot(aes(x=combo_class, y=MFI)) + facet_grid(~roost_site)

ggplot(data= subset(merge.sero.dat, antigen=="GhV" & seropos==1 & bat_species=="Eidolon dupreanum" & combo_class!="juvenile")) + 
  geom_point(aes(x=mass_forearm_residual, y=log10MFI, color=combo_class, shape=bat_sex))

m1 <- glm(MFI ~mass_forearm_residual:roost_site, 
              data = subset(merge.sero.dat, antigen=="GhV" & seropos_uci==1 & 
                              bat_species=="Eidolon dupreanum" & 
                              combo_class!="juvenile" & combo_class!="female_P"))
summary(m1)


merge.sero.dat$mass_cat <- NA
merge.sero.dat$mass_cat[merge.sero.dat$mass_forearm_residual>0] = "pos"
merge.sero.dat$mass_cat[merge.sero.dat$mass_forearm_residual<0] = "neg"

ggplot(data= subset(merge.sero.dat, antigen=="GhV" & seropos==1 & bat_species=="Eidolon dupreanum" & combo_class!="juvenile" & combo_class!="female_P" & !is.na(mass_cat))) + 
  geom_boxplot(aes(x=mass_cat, y=MFI)) 

ggplot(data= subset(merge.sero.dat, antigen=="GhV" & seropos==1 & bat_species=="Eidolon dupreanum" & combo_class!="juvenile" & combo_class!="female_P" & !is.na(mass_cat))) + 
  geom_boxplot(aes(x=mass_cat, y=log10MFI)) 

ggplot(data= subset(merge.sero.dat, antigen=="GhV" & bat_species=="Eidolon dupreanum" & combo_class!="juvenile" & combo_class!="female_P")) + 
  geom_boxplot(aes(x=as.factor(seropos), y=mass_forearm_residual)) 


ggplot(data= subset(merge.sero.dat, antigen=="GhV" & bat_species=="Eidolon dupreanum")) + 
  geom_boxplot(aes(x=mid_date, y=MFI, group=mid_date)) 


#and look at correlation Ct and MFR
ggplot(data= subset(tot.merge, GhV_ct>0 & bat_species=="Eidolon dupreanum"  & GhV_pos==1)) + 
  geom_point(aes(x=mass_forearm_residual, y=GhV_ct, color=bat_age_class, shape=roost_site), size=3) + facet_grid(~roost_site)

m2 <- glm(GhV_ct ~mass_forearm_residual, 
          data = subset(tot.merge, GhV_ct>0  & 
                          bat_species=="Eidolon dupreanum" ))
summary(m2)

merge.sero.dat$mass_cat <- NA
merge.sero.dat$mass_cat[merge.sero.dat$mass_forearm_residual>0] = "pos"
merge.sero.dat$mass_cat[merge.sero.dat$mass_forearm_residual<0] = "neg"

ggplot(data= subset(merge.sero.dat, antigen=="GhV" & seropos==1 & bat_species=="Eidolon dupreanum" & combo_class!="juvenile" & !is.na(mass_cat))) + 
  geom_boxplot(aes(x=mass_cat, y=MFI)) + facet_grid(~bat_sex)

ggplot(data= subset(tot.merge, bat_species=="Eidolon dupreanum")) + 
  geom_boxplot(aes(x=as.factor(GhV_pos), y=mass_forearm_residual)) 

ggplot(data= subset(tot.merge, bat_species=="Eidolon dupreanum" & GhV_ct>0)) + 
  geom_boxplot(aes(x=bat_age_class, y=GhV_ct)) + facet_grid(~roost_site)  + scale_y_reverse()


ggplot(data= subset(tot.merge, bat_species=="Eidolon dupreanum" & GhV_ct>0)) + 
  geom_boxplot(aes(x=mid_date, y=GhV_ct, group=mid_date)) + 
  facet_grid(~roost_site, space="free")  + 
  scale_y_reverse()


#Ct and MFI
head(tot.merge)
head(merge.sero.dat)

sero.dat <- dplyr::select(merge.sero.dat, sampleid, antigen, MFI, seropos, seropos_lci, seropos_uci)
sero.dat = subset(sero.dat, antigen=="GhV")

tot.new <- merge(tot.merge, sero.dat, by = "sampleid")

ggplot(subset(tot.new, bat_species=="Eidolon dupreanum" &GhV_ct>0 & seropos==1)) + geom_point(aes(x=MFI, y= GhV_ct, color=bat_age_class, shape=roost_site), size=3)


ggplot(data= subset(tot.new, bat_species=="Eidolon dupreanum" & GhV_ct>0)) + 
  geom_violin(aes(x=as.factor(seropos_lci), y=GhV_ct), draw_quantiles = c(0,.25,.5,.75)) + facet_grid(~bat_sex) +
  geom_jitter(aes(x = as.factor(seropos_lci), y=GhV_ct)) +
  scale_y_reverse()


