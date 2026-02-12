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




################################################################################
################################################################################

homewd = "/Users/carabrook/Developer/luminex-serology-madabat"
batwd = "/Users/carabrook/Developer/Madagascar-Bat-Data"

# plot of seasonality of antigens for analysis
# also GAMs of predictors of seropositivity and/or seroprevalence
################################################################################
################################################################################
## Section 1. Prep the data

# First, load the luminex data after determining who is seropositive/
# seronegative in the other scripts
lum.dat <- read.csv(file = paste0(homewd, "/data/mix_sub_data_75_90_96.csv"))
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

# now, attach the corrected doy that corresponds to the birth pulse for each region
merge.dat <- subset(merge.dat, mid_date>"2018-01-01")
merge.dat <- subset(merge.dat, mid_date>"2018-01-15")


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
## Fig. SXX: Intra-annual seasonality and other correlates
head(merge.dat)

## check seasonality intra-annual seroprevalence
merge.sum <- ddply(merge.dat, .(bat_species, roost_site, antigen, mid_date), summarise, N=length(seropos), seropos=sum(seropos))
merge.sum$doy <- yday(merge.sum$mid_date)
merge.sum$year <- year(merge.sum$mid_date)


GhVintra <- gam(cbind(seropos, seroneg) ~ 
              roost_site + 
              s(as.numeric(doy), bs = "cc", k =-1, by = roost_site) +
              s(year_f,  bs="fs", k=3, by= roost_site),
            data= subset(merge.sum, bat_species=="Eidolon dupreanum" & antigen=="GhV"), 
            family = "binomial", 
            method = "REML", 
            gamma = 1,
            select = TRUE) #

## add columns for day of year and year -- this time we use the processing date for the 
## specific bat rather than the mid_date of the sampling session
## and let's also test the doy_corrected that corresponds to 
## a bat's age with respect to the reproductive calendar which will
## allow for comparisons across sites

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

#GhV
# no seasonality by sex
m1GhVnoSeX <- gam(seropos ~ 
                  s(doy,k=7, bs="cc")+
                  s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex) +
                  s(combo_class, bs="re") +
                  s(roost_site, bs="re") +
                  s(year, bs = "re"), 
                dat = subset(merge.dat, bat_species=="Eidolon dupreanum"  & antigen=="GhV"),
                select = TRUE,
                family = "binomial") 
summary(m1GhVnoSeX)
#significant seasonality by doy and MFR in males and year at the periphery
plot(m1GhVnoSeX, pages = 1)

#with sex in seasonality
m1GhVSeX <- gam(seropos ~
             s(doy_corrected,k=7, bs="cc", by=bat_sex)+
             s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
             s(combo_class, bs="re") +
             #s(age, bs="cr") +
             s(roost_site, bs="re") +
             s(year, bs = "re"), 
           dat = subset(merge.dat, bat_species=="Eidolon dupreanum"  & antigen=="GhV"),
           select = TRUE,
           family = "binomial") 
summary(m1GhVSeX)
#significant seasonality by doy in males and MFR in males and year at the periphery
plot(m1GhVSeX, pages = 1)


AIC(m1GhVnoSeX, m1GhVSeX)
#simpler no sex model is favored though they are close

#AngV
m1AngVnoSeX <- gam(seropos ~
                    s(doy_corrected,k=7, bs="cc")+
                    #s(doy,k=-1, bs="cc")+
                    s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                    s(combo_class, bs="re") +
                    #s(age, bs="cr") +
                    s(roost_site, bs="re") +
                    s(year, bs = "re"), 
                    dat = subset(merge.dat, bat_species=="Eidolon dupreanum"  & antigen=="AngV"),
                   select = TRUE,
                    family = "binomial") 
summary(m1AngVnoSeX) #combo class borderline sig
plot(m1AngVnoSeX, pages = 1)

m1AngVSeX <- gam(seropos ~
                   s(doy_corrected,k=7, bs="cc", by=bat_sex)+
                   #s(doy,k=-1, bs="cc")+
                   s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                   s(combo_class, bs="re") +
                   #s(age, bs="cr") +
                   s(roost_site, bs="re") +
                   s(year, bs = "re"), 
                 dat = subset(merge.dat, bat_species=="Eidolon dupreanum"  & antigen=="AngV"),
                 select = TRUE,
                 family = "binomial") 
summary(m1AngVSeX)
#significant seasonality by doy in females and by combo class
plot(m1AngVSeX, pages = 1)


AIC(m1AngVnoSeX, m1AngVSeX)
#including sex improves model fit very slightly 


#NiV
m1NiVnoSeX <- gam(seropos ~
                     s(doy_corrected,k=7, bs="cc")+
                     #s(doy,k=-1, bs="cc")+
                     s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                     s(combo_class, bs="re") +
                     #s(age, bs="cr") +
                     s(roost_site, bs="re") +
                     s(year, bs = "re"), 
                    dat = subset(merge.dat, bat_species=="Eidolon dupreanum"  & antigen=="NiV"),
                    select = TRUE,
                   family = "binomial") 
summary(m1NiVnoSeX) #nothing sig
plot(m1NiVnoSeX, pages = 1)

m1NiVSeX <- gam(seropos ~
                   s(doy_corrected,k=7, bs="cc", by=bat_sex)+
                   #s(doy,k=-1, bs="cc")+
                   s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                   s(combo_class, bs="re") +
                   #s(age, bs="cr") +
                   s(roost_site, bs="re") +
                   s(year, bs = "re"), 
                 dat = subset(merge.dat, bat_species=="Eidolon dupreanum"  & antigen=="NiV"),
                select = TRUE,
                 family = "binomial")  
summary(m1NiVSeX)
#nothing significant.
plot(m1NiVSeX, pages = 1)


AIC(m1NiVnoSeX, m1NiVSeX)
# same

#HeV
m1HeVnoSeX <- gam(seropos ~
                     s(doy_corrected,k=7, bs="cc")+
                     #s(doy,k=-1, bs="cc")+
                     s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                     s(combo_class, bs="re") +
                     #s(age, bs="cr") +
                     s(roost_site, bs="re") +
                     s(year, bs = "re"), 
                   dat = subset(merge.dat, bat_species=="Eidolon dupreanum"  & antigen=="HeV"),
                  select = TRUE,
                   family = "binomial") 
summary(m1HeVnoSeX) #nothing sig
plot(m1HeVnoSeX, pages = 1)

m1HeVSeX <- gam(seropos ~
                   s(doy_corrected,k=7, bs="cc", by=bat_sex)+
                   #s(doy,k=-1, bs="cc")+
                   s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                   s(combo_class, bs="re") +
                   #s(age, bs="cr") +
                   s(roost_site, bs="re") +
                   s(year, bs = "re"), 
                 dat = subset(merge.dat, bat_species=="Eidolon dupreanum"  & antigen=="HeV"),
                select = TRUE,
                 family = "binomial")  
summary(m1HeVSeX)
#nothing sig except year
plot(m1HeVSeX, pages = 1)


AIC(m1HeVnoSeX, m1HeVSeX)
#models both bad. No sex a little better


# and Rousettus models

#GhV
# no seasonality by sex
m2GhVnoSeX <- gam(seropos ~ 
                    s(doy_corrected,k=7, bs="cc")+
                    s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                    s(combo_class, bs="re") +
                    s(roost_site, bs="re") +
                    s(year, bs = "re"), 
                  dat = subset(merge.dat, bat_species=="Rousettus madagascariensis"  & antigen=="GhV"),
                  select = TRUE,
                  family = "binomial") 
summary(m2GhVnoSeX)
#Nothing significant
plot(m2GhVnoSeX, pages = 1)

#with sex in seasonality
m2GhVSeX <- gam(seropos ~
                  s(doy_corrected,k=7, bs="cc", by=bat_sex)+
                  s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                  s(combo_class, bs="re") +
                  #s(age, bs="cr") +
                  s(roost_site, bs="re") +
                  s(year, bs = "re"), 
                dat = subset(merge.dat, bat_species=="Rousettus madagascariensis"  & antigen=="GhV"),
                select = TRUE,
                family = "binomial") 
summary(m2GhVSeX)
#nothing significant
plot(m2GhVSeX, pages = 1)

AIC(m2GhVnoSeX, m2GhVSeX) #same

#RavV
# no seasonality by sex
m2RavVnoSeX <- gam(seropos ~ 
                    s(doy_corrected,k=7, bs="cc")+
                    s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                    s(combo_class, bs="re") +
                    s(roost_site, bs="re") +
                    s(year, bs = "re"), 
                  dat = subset(merge.dat, bat_species=="Rousettus madagascariensis"  & antigen=="RavV"),
                  select = TRUE,
                  family = "binomial") 
summary(m2RavVnoSeX)
#female MFR a little significant. also, year
plot(m2RavVnoSeX, pages = 1)

#with sex in seasonality
m2RavVSeX <- gam(seropos ~
                  s(doy_corrected,k=7, bs="cc", by=bat_sex)+
                  s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                  s(combo_class, bs="re") +
                  #s(age, bs="cr") +
                  s(roost_site, bs="re") +
                  s(year, bs = "re"), 
                dat = subset(merge.dat, bat_species=="Rousettus madagascariensis"  & antigen=="RavV"),
                select = TRUE,
                family = "binomial") 
summary(m2RavVSeX)
#male season significant and year. female season borderline
plot(m2RavVSeX, pages = 1)


AIC(m2RavVnoSeX, m2RavVSeX) #including sex improves model fit



#BdbV
# no seasonality by sex
m2BdbVnoSeX <- gam(seropos ~ 
                     s(doy_corrected,k=7, bs="cc")+
                     s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                     s(combo_class, bs="re") +
                     s(roost_site, bs="re") +
                     s(year, bs = "re"), 
                   dat = subset(merge.dat, bat_species=="Rousettus madagascariensis"  & antigen=="BdbV"),
                   select = TRUE,
                   family = "binomial") 
summary(m2BdbVnoSeX)
#male MFR is sig
plot(m2BdbVnoSeX, pages = 1)

#with sex in seasonality
m2BdbVSeX <- gam(seropos ~
                   s(doy_corrected,k=7, bs="cc", by=bat_sex)+
                   s(mass_forearm_residual, k=7, fx=FALSE, bs="cr", by=bat_sex)+
                   s(combo_class, bs="re") +
                   #s(age, bs="cr") +
                   s(roost_site, bs="re") +
                   s(year, bs = "re"), 
                 dat = subset(merge.dat, bat_species=="Rousettus madagascariensis"  & antigen=="BdbV"),
                 select = TRUE,
                 family = "binomial") 
summary(m2BdbVSeX)
# roost site slightly sig and male MFR
plot(m2BdbVSeX, pages = 1)


AIC(m2BdbVnoSeX, m2BdbVSeX) #including sex improves model fit a tiny amount




# seasonal plotting

pred.df.Eid <- cbind.data.frame(doy_corrected=rep(1:365,3*4), year=rep(c(2018,2019,2020), each=(365*4)), antigen = rep(rep(c("GhV", "HeV", "AngV", "NiV"), each = 365), 3) )
#pred.df.Eid <- cbind.data.frame(doy_corrected=rep(1:365,3*4), year=rep(c(2018,2019,2020), each=(365*4)), antigen = rep(rep(c("GhV", "HeV", "AngV", "NiV"), each = 365), 3) )
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
pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="GhV"] <- predict.gam(m1GhVSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="GhV",], type="response", exclude = "s(year)")
pred.df.Eid$predicted_seropos_lci[pred.df.Eid$antigen=="GhV"] <- pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="GhV"] - 1.96*predict.gam(m1GhVSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="GhV",], type="response", se.fit = T, exclude = "s(year)")$se
pred.df.Eid$predicted_seropos_uci[pred.df.Eid$antigen=="GhV"] <- pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="GhV"] + 1.96*predict.gam(m1GhVSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="GhV",], type="response", se.fit = T, exclude = "s(year)")$se

#AngV
pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="AngV"] <- predict.gam(m1AngVSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="AngV",], type="response", exclude = "s(year)")
pred.df.Eid$predicted_seropos_lci[pred.df.Eid$antigen=="AngV"] <- pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="AngV"] - 1.96*predict.gam(m1AngVSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="AngV",], type="response", se.fit = T, exclude = "s(year)")$se
pred.df.Eid$predicted_seropos_uci[pred.df.Eid$antigen=="AngV"] <- pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="AngV"] + 1.96*predict.gam(m1AngVSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="AngV",], type="response", se.fit = T, exclude = "s(year)")$se

#HeV
pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="HeV"] <- predict.gam(m1HeVSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="HeV",], type="response", exclude = "s(year)")
pred.df.Eid$predicted_seropos_lci[pred.df.Eid$antigen=="HeV"] <- pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="HeV"] - 1.96*predict.gam(m1HeVSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="HeV",], type="response", se.fit = T, exclude = "s(year)")$se
pred.df.Eid$predicted_seropos_uci[pred.df.Eid$antigen=="HeV"] <- pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="HeV"] + 1.96*predict.gam(m1HeVSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="HeV",], type="response", se.fit = T, exclude = "s(year)")$se

#NiV
pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="NiV"] <- predict.gam(m1NiVSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="NiV",], type="response", exclude = "s(year)")
pred.df.Eid$predicted_seropos_lci[pred.df.Eid$antigen=="NiV"] <- pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="NiV"] - 1.96*predict.gam(m1NiVSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="NiV",], type="response", se.fit = T, exclude = "s(year)")$se
pred.df.Eid$predicted_seropos_uci[pred.df.Eid$antigen=="NiV"] <- pred.df.Eid$predicted_seropos[pred.df.Eid$antigen=="NiV"] + 1.96*predict.gam(m1NiVSeX, newdata = pred.df.Eid[pred.df.Eid$antigen=="NiV",], type="response", se.fit = T, exclude = "s(year)")$se


pred.df.Eid$predicted_seropos_lci[pred.df.Eid$predicted_seropos_lci<0] <- 0
pred.df.Eid$predicted_seropos_uci[pred.df.Eid$predicted_seropos_uci>1] <- 1
pred.df.Eid$bat_species="Eidolon dupreanum"



#rousettus seasonal plotting


pred.df.Rou <- cbind.data.frame(doy_corrected=rep(1:365,3*3), year=rep(c(2018,2019,2020), each=(365*3)), antigen = rep(rep(c("GhV", "RavV", "BdbV"), each = 365), 3) )

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
pred.df.Rou$predicted_seropos[pred.df.Rou$antigen=="RavV"] <- predict.gam(m2RavVSeX, newdata = pred.df.Rou[pred.df.Rou$antigen=="RavV",], type="response", exclude = "s(year)")
pred.df.Rou$predicted_seropos_lci[pred.df.Rou$antigen=="RavV"] <- pred.df.Rou$predicted_seropos[pred.df.Rou$antigen=="RavV"] - 1.96*predict.gam(m2RavVSeX, newdata = pred.df.Rou[pred.df.Rou$antigen=="RavV",], type="response", se.fit = T, exclude = "s(year)")$se
pred.df.Rou$predicted_seropos_uci[pred.df.Rou$antigen=="RavV"] <- pred.df.Rou$predicted_seropos[pred.df.Rou$antigen=="RavV"] + 1.96*predict.gam(m2RavVSeX, newdata = pred.df.Rou[pred.df.Rou$antigen=="RavV",], type="response", se.fit = T, exclude = "s(year)")$se

pred.df.Rou$predicted_seropos[pred.df.Rou$antigen=="BdbV"] <- predict.gam(m2BdbVSeX, newdata = pred.df.Rou[pred.df.Rou$antigen=="BdbV",], type="response", exclude = "s(year)")
pred.df.Rou$predicted_seropos_lci[pred.df.Rou$antigen=="BdbV"] <- pred.df.Rou$predicted_seropos[pred.df.Rou$antigen=="BdbV"] - 1.96*predict.gam(m2BdbVSeX, newdata = pred.df.Rou[pred.df.Rou$antigen=="BdbV",], type="response", se.fit = T, exclude = "s(year)")$se
pred.df.Rou$predicted_seropos_uci[pred.df.Rou$antigen=="BdbV"] <- pred.df.Rou$predicted_seropos[pred.df.Rou$antigen=="BdbV"] + 1.96*predict.gam(m2BdbVSeX, newdata = pred.df.Rou[pred.df.Rou$antigen=="BdbV",], type="response", se.fit = T, exclude = "s(year)")$se



pred.df.Rou$bat_species="Rousettus madagascariensis"


pred.df.supp <- rbind(pred.df.Eid, pred.df.Rou)# add in the timing of the reproductive season in Moramanga for both species
pred.df.supp$sig = 0
pred.df.supp$sig[pred.df.supp$bat_species=="Rousettus madagascariensis" &  pred.df.supp$antigen=="RavV"] <- 1
pred.df.supp$sig[pred.df.supp$bat_species=="Rousettus madagascariensis" &  pred.df.supp$antigen=="BdbV"] <- 2
pred.df.supp$sig[pred.df.supp$bat_species=="Eidolon dupreanum" & pred.df.supp$antigen=="GhV"] <- 1
pred.df.supp$sig[pred.df.supp$bat_species=="Eidolon dupreanum" & pred.df.supp$antigen=="NiV"] <- 0
pred.df.supp$sig[pred.df.supp$bat_species=="Eidolon dupreanum" & pred.df.supp$antigen=="AngV"] <- 2
pred.df.supp$sig[pred.df.supp$bat_species=="Eidolon dupreanum" & pred.df.supp$antigen=="HeV"] <- 0



# add in the timing of the reproductive season in Moramanga for both species
# from angelo and santino's paper
# eidolon, gestation is from: Aug-3 - Nov 15 (104 days); lactation Nov 16 - Feb 2 (78 days)
# rousettus, gestation is from: Sept 11 - Dec 11 (91 days); lactation Dec 12 - Feb 19 (69 days)
# dry season = May 15 - October 10

#here, we assume that day 1 is day 1 of lactation

#gest.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis"), start = c(215,  254), stop=c(319, 345))
gest.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis"), start = c(261,  274), stop=c(365, 365))

gest.dat$season = "gestation"
gest.dat$bat_sex="female"

#lact.dat <- cbind.data.frame(bat_species=rep(c("Eidolon dupreanum", "Rousettus madagascariensis"), each=2), start = c(320, 1, 346,1), stop=c(365, 50, 365, 33))
lact.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis"), start = c(1, 1), stop=c(78,69))
lact.dat$season = "lactation"
lact.dat$bat_sex="female"

#and winter
# let's say a moramanga winter (may 15-Oct 10) (148 days... starts 102 days after lactation for Eid and 85 days after for Rou)
cal.dat <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis"), start = c(180,154), stop=c(328,302))

cal.dat$season <- "winter"
cal.dat$bat_sex="male"

cal.dat <- rbind(cal.dat, gest.dat, lact.dat)

#collect data on prevalence just for plotting
head(merge.dat)
merge.dat$seropos <- as.numeric(as.character(merge.dat$seropos))
merge.dat$day_of_year_session <- yday(merge.dat$mid_date)
merge.plot <- ddply(merge.dat,.(bat_species,bat_sex, year, day_of_year_session, antigen), summarise, N=length(seropos), seropos=sum(seropos))
merge.plot$seroprevalence <-   merge.plot$seropos/merge.plot$N
merge.plot$seroprevalence_lci <-  binom.exact(merge.plot$seropos,merge.plot$N)$lower
merge.plot$seroprevalence_uci <-  binom.exact(merge.plot$seropos,merge.plot$N)$upper


merge.plot$sig = 0
merge.plot$sig[merge.plot$bat_species=="Rousettus madagascariensis" & merge.plot$antigen=="RavV"] <- 1
merge.plot$sig[merge.plot$bat_species=="Eidolon dupreanum" & merge.plot$antigen=="GhV"] <- 1

alphaz = c('0' = .3, '1'=.8)
alphaz.new <- c('0'=.4, '1'=.8)
pred.df.supp$sig <- as.factor(pred.df.supp$sig)
merge.plot$sig <- as.factor(merge.plot$sig)

pred.df.supp$predicted_seropos_lci[pred.df.supp$predicted_seropos_lci<0] <- 0


colorz = c('AngV' = "#F8766D", 'HeV' ="#A58AFF", 
           'NiV' = "#DB72FB", 'GhV' = "#00C094", 
           'CedV'= "#53B400", "MojV" = "#619CFF",
           'EboV' = "#06A4FF", 'BdbV' = "#00BFC4",
           'BomV' = "#00C094", 'SudV' = "#00BC56", 
           'RestV' = "#53B400", 'LloV'= "#A58AFF",
           'RavV' = "#E38900",'MlaV' = "#FF66A8")

filz <- c(colorz, 'gestation' = "yellow", 'lactation' = "greenyellow", 'winter' = "skyblue")

filz1 <- filz[1:(length(filz)-3)]
filz2 <- filz[(length(filz)-2):length(filz)]

#slim to just 2019 for plotting purposes
pIntra <- ggplot(subset(pred.df.supp, sig==1 )) + 
  geom_point(data=merge.plot, aes(x=day_of_year_session, y=seroprevalence, size=N, color=antigen, shape=year, alpha = sig)) +
  geom_linerange(data=merge.plot, aes(x=day_of_year_session, ymin=seroprevalence_lci, ymax=seroprevalence_uci, color=antigen, alpha = sig, group=year), linetype=2, show.legend = F) +
  geom_rect(data=cal.dat, aes(xmin=start, xmax=stop, ymin=0, ymax=1, fill=season), alpha=.3) +
  scale_fill_manual(values=filz2) +
  ggnewscale::new_scale_fill() +
  scale_alpha_manual(values=alphaz.new) +
  geom_ribbon(data = subset(pred.df.supp, sig==2), aes(x=doy_corrected, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.15, show.legend = F) +
  geom_line(data = subset(pred.df.supp, sig==2),aes(x=doy_corrected, y=predicted_seropos, color=antigen), alpha=.3) +
  geom_ribbon(aes(x=doy_corrected, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.3, show.legend = F) +
  geom_line(aes(x=doy_corrected, y=predicted_seropos, color=antigen)) + facet_grid(bat_sex~bat_species) + theme_bw() +
  scale_color_manual(values=colorz) + scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("1", "3", "5", "7", "9", "11")) +
  scale_fill_manual(values=filz1) +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  theme(axis.text = element_text(size=14), #axis.title.x = element_blank(), 
        axis.title = element_text(size=16), strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16),panel.spacing = unit(1.5, "lines"),
        strip.text.x = element_text(size=16, face = "italic"), panel.grid = element_blank()) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + xlab("months post-birth pulse")  + guides(alpha="none");pIntra


legendAB <- cowplot::get_legend(pIntra)


pB <- ggplot(subset(pred.df.supp, sig==1 & bat_species =="Eidolon dupreanum")) + 
  geom_point(data=subset(merge.plot, bat_species =="Eidolon dupreanum"), aes(x=day_of_year_session, y=seroprevalence, size=N, color=antigen, shape=year, alpha = sig), show.legend = F) +
  geom_linerange(data=subset(merge.plot,bat_species =="Eidolon dupreanum"), aes(x=day_of_year_session, ymin=seroprevalence_lci, ymax=seroprevalence_uci, color=antigen, alpha = sig, group=year), linetype=2, show.legend = F) +
  geom_rect(data=subset(cal.dat, bat_species =="Eidolon dupreanum"), aes(xmin=start, xmax=stop, ymin=0, ymax=1, fill=season), alpha=.3, show.legend = F) +
  scale_fill_manual(values=filz2) +
  ggnewscale::new_scale_fill() +
  scale_alpha_manual(values=alphaz.new) +
  geom_ribbon(data = subset(pred.df.supp, sig==2 & bat_species =="Eidolon dupreanum"), aes(x=doy_corrected, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.15, show.legend = F) +
  geom_line(data = subset(pred.df.supp, sig==2 & bat_species =="Eidolon dupreanum"),aes(x=doy_corrected, y=predicted_seropos, color=antigen), alpha=.3, show.legend = F) +
  geom_ribbon(aes(x=doy_corrected, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.3, show.legend = F) +
  geom_line(aes(x=doy_corrected, y=predicted_seropos, color=antigen), show.legend = F) + facet_grid(bat_sex~bat_species) + theme_bw() +
  scale_color_manual(values=colorz) + scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("1", "3", "5", "7", "9", "11")) +
  scale_fill_manual(values=filz1) +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  theme(axis.text = element_text(size=12), #axis.title.x = element_blank(), 
        axis.title = element_text(size=14), 
        strip.background.x = element_blank(),
        plot.margin = unit(c(.1,.5,.1,.1), "cm"),
        strip.background.y = element_rect(fill="white"),
        strip.text.x = element_blank(),panel.spacing = unit(.5, "lines"),
        strip.text.y = element_text(size=16), panel.grid = element_blank()) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + xlab("months post-birth pulse")  + guides(alpha="none");pB



pD <- ggplot(subset(pred.df.supp, sig==1 & bat_species =="Rousettus madagascariensis")) + 
  geom_point(data=subset(merge.plot, bat_species =="Rousettus madagascariensis"), aes(x=day_of_year_session, y=seroprevalence, size=N, color=antigen, shape=year, alpha = sig), show.legend = F) +
  geom_linerange(data=subset(merge.plot,bat_species =="Rousettus madagascariensis"), aes(x=day_of_year_session, ymin=seroprevalence_lci, ymax=seroprevalence_uci, color=antigen, alpha = sig, group=year), linetype=2, show.legend = F) +
  geom_rect(data=subset(cal.dat, bat_species =="Rousettus madagascariensis"), aes(xmin=start, xmax=stop, ymin=0, ymax=1, fill=season), alpha=.3, show.legend = F) +
  scale_fill_manual(values=filz2) +
  ggnewscale::new_scale_fill() +
  scale_alpha_manual(values=alphaz.new) +
  geom_ribbon(data = subset(pred.df.supp, sig==2 & bat_species =="Rousettus madagascariensis"), aes(x=doy_corrected, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.15, show.legend = F) +
  geom_line(data = subset(pred.df.supp, sig==2 & bat_species =="Rousettus madagascariensis"),aes(x=doy_corrected, y=predicted_seropos, color=antigen), alpha=.3, show.legend = F) +
  geom_ribbon(aes(x=doy_corrected, ymin=predicted_seropos_lci, ymax= predicted_seropos_uci, fill = antigen), alpha=.3, show.legend = F) +
  geom_line(aes(x=doy_corrected, y=predicted_seropos, color=antigen), show.legend = F) + facet_grid(bat_sex~bat_species) + theme_bw() +
  scale_color_manual(values=colorz) + scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("1", "3", "5", "7", "9", "11")) +
  scale_fill_manual(values=filz1) +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  theme(axis.text = element_text(size=12), #axis.title.x = element_blank(), 
        axis.title = element_text(size=14), 
        strip.background.x = element_blank(),
        plot.background = element_rect(fill="transparent"),
        plot.margin = unit(c(.1,.1,.1,.1), "cm"),
        strip.background.y = element_rect(fill="white"),
        strip.text.x = element_blank(),panel.spacing = unit(.5, "lines"),
        strip.text.y = element_text(size=16), panel.grid = element_blank()) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + xlab("months post-birth pulse")  + guides(alpha="none");pD



#now add in the other correlates to this plot


# MFR

#and pair with the MFR relationship

#and to get the interacting effects of sex and MFR:

#GhV
Eid.f.MFR.GhV <- plot(m1GhVSeX, pages=1)[[3]]
Eid.m.MFR.GhV <- plot(m1GhVSeX, pages=1)[[4]]

Eid.f.dat.GhV <- cbind.data.frame(MFR=Eid.f.MFR.GhV$x, effect=Eid.f.MFR.GhV$fit, se = Eid.f.MFR.GhV$se)
Eid.f.dat.GhV$yupper <- Eid.f.dat.GhV$effect+1.96*Eid.f.dat.GhV$se
Eid.f.dat.GhV$ylower <- Eid.f.dat.GhV$effect-1.96*Eid.f.dat.GhV$se
Eid.f.dat.GhV$bat_sex = "female"
#raw data included in fem.MFR$raw


Eid.m.dat.GhV <- cbind.data.frame(MFR=Eid.m.MFR.GhV$x, effect=Eid.m.MFR.GhV$fit, se = Eid.m.MFR.GhV$se)
Eid.m.dat.GhV$yupper <- Eid.m.dat.GhV$effect+1.96*Eid.m.dat.GhV$se
Eid.m.dat.GhV$ylower <- Eid.m.dat.GhV$effect-1.96*Eid.m.dat.GhV$se
Eid.m.dat.GhV$bat_sex = "male"



Eid.MFR.dat.GhV <- rbind(Eid.f.dat.GhV, Eid.m.dat.GhV)
Eid.MFR.dat.GhV$IsSignificant <- "NotSig"
Eid.MFR.dat.GhV$IsSignificant[Eid.MFR.dat.GhV $bat_sex=="male"] <- "Neg"
Eid.MFR.dat.GhV$antigen <- "GhV"


#AngV
Eid.f.MFR.AngV <- plot(m1AngVSeX, pages=1)[[3]]
Eid.m.MFR.AngV <- plot(m1AngVSeX, pages=1)[[4]]

Eid.f.dat.AngV <- cbind.data.frame(MFR=Eid.f.MFR.AngV$x, effect=Eid.f.MFR.AngV$fit, se = Eid.f.MFR.AngV$se)
Eid.f.dat.AngV$yupper <- Eid.f.dat.AngV$effect+1.96*Eid.f.dat.AngV$se
Eid.f.dat.AngV$ylower <- Eid.f.dat.AngV$effect-1.96*Eid.f.dat.AngV$se
Eid.f.dat.AngV$bat_sex = "female"
#raw data included in fem.MFR$raw


Eid.m.dat.AngV <- cbind.data.frame(MFR=Eid.m.MFR.AngV$x, effect=Eid.m.MFR.AngV$fit, se = Eid.m.MFR.AngV$se)
Eid.m.dat.AngV$yupper <- Eid.m.dat.AngV$effect+1.96*Eid.m.dat.AngV$se
Eid.m.dat.AngV$ylower <- Eid.m.dat.AngV$effect-1.96*Eid.m.dat.AngV$se
Eid.m.dat.AngV$bat_sex = "male"



Eid.MFR.dat.AngV <- rbind(Eid.f.dat.AngV, Eid.m.dat.AngV)
Eid.MFR.dat.AngV$IsSignificant <- "NotSig"
#Eid.MFR.dat.AngV$IsSignificant[Eid.MFR.dat.AngV $bat_sex=="female"] <- "Neg"
Eid.MFR.dat.AngV$antigen <- "AngV"



#NiV
Eid.f.MFR.NiV <- plot(m1NiVSeX, pages=1)[[3]]
Eid.m.MFR.NiV <- plot(m1NiVSeX, pages=1)[[4]]

Eid.f.dat.NiV <- cbind.data.frame(MFR=Eid.f.MFR.NiV$x, effect=Eid.f.MFR.NiV$fit, se = Eid.f.MFR.NiV$se)
Eid.f.dat.NiV$yupper <- Eid.f.dat.NiV$effect+1.96*Eid.f.dat.NiV$se
Eid.f.dat.NiV$ylower <- Eid.f.dat.NiV$effect-1.96*Eid.f.dat.NiV$se
Eid.f.dat.NiV$bat_sex = "female"
#raw data included in fem.MFR$raw


Eid.m.dat.NiV <- cbind.data.frame(MFR=Eid.m.MFR.NiV$x, effect=Eid.m.MFR.NiV$fit, se = Eid.m.MFR.NiV$se)
Eid.m.dat.NiV$yupper <- Eid.m.dat.NiV$effect+1.96*Eid.m.dat.NiV$se
Eid.m.dat.NiV$ylower <- Eid.m.dat.NiV$effect-1.96*Eid.m.dat.NiV$se
Eid.m.dat.NiV$bat_sex = "male"



Eid.MFR.dat.NiV <- rbind(Eid.f.dat.NiV, Eid.m.dat.NiV)
Eid.MFR.dat.NiV$IsSignificant <- "NotSig"
#Eid.MFR.dat.NiV$IsSignificant[Eid.MFR.dat.NiV $bat_sex=="female"] <- "Neg"
Eid.MFR.dat.NiV$antigen <- "NiV"

#HeV
Eid.f.MFR.HeV <- plot(m1HeVSeX, pages=1)[[3]]
Eid.m.MFR.HeV <- plot(m1HeVSeX, pages=1)[[4]]

Eid.f.dat.HeV <- cbind.data.frame(MFR=Eid.f.MFR.HeV$x, effect=Eid.f.MFR.HeV$fit, se = Eid.f.MFR.HeV$se)
Eid.f.dat.HeV$yupper <- Eid.f.dat.HeV$effect+1.96*Eid.f.dat.HeV$se
Eid.f.dat.HeV$ylower <- Eid.f.dat.HeV$effect-1.96*Eid.f.dat.HeV$se
Eid.f.dat.HeV$bat_sex = "female"
#raw data included in fem.MFR$raw


Eid.m.dat.HeV <- cbind.data.frame(MFR=Eid.m.MFR.HeV$x, effect=Eid.m.MFR.HeV$fit, se = Eid.m.MFR.HeV$se)
Eid.m.dat.HeV$yupper <- Eid.m.dat.HeV$effect+1.96*Eid.m.dat.HeV$se
Eid.m.dat.HeV$ylower <- Eid.m.dat.HeV$effect-1.96*Eid.m.dat.HeV$se
Eid.m.dat.HeV$bat_sex = "male"



Eid.MFR.dat.HeV <- rbind(Eid.f.dat.HeV, Eid.m.dat.HeV)
Eid.MFR.dat.HeV$IsSignificant <- "NotSig"
#Eid.MFR.dat.HeV$IsSignificant[Eid.MFR.dat.HeV $bat_sex=="female"] <- "Neg"
Eid.MFR.dat.HeV$antigen <- "HeV"


#consolidate MFR dat for eidolon
Eid.MFR.dat <- rbind(Eid.MFR.dat.GhV, Eid.MFR.dat.AngV, Eid.MFR.dat.NiV, Eid.MFR.dat.HeV)
Eid.MFR.dat$bat_species="Eidolon dupreanum"


# and add in Rousettus


#GhV
Rou.f.MFR.GhV <- plot(m2GhVSeX, pages=1)[[3]]
Rou.m.MFR.GhV <- plot(m2GhVSeX, pages=1)[[4]]

Rou.f.dat.GhV <- cbind.data.frame(MFR=Rou.f.MFR.GhV$x, effect=Rou.f.MFR.GhV$fit, se = Rou.f.MFR.GhV$se)
Rou.f.dat.GhV$yupper <- Rou.f.dat.GhV$effect+1.96*Rou.f.dat.GhV$se
Rou.f.dat.GhV$ylower <- Rou.f.dat.GhV$effect-1.96*Rou.f.dat.GhV$se
Rou.f.dat.GhV$bat_sex = "female"
#raw data included in fem.MFR$raw


Rou.m.dat.GhV <- cbind.data.frame(MFR=Rou.m.MFR.GhV$x, effect=Rou.m.MFR.GhV$fit, se = Rou.m.MFR.GhV$se)
Rou.m.dat.GhV$yupper <- Rou.m.dat.GhV$effect+1.96*Rou.m.dat.GhV$se
Rou.m.dat.GhV$ylower <- Rou.m.dat.GhV$effect-1.96*Rou.m.dat.GhV$se
Rou.m.dat.GhV$bat_sex = "male"



Rou.MFR.dat.GhV <- rbind(Rou.f.dat.GhV, Rou.m.dat.GhV)
Rou.MFR.dat.GhV$IsSignificant <- "NotSig"
#Rou.MFR.dat.GhV$IsSignificant[Rou.MFR.dat.GhV $bat_sex=="male"] <- "Neg"
Rou.MFR.dat.GhV$antigen <- "GhV"



#RavV
Rou.f.MFR.RavV <- plot(m2RavVSeX, pages=1)[[3]]
Rou.m.MFR.RavV <- plot(m2RavVSeX, pages=1)[[4]]

Rou.f.dat.RavV <- cbind.data.frame(MFR=Rou.f.MFR.RavV$x, effect=Rou.f.MFR.RavV$fit, se = Rou.f.MFR.RavV$se)
Rou.f.dat.RavV$yupper <- Rou.f.dat.RavV$effect+1.96*Rou.f.dat.RavV$se
Rou.f.dat.RavV$ylower <- Rou.f.dat.RavV$effect-1.96*Rou.f.dat.RavV$se
Rou.f.dat.RavV$bat_sex = "female"
#raw data included in fem.MFR$raw


Rou.m.dat.RavV <- cbind.data.frame(MFR=Rou.m.MFR.RavV$x, effect=Rou.m.MFR.RavV$fit, se = Rou.m.MFR.RavV$se)
Rou.m.dat.RavV$yupper <- Rou.m.dat.RavV$effect+1.96*Rou.m.dat.RavV$se
Rou.m.dat.RavV$ylower <- Rou.m.dat.RavV$effect-1.96*Rou.m.dat.RavV$se
Rou.m.dat.RavV$bat_sex = "male"



Rou.MFR.dat.RavV <- rbind(Rou.f.dat.RavV, Rou.m.dat.RavV)
Rou.MFR.dat.RavV$IsSignificant <- "NotSig"
#Rou.MFR.dat.RavV$IsSignificant[Rou.MFR.dat.RavV $bat_sex=="male"] <- "Neg"
Rou.MFR.dat.RavV$antigen <- "RavV"


#BdbV
Rou.f.MFR.BdbV <- plot(m2BdbVSeX, pages=1)[[3]]
Rou.m.MFR.BdbV <- plot(m2BdbVSeX, pages=1)[[4]]

Rou.f.dat.BdbV <- cbind.data.frame(MFR=Rou.f.MFR.BdbV$x, effect=Rou.f.MFR.BdbV$fit, se = Rou.f.MFR.BdbV$se)
Rou.f.dat.BdbV$yupper <- Rou.f.dat.BdbV$effect+1.96*Rou.f.dat.BdbV$se
Rou.f.dat.BdbV$ylower <- Rou.f.dat.BdbV$effect-1.96*Rou.f.dat.BdbV$se
Rou.f.dat.BdbV$bat_sex = "female"
#raw data included in fem.MFR$raw


Rou.m.dat.BdbV <- cbind.data.frame(MFR=Rou.m.MFR.BdbV$x, effect=Rou.m.MFR.BdbV$fit, se = Rou.m.MFR.BdbV$se)
Rou.m.dat.BdbV$yupper <- Rou.m.dat.BdbV$effect+1.96*Rou.m.dat.BdbV$se
Rou.m.dat.BdbV$ylower <- Rou.m.dat.BdbV$effect-1.96*Rou.m.dat.BdbV$se
Rou.m.dat.BdbV$bat_sex = "male"



Rou.MFR.dat.BdbV <- rbind(Rou.f.dat.BdbV, Rou.m.dat.BdbV)
Rou.MFR.dat.BdbV$IsSignificant <- "NotSig"
Rou.MFR.dat.BdbV$IsSignificant[Rou.MFR.dat.BdbV $bat_sex=="male"] <- "Neg"
Rou.MFR.dat.BdbV$antigen <- "BdbV"



#consolidate MFR dat for eidolon
Rou.MFR.dat <- rbind(Rou.MFR.dat.GhV, Rou.MFR.dat.RavV, Rou.MFR.dat.BdbV)
Rou.MFR.dat$bat_species="Rousettus madagascariensis"


Rou.MFR.dat$plot <- 0
Rou.MFR.dat$plot[Rou.MFR.dat$antigen=="RavV" | Rou.MFR.dat$antigen=="BdbV" ] <- 1

Eid.MFR.dat$plot <- 0
Eid.MFR.dat$plot[Eid.MFR.dat$antigen=="GhV" | Eid.MFR.dat$antigen=="AngV" ] <- 1

# and have a combined plot for just the antigens of interest:
# GhV/AngV in Eidolon and BdbV/RavV in Rousettus
MFR.dat <- rbind(Eid.MFR.dat, Rou.MFR.dat)
MFR.dat$antigen <- as.factor(MFR.dat$antigen)
MFR.dat$IsSignificant <- as.factor(MFR.dat$IsSignificant)


alphaz.ribbon = c('Neg'=.5, 'NotSig'=.1)
alphaz.line = c('Neg'=1, 'NotSig'=.1)


p8 <- ggplot(data=MFR.dat) + 
  facet_grid(bat_sex~bat_species, scales = "free_x") +
  geom_ribbon(aes(x=MFR, ymin=ylower, ymax=yupper, fill=antigen, alpha = IsSignificant),  show.legend = F) +
  scale_color_manual(values = colorz) +
  scale_alpha_manual(values = alphaz.ribbon) +
  scale_fill_manual(values = colorz) +
  ggnewscale::new_scale("alpha")+
  geom_line(data=MFR.dat, aes(x=MFR, y=effect, color=antigen, alpha = IsSignificant), show.legend = F) +
  scale_alpha_manual(values = alphaz.line) +
  theme_bw() +
  ylab("effect on seropositivity") +
  xlab("mass : forearm residual") +
  theme(strip.background = element_rect(fill="white"), 
        strip.text.x = element_text(size=16, face="italic"),
        strip.text.y = element_text(size=16),
        panel.spacing = unit(1.4, "cm"),
        panel.grid = element_blank(),  
        axis.title = element_text(size=16), axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle=360),  
        plot.margin = unit(c(.2,.1,.1,.1), "cm"));p8



pA <- ggplot(data=subset(MFR.dat, bat_species=="Eidolon dupreanum")) + 
  facet_grid(bat_sex~bat_species) + coord_cartesian(ylim=c(-45,20)) +
  geom_ribbon(aes(x=MFR, ymin=ylower, ymax=yupper, fill=antigen, alpha = IsSignificant),  show.legend = F) +
  scale_color_manual(values = colorz) +
  scale_alpha_manual(values = alphaz.ribbon) +
  scale_fill_manual(values = colorz) +
  ggnewscale::new_scale("alpha")+
  geom_line(data=subset(MFR.dat, bat_species=="Eidolon dupreanum"), aes(x=MFR, y=effect, color=antigen, alpha = IsSignificant), show.legend = F) +
  scale_alpha_manual(values = alphaz.line) +
  theme_bw() +
  ylab("effect on seropositivity") +
  xlab("mass : forearm residual") +
  theme(#strip.background = element_rect(fill="white"), 
        #strip.text.x = element_text(size=16, face="italic"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(.5, "cm"),
        panel.grid = element_blank(),  
        axis.title = element_text(size=14), axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle=360),  
        plot.margin = unit(c(.1,.1,.1,.1), "cm"));pA

pC <- ggplot(data=subset(MFR.dat, bat_species=="Rousettus madagascariensis")) + 
  facet_grid(bat_sex~bat_species) + coord_cartesian(ylim=c(-45,20)) +
  geom_ribbon(aes(x=MFR, ymin=ylower, ymax=yupper, fill=antigen, alpha = IsSignificant),  show.legend = F) +
  scale_color_manual(values = colorz) +
  scale_alpha_manual(values = alphaz.ribbon) +
  scale_fill_manual(values = colorz) +
  ggnewscale::new_scale("alpha")+
  geom_line(data=subset(MFR.dat,bat_species=="Rousettus madagascariensis"), aes(x=MFR, y=effect, color=antigen, alpha = IsSignificant), show.legend = F) +
  scale_alpha_manual(values = alphaz.line) +
  theme_bw() +
  ylab("effect on seropositivity") +
  xlab("mass : forearm residual") +
  theme(#strip.background = element_rect(fill="white"), 
        strip.text = element_blank(),
        strip.background = element_blank(),
        strip.background.y = element_blank(),
        plot.background = element_rect(fill="transparent"),
        panel.spacing = unit(.5, "cm"),
        panel.grid = element_blank(),  
        axis.title = element_text(size=14), 
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle=360),  
        plot.margin = unit(c(.1,.1,.1,.5), "cm"));pC



titleEid <- ggdraw() + 
  draw_label(
    "Eidolon dupreanum",
    fontface = 'italic',
    size=16,
    x = 0,
    vjust = .65,
    hjust = -.1
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(.1, .1, .1, 20)
  )


titleRou <- ggdraw() + 
  draw_label(
    "Rousettus madagascariensis",
    fontface = 'italic',
    size=16,
    x = 0,
    vjust = .65,
    hjust = -.1
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(.1, .1, .1, 20)
  )


pEid1 <- cowplot::plot_grid(pA,pB, ncol = 2,  rel_widths = c(.75,1), align = "h");pEid1
pEid <- cowplot::plot_grid(titleEid, pEid1, nrow=2,  rel_heights = c(.1,1));pEid


pRou1 <- cowplot::plot_grid(pC,pD, ncol = 2,  rel_widths = c(.8,1), align = "h");pRou1
pRou <- cowplot::plot_grid(titleRou, pRou1, nrow=2,  rel_heights = c(.1,1));pRou


pFigSX <- cowplot::plot_grid(pEid,pRou,legendAB, ncol = 3, align = "h", labels = c("A", "B"), label_y = .99, rel_widths = c(1,1,.15), label_size = 20);pFigSX



ggsave(file = paste0(homewd,"/figures/FigSX-intraannual-MFR.png"),
       plot = pFigSX,
       bg="white",
       units="mm",  
       width=110, 
       height=45, 
       scale=3.2, 
       dpi=300)







# let's look at the partial effects
# because of the interaction term, the doy function
# won't run, but you could pull the partial effects from 
# the plot() output above and visualize it if desired

# source script from Guth/Streicker/Mollentze
source(paste0(homewd, "/scripts/partial-gam-functions.R"))

# combo class first
Eid.combo.df.GhV <- get_partial_effects(m1GhVSeX, var = "combo_class")
plot.partial(Eid.combo.df.GhV, var="combo_class", response_var = "prob. seropositive") 
Eid.combo.df.AngV <- get_partial_effects(m1AngVSeX, var = "combo_class") # this is the only significant addition in the model
plot.partial(Eid.combo.df.AngV, var="combo_class", response_var = "prob. seropositive") # no sig differences though allowing flexibility here improves the whole model
Eid.combo.df.NiV <- get_partial_effects(m1NiVSeX, var = "combo_class")
plot.partial(Eid.combo.df.NiV, var="combo_class", response_var = "prob. seropositive") 
Eid.combo.df.HeV <- get_partial_effects(m1HeVSeX, var = "combo_class")
plot.partial(Eid.combo.df.HeV, var="combo_class", response_var = "prob. seropositive") 

Rou.combo.df.GhV <- get_partial_effects(m2GhVSeX, var = "combo_class")
plot.partial(Rou.combo.df.GhV, var="combo_class", response_var = "prob. seropositive") 
Rou.combo.df.RavV <- get_partial_effects(m2RavVSeX, var = "combo_class")
plot.partial(Rou.combo.df.RavV, var="combo_class", response_var = "prob. seropositive") 
Rou.combo.df.BdbV <- get_partial_effects(m2BdbVSeX, var = "combo_class") 
plot.partial(Rou.combo.df.BdbV, var="combo_class", response_var = "prob. seropositive") 
Rou.roost.df.BdbV <- get_partial_effects(m2BdbVSeX, var = "roost_site") # this is the only significant addition in the model
plot.partial(Rou.roost.df.BdbV, var="roost_site", response_var = "prob. seropositive") # no sig differences though allowing flexibility here improves the whole model


# because there were so few significant effects, we chose not to include these in the final model


