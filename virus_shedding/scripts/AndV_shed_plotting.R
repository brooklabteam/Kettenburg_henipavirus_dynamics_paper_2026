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
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/virus_shedding/data"
datawd = "/Users/gwenddolenkettenburg/Desktop/brook_lab/malagasy_fruit_bats_population_level_work/molecular_results"

# plot of seasonality of viruss for analysis
# also GAMs of predictors of positivity and/or prevalence
################################################################################
################################################################################
## Section 1. Prep the data

# First, load the viral shedding data - there is also a file for all molecular pos and all mNGS but we won't use them here
shed.dat <- read.csv(file = paste0(homewd, "/henipavirus_shedding_data.csv"))
head(shed.dat)

#slim to just AndV, AngV, GhV, and GhV-like data: one for just being positive and one for the ct values
shed.dat.pos <- dplyr::select(shed.dat, sampleid, AndV_shedding_ur,
                              AndV_shedding_ct_ur)

## merge the shedding data with the catching data

# and make into long format for prevalence
molecular.long.pos <- melt(shed.dat.pos, id.vars = c("sampleid"), value.name = "pos", variable.name = "virus")
head(molecular.long.pos)
molecular.long.pos$virus <- as.character(molecular.long.pos$virus) 
molecular.long.pos$virus <- gsub(pattern="qPCR_", replacement = "", x=molecular.long.pos$virus, fixed = T)
names(molecular.long.pos)[1] <- "sampleid"
#complete cases
molecular.long.pos <- molecular.long.pos[complete.cases(molecular.long.pos),]

# link this to the catching data to plot by sampling session, add the aging data as well
catch.dat <- read.csv(file = paste0(datawd, "/catching_data_clean_asofApril2025_edited.csv"), header = T, stringsAsFactors = F)
barcode.dat<- read.csv(file = paste0(datawd, "/barcode_ids.csv"), header = T, stringsAsFactors = F)
age.dat<-read.csv(file = paste0(datawd, "/aged_bats.csv"), header = T, stringsAsFactors = F)

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, region, roost_site, bat_species, processing_date, sampleid, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new,
                          year, quarter, month, doy, dataset, combo_class)
age.dat <-dplyr::select(age.dat, sampleid, age)

#and merge the catching data with the barcode and age data first
merge.dat <- merge(catch.dat, barcode.dat, by="sampleid", all.x = T)
head(merge.dat)

merge.dat <- merge(merge.dat, age.dat, by="sampleid", all.x = T)
head(merge.dat)

#now add the viral shedding data to the merged catching data
merge.dat.pos <- merge(molecular.long.pos, merge.dat, by="sampleid", all.x = T)
head(merge.dat.pos)

#get the midpoint processing date per sampling session
merge.dat.pos$processing_date <- as.Date(merge.dat.pos$processing_date)
sampling.date <- ddply(merge.dat.pos, .(bat_species, sampling_session), summarise, min_date = min(processing_date), max_date = max(processing_date))
sampling.date$mid_date <- (sampling.date$max_date-sampling.date$min_date)/2 + sampling.date$min_date

#and slim to merge
sampling.merge <- dplyr::select(sampling.date, bat_species, sampling_session, mid_date)

merge.dat.pos <- merge(merge.dat.pos, sampling.merge, by = c("bat_species", "sampling_session"))
merge.dat.pos$year <- year(merge.dat.pos$processing_date)
merge.dat.pos$month <- as.character(month(merge.dat.pos$processing_date))
head(merge.dat.pos)

#what species?
unique(merge.dat.pos$bat_species) #NA is the owl, get rid of that

merge.dat.pos<-subset(merge.dat.pos, bat_species=="Eidolon dupreanum" | bat_species=="Rousettus madagascariensis" | bat_species=="Pteropus rufus")
unique(merge.dat.pos$bat_species)

####################################################################################################
####################################################################################################
## SECTION 2: Seasonality of prevalence 

# now look at prevalence by sampling session by species and virus
# plot data first

#slim to those that are positive with any target for now, can come back and change

#rename virus for simplicity
merge.dat.pos$virus[merge.dat.pos$virus=="AndV_shedding_ur"] <- "AndV"
merge.dat.pos$virus[merge.dat.pos$virus=="AndV_shedding_ct_ur"] <- "AndV_ct"

#subset into the three time series datasets
merge.dat.pos_old <- merge.dat.pos[merge.dat.pos$mid_date > "2014-10-13" & merge.dat.pos$mid_date < "2016-01-19", ]
merge.dat.pos_mid <- merge.dat.pos[merge.dat.pos$mid_date > "2018-02-15" & merge.dat.pos$mid_date < "2019-09-21", ]
merge.dat.pos_new <- merge.dat.pos[merge.dat.pos$mid_date > "2022-07-19" & merge.dat.pos$mid_date < "2023-05-12", ]

prev.sum.old <- ddply(merge.dat.pos_old, .(bat_species, region,virus, mid_date), summarise, N=length(pos), pos=sum(pos))
prev.sum.old$neg <- prev.sum.old$N - prev.sum.old$pos
prev.sum.old$prev <- prev.sum.old$pos/prev.sum.old$N
prev.sum.old$prev_lci <-binom.exact(x=prev.sum.old$pos, n=prev.sum.old$N)$lower
prev.sum.old$prev_uci <-binom.exact(x=prev.sum.old$pos, n=prev.sum.old$N)$upper
head(prev.sum.old)

prev.sum.mid <- ddply(merge.dat.pos_mid, .(bat_species, region,virus, mid_date), summarise, N=length(pos), pos=sum(pos))
prev.sum.mid$neg <- prev.sum.mid$N - prev.sum.mid$pos
prev.sum.mid$prev <- prev.sum.mid$pos/prev.sum.mid$N
prev.sum.mid$prev_lci <-binom.exact(x=prev.sum.mid$pos, n=prev.sum.mid$N)$lower
prev.sum.mid$prev_uci <-binom.exact(x=prev.sum.mid$pos, n=prev.sum.mid$N)$upper
head(prev.sum.mid)

prev.sum.new <- ddply(merge.dat.pos_new, .(bat_species, region,virus, mid_date), summarise, N=length(pos), pos=sum(pos))
prev.sum.new$neg <- prev.sum.new$N - prev.sum.new$pos
prev.sum.new$prev <- prev.sum.new$pos/prev.sum.new$N
prev.sum.new$prev_lci <-binom.exact(x=prev.sum.new$pos, n=prev.sum.new$N)$lower
prev.sum.new$prev_uci <-binom.exact(x=prev.sum.new$pos, n=prev.sum.new$N)$upper
head(prev.sum.new)

#define colors per virus

# colz1 = c('winter' = "#678A76", "gestation"="#0A3C79", "lactation"="black",
#           'GhV' = "deepskyblue2", 'GhV-like' = "tomato1",
#           'AngV' = "#C49A00", 'AndV' = "#525")
# filz = c('winter' = "#678A76", "gestation"="#0A3C79", "lactation"="black",
#          'GhV' = "deepskyblue2", 'GhV-like' = "tomato1",
#          'AngV' = "#C49A00", 'AndV' = "#525")

colz1 = c('winter' = "#678A76", "gestation"="#0A3C79", "lactation"="black",
           
           'AngV' = "#F8766D", 'HeV' ="#E38900", 
           'NiV' = "goldenrod1", 'GhV' = "#7CAE00", 'GhV-like' = "#525", 
           'AndV'= "darkgreen",
           
           'Pruf_CoV'="#C77CFF", 'Rmad_CoV'="#FF6A98",
           'HKU9'="#35A2FF")

filz = c('winter' = "#678A76", "gestation"="#0A3C79", "lactation"="black",
         
         'AngV' = "#F8766D", 'HeV' ="#E38900", 
         'NiV' = "goldenrod1", 'GhV' = "#7CAE00", 'GhV-like' = "#525", 
         'AndV'= "darkgreen",
         
         'Pruf_CoV'="#C77CFF", 'Rmad_CoV'="#FF6A98",
         'HKU9'="#35A2FF")

p1 <- ggplot(data=prev.sum.mid) + geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) +
  geom_line(aes(x=mid_date, y=prev, color=virus)) + scale_color_manual(values = colz1) +
  facet_grid(bat_species~.) + theme_bw() 
p1

# add in the timing of the reproductive season for both species
# from angelo and santino's paper
# eidolon, gestation is from: Aug-3 - Nov 15; lactation Nov 16- 2 feb
# rousettus, gestation is from: sept 11 - Dec 11; lactation Dec 12- 19 feb
# pteropus, gestation is from: July 7 - Aug 2; lactation Aug 3- Jan 21
# dry season = may 15 - october 10

#all dates together, keep in case I need it
#old dataset 2014-2016
lact.dat.old <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum",
                                               "Rousettus madagascariensis","Rousettus madagascariensis", "Rousettus madagascariensis"), 
                                 start = c(as.Date("2013-11-16"), as.Date("2014-11-16"), as.Date("2015-11-16"), 
                                           as.Date("2013-12-12"), as.Date("2014-12-12"), as.Date("2015-12-12")), 
                                 stop=c(as.Date("2014-02-02"), as.Date("2015-02-02"), as.Date("2016-02-02"),
                                        as.Date("2014-02-19"), as.Date("2015-02-19"), as.Date("2016-02-19")))
lact.dat.old$type = "lactation"
lact.dat.old$bat_sex ="female"
gest.dat.old <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum",
                                               "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis"), 
                                 start = c(as.Date("2013-08-03"), as.Date("2014-08-03"), as.Date("2015-08-03"), 
                                           as.Date("2013-09-11"), as.Date("2014-09-11"), as.Date("2015-09-11")), 
                                 stop=c(as.Date("2013-11-15"), as.Date("2014-11-15"), as.Date("2015-11-15"),
                                        as.Date("2013-12-11"), as.Date("2014-12-11"), as.Date("2015-12-11")))
gest.dat.old$type = "gestation"
gest.dat.old$bat_sex="female"

#and winter
cal.dat.old <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", 
                                               "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis"), 
                                 start = c(as.Date("2013-05-15"), as.Date("2014-05-15"), as.Date("2015-05-15"), 
                                           as.Date("2013-05-15"), as.Date("2014-05-15"), as.Date("2015-05-15")), 
                                 stop=c(as.Date("2013-10-10"), as.Date("2014-10-10"), as.Date("2015-10-10"), 
                                        as.Date("2013-10-10"), as.Date("2014-10-10"), as.Date("2015-10-10")))
cal.dat.old$type <- "winter"
cal.dat.old$bat_sex = "male"

library(reshape2)

repro.dat.old <- rbind(lact.dat.old, gest.dat.old)
demo.dat.old <- rbind(gest.dat.old, lact.dat.old, cal.dat.old)


#mid dataset 2018-2020
lact.dat.mid <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum",  
                                               "Rousettus madagascariensis","Rousettus madagascariensis"), 
                                 start = c(as.Date("2018-11-16"), as.Date("2019-11-16"), 
                                           as.Date("2018-12-12"), as.Date("2019-12-12")), 
                                 stop=c(as.Date("2019-02-02"), as.Date("2020-02-02"),
                                        as.Date("2019-02-19"), as.Date("2020-02-19")))
lact.dat.mid$type = "lactation"
lact.dat.mid$bat_sex ="female"
gest.dat.mid <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum",
                                               "Rousettus madagascariensis", "Rousettus madagascariensis"), 
                                 start = c(as.Date("2018-08-03"), as.Date("2019-08-03"), 
                                           as.Date("2018-09-11"), as.Date("2019-09-11")), 
                                 stop=c(as.Date("2018-11-15"), as.Date("2019-11-15"),
                                        as.Date("2018-12-11"), as.Date("2019-12-11")))
gest.dat.mid$type = "gestation"
gest.dat.mid$bat_sex= "female"

#and winter
cal.dat.mid <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", 
                                               "Rousettus madagascariensis", "Rousettus madagascariensis"), 
                                 start = c(as.Date("2018-05-15"), as.Date("2019-05-15"), 
                                           as.Date("2018-05-15"), as.Date("2019-05-15")), 
                                 stop=c(as.Date("2018-10-10"), as.Date("2019-10-10"),
                                        as.Date("2018-10-10"), as.Date("2019-10-10")))
cal.dat.mid$type <- "winter"
cal.dat.mid$bat_sex ="male"

library(reshape2)

repro.dat.mid <- rbind(lact.dat.mid, gest.dat.mid)
demo.dat.mid <- rbind(gest.dat.mid, lact.dat.mid, cal.dat.mid)

demo.dat.old.ed = subset(demo.dat.old, bat_species=="Eidolon dupreanum")
demo.dat.mid.ed = subset(demo.dat.mid, bat_species=="Eidolon dupreanum")
demo.dat.old.rm = subset(demo.dat.old, bat_species=="Rousettus madagascariensis")
demo.dat.mid.rm = subset(demo.dat.mid, bat_species=="Rousettus madagascariensis")


#now add to plot

p2 <- ggplot(data=subset(prev.sum.mid)) + 
  geom_rect(data=demo.dat.mid, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) +
  geom_line(aes(x=mid_date, y=prev, color=virus)) + scale_color_manual(values = colz1) +
  facet_grid(bat_species~.) + theme_bw()

p2


#what if you did males and females separately?
prev.sum.sex.old <- ddply(merge.dat.pos_old, .(bat_species, region, bat_sex, virus, sampling_session,  mid_date), summarise, N=length(pos), pos=sum(pos))
prev.sum.sex.old$neg <- prev.sum.sex.old$N - prev.sum.sex.old$pos
prev.sum.sex.old$prev <- prev.sum.sex.old$pos/prev.sum.sex.old$N
prev.sum.sex.old$prev_lci <-binom.exact(x=prev.sum.sex.old$pos, n=prev.sum.sex.old$N)$lower
prev.sum.sex.old$prev_uci <-binom.exact(x=prev.sum.sex.old$pos, n=prev.sum.sex.old$N)$upper
head(prev.sum.sex.old)
prev.sum.sex.old = subset(prev.sum.sex.old, bat_sex=="female" | bat_sex=="male")

prev.sum.sex.mid <- ddply(merge.dat.pos_mid, .(bat_species, region, bat_sex, virus, sampling_session, mid_date), summarise, N=length(pos), pos=sum(pos))
prev.sum.sex.mid$neg <- prev.sum.sex.mid$N - prev.sum.sex.mid$pos
prev.sum.sex.mid$prev <- prev.sum.sex.mid$pos/prev.sum.sex.mid$N
prev.sum.sex.mid$prev_lci <-binom.exact(x=prev.sum.sex.mid$pos, n=prev.sum.sex.mid$N)$lower
prev.sum.sex.mid$prev_uci <-binom.exact(x=prev.sum.sex.mid$pos, n=prev.sum.sex.mid$N)$upper
head(prev.sum.sex.mid)
prev.sum.sex.mid = subset(prev.sum.sex.mid, bat_sex=="female" | bat_sex=="male")

prev.sum.sex.new <- ddply(merge.dat.pos_new, .(bat_species, region, bat_sex, virus, sampling_session, mid_date), summarise, N=length(pos), pos=sum(pos))
prev.sum.sex.new$neg <- prev.sum.sex.new$N - prev.sum.sex.new$pos
prev.sum.sex.new$prev <- prev.sum.sex.new$pos/prev.sum.sex.new$N
prev.sum.sex.new$prev_lci <-binom.exact(x=prev.sum.sex.new$pos, n=prev.sum.sex.new$N)$lower
prev.sum.sex.new$prev_uci <-binom.exact(x=prev.sum.sex.new$pos, n=prev.sum.sex.new$N)$upper
head(prev.sum.sex.new)
prev.sum.sex.new = subset(prev.sum.sex.new, bat_sex=="female" | bat_sex=="male")

p3 <- ggplot(data=subset(prev.sum.sex.mid)) + 
  geom_rect(data=demo.dat.mid, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values = colz1) +
  geom_line(aes(x=mid_date, y=prev, color=virus)) +
  facet_grid(bat_species~bat_sex) + theme_bw()

p3

# time series is too short here, so we will focus seasonality analyses on the 

# you want to test whether including sex or virus matters in driving prevalence
prev.sum.sex.old$virus <- as.factor(prev.sum.sex.old$virus)
prev.sum.sex.old$bat_sex <- as.factor(prev.sum.sex.old$bat_sex)
prev.sum.sex.old$region <- as.factor(prev.sum.sex.old$region)

prev.sum.sex.mid$virus <- as.factor(prev.sum.sex.mid$virus)
prev.sum.sex.mid$bat_sex <- as.factor(prev.sum.sex.mid$bat_sex)
prev.sum.sex.mid$region <- as.factor(prev.sum.sex.mid$region)

# prior work suggests important sex-based differences, so we will 
# model sex interactions using the "by" term and as a fixed effect
# but model each virus separately

# research question to be answered: is there sex-specific seasonality in
# the qPCR data for each virus?

# for each species and virus, we first test: 
# is there seasonality in prevalence by virus without accounting for sex?
# we restrict our analyses to just the Moramanga data


##old dataset
# GhV
gamEidNoSexGhV.old <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.old, bat_species=="Eidolon dupreanum"  & virus=="GhV"), family = "binomial", select = TRUE)
summary(gamEidNoSexGhV.old) #significant
plot(gamEidNoSexGhV.old, pages = 1)

# AngV
gamEidNoSexAngV.old <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.old, bat_species=="Eidolon dupreanum"   & virus=="AngV"), family = "binomial", select = TRUE)
summary(gamEidNoSexAngV.old) #significant
plot(gamEidNoSexAngV.old, pages = 1)

# GhV-like
gamEidNoSexGhVlike.old <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.old, bat_species=="Eidolon dupreanum"  & virus=="GhV-like"), family = "binomial", select = TRUE)
summary(gamEidNoSexGhVlike.old) #sig
plot(gamEidNoSexGhVlike.old, pages = 1)

# AndV
gamEidNoSexAndV.old <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.old, bat_species=="Eidolon dupreanum"  & virus=="AndV"), family = "binomial", select = TRUE)
summary(gamEidNoSexAndV.old) #sig
plot(gamEidNoSexAndV.old, pages = 1)

##mid dataset
# GhV
prev.sum.mid$virus <- as.factor(prev.sum.mid$virus)
gamEidNoSexGhV.mid <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.mid, bat_species=="Eidolon dupreanum"  & virus=="GhV"), family = "binomial", select = TRUE)
summary(gamEidNoSexGhV.mid) #almost sig
plot(gamEidNoSexGhV.mid, pages = 1)

# AngV
gamEidNoSexAngV.mid <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.mid, bat_species=="Eidolon dupreanum"  & virus=="AngV"), family = "binomial", select = TRUE)
summary(gamEidNoSexAngV.mid) #not sig
plot(gamEidNoSexAngV.mid, pages = 1)

# GhV-like
gamEidNoSexGhVlike.mid <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.mid, bat_species=="Eidolon dupreanum"  & virus=="GhV-like"), family = "binomial", select = TRUE)
summary(gamEidNoSexGhVlike.mid) #not sig
plot(gamEidNoSexGhVlike.mid, pages = 1)

# AndV
gamEidNoSexAndV.mid <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.mid, bat_species=="Eidolon dupreanum"  & virus=="AndV"), family = "binomial", select = TRUE)
summary(gamEidNoSexAndV.mid) #not sig
plot(gamEidNoSexAndV.mid, pages = 1)

##new dataset
# GhV
# prev.sum.new$virus <- as.factor(prev.sum.new$virus)
# gamEidNoSexGhV.new <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.new, bat_species=="Eidolon dupreanum"  & virus=="GhV"), family = "binomial", select = TRUE)
# summary(gamEidNoSexGhV.new) #
# plot(gamEidNoSexGhV.new, pages = 1)

# AngV
# gamEidNoSexAngV.new <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.new, bat_species=="Eidolon dupreanum"  & virus=="AngV"), family = "binomial", select = TRUE)
# summary(gamEidNoSexAngV.new) #
# plot(gamEidNoSexAngV.new, pages = 1)
# 
# # GhV-like
# gamEidNoSexGhVlike.new <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.new, bat_species=="Eidolon dupreanum"  & virus=="GhV-like"), family = "binomial", select = TRUE)
# summary(gamEidNoSexGhVlike.new) #
# plot(gamEidNoSexGhVlike.new, pages = 1)
# 
# # AndV
# gamEidNoSexAndV.new <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.new, bat_species=="Eidolon dupreanum"  & virus=="AndV"), family = "binomial", select = TRUE)
# summary(gamEidNoSexAndV.new) #
# plot(gamEidNoSexAndV.new, pages = 1)



# Now, model by sex

##old dataset
#GhV
gamEidSexGhV.old <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=-1, bs="cr", by=bat_sex), data=subset(prev.sum.sex.old, bat_species=="Eidolon dupreanum" & virus=="GhV"), family = "binomial", select = TRUE)
plot(gamEidSexGhV.old, pages = 1)
summary(gamEidSexGhV.old) #female sig

# compare model with and without sex
AIC(gamEidNoSexGhV.old, gamEidSexGhV.old) #model with sex better

#AngV
gamEidSexAngV.old <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=-1, bs="cr", by=bat_sex), data=subset(prev.sum.sex.old, bat_species=="Eidolon dupreanum" & virus=="AngV"),family = "binomial", select = TRUE)
plot(gamEidSexAngV.old, pages = 1)
summary(gamEidSexAngV.old) #no sig

# compare model with and without sex
AIC(gamEidNoSexAngV.old, gamEidSexAngV.old) # model with no sex better

#AndV
gamEidSexAndV.old <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=-1, bs="cr", by=bat_sex), data=subset(prev.sum.sex.old, bat_species=="Eidolon dupreanum"  & virus=="AndV"),family = "binomial", select = TRUE)
plot(gamEidSexAndV.old, pages = 1)
summary(gamEidSexAndV.old) #neither sex sig

# compare model with and without sex
AIC(gamEidNoSexAndV.old, gamEidSexAndV.old) #model with no sex better

#GhV like
gamEidSexGhVlike.old <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=-1, bs="cr", by=bat_sex), data=subset(prev.sum.sex.old, bat_species=="Eidolon dupreanum"  & virus=="GhV-like"),family = "binomial", select = TRUE)
plot(gamEidSexGhVlike.old, pages = 1)
summary(gamEidSexGhVlike.old) #female almost sig

# compare model with and without sex
AIC(gamEidNoSexGhVlike.old, gamEidSexGhVlike.old) #model with no sex better


##mid dataset
#GhV
gamEidSexGhV.mid <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=-1, bs="cr", by=bat_sex), data=subset(prev.sum.sex.mid, bat_species=="Eidolon dupreanum"  & virus=="GhV"), family = "binomial", select = TRUE)
plot(gamEidSexGhV.mid, pages = 1)
summary(gamEidSexGhV.mid) #males sig

# compare model with and without sex
AIC(gamEidNoSexGhV.mid, gamEidSexGhV.mid) # model with no sex better

#AngV
gamEidSexAngV.mid <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=-1, bs="cr", by=bat_sex), data=subset(prev.sum.sex.mid, bat_species=="Eidolon dupreanum"  & virus=="AngV"),family = "binomial", select = TRUE)
plot(gamEidSexAngV.mid, pages = 1)
summary(gamEidSexAngV.mid) #neither sex sig

# compare model with and without sex
AIC(gamEidNoSexAngV.mid, gamEidSexAngV.mid) # model with no sex better

#AndV
gamEidSexAndV.mid <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=-1, bs="cr", by=bat_sex), data=subset(prev.sum.sex.mid, bat_species=="Eidolon dupreanum"  & virus=="AndV"),family = "binomial", select = TRUE)
plot(gamEidSexAndV.mid, pages = 1)
summary(gamEidSexAndV.mid) #neither sex sig

# compare model with and without sex
AIC(gamEidNoSexAndV.mid, gamEidSexAndV.mid) # model with no sex better

#GhV like
gamEidSexGhVlike.mid <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=-1, bs="cr", by=bat_sex), data=subset(prev.sum.sex.mid, bat_species=="Eidolon dupreanum"  & virus=="GhV-like"),family = "binomial", select = TRUE)
plot(gamEidSexGhVlike.mid, pages = 1)
summary(gamEidSexGhVlike.mid) #neither sex sig

# compare model with and without sex
AIC(gamEidNoSexGhVlike.mid, gamEidSexGhVlike.mid) # model with no sex better



#use above to predict and overlay on the data
#make dataset to capture predictions


##old dataset first
pred.df.old.ed <- cbind.data.frame(mid_date = rep(seq(as.Date("2013-11-25"), as.Date("2016-01-19"), by = 1), 
                                               (1*length(unique(subset(prev.sum.sex.old, bat_species=="Eidolon dupreanum")$virus)))), 
                                bat_sex=rep(c("female", "male"), each=length(seq(as.Date("2013-11-25"), as.Date("2016-01-19"), by = 1))*length(unique(subset(prev.sum.sex.old, bat_species=="Eidolon dupreanum")$virus))), 
                                virus= rep(rep(unique(subset(prev.sum.sex.old, bat_species=="Eidolon dupreanum")$virus), each = length(seq(as.Date("2013-11-25"), as.Date("2016-01-19"), by = 1))),1))
head(pred.df.old.ed)

#GhV
pred.df.old.ed$predicted_prevalence <- NA
pred.df.old.ed$predicted_prevalence_lci <- NA
pred.df.old.ed$predicted_prevalence_uci <- NA
pred.df.old.ed$predicted_prevalence[pred.df.old.ed$virus=="GhV"] <- predict.gam(gamEidSexGhV.old, newdata=subset(pred.df.old.ed,virus=="GhV"), type = "response", se.fit = T)$fit
pred.df.old.ed$predicted_prevalence_lci[pred.df.old.ed$virus=="GhV"] <-pred.df.old.ed$predicted_prevalence[pred.df.old.ed$virus=="GhV"] - 1.96*predict.gam(gamEidSexGhV.old, newdata=subset(pred.df.old.ed,virus=="GhV"), type = "response", se.fit = T)$se
pred.df.old.ed$predicted_prevalence_uci[pred.df.old.ed$virus=="GhV"] <-pred.df.old.ed$predicted_prevalence[pred.df.old.ed$virus=="GhV"] + 1.96*predict.gam(gamEidSexGhV.old, newdata=subset(pred.df.old.ed,virus=="GhV"), type = "response", se.fit = T)$se


#AngV
pred.df.old.ed$predicted_prevalence[pred.df.old.ed$virus=="AngV"] <- predict.gam(gamEidSexAngV.old, newdata=subset(pred.df.old.ed,virus=="AngV"), type = "response", se.fit = T)$fit
pred.df.old.ed$predicted_prevalence_lci[pred.df.old.ed$virus=="AngV"] <-pred.df.old.ed$predicted_prevalence[pred.df.old.ed$virus=="AngV"] - 1.96*predict.gam(gamEidSexAngV.old, newdata=subset(pred.df.old.ed,virus=="AngV"), type = "response", se.fit = T)$se
pred.df.old.ed$predicted_prevalence_uci[pred.df.old.ed$virus=="AngV"] <-pred.df.old.ed$predicted_prevalence[pred.df.old.ed$virus=="AngV"] + 1.96*predict.gam(gamEidSexAngV.old, newdata=subset(pred.df.old.ed,virus=="AngV"), type = "response", se.fit = T)$se

#GhVlike
pred.df.old.ed$predicted_prevalence[pred.df.old.ed$virus=="GhV-like"] <- predict.gam(gamEidSexGhVlike.old, newdata=subset(pred.df.old.ed,virus=="GhV-like"), type = "response", se.fit = T)$fit
pred.df.old.ed$predicted_prevalence_lci[pred.df.old.ed$virus=="GhV-like"] <-pred.df.old.ed$predicted_prevalence[pred.df.old.ed$virus=="GhV-like"] - 1.96*predict.gam(gamEidSexGhVlike.old, newdata=subset(pred.df.old.ed,virus=="GhV-like"), type = "response", se.fit = T)$se
pred.df.old.ed$predicted_prevalence_uci[pred.df.old.ed$virus=="GhV-like"] <-pred.df.old.ed$predicted_prevalence[pred.df.old.ed$virus=="GhV-like"] + 1.96*predict.gam(gamEidSexGhVlike.old, newdata=subset(pred.df.old.ed,virus=="GhV-like"), type = "response", se.fit = T)$se

#AndV
pred.df.old.ed$predicted_prevalence[pred.df.old.ed$virus=="AndV"] <- predict.gam(gamEidSexAndV.old, newdata=subset(pred.df.old.ed,virus=="AndV"), type = "response", se.fit = T)$fit
pred.df.old.ed$predicted_prevalence_lci[pred.df.old.ed$virus=="AndV"] <-pred.df.old.ed$predicted_prevalence[pred.df.old.ed$virus=="AndV"] - 1.96*predict.gam(gamEidSexAndV.old, newdata=subset(pred.df.old.ed,virus=="AndV"), type = "response", se.fit = T)$se
pred.df.old.ed$predicted_prevalence_uci[pred.df.old.ed$virus=="AndV"] <-pred.df.old.ed$predicted_prevalence[pred.df.old.ed$virus=="AndV"] + 1.96*predict.gam(gamEidSexAndV.old, newdata=subset(pred.df.old.ed,virus=="AndV"), type = "response", se.fit = T)$se

pred.df.old.ed$bat_species = "Eidolon dupreanum"
pred.df.old.ed$dataset = "2013-2016"
pred.df.old.ed$predicted_prevalence_lci[pred.df.old.ed$predicted_prevalence_lci <0] <- 0
pred.df.old.ed$predicted_prevalence_uci[pred.df.old.ed$predicted_prevalence_uci >1] <- 1


##mid dataset
pred.df.mid.ed <- cbind.data.frame(mid_date = rep(seq(as.Date("2018-02-14"), as.Date("2020-03-12"), by = 1), 
                                               (1*length(unique(subset(prev.sum.sex.mid, bat_species=="Eidolon dupreanum")$virus)))), 
                                bat_sex=rep(c("female", "male"), each=length(seq(as.Date("2018-02-14"), as.Date("2020-03-12"), by = 1))*length(unique(subset(prev.sum.sex.mid, bat_species=="Eidolon dupreanum")$virus))), 
                                virus= rep(rep(unique(subset(prev.sum.sex.mid, bat_species=="Eidolon dupreanum")$virus), each = length(seq(as.Date("2018-02-14"), as.Date("2020-03-12"), by = 1))),1))
head(pred.df.mid.ed)

#GhV
pred.df.mid.ed$predicted_prevalence <- NA
pred.df.mid.ed$predicted_prevalence_lci <- NA
pred.df.mid.ed$predicted_prevalence_uci <- NA
pred.df.mid.ed$predicted_prevalence[pred.df.mid.ed$virus=="GhV"] <- predict.gam(gamEidSexGhV.mid, newdata=subset(pred.df.mid.ed,virus=="GhV"), type = "response", se.fit = T)$fit
pred.df.mid.ed$predicted_prevalence_lci[pred.df.mid.ed$virus=="GhV"] <-pred.df.mid.ed$predicted_prevalence[pred.df.mid.ed$virus=="GhV"] - 1.96*predict.gam(gamEidSexGhV.mid, newdata=subset(pred.df.mid.ed,virus=="GhV"), type = "response", se.fit = T)$se
pred.df.mid.ed$predicted_prevalence_uci[pred.df.mid.ed$virus=="GhV"] <-pred.df.mid.ed$predicted_prevalence[pred.df.mid.ed$virus=="GhV"] + 1.96*predict.gam(gamEidSexGhV.mid, newdata=subset(pred.df.mid.ed,virus=="GhV"), type = "response", se.fit = T)$se


#AngV
pred.df.mid.ed$predicted_prevalence[pred.df.mid.ed$virus=="AngV"] <- predict.gam(gamEidSexAngV.mid, newdata=subset(pred.df.mid.ed,virus=="AngV"), type = "response", se.fit = T)$fit
pred.df.mid.ed$predicted_prevalence_lci[pred.df.mid.ed$virus=="AngV"] <-pred.df.mid.ed$predicted_prevalence[pred.df.mid.ed$virus=="AngV"] - 1.96*predict.gam(gamEidSexAngV.mid, newdata=subset(pred.df.mid.ed,virus=="AngV"), type = "response", se.fit = T)$se
pred.df.mid.ed$predicted_prevalence_uci[pred.df.mid.ed$virus=="AngV"] <-pred.df.mid.ed$predicted_prevalence[pred.df.mid.ed$virus=="AngV"] + 1.96*predict.gam(gamEidSexAngV.mid, newdata=subset(pred.df.mid.ed,virus=="AngV"), type = "response", se.fit = T)$se

#GhVlike
pred.df.mid.ed$predicted_prevalence[pred.df.mid.ed$virus=="GhV-like"] <- predict.gam(gamEidSexGhVlike.mid, newdata=subset(pred.df.mid.ed,virus=="GhV-like"), type = "response", se.fit = T)$fit
pred.df.mid.ed$predicted_prevalence_lci[pred.df.mid.ed$virus=="GhV-like"] <-pred.df.mid.ed$predicted_prevalence[pred.df.mid.ed$virus=="GhV-like"] - 1.96*predict.gam(gamEidSexGhVlike.mid, newdata=subset(pred.df.mid.ed,virus=="GhV-like"), type = "response", se.fit = T)$se
pred.df.mid.ed$predicted_prevalence_uci[pred.df.mid.ed$virus=="GhV-like"] <-pred.df.mid.ed$predicted_prevalence[pred.df.mid.ed$virus=="GhV-like"] + 1.96*predict.gam(gamEidSexGhVlike.mid, newdata=subset(pred.df.mid.ed,virus=="GhV-like"), type = "response", se.fit = T)$se

#AndV
pred.df.mid.ed$predicted_prevalence[pred.df.mid.ed$virus=="AndV"] <- predict.gam(gamEidSexAndV.mid, newdata=subset(pred.df.mid.ed,virus=="AndV"), type = "response", se.fit = T)$fit
pred.df.mid.ed$predicted_prevalence_lci[pred.df.mid.ed$virus=="AndV"] <-pred.df.mid.ed$predicted_prevalence[pred.df.mid.ed$virus=="AndV"] - 1.96*predict.gam(gamEidSexAndV.mid, newdata=subset(pred.df.mid.ed,virus=="AndV"), type = "response", se.fit = T)$se
pred.df.mid.ed$predicted_prevalence_uci[pred.df.mid.ed$virus=="AndV"] <-pred.df.mid.ed$predicted_prevalence[pred.df.mid.ed$virus=="AndV"] + 1.96*predict.gam(gamEidSexAndV.mid, newdata=subset(pred.df.mid.ed,virus=="AndV"), type = "response", se.fit = T)$se

pred.df.mid.ed$bat_species = "Eidolon dupreanum"
pred.df.mid.ed$dataset = "2018-2020"
pred.df.mid.ed$predicted_prevalence_lci[pred.df.mid.ed$predicted_prevalence_lci <0] <- 0
pred.df.mid.ed$predicted_prevalence_uci[pred.df.mid.ed$predicted_prevalence_uci >1] <- 1


#then, mark those projections that were from significant model fits only for plotting

##old dataset
#Overall: GhV, AngV, AndV, GhV-like sig
#By sex: GhV-female sig, GhVlike female sig, AngV AndV neither sex sig
# pred.df.old.ed$sig <- 0
# pred.df.old.ed$sig[pred.df.old.ed$virus=="GhV"& pred.df.old.ed$bat_sex=="male"] <- 2
# pred.df.old.ed$sig[pred.df.old.ed$virus=="GhV"& pred.df.old.ed$bat_sex=="female"] <- 1
# pred.df.old.ed$sig[pred.df.old.ed$virus=="GhV-like"& pred.df.old.ed$bat_sex=="male"] <- 2
# pred.df.old.ed$sig[pred.df.old.ed$virus=="GhV-like"& pred.df.old.ed$bat_sex=="female"] <- 1
# pred.df.old.ed$sig[pred.df.old.ed$virus=="AngV"] <- 2
# pred.df.old.ed$sig[pred.df.old.ed$virus=="AndV"] <- 2
# 
# 
# ##mid dataset
# #Overall: GhV sig, AngV, GhV-like, AndV not sig
# #By sex: GhV-male sig, AndV, AngV, GhVlike neither sex sig
# pred.df.mid.ed$sig <- 0
# pred.df.mid.ed$sig[pred.df.mid.ed$virus=="GhV"] <- 1
# pred.df.mid.ed$sig[pred.df.mid.ed$virus=="AngV"] <- 1
# pred.df.mid.ed$sig[pred.df.mid.ed$virus=="AndV"] <- 1
# pred.df.mid.ed$sig[pred.df.mid.ed$virus=="GhV-like"] <- 1
# pred.df.mid.ed$sig[pred.df.mid.ed$virus=="GhV-like" &pred.df.mid.ed$bat_sex=="male"] <- 1
# pred.df.mid.ed$sig[pred.df.mid.ed$virus=="GhV-like" &pred.df.mid.ed$bat_sex=="female"] <- 1

#and plot with data:

old.ed <- ggplot(data=subset(prev.sum.sex.old, region=="Moramanga" & bat_species=="Eidolon dupreanum")) +
  geom_rect(data=demo.dat.old.ed, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values=colz1) + scale_fill_manual(values=filz) +
  #geom_line(aes(x=mid_date, y=prev, color=virus), linetype=2) +
  #geom_linerange(aes(x=mid_date, ymin=prev_lci, ymax=prev_uci, color=virus), linetype=2, alpha=.3, show.legend = F) +
  #geom_ribbon(data=subset(pred.df.old.ed, sig==1), aes(x=mid_date, ymin= predicted_prevalence_lci, ymax=predicted_prevalence_uci, fill=virus), alpha=.3) +
  #geom_line(data=subset(pred.df.old.ed, sig==1), aes(x=mid_date, y= predicted_prevalence, color=virus)) +
  geom_ribbon(data=subset(pred.df.old.ed), aes(x=mid_date, ymin= predicted_prevalence_lci, ymax=predicted_prevalence_uci, fill=virus), alpha=.3) +
  geom_line(data=subset(pred.df.old.ed), aes(x=mid_date, y= predicted_prevalence, color=virus)) +
  facet_grid(bat_sex~bat_species) + 
  theme_bw()+
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(), 
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=9), panel.grid = element_blank(), 
        legend.title = element_blank())
old.ed

mid.ed <- ggplot(data=subset(prev.sum.sex.mid, region=="Moramanga" & bat_species=="Eidolon dupreanum")) +
  geom_rect(data=demo.dat.mid.ed, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values=colz1) + scale_fill_manual(values=filz) +
  #geom_line(aes(x=mid_date, y=prev, color=virus), linetype=2) +
  #geom_linerange(aes(x=mid_date, ymin=prev_lci, ymax=prev_uci, color=virus), linetype=2, alpha=.3, show.legend = F) +
  #geom_ribbon(data=subset(pred.df.mid.ed, sig==1), aes(x=mid_date, ymin= predicted_prevalence_lci, ymax=predicted_prevalence_uci, fill=virus), alpha=.3) +
  #geom_line(data=subset(pred.df.mid.ed, sig==1), aes(x=mid_date, y= predicted_prevalence, color=virus)) +
  geom_ribbon(data=subset(pred.df.mid.ed), aes(x=mid_date, ymin= predicted_prevalence_lci, ymax=predicted_prevalence_uci, fill=virus), alpha=.3) +
  geom_line(data=subset(pred.df.mid.ed), aes(x=mid_date, y= predicted_prevalence, color=virus)) +
  facet_grid(bat_species+bat_sex~virus) + 
  scale_x_date(limits = as.Date(c("2018-02-18", "2019-09-20")))+
  theme_bw()+
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(), 
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=9), panel.grid = element_blank(), 
        legend.title = element_blank())
mid.ed


##Now do rousettus
##old dataset
#GhV
prev.sum.old$virus <- as.factor(prev.sum.old$virus)
gamRouNoSexGhV.old <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=5, bs="cr"), data=subset(prev.sum.old, bat_species=="Rousettus madagascariensis"  & virus=="GhV"), family = "binomial", select = TRUE)
summary(gamRouNoSexGhV.old) #no sig
plot(gamRouNoSexGhV.old, pages = 1)

# AngV
gamRouNoSexAngV.old <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=5, bs="cr"), data=subset(prev.sum.old, bat_species=="Rousettus madagascariensis"  & virus=="AngV"), family = "binomial", select = TRUE)
summary(gamRouNoSexAngV.old) #no sig
plot(gamRouNoSexAngV.old, pages = 1)

# GhV-like
gamRouNoSexGhVlike.old <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=5, bs="cr"), data=subset(prev.sum.old, bat_species=="Rousettus madagascariensis"  & virus=="GhV-like"), family = "binomial", select = TRUE)
summary(gamRouNoSexGhVlike.old) #no sig
plot(gamRouNoSexGhVlike.old, pages = 1)
 
# AndV
gamRouNoSexAndV.old <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=5, bs="cr"), data=subset(prev.sum.old, bat_species=="Rousettus madagascariensis"  & virus=="AndV"), family = "binomial", select = TRUE)
summary(gamRouNoSexAndV.old) #no sig
plot(gamRouNoSexAndV.old, pages = 1)

##mid dataset
# GhV
prev.sum.mid$virus <- as.factor(prev.sum.mid$virus)
gamRouNoSexGhV.mid <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.mid, bat_species=="Rousettus madagascariensis"  & virus=="GhV"), family = "binomial", select = TRUE)
summary(gamRouNoSexGhV.mid) #no sig
plot(gamRouNoSexGhV.mid, pages = 1)

# AngV
gamRouNoSexAngV.mid <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.mid, bat_species=="Rousettus madagascariensis"  & virus=="AngV"), family = "binomial", select = TRUE)
summary(gamRouNoSexAngV.mid) #not sig
plot(gamRouNoSexAngV.mid, pages = 1)

# GhV-like
gamRouNoSexGhVlike.mid <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.mid, bat_species=="Rousettus madagascariensis"  & virus=="GhV-like"), family = "binomial", select = TRUE)
summary(gamRouNoSexGhVlike.mid) #not sig
plot(gamRouNoSexGhVlike.mid, pages = 1)

# AndV
gamRouNoSexAndV.mid <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.mid, bat_species=="Rousettus madagascariensis"  & virus=="AndV"), family = "binomial", select = TRUE)
summary(gamRouNoSexAndV.mid) #not sig
plot(gamRouNoSexAndV.mid, pages = 1)

##new dataset
# GhV
# prev.sum.new$virus <- as.factor(prev.sum.new$virus)
# gamRouNoSexGhV.new <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.new, bat_species=="Rousettus madagascariensis"  & virus=="GhV"), family = "binomial", select = TRUE)
# summary(gamRouNoSexGhV.new) #
# plot(gamRouNoSexGhV.new, pages = 1)

# AngV
# gamRouNoSexAngV.new <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.new, bat_species=="Rousettus madagascariensis"  & virus=="AngV"), family = "binomial", select = TRUE)
# summary(gamRouNoSexAngV.new) #
# plot(gamRouNoSexAngV.new, pages = 1)
# 
# # GhV-like
# gamRouNoSexGhVlike.new <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.new, bat_species=="Rousettus madagascariensis"  & virus=="GhV-like"), family = "binomial", select = TRUE)
# summary(gamRouNoSexGhVlike.new) #
# plot(gamRouNoSexGhVlike.new, pages = 1)
# 
# # AndV
# gamRouNoSexAndV.new <- gam(cbind(pos, neg)~s(as.numeric(mid_date), k=-1, bs="cr"), data=subset(prev.sum.new, bat_species=="Rousettus madagascariensis"  & virus=="AndV"), family = "binomial", select = TRUE)
# summary(gamRouNoSexAndV.new) #
# plot(gamRouNoSexAndV.new, pages = 1)



# Now, model by sex

##old dataset
#GhV
gamRouSexGhV.old <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=5, bs="cr", by=bat_sex), data=subset(prev.sum.sex.old, bat_species=="Rousettus madagascariensis"  & virus=="GhV"), family = "binomial", select = TRUE)
plot(gamRouSexGhV.old, pages = 1)
summary(gamRouSexGhV.old) #no sig

# compare model with and without sex
AIC(gamRouNoSexGhV.old, gamRouSexGhV.old) #model with no sex better

#AngV
gamRouSexAngV.old <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=5, bs="cr", by=bat_sex), data=subset(prev.sum.sex.old, bat_species=="Rousettus madagascariensis"  & virus=="AngV"),family = "binomial", select = TRUE)
plot(gamRouSexAngV.old, pages = 1)
summary(gamRouSexAngV.old) #neither sex sig

# compare model with and without sex
AIC(gamRouNoSexAngV.old, gamRouSexAngV.old) # model with no sex better

#AndV
gamRouSexAndV.old <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=5, bs="cr", by=bat_sex), data=subset(prev.sum.sex.old, bat_species=="Rousettus madagascariensis"  & virus=="AndV"),family = "binomial", select = TRUE)
plot(gamRouSexAndV.old, pages = 1)
summary(gamRouSexAndV.old) #neither sex sig

# compare model with and without sex
AIC(gamRouNoSexAndV.old, gamRouSexAndV.old) #model with no sex better

#GhV like
gamRouSexGhVlike.old <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=5, bs="cr", by=bat_sex), data=subset(prev.sum.sex.old, bat_species=="Rousettus madagascariensis"  & virus=="GhV-like"),family = "binomial", select = TRUE)
plot(gamRouSexGhVlike.old, pages = 1)
summary(gamRouSexGhVlike.old) #neither sex sig

# compare model with and without sex
AIC(gamRouNoSexGhVlike.old, gamRouSexGhVlike.old) #model with no sex better


##mid dataset
#GhV
gamRouSexGhV.mid <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=-1, bs="cr", by=bat_sex), data=subset(prev.sum.sex.mid, bat_species=="Rousettus madagascariensis"  & virus=="GhV"), family = "binomial", select = TRUE)
plot(gamRouSexGhV.mid, pages = 1)
summary(gamRouSexGhV.mid) #no sig

# compare model with and without sex
AIC(gamRouNoSexGhV.mid, gamRouSexGhV.mid) # model with no sex better

#AngV
gamRouSexAngV.mid <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=-1, bs="cr", by=bat_sex), data=subset(prev.sum.sex.mid, bat_species=="Rousettus madagascariensis"  & virus=="AngV"),family = "binomial", select = TRUE)
plot(gamRouSexAngV.mid, pages = 1)
summary(gamRouSexAngV.mid) #no sig

# compare model with and without sex
AIC(gamRouNoSexAngV.mid, gamRouSexAngV.mid) # model with no sex better

#AndV
gamRouSexAndV.mid <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=-1, bs="cr", by=bat_sex), data=subset(prev.sum.sex.mid, bat_species=="Rousettus madagascariensis"  & virus=="AndV"),family = "binomial", select = TRUE)
plot(gamRouSexAndV.mid, pages = 1)
summary(gamRouSexAndV.mid) #neither sex sig

# compare model with and without sex
AIC(gamRouNoSexAndV.mid, gamRouSexAndV.mid) # model with no sex better

#GhV like
gamRouSexGhVlike.mid <- gam(cbind(pos, neg)~ bat_sex + s(as.numeric(mid_date), k=-1, bs="cr", by=bat_sex), data=subset(prev.sum.sex.mid, bat_species=="Rousettus madagascariensis"  & virus=="GhV-like"),family = "binomial", select = TRUE)
plot(gamRouSexGhVlike.mid, pages = 1)
summary(gamRouSexGhVlike.mid) #neither sex sig

# compare model with and without sex
AIC(gamRouNoSexGhVlike.mid, gamRouSexGhVlike.mid) # model with no sex better



#use above to predict and overlay on the data
#make dataset to capture predictions


##old dataset first
# pred.df.old.rm <- cbind.data.frame(mid_date = rep(seq(as.Date("2013-11-25"), as.Date("2016-01-19"), by = 1), 
#                                                (2*length(unique(subset(prev.sum.sex.old, bat_species=="Rousettus madagascariensis")$virus)))), 
#                                 bat_sex=rep(c("female", "male"), each=length(seq(as.Date("2013-11-25"), as.Date("2016-01-19"), by = 1))*length(unique(subset(prev.sum.sex.old, bat_species=="Rousettus madagascariensis")$virus))), 
#                                 virus= rep(rep(unique(subset(prev.sum.sex.old, bat_species=="Rousettus madagascariensis")$virus), each = length(seq(as.Date("2013-11-25"), as.Date("2016-01-19"), by = 1))),2))
# head(pred.df.old.rm)
# 
# #GhV
# pred.df.old.rm$predicted_prevalence <- NA
# pred.df.old.rm$predicted_prevalence_lci <- NA
# pred.df.old.rm$predicted_prevalence_uci <- NA
# pred.df.old.rm$predicted_prevalence[pred.df.old.rm$virus=="GhV"] <- predict.gam(gamRouSexGhV.old, newdata=subset(pred.df.old.rm,virus=="GhV"), type = "response", se.fit = T)$fit
# pred.df.old.rm$predicted_prevalence_lci[pred.df.old.rm$virus=="GhV"] <-pred.df.old.rm$predicted_prevalence[pred.df.old.rm$virus=="GhV"] - 1.96*predict.gam(gamRouSexGhV.old, newdata=subset(pred.df.old.rm,virus=="GhV"), type = "response", se.fit = T)$se
# pred.df.old.rm$predicted_prevalence_uci[pred.df.old.rm$virus=="GhV"] <-pred.df.old.rm$predicted_prevalence[pred.df.old.rm$virus=="GhV"] + 1.96*predict.gam(gamRouSexGhV.old, newdata=subset(pred.df.old.rm,virus=="GhV"), type = "response", se.fit = T)$se
# 
# 
# #AngV
# pred.df.old.rm$predicted_prevalence[pred.df.old.rm$virus=="AngV"] <- predict.gam(gamRouSexAngV.old, newdata=subset(pred.df.old.rm,virus=="AngV"), type = "response", se.fit = T)$fit
# pred.df.old.rm$predicted_prevalence_lci[pred.df.old.rm$virus=="AngV"] <-pred.df.old.rm$predicted_prevalence[pred.df.old.rm$virus=="AngV"] - 1.96*predict.gam(gamRouSexAngV.old, newdata=subset(pred.df.old.rm,virus=="AngV"), type = "response", se.fit = T)$se
# pred.df.old.rm$predicted_prevalence_uci[pred.df.old.rm$virus=="AngV"] <-pred.df.old.rm$predicted_prevalence[pred.df.old.rm$virus=="AngV"] + 1.96*predict.gam(gamRouSexAngV.old, newdata=subset(pred.df.old.rm,virus=="AngV"), type = "response", se.fit = T)$se
# 
# #GhVlike
# pred.df.old.rm$predicted_prevalence[pred.df.old.rm$virus=="GhV-like"] <- predict.gam(gamRouSexGhVlike.old, newdata=subset(pred.df.old.rm,virus=="GhV-like"), type = "response", se.fit = T)$fit
# pred.df.old.rm$predicted_prevalence_lci[pred.df.old.rm$virus=="GhV-like"] <-pred.df.old.rm$predicted_prevalence[pred.df.old.rm$virus=="GhV-like"] - 1.96*predict.gam(gamRouSexGhVlike.old, newdata=subset(pred.df.old.rm,virus=="GhV-like"), type = "response", se.fit = T)$se
# pred.df.old.rm$predicted_prevalence_uci[pred.df.old.rm$virus=="GhV-like"] <-pred.df.old.rm$predicted_prevalence[pred.df.old.rm$virus=="GhV-like"] + 1.96*predict.gam(gamRouSexGhVlike.old, newdata=subset(pred.df.old.rm,virus=="GhV-like"), type = "response", se.fit = T)$se
# 
# #AndV
# pred.df.old.rm$predicted_prevalence[pred.df.old.rm$virus=="AndV"] <- predict.gam(gamRouSexAndV.old, newdata=subset(pred.df.old.rm,virus=="AndV"), type = "response", se.fit = T)$fit
# pred.df.old.rm$predicted_prevalence_lci[pred.df.old.rm$virus=="AndV"] <-pred.df.old.rm$predicted_prevalence[pred.df.old.rm$virus=="AndV"] - 1.96*predict.gam(gamRouSexAndV.old, newdata=subset(pred.df.old.rm,virus=="AndV"), type = "response", se.fit = T)$se
# pred.df.old.rm$predicted_prevalence_uci[pred.df.old.rm$virus=="AndV"] <-pred.df.old.rm$predicted_prevalence[pred.df.old.rm$virus=="AndV"] + 1.96*predict.gam(gamRouSexAndV.old, newdata=subset(pred.df.old.rm,virus=="AndV"), type = "response", se.fit = T)$se
# 
# pred.df.old.rm$bat_species = "Rousettus madagascariensis"
# pred.df.old.rm$dataset = "2013-2016"
# pred.df.old.rm$predicted_prevalence_lci[pred.df.old.rm$predicted_prevalence_lci <0] <- 0
# pred.df.old.rm$predicted_prevalence_uci[pred.df.old.rm$predicted_prevalence_uci >1] <- 1


##mid dataset
pred.df.mid.rm <- cbind.data.frame(mid_date = rep(seq(as.Date("2018-02-14"), as.Date("2020-03-12"), by = 1), 
                                               (1*length(unique(subset(prev.sum.sex.mid, bat_species=="Rousettus madagascariensis")$virus)))), 
                                bat_sex=rep(c("female", "male"), each=length(seq(as.Date("2018-02-14"), as.Date("2020-03-12"), by = 1))*length(unique(subset(prev.sum.sex.mid, bat_species=="Rousettus madagascariensis")$virus))), 
                                virus= rep(rep(unique(subset(prev.sum.sex.mid, bat_species=="Rousettus madagascariensis")$virus), each = length(seq(as.Date("2018-02-14"), as.Date("2020-03-12"), by = 1))),1))
head(pred.df.mid.rm)

#GhV
pred.df.mid.rm$predicted_prevalence <- NA
pred.df.mid.rm$predicted_prevalence_lci <- NA
pred.df.mid.rm$predicted_prevalence_uci <- NA
pred.df.mid.rm$predicted_prevalence[pred.df.mid.rm$virus=="GhV"] <- predict.gam(gamRouSexGhV.mid, newdata=subset(pred.df.mid.rm,virus=="GhV"), type = "response", se.fit = T)$fit
pred.df.mid.rm$predicted_prevalence_lci[pred.df.mid.rm$virus=="GhV"] <-pred.df.mid.rm$predicted_prevalence[pred.df.mid.rm$virus=="GhV"] - 1.96*predict.gam(gamRouSexGhV.mid, newdata=subset(pred.df.mid.rm,virus=="GhV"), type = "response", se.fit = T)$se
pred.df.mid.rm$predicted_prevalence_uci[pred.df.mid.rm$virus=="GhV"] <-pred.df.mid.rm$predicted_prevalence[pred.df.mid.rm$virus=="GhV"] + 1.96*predict.gam(gamRouSexGhV.mid, newdata=subset(pred.df.mid.rm,virus=="GhV"), type = "response", se.fit = T)$se


#AngV
pred.df.mid.rm$predicted_prevalence[pred.df.mid.rm$virus=="AngV"] <- predict.gam(gamRouSexAngV.mid, newdata=subset(pred.df.mid.rm,virus=="AngV"), type = "response", se.fit = T)$fit
pred.df.mid.rm$predicted_prevalence_lci[pred.df.mid.rm$virus=="AngV"] <-pred.df.mid.rm$predicted_prevalence[pred.df.mid.rm$virus=="AngV"] - 1.96*predict.gam(gamRouSexAngV.mid, newdata=subset(pred.df.mid.rm,virus=="AngV"), type = "response", se.fit = T)$se
pred.df.mid.rm$predicted_prevalence_uci[pred.df.mid.rm$virus=="AngV"] <-pred.df.mid.rm$predicted_prevalence[pred.df.mid.rm$virus=="AngV"] + 1.96*predict.gam(gamRouSexAngV.mid, newdata=subset(pred.df.mid.rm,virus=="AngV"), type = "response", se.fit = T)$se

#GhVlike
pred.df.mid.rm$predicted_prevalence[pred.df.mid.rm$virus=="GhV-like"] <- predict.gam(gamRouSexGhVlike.mid, newdata=subset(pred.df.mid.rm,virus=="GhV-like"), type = "response", se.fit = T)$fit
pred.df.mid.rm$predicted_prevalence_lci[pred.df.mid.rm$virus=="GhV-like"] <-pred.df.mid.rm$predicted_prevalence[pred.df.mid.rm$virus=="GhV-like"] - 1.96*predict.gam(gamRouSexGhVlike.mid, newdata=subset(pred.df.mid.rm,virus=="GhV-like"), type = "response", se.fit = T)$se
pred.df.mid.rm$predicted_prevalence_uci[pred.df.mid.rm$virus=="GhV-like"] <-pred.df.mid.rm$predicted_prevalence[pred.df.mid.rm$virus=="GhV-like"] + 1.96*predict.gam(gamRouSexGhVlike.mid, newdata=subset(pred.df.mid.rm,virus=="GhV-like"), type = "response", se.fit = T)$se

#AndV
pred.df.mid.rm$predicted_prevalence[pred.df.mid.rm$virus=="AndV"] <- predict.gam(gamRouSexAndV.mid, newdata=subset(pred.df.mid.rm,virus=="AndV"), type = "response", se.fit = T)$fit
pred.df.mid.rm$predicted_prevalence_lci[pred.df.mid.rm$virus=="AndV"] <-pred.df.mid.rm$predicted_prevalence[pred.df.mid.rm$virus=="AndV"] - 1.96*predict.gam(gamRouSexAndV.mid, newdata=subset(pred.df.mid.rm,virus=="AndV"), type = "response", se.fit = T)$se
pred.df.mid.rm$predicted_prevalence_uci[pred.df.mid.rm$virus=="AndV"] <-pred.df.mid.rm$predicted_prevalence[pred.df.mid.rm$virus=="AndV"] + 1.96*predict.gam(gamRouSexAndV.mid, newdata=subset(pred.df.mid.rm,virus=="AndV"), type = "response", se.fit = T)$se

pred.df.mid.rm$bat_species = "Rousettus madagascariensis"
pred.df.mid.rm$dataset = "2018-2020"
pred.df.mid.rm$predicted_prevalence_lci[pred.df.mid.rm$predicted_prevalence_lci <0] <- 0
pred.df.mid.rm$predicted_prevalence_uci[pred.df.mid.rm$predicted_prevalence_uci >1] <- 1


#then, mark those projections that were from significant model fits only for plotting

##old dataset
#Overall: GhV, AngV, GhV-like sig, AndV not sig
#By sex: GhV-female sig, AndV, AngV, GhVlike neither sex sig
# pred.df.old.rm$sig <- 0
# pred.df.old.rm$sig[pred.df.old.rm$virus=="GhV"& pred.df.old.rm$bat_sex=="female"] <- 1
# pred.df.old.rm$sig[pred.df.old.rm$virus=="GhV"& pred.df.old.rm$bat_sex=="male"] <- 2
# pred.df.old.rm$sig[pred.df.old.rm$virus=="AngV"] <- 2
# pred.df.old.rm$sig[pred.df.old.rm$virus=="AndV"] <- 2
# pred.df.old.rm$sig[pred.df.old.rm$virus=="GhV-like"] <- 2


##mid dataset
#Overall: GhV sig, AngV, GhV-like, AndV not sig
#By sex: GhV-female almost sig, AndV, AngV, GhVlike neither sex sig
# pred.df.mid.rm$sig <- 0
# pred.df.mid.rm$sig[pred.df.mid.rm$virus=="GhV"] <- 1
# pred.df.mid.rm$sig[pred.df.mid.rm$virus=="AngV"] <- 1
# pred.df.mid.rm$sig[pred.df.mid.rm$virus=="AndV"] <- 1
# pred.df.mid.rm$sig[pred.df.mid.rm$virus=="GhV-like"] <- 1

#and plot with data:

#filz <- c(colz1, 'gestation' = "yellow", 'lactation' = "greenyellow", 'winter' = "skyblue")

# old.rm <- ggplot(data=subset(prev.sum.old, region=="Moramanga" & bat_species=="Rousettus madagascariensis")) +
#   geom_rect(data=demo.dat.old.ed, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
#   geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values=colz1) + scale_fill_manual(values=filz) +
#   geom_line(aes(x=mid_date, y=prev, color=virus), linetype=2) +
#   geom_linerange(aes(x=mid_date, ymin=prev_lci, ymax=prev_uci, color=virus), linetype=2, alpha=.3, show.legend = F) +
#   geom_ribbon(data=subset(pred.df.old, sig==1), aes(x=mid_date, ymin= predicted_prevalence_lci, ymax=predicted_prevalence_uci, fill=virus), alpha=.3) +
#   geom_line(data=subset(pred.df.old, sig==1), aes(x=mid_date, y= predicted_prevalence, color=virus)) +
#   facet_grid(bat_sex~bat_species) + 
#   theme_bw()
# old.rm

mid.rm <- ggplot(data=subset(prev.sum.sex.mid, region=="Moramanga" & bat_species=="Rousettus madagascariensis")) +
  geom_rect(data=demo.dat.mid.rm, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values=colz1) + scale_fill_manual(values=filz) +
  #geom_line(aes(x=mid_date, y=prev, color=virus), linetype=2) +
  geom_linerange(aes(x=mid_date, ymin=prev_lci, ymax=prev_uci, color=virus), linetype=2, alpha=.3, show.legend = F) +
  #geom_ribbon(data=subset(pred.df.mid, sig==1), aes(x=mid_date, ymin= predicted_prevalence_lci, ymax=predicted_prevalence_uci, fill=virus), alpha=.3) +
  #geom_line(data=subset(pred.df.mid, sig==1), aes(x=mid_date, y= predicted_prevalence, color=virus)) +
  geom_ribbon(data=subset(pred.df.mid.rm), aes(x=mid_date, ymin= predicted_prevalence_lci, ymax=predicted_prevalence_uci, fill=virus), alpha=.3) +
  geom_line(data=subset(pred.df.mid.rm), aes(x=mid_date, y= predicted_prevalence, color=virus)) +
  facet_grid(bat_sex~bat_species) + 
  theme_bw()+
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(), 
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=9), panel.grid = element_blank(), 
        legend.title = element_blank())
mid.rm ##one plot for both e dupreanum and r madagascariensis

pred.df.mid.all <- rbind(pred.df.mid.ed, pred.df.mid.rm)
demo.dat.mid.ed.rm <- rbind(demo.dat.mid.ed, demo.dat.mid.rm)

#both together
prev.sum.sex.mid.plot<-prev.sum.sex.mid[prev.sum.sex.mid$bat_species %in% c("Eidolon dupreanum", "Rousettus madagascariensis"),]

mid <- ggplot(data=subset(prev.sum.sex.mid.plot, region=="Moramanga")) +
  geom_rect(data=demo.dat.mid, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values=colz1) + scale_fill_manual(values=filz) +
  #geom_line(aes(x=mid_date, y=prev, color=virus), linetype=2) +
  geom_linerange(aes(x=mid_date, ymin=prev_lci, ymax=prev_uci, color=virus), linetype=2, alpha=.3, show.legend = F) +
  geom_ribbon(data=subset(pred.df.mid.all), aes(x=mid_date, ymin= predicted_prevalence_lci, ymax=predicted_prevalence_uci, fill=virus), alpha=.3) +
  geom_line(data=subset(pred.df.mid.all), aes(x=mid_date, y= predicted_prevalence, color=virus)) +
  facet_grid(bat_species~bat_sex+virus) + 
  scale_x_date(limits = as.Date(c("2018-02-18","2019-09-20")))+
  theme_bw()+
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(), 
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=9), panel.grid = element_blank(), 
        legend.title = element_blank())
mid



####################################################################################################
####################################################################################################
## SECTION 3: Correlates of positivity by infection status and nutrition
## nutrition bit can only be computed for adults 

## what are the "risk factors" for being positive ?
## look at individuals this time instead of prevalence

# # remove bats from the dataset with no age class listed
molecular.dat = subset(merge.dat.pos, !is.na(bat_age_class))

#factor the categorical data
molecular.dat$year <- as.factor(molecular.dat$year)
molecular.dat$bat_sex <- as.factor(molecular.dat$bat_sex)
molecular.dat$roost_site <- as.factor(molecular.dat$roost_site)
molecular.dat$region <- as.factor(molecular.dat$region)
molecular.dat$virus <- as.factor(molecular.dat$virus)
molecular.dat$bat_species <- as.factor(molecular.dat$bat_species)
molecular.dat$bat_age_class <- as.factor(molecular.dat$bat_age_class)
molecular.dat$young_of_year <- as.factor(molecular.dat$young_of_year)
molecular.dat$pupmom_pair <- as.factor(molecular.dat$pupmom_pair)
molecular.dat$region <- as.factor(molecular.dat$region)
molecular.dat$combo_class<-as.factor(molecular.dat$combo_class)
molecular.dat$month<-as.numeric(molecular.dat$month)

##New seasonality timing
#day
lact.dat.day <- cbind.data.frame(bat_species=rep(c("Eidolon dupreanum", "Rousettus madagascariensis"), each=2), 
                                 start = c(320,1, 346,1), 
                                 stop=c(365,33,365,50))
lact.dat.day$type = "lactation"
lact.dat.day$bat_sex="female"
gest.dat.day <- cbind.data.frame(bat_species=c("Eidolon dupreanum","Rousettus madagascariensis"), 
                                 start = c(215, 254),
                                 stop=c(320, 345))
gest.dat.day$type = "gestation"
gest.dat.day$bat_sex = "female"
#and winter
cal.dat.day <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis"), 
                                 start = c(135,135),
                                 stop=c(283,283))
cal.dat.day$type <- "winter"
cal.dat.day$bat_sex="male"

library(reshape2)

repro.dat.day <- rbind(lact.dat.day, gest.dat.day)
demo.dat.day <- rbind(gest.dat.day, lact.dat.day, cal.dat.day)

#month
lact.dat.month <- cbind.data.frame(bat_species=rep(c("Eidolon dupreanum", "Rousettus madagascariensis"), each=2), 
                                   start = c(11.5,1, 11.4,1), 
                                   stop=c(12,2.05,12,2.7))
lact.dat.month$type = "lactation"
lact.dat.month$bat_sex="female"

gest.dat.month <- cbind.data.frame(bat_species=c("Eidolon dupreanum","Rousettus madagascariensis"), 
                                   start = c(8.15, 9.45),
                                   stop=c(11.5, 11.39))
gest.dat.month$type = "gestation"
gest.dat.month$bat_sex = "female"

#and winter
cal.dat.month <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis"), 
                                   start = c(5.5,5.5),
                                   stop=c(10.3,10.3))
cal.dat.month$type <- "winter"
cal.dat.month$bat_sex="male"

library(reshape2)

repro.dat.month <- rbind(lact.dat.month, gest.dat.month)
demo.dat.month <- rbind(gest.dat.month, lact.dat.month, cal.dat.month)

### first, confirm that seasonality matters

####Eidolon dupreanum 
ed1_molecular <- gam(pos~ 
                         s(month,k=4, bs="cc", by=virus)+
                         s(virus, bs="re") +
                         s(bat_sex, bs="re")+
                         s(year, bs ="re"), 
                       dat = subset(molecular.dat, bat_species=="Eidolon dupreanum"), 
                       family = "binomial")
summary(ed1_molecular) #angv, ghv, virus, sex, and year sig
plot(ed1_molecular, pages=1) 

ed1_molecular.f <- gam(pos~ 
                       s(month,k=4, bs="cc", by=virus)+
                       s(virus, bs="re") +
                       s(year, bs ="re"), 
                     dat = subset(molecular.dat, bat_species=="Eidolon dupreanum" & bat_sex =="female"), 
                     family = "binomial")
summary(ed1_molecular.f) #ghv, virus, year sig
plot(ed1_molecular.f, pages=1) 

ed1_molecular.m <- gam(pos~ 
                         s(month,k=4, bs="cc", by=virus)+
                         s(virus, bs="re") +
                         s(year, bs ="re"), 
                       dat = subset(molecular.dat, bat_species=="Eidolon dupreanum" & bat_sex =="male"), 
                       family = "binomial")
summary(ed1_molecular.m) #angv, ghv, virus, year sig
plot(ed1_molecular.m, pages=1) 


###Rousettus madagascariensis
rm1_molecular <- gam(pos~ 
                       s(month,k=4, bs="cc", by=virus)+
                       s(virus, bs="re") +
                       s(bat_sex, bs="re")+
                       s(year, bs ="re"), 
                     dat = subset(molecular.dat, bat_species=="Rousettus madagascariensis"), 
                     family = "binomial")
summary(rm1_molecular) #angv, ghvlike, year sig, sex almost sig
plot(rm1_molecular, pages=1) 

rm1_molecular.f <- gam(pos~ 
                       s(month,k=4, bs="cc", by=virus)+
                       s(virus, bs="re") +
                       s(year, bs ="re"), 
                     dat = subset(molecular.dat, bat_species=="Rousettus madagascariensis" & bat_sex =="female"), 
                     family = "binomial")
summary(rm1_molecular.f) #ghvlike abd virus sig
plot(rm1_molecular.f, pages=1) 

rm1_molecular.m <- gam(pos~ 
                         s(month,k=4, bs="cc", by=virus)+
                         s(virus, bs="re") +
                         s(year, bs ="re"), 
                       dat = subset(molecular.dat, bat_species=="Rousettus madagascariensis" & bat_sex =="male"), 
                       family = "binomial")
summary(rm1_molecular.m) #angv ghv virus and year sig
plot(rm1_molecular.m, pages=1) 



####Eidolon dupreanum
#predict without effect of year for plotting
pred.df.ED_molecular <- cbind.data.frame(month=rep(seq(1,12, by=.1),8*8), year=rep(c(2013,2014,2015,2016, 2018, 2019, 2020, 2022), each=(length(seq(1,12,by=.1)*8))),
                                         bat_sex = rep(rep(c("female","male"), each=length(seq(1,12,by=.1))),8),
                                         virus = rep(rep(c("AndV", "AngV", "GhV", "GhV-like"), each = length(seq(1,12,by=.1))),8))
#predict without effect of year for plotting
pred.df.ED_molecular$predicted_pos <- predict.gam(ed1_molecular.f, newdata = pred.df.ED_molecular, type="response", exclude = "s(year)")
pred.df.ED_molecular$predicted_pos_lci <- pred.df.ED_molecular$predicted_pos - 1.96*predict.gam(ed1_molecular.f, newdata = pred.df.ED_molecular, type="response", se.fit = T, exclude = "s(year)")$se
pred.df.ED_molecular$predicted_pos_uci <- pred.df.ED_molecular$predicted_pos + 1.96*predict.gam(ed1_molecular.f, newdata = pred.df.ED_molecular, type="response", se.fit = T, exclude = "s(year)")$se
pred.df.ED_molecular$bat_species="Eidolon dupreanum"
#correct any predictions < 0:
pred.df.ED_molecular$predicted_pos_lci[pred.df.ED_molecular$predicted_pos_lci<0] <- 0
#correct any predictions > 1:
pred.df.ED_molecular$predicted_pos_uci[pred.df.ED_molecular$predicted_pos_uci>1] <- 1

pred.df.EDF_molecular <- cbind.data.frame(month=rep(seq(1,12, by=.1),7*4), year=rep(c(2013,2014,2015,2016, 2018, 2019, 2020), each=(length(seq(1,12,by=.1)*4))), 
                                          virus = rep(rep(c("AndV", "AngV", "GhV", "GhV-like"), each = length(seq(1,12,by=.1))),7))
#predict without effect of year for plotting
pred.df.EDF_molecular$predicted_pos <- predict.gam(ed1_molecular.f, newdata = pred.df.EDF_molecular, type="response", exclude = "s(year)")
pred.df.EDF_molecular$predicted_pos_lci <- pred.df.EDF_molecular$predicted_pos - 1.96*predict.gam(ed1_molecular.f, newdata = pred.df.EDF_molecular, type="response", se.fit = T, exclude = "s(year)")$se
pred.df.EDF_molecular$predicted_pos_uci <- pred.df.EDF_molecular$predicted_pos + 1.96*predict.gam(ed1_molecular.f, newdata = pred.df.EDF_molecular, type="response", se.fit = T, exclude = "s(year)")$se
pred.df.EDF_molecular$bat_species="Eidolon dupreanum"
pred.df.EDF_molecular$bat_sex="female"
#correct any predictions < 0:
pred.df.EDF_molecular$predicted_pos_lci[pred.df.EDF_molecular$predicted_pos_lci<0] <- 0
#correct any predictions > 1:
pred.df.EDF_molecular$predicted_pos_uci[pred.df.EDF_molecular$predicted_pos_uci>1] <- 1

pred.df.EDM_molecular <- cbind.data.frame(month=rep(seq(1,12, by=.1),7*4), year=rep(c(2013,2014,2015,2016, 2018, 2019, 2020), each=(length(seq(1,12,by=.1)*4))), 
                                         virus = rep(rep(c("AndV", "AngV", "GhV", "GhV-like"), each = length(seq(1,12,by=.1))),7))

#predict without effect of year for plotting
pred.df.EDM_molecular$predicted_pos <- predict.gam(ed1_molecular.m, newdata = pred.df.EDM_molecular, type="response", exclude = "s(year)")
pred.df.EDM_molecular$predicted_pos_lci <- pred.df.EDM_molecular$predicted_pos - 1.96*predict.gam(ed1_molecular.m, newdata = pred.df.EDM_molecular, type="response", se.fit = T, exclude = "s(year)")$se
pred.df.EDM_molecular$predicted_pos_uci <- pred.df.EDM_molecular$predicted_pos + 1.96*predict.gam(ed1_molecular.m, newdata = pred.df.EDM_molecular, type="response", se.fit = T, exclude = "s(year)")$se
pred.df.EDM_molecular$bat_species="Eidolon dupreanum"
pred.df.EDM_molecular$bat_sex="male"
#correct any predictions < 0:
pred.df.EDM_molecular$predicted_pos_lci[pred.df.EDM_molecular$predicted_pos_lci<0] <- 0
#correct any predictions > 1:
pred.df.EDM_molecular$predicted_pos_uci[pred.df.EDM_molecular$predicted_pos_uci>1] <- 1


####Rousettus madagascariesis
pred.df.RM_molecular <- cbind.data.frame(month=rep(seq(1,12, by=.1),7*8), year=rep(c(2014,2015,2018, 2019, 2020, 2022, 2023), each=(length(seq(1,12,by=.1)*8))),
                                         bat_sex = rep(rep(c("female","male"), each=length(seq(1,12,by=.1))),7),
                                         virus = rep(rep(c("AndV", "AngV", "GhV", "GhV-like"), each = length(seq(1,12,by=.1))),7))
#predict without effect of year for plotting
pred.df.RM_molecular$predicted_pos <- predict.gam(rm1_molecular.f, newdata = pred.df.RM_molecular, type="response", exclude = "s(year)")
pred.df.RM_molecular$predicted_pos_lci <- pred.df.RM_molecular$predicted_pos - 1.96*predict.gam(rm1_molecular.f, newdata = pred.df.RM_molecular, type="response", se.fit = T, exclude = "s(year)")$se
pred.df.RM_molecular$predicted_pos_uci <- pred.df.RM_molecular$predicted_pos + 1.96*predict.gam(rm1_molecular.f, newdata = pred.df.RM_molecular, type="response", se.fit = T, exclude = "s(year)")$se
pred.df.RM_molecular$bat_species="Rousettus madagascariensis"
#correct any predictions < 0:
pred.df.RM_molecular$predicted_pos_lci[pred.df.RM_molecular$predicted_pos_lci<0] <- 0
#correct any predictions > 1:
pred.df.RM_molecular$predicted_pos_uci[pred.df.RM_molecular$predicted_pos_uci>1] <- 1

pred.df.RMF_molecular <- cbind.data.frame(month=rep(seq(1,12, by=.1),7*4), year=rep(c(2014,2015,2018, 2019, 2020, 2022, 2023), each=(length(seq(1,12,by=.1)*4))), 
                                          virus = rep(rep(c("AndV", "AngV", "GhV", "GhV-like"), each = length(seq(1,12,by=.1))),7))
#predict without effect of year for plotting
pred.df.RMF_molecular$predicted_pos <- predict.gam(rm1_molecular.f, newdata = pred.df.RMF_molecular, type="response", exclude = "s(year)")
pred.df.RMF_molecular$predicted_pos_lci <- pred.df.RMF_molecular$predicted_pos - 1.96*predict.gam(rm1_molecular.f, newdata = pred.df.RMF_molecular, type="response", se.fit = T, exclude = "s(year)")$se
pred.df.RMF_molecular$predicted_pos_uci <- pred.df.RMF_molecular$predicted_pos + 1.96*predict.gam(rm1_molecular.f, newdata = pred.df.RMF_molecular, type="response", se.fit = T, exclude = "s(year)")$se
pred.df.RMF_molecular$bat_species="Rousettus madagascariensis"
pred.df.RMF_molecular$bat_sex="female"
#correct any predictions < 0:
pred.df.RMF_molecular$predicted_pos_lci[pred.df.RMF_molecular$predicted_pos_lci<0] <- 0
#correct any predictions > 1:
pred.df.RMF_molecular$predicted_pos_uci[pred.df.RMF_molecular$predicted_pos_uci>1] <- 1

pred.df.RMM_molecular <- cbind.data.frame(month=rep(seq(1,12, by=.1),7*4), year=rep(c(2014,2015,2018, 2019, 2020, 2022, 2023), each=(length(seq(1,12,by=.1)*4))), 
                                          virus = rep(rep(c("AndV", "AngV", "GhV", "GhV-like"), each = length(seq(1,12,by=.1))),7))
#predict without effect of year for plotting
pred.df.RMM_molecular$predicted_pos <- predict.gam(rm1_molecular.m, newdata = pred.df.RMM_molecular, type="response", exclude = "s(year)")
pred.df.RMM_molecular$predicted_pos_lci <- pred.df.RMM_molecular$predicted_pos - 1.96*predict.gam(rm1_molecular.m, newdata = pred.df.RMM_molecular, type="response", se.fit = T, exclude = "s(year)")$se
pred.df.RMM_molecular$predicted_pos_uci <- pred.df.RMM_molecular$predicted_pos + 1.96*predict.gam(rm1_molecular.m, newdata = pred.df.RMM_molecular, type="response", se.fit = T, exclude = "s(year)")$se
pred.df.RMM_molecular$bat_species="Rousettus madagascariensis"
pred.df.RMM_molecular$bat_sex="male"
#correct any predictions < 0:
pred.df.RMM_molecular$predicted_pos_lci[pred.df.RMM_molecular$predicted_pos_lci<0] <- 0
#correct any predictions > 1:
pred.df.RMM_molecular$predicted_pos_uci[pred.df.RMM_molecular$predicted_pos_uci>1] <- 1

#mark which ones are sig

#Eidolon
#Overall: GhV sig, AngV, GhV-like, AndV not sig
#By sex: GhV-female almost sig, AndV, AngV, GhVlike neither sex sig
# pred.df.ED_molecular$sig <- 0
# pred.df.ED_molecular$sig[pred.df.ED_molecular$virus=="GhV"] <- 1
# pred.df.ED_molecular$sig[pred.df.ED_molecular$virus=="AngV"] <- 1
# pred.df.ED_molecular$sig[pred.df.ED_molecular$virus=="AndV"] <- 2
# pred.df.ED_molecular$sig[pred.df.ED_molecular$virus=="GhV-like"] <- 2
# 
# pred.df.EDM_molecular$sig <- 0
# pred.df.EDM_molecular$sig[pred.df.EDM_molecular$virus=="GhV"] <- 1
# pred.df.EDM_molecular$sig[pred.df.EDM_molecular$virus=="AngV"] <- 1
# pred.df.EDM_molecular$sig[pred.df.EDM_molecular$virus=="AndV"] <- 2
# pred.df.EDM_molecular$sig[pred.df.EDM_molecular$virus=="GhV-like"] <- 2
# 
# pred.df.EDF_molecular$sig <- 0
# pred.df.EDF_molecular$sig[pred.df.EDF_molecular$virus=="GhV"] <- 1
# pred.df.EDF_molecular$sig[pred.df.EDF_molecular$virus=="AngV"] <- 2
# pred.df.EDF_molecular$sig[pred.df.EDF_molecular$virus=="AndV"] <- 2
# pred.df.EDF_molecular$sig[pred.df.EDF_molecular$virus=="GhV-like"] <- 2

#Rousettus
#Overall: GhV sig, AngV, GhV-like, AndV not sig
#By sex: GhV-female almost sig, AndV, AngV, GhVlike neither sex sig
# pred.df.RM_molecular$sig <- 0
# pred.df.RM_molecular$sig[pred.df.RM_molecular$virus=="GhV"] <- 2
# pred.df.RM_molecular$sig[pred.df.RM_molecular$virus=="AngV"] <- 1
# pred.df.RM_molecular$sig[pred.df.RM_molecular$virus=="AndV"] <- 2
# pred.df.RM_molecular$sig[pred.df.RM_molecular$virus=="GhV-like"] <- 1
# 
# #males
# pred.df.RMM_molecular$sig <- 0
# pred.df.RMM_molecular$sig[pred.df.RMM_molecular$virus=="GhV"] <- 1
# pred.df.RMM_molecular$sig[pred.df.RMM_molecular$virus=="AngV"] <- 1
# pred.df.RMM_molecular$sig[pred.df.RMM_molecular$virus=="AndV"] <- 2
# pred.df.RMM_molecular$sig[pred.df.RMM_molecular$virus=="GhV-like"] <- 2
# 
# #females
# pred.df.RMF_molecular$sig <- 0
# pred.df.RMF_molecular$sig[pred.df.RMF_molecular$virus=="GhV"] <- 2
# pred.df.RMF_molecular$sig[pred.df.RMF_molecular$virus=="AngV"] <- 2
# pred.df.RMF_molecular$sig[pred.df.RMF_molecular$virus=="AndV"] <- 2
# pred.df.RMF_molecular$sig[pred.df.RMF_molecular$virus=="GhV-like"] <- 1


#add together
pred.df_molecular <- rbind(pred.df.ED_molecular, pred.df.RM_molecular)
pred.df.sex_molecular <- rbind(pred.df.EDF_molecular, pred.df.RMF_molecular, pred.df.EDM_molecular, pred.df.RMM_molecular)

# let's look at the partial effects
# because of the interaction term, the day of year function
# won't run, but you could pull the partial effects from 
# the plot() output above and visualize it if desired

# source script from Guth/Streicker/Mollentze
# generalcodewd = "/Users/gwenddolenkettenburg/Desktop/brook_lab/general_code_datasets"
# source(paste0(generalcodewd, "/partial-gam-functions.R"))
# 
# ###Eidolon 
# ed.molecular.virus.df <- get_partial_effects(ed1_molecular, var = "virus")
# ed.molecular.year.df <- get_partial_effects(ed1_molecular, var = "year")
# plot.partial(ed.molecular.virus.df, var="virus", response_var = "prob.positive") # no sig
# plot.partial(ed.molecular.year.df, var="year", response_var = "prob.positive") #2015 and 2016 pos sig, 2018 2019 neg sig
# 
# ed.molecular.virus.df.f <- get_partial_effects(ed1_molecular.f, var = "virus")
# ed.molecular.year.df.f <- get_partial_effects(ed1_molecular.f, var = "year")
# plot.partial(ed.molecular.virus.df.f, var="virus", response_var = "prob.positive") # no sig
# plot.partial(ed.molecular.year.df.f, var="year", response_var = "prob.positive") #2016 pos sig
# 
# ed.molecular.virus.df.m <- get_partial_effects(ed1_molecular.m, var = "virus")
# ed.molecular.year.df.m <- get_partial_effects(ed1_molecular.m, var = "year")
# plot.partial(ed.molecular.virus.df.m, var="virus", response_var = "prob.positive") # no sig
# plot.partial(ed.molecular.year.df.m, var="year", response_var = "prob.positive") #2015 pos sig, 2018 neg sig
# 
# ###Rousettus 
# rm.molecular.virus.df <- get_partial_effects(rm1_molecular, var = "virus")
# rm.molecular.year.df <- get_partial_effects(rm1_molecular, var = "year")
# plot.partial(rm.molecular.virus.df, var="virus", response_var = "prob.positive") # no sig
# plot.partial(rm.molecular.year.df, var="year", response_var = "prob.positive") #2022 pos sig
# 
# rm.molecular.virus.df.f <- get_partial_effects(rm1_molecular.f, var = "virus")
# rm.molecular.year.df.f <- get_partial_effects(rm1_molecular.f, var = "year")
# plot.partial(rm.molecular.virus.df.f, var="virus", response_var = "prob.positive") #no sig
# plot.partial(rm.molecular.year.df.f, var="year", response_var = "prob.positive") #no sig
# 
# rm.molecular.virus.df.m <- get_partial_effects(rm1_molecular.m, var = "virus")
# rm.molecular.year.df.m <- get_partial_effects(rm1_molecular.m, var = "year")
# plot.partial(rm.molecular.virus.df.m, var="virus", response_var = "prob.positive") #no sig
# plot.partial(rm.molecular.year.df.m, var="year", response_var = "prob.positive") #2022 pos sig

#for plotting
merge.plot <- ddply(molecular.dat,.(bat_species, year, month, virus), summarise, N=length(pos), pos=sum(pos))
merge.plot$molecular <-   merge.plot$pos/merge.plot$N

merge.plot.sex <- ddply(molecular.dat,.(bat_species, bat_sex, year, month, virus), summarise, N=length(pos), pos=sum(pos))
merge.plot.sex$molecular <-   merge.plot.sex$pos/merge.plot.sex$N

shapez=c("2013"=16,"2014"=17,"2015"=15, "2016"=3, "2018"=7, "2019"=8, "2020"=18, "2022"=9, "2023"=10)

###plotting by day of year
#all datasets
intra <- ggplot(pred.df_molecular) + 
  geom_rect(data=demo.dat.month, aes(xmin=start, xmax=stop, ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(data=merge.plot, aes(x=month, y=molecular, size=N, color=virus, shape=year)) +
  #geom_line(aes(x=month, y=predicted_pos, color=virus)) + 
  #geom_ribbon(data=subset(pred.df_molecular, sig==1), aes(x=month, ymin= predicted_pos_lci, ymax=predicted_pos_uci, fill=virus), alpha=.3) +
  #geom_line(data=subset(pred.df_molecular, sig==1), aes(x=month, y=predicted_pos, color=virus)) +
  geom_ribbon(data=subset(pred.df_molecular), aes(x=month, ymin= predicted_pos_lci, ymax=predicted_pos_uci, fill=virus), alpha=.3) +
  geom_line(data=subset(pred.df_molecular), aes(x=month, y=predicted_pos, color=virus)) +
  facet_grid(bat_species~virus) + 
  theme_bw() +
  scale_color_manual(values=colz1) +
  scale_fill_manual(values=filz) +
  scale_shape_manual(values=shapez) +
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(), 
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=9), panel.grid = element_blank(), legend.title = element_blank()) +
  guides(fill="none")+
  coord_cartesian(expand = F) + ylab("Predicted prevalence") + xlab("Month")
intra

intra.sex <- ggplot(pred.df.sex_molecular) + 
  geom_rect(data=demo.dat.month, aes(xmin=start, xmax=stop, ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(data=merge.plot.sex, aes(x=month, y=molecular, size=N, color=virus, shape=year)) +
  #geom_line(aes(x=month, y=predicted_pos, color=virus)) + 
  #geom_ribbon(data=subset(pred.df.sex_molecular, sig==1), aes(x=month, ymin= predicted_pos_lci, ymax=predicted_pos_uci, fill=virus), alpha=.3) +
  #geom_line(data=subset(pred.df.sex_molecular, sig==1), aes(x=month, y=predicted_pos, color=virus)) +
  geom_ribbon(data=subset(pred.df.sex_molecular), aes(x=month, ymin= predicted_pos_lci, ymax=predicted_pos_uci, fill=virus), alpha=.3) +
  geom_line(data=subset(pred.df.sex_molecular), aes(x=month, y=predicted_pos, color=virus)) +
  facet_grid(bat_species~virus+bat_sex) + 
  theme_bw() +
  scale_color_manual(values=colz1) +
  scale_fill_manual(values=filz) +
  scale_shape_manual(values=shapez) +
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(),
        strip.text.x=element_text(face="italic"),
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text = element_text(size=9), panel.grid = element_blank(), 
        legend.title = element_blank()) +
  guides(fill="none")+
  coord_cartesian(expand = F) + ylab("Predicted prevalence") + xlab("Month") 
intra.sex



####################################################################################################
####################################################################################################
## SECTION 4: look at age class and recaps

# could plot age-molecular?
age.Eid <- subset(molecular.dat, bat_species=="Eidolon dupreanum")
age.Eid <- ddply(age.Eid, .(bat_species, virus, age), summarise, N=length(pos), pos=sum(pos))
age.Eid$molecular <- age.Eid$pos/age.Eid$N

age.Rou <- subset(molecular.dat, bat_species=="Rousettus madagascariensis")
age.Rou <- ddply(age.Rou, .(bat_species, virus, age), summarise, N=length(pos), pos=sum(pos))
age.Rou$molecular <- age.Rou$pos/age.Rou$N

age.Pruf <- subset(molecular.dat, bat_species=="Pteropus rufus")
age.Pruf <- ddply(age.Pruf, .(bat_species, virus, age), summarise, N=length(pos), pos=sum(pos))
age.Pruf$molecular <- age.Pruf$pos/age.Pruf$N

ggplot(age.Eid) + geom_point(aes(x=age, y=molecular, size=N)) 
ggplot(age.Rou) + geom_point(aes(x=age, y=molecular, size=N)) 
ggplot(age.Pruf) + geom_point(aes(x=age, y=molecular, size=N)) 


####################################################################################################
####################################################################################################
## SECTION 5: #and look at data for the recaptured individuals

#collect the recapture date 
recap.dat = subset(molecular.dat, recap_new=="yes")
recap.dat$barcode_id <- as.factor(recap.dat$barcode_id)

#subset to the MFIs of bats that were positive for one antigen once
pos.once <- ddply(recap.dat, .(barcode_id, virus), summarise, sum_pos=sum(pos),  Ncap=length(unique(processing_date)))
#retain those pos once
pos.once = subset(pos.once, sum_pos>0 & Ncap>1 )
pos.once <- dplyr::select(pos.once, -(virus))

recap.plot <- merge(recap.dat, pos.once, by="barcode_id")
head(recap.plot)
#recap.plot = subset(recap.plot, virus == "pmx")
# recap.sum = ddply(recap.plot, .(barcode_id), summarise, N= sum(pos))
# recap.plot <- merge(recap.plot, recap.sum, by="barcode_id")
# recap.plot <- subset(recap.plot, N>1)

# Could repeat across the other sample types
Eid.cal = subset(demo.dat.day, bat_species=="Eidolon dupreanum")
Rou.cal = subset(demo.dat.day, bat_species=="Rousettus madagascariensis")
Pruf.cal = subset(demo.dat.day, bat_species=="Pteropus rufus")

shapez = c("female"=21, "male"=24)
filz = c('winter' = "#00A9FF", 'gestation' = "#FF61CC", 'lactation' = "#C77CFF", '0'="white", '1'="black")
colz = c('0'="white", '1'="black")
filz = c("winter"  = "yellow", "gestation" = "darksalmon", "lactation" = "cornflowerblue", '0'="white", '1'="black")

recap.plot$pos <- as.factor(recap.plot$pos)


recap_Eid <- ggplot(data=subset(recap.plot, bat_species=="Eidolon dupreanum")) +theme_bw() +
  #geom_rect(data=Eid.cal, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=-Inf, ymax=Inf, fill=type), alpha=.3) +
  geom_line(aes(x=processing_date, y=barcode_id), show.legend = F, size=1) +
  geom_point(aes(x=processing_date, y=barcode_id, fill=pos, shape=bat_sex),  size=3)+
  scale_shape_manual(values=shapez) +scale_fill_manual(values=filz) + scale_color_manual(values=colz) +
  theme(strip.background = element_rect(fill="white"), strip.text = element_text(size=16),
        #panel.spacing = unit(1.4, "cm"),
        panel.grid = element_blank(),  axis.title.x = element_blank(),
        axis.title.y = element_text(size=16), axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),  
        plot.margin = unit(c(.2,.1,.1,.1), "cm"))
recap_Eid

#and try rousettus
recap_Rou <- ggplot(data=subset(recap.plot, bat_species=="Rousettus madagascariensis")) +theme_bw() +
  #geom_rect(data=Rou.cal, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=-Inf, ymax=Inf, fill=type), alpha=.3) +
  geom_line(aes(x=processing_date, y=barcode_id), show.legend = F, size=1) +
  geom_point(aes(x=processing_date, y=barcode_id, fill=pos, shape=bat_sex),  size=3)+
  scale_shape_manual(values=shapez) +scale_fill_manual(values=filz) + scale_color_manual(values=colz) +
  theme(strip.background = element_rect(fill="white"), strip.text = element_text(size=16),
        #panel.spacing = unit(1.4, "cm"),
        panel.grid = element_blank(),  axis.title.x = element_blank(),
        axis.title.y = element_text(size=16), axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),  
        plot.margin = unit(c(.2,.1,.1,.1), "cm"))
recap_Rou

#and try pteropus
recap_Pruf <- ggplot(data=subset(recap.plot, bat_species=="Pteropus rufus")) +theme_bw() +
  #geom_rect(data=Pruf.cal, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=-Inf, ymax=Inf, fill=type), alpha=.3) +
  geom_line(aes(x=processing_date, y=barcode_id), show.legend = F, size=1) +
  geom_point(aes(x=processing_date, y=barcode_id, fill=pos, shape=bat_sex),  size=3)+
  scale_shape_manual(values=shapez) +scale_fill_manual(values=filz) + scale_color_manual(values=colz) +
  theme(strip.background = element_rect(fill="white"), strip.text = element_text(size=16),
        #panel.spacing = unit(1.4, "cm"),
        panel.grid = element_blank(),  axis.title.x = element_blank(),
        axis.title.y = element_text(size=16), axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),  
        plot.margin = unit(c(.2,.1,.1,.1), "cm"))
recap_Pruf
