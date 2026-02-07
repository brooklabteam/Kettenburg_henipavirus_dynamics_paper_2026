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
shed.dat.pos <- dplyr::select(shed.dat, sampleid, AndV_shedding_ur)

## merge the shedding data with the catching data

#and make into long format for prevalence
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
#merge.dat.pos$virus[merge.dat.pos$virus=="AndV_shedding_ct_ur"] <- "AndV_ct"

# merge.dat.pos <- rename(merge.dat.pos, pos = AndV_shedding_ur)
# merge.dat.pos <- rename(merge.dat.pos, ct = AndV_shedding_ct_ur)

#subset into the three time series datasets
merge.dat.pos_old <- merge.dat.pos[merge.dat.pos$mid_date > "2014-10-13" & merge.dat.pos$mid_date < "2016-01-19", ]
merge.dat.pos_mid <- merge.dat.pos[merge.dat.pos$mid_date > "2018-02-15" & merge.dat.pos$mid_date < "2019-09-21", ]
merge.dat.pos_new <- merge.dat.pos[merge.dat.pos$mid_date > "2022-07-19" & merge.dat.pos$mid_date < "2023-05-12", ]

prev.sum.old <- ddply(merge.dat.pos_old, .(bat_species, region, virus, mid_date), summarise, N=length(pos), pos=sum(pos))
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

colz1 = c('winter' = "#35A2FF", "gestation"="#E38900", "lactation"="#7CAE00",
           
           'AndV' = "black")

filz = c('winter' = "#35A2FF", "gestation"="#E38900", "lactation"="#7CAE00",
         
         'AndV' = "black")

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
lact.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", 
                                           "Rousettus madagascariensis","Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis",
                                           "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus"), 
                             start = c(as.Date("2013-11-16"), as.Date("2014-11-16"), as.Date("2015-11-16"), as.Date("2016-11-16"), as.Date("2017-11-16"), as.Date("2018-11-16"), as.Date("2019-11-16"), as.Date("2020-11-16"),as.Date("2022-11-16"),
                                       as.Date("2013-12-12"), as.Date("2014-12-12"), as.Date("2015-12-12"), as.Date("2016-12-12"), as.Date("2017-12-12"), as.Date("2018-12-12"), as.Date("2019-12-12"), as.Date("2020-12-12"), as.Date("2022-12-12"),
                                       as.Date("2013-08-03"), as.Date("2014-08-03"), as.Date("2015-08-03"), as.Date("2016-08-03"), as.Date("2017-08-03"), as.Date("2018-08-03"), as.Date("2019-08-03"), as.Date("2020-08-03"), as.Date("2022-08-03")), 
                             stop=c(as.Date("2014-02-02"), as.Date("2015-02-02"), as.Date("2016-02-02"), as.Date("2017-02-02"), as.Date("2018-02-02"), as.Date("2019-02-02"), as.Date("2020-02-02"), as.Date("2021-02-02"),  as.Date("2023-02-02"),
                                    as.Date("2014-02-19"), as.Date("2015-02-19"), as.Date("2016-02-19"), as.Date("2017-02-19"), as.Date("2018-02-19"), as.Date("2019-02-19"), as.Date("2020-02-19"), as.Date("2021-02-19"),  as.Date("2023-02-19"),
                                    as.Date("2014-01-21"), as.Date("2015-01-21"), as.Date("2016-01-21"), as.Date("2017-01-21"), as.Date("2018-01-21"), as.Date("2019-01-21"), as.Date("2020-01-21"), as.Date("2021-01-21"),  as.Date("2023-01-21")))
lact.dat$type = "lactation"
gest.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum",
                                           "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis",
                                           "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus"), 
                             start = c(as.Date("2013-08-03"), as.Date("2014-08-03"), as.Date("2015-08-03"), as.Date("2016-08-03"),as.Date("2017-08-03"), as.Date("2018-08-03"), as.Date("2019-08-03"), as.Date("2020-08-03"),  as.Date("2022-08-03"),
                                       as.Date("2013-09-11"), as.Date("2014-09-11"), as.Date("2015-09-11"), as.Date("2016-09-11"), as.Date("2017-09-11"), as.Date("2018-09-11"), as.Date("2019-09-11"), as.Date("2020-09-11"),as.Date("2022-09-11"),
                                       as.Date("2013-07-07"), as.Date("2014-07-07"), as.Date("2015-07-07"), as.Date("2016-07-07"), as.Date("2017-07-07"), as.Date("2018-07-07"), as.Date("2019-07-07"), as.Date("2020-07-07"),  as.Date("2022-07-07")), 
                             stop=c(as.Date("2013-11-15"), as.Date("2014-11-15"), as.Date("2015-11-15"), as.Date("2016-11-15"), as.Date("2017-11-15"), as.Date("2018-11-15"), as.Date("2019-11-15"), as.Date("2020-11-15"),  as.Date("2022-11-15"),
                                    as.Date("2013-12-11"), as.Date("2014-12-11"), as.Date("2015-12-11"), as.Date("2016-12-11"), as.Date("2017-12-11"), as.Date("2018-12-11"), as.Date("2019-12-11"), as.Date("2020-12-11"),  as.Date("2022-12-11"),
                                    as.Date("2013-08-02"), as.Date("2014-08-02"), as.Date("2015-08-02"), as.Date("2016-08-02"), as.Date("2017-08-02"), as.Date("2018-08-02"), as.Date("2019-08-02"), as.Date("2020-08-02"),  as.Date("2022-08-02")))
gest.dat$type = "gestation"


#and winter
cal.dat <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum",
                                           "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis","Rousettus madagascariensis", "Rousettus madagascariensis",
                                           "Pteropus rufus","Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus", "Pteropus rufus"), 
                             start = c(as.Date("2013-05-15"), as.Date("2014-05-15"), as.Date("2015-05-15"), as.Date("2016-05-15"),as.Date("2017-05-15"), as.Date("2018-05-15"), as.Date("2019-05-15"), as.Date("2020-05-15"), as.Date("2021-05-15"), 
                                       as.Date("2013-05-15"), as.Date("2014-05-15"), as.Date("2015-05-15"), as.Date("2016-05-15"), as.Date("2017-05-15"), as.Date("2018-05-15"), as.Date("2019-05-15"), as.Date("2020-05-15"),  as.Date("2021-05-15"),
                                       as.Date("2013-05-15"), as.Date("2014-05-15"), as.Date("2015-05-15"), as.Date("2016-05-15"), as.Date("2017-05-15"), as.Date("2018-05-15"), as.Date("2019-05-15"), as.Date("2020-05-15"),  as.Date("2021-05-15")), 
                             stop=c(as.Date("2013-10-10"), as.Date("2014-10-10"), as.Date("2015-10-10"), as.Date("2016-10-10"), as.Date("2017-10-10"), as.Date("2018-10-10"), as.Date("2019-10-10"), as.Date("2020-10-10"),  as.Date("2021-10-10"),
                                    as.Date("2013-10-10"), as.Date("2014-10-10"), as.Date("2015-10-10"), as.Date("2016-10-10"), as.Date("2017-10-10"), as.Date("2018-10-10"), as.Date("2019-10-10"), as.Date("2020-10-10"),  as.Date("2021-10-10"),
                                    as.Date("2013-10-10"), as.Date("2014-10-10"), as.Date("2015-10-10"), as.Date("2016-10-10"), as.Date("2017-10-10"), as.Date("2018-10-10"), as.Date("2019-10-10"), as.Date("2020-10-10"), as.Date("2021-10-10")))
cal.dat$type <- "winter"

library(reshape2)

repro.dat <- rbind(lact.dat, gest.dat)
demo.dat <- rbind(gest.dat, lact.dat, cal.dat)


#old dataset 2013-2016
lact.dat.old <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum",
                                               "Rousettus madagascariensis", "Rousettus madagascariensis",
                                               "Pteropus rufus", "Pteropus rufus"), 
                                 start = c(as.Date("2014-11-16"), as.Date("2015-11-16"), 
                                            as.Date("2014-12-12"), as.Date("2015-12-12"),
                                            as.Date("2014-08-13"), as.Date("2015-08-03")), 
                                 stop=c(as.Date("2015-02-02"), as.Date("2016-02-02"),
                                        as.Date("2015-02-19"), as.Date("2016-02-19"),
                                        as.Date("2015-01-21"), as.Date("2016-01-21")))
lact.dat.old$type = "lactation"
lact.dat.old$bat_sex = "female"

gest.dat.old <- cbind.data.frame(bat_species=c("Eidolon dupreanum",
                                               "Rousettus madagascariensis", 
                                               "Pteropus rufus"), 
                                 start = c(as.Date("2015-08-03"), 
                                           as.Date("2015-09-11"),
                                           as.Date("2015-07-07")), 
                                 stop=c(as.Date("2015-11-15"),
                                        as.Date("2015-12-11"),
                                        as.Date("2015-08-02")))
gest.dat.old$type = "gestation"
gest.dat.old$bat_sex="female"

#and winter
cal.dat.old <-  cbind.data.frame(bat_species=c("Eidolon dupreanum",
                                               "Rousettus madagascariensis",
                                               "Pteropus rufus"), 
                                 start = c(as.Date("2015-05-15"), 
                                            as.Date("2015-05-15"),
                                            as.Date("2015-05-15")), 
                                 stop=c(as.Date("2015-10-10"), 
                                        as.Date("2015-10-10"),
                                        as.Date("2015-10-10")))
cal.dat.old$type <- "winter"
cal.dat.old$bat_sex <-"male"

library(reshape2)

repro.dat.old <- rbind(lact.dat.old, gest.dat.old)
demo.dat.old <- rbind(gest.dat.old, lact.dat.old, cal.dat.old)




#mid dataset 2018-2020
lact.dat.mid <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum",
                                               "Rousettus madagascariensis","Rousettus madagascariensis", 
                                               "Pteropus rufus"), 
                                 start = c(as.Date("2018-11-16"), as.Date("2019-11-16"),
                                           as.Date("2018-12-12"), as.Date("2019-12-12"),
                                           as.Date("2019-08-03")), 
                                 stop=c(as.Date("2019-02-02"), as.Date("2020-02-02") ,
                                        as.Date("2019-02-19"), as.Date("2020-02-19") ,
                                        as.Date("2020-01-21")))
lact.dat.mid$type = "lactation"
lact.dat.mid$bat_sex ="female"

gest.dat.mid <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", 
                                               "Rousettus madagascariensis", "Rousettus madagascariensis", 
                                               "Pteropus rufus"), 
                                 start = c(as.Date("2018-08-03"), as.Date("2019-08-03"),
                                           as.Date("2018-09-11"), as.Date("2019-09-11"),
                                           as.Date("2019-07-07")), 
                                 stop=c(as.Date("2018-11-15"), as.Date("2019-11-15"),
                                        as.Date("2018-12-11"), as.Date("2019-12-11"),
                                        as.Date("2019-08-02")))
gest.dat.mid$type = "gestation"
gest.dat.mid$bat_sex ="female"

#and winter
cal.dat.mid <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", 
                                               "Rousettus madagascariensis", "Rousettus madagascariensis",
                                               "Pteropus rufus","Pteropus rufus"), 
                                 start = c(as.Date("2018-05-15"), as.Date("2019-05-15"),
                                           as.Date("2018-05-15"), as.Date("2019-05-15"),
                                           as.Date("2018-05-15"), as.Date("2019-05-15")), 
                                 stop=c(as.Date("2018-10-10"), as.Date("2019-10-10"),
                                        as.Date("2018-10-10"), as.Date("2019-10-10"),
                                        as.Date("2018-10-10"), as.Date("2019-10-10")))
cal.dat.mid$type <- "winter"
cal.dat.mid$bat_sex="male"

library(reshape2)

repro.dat.mid <- rbind(lact.dat.mid, gest.dat.mid)
demo.dat.mid <- rbind(gest.dat.mid, lact.dat.mid, cal.dat.mid)

demo.dat.old.ed = subset(demo.dat.old, bat_species=="Eidolon dupreanum")
demo.dat.mid.ed = subset(demo.dat.mid, bat_species=="Eidolon dupreanum")
demo.dat.old.rm = subset(demo.dat.old, bat_species=="Rousettus madagascariensis")
demo.dat.mid.rm = subset(demo.dat.mid, bat_species=="Rousettus madagascariensis")
demo.dat.old.pr = subset(demo.dat.old, bat_species=="Pteropus rufus")
demo.dat.mid.pr = subset(demo.dat.mid, bat_species=="Pteropus rufus")

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


#and plot with data:

old.ed <- ggplot(data=subset(prev.sum.sex.old, bat_species=="Eidolon dupreanum")) +
  geom_rect(data=demo.dat.old.ed, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values=colz1) + scale_fill_manual(values=filz) +
  facet_grid(bat_sex~bat_species+region) + 
  theme_bw()+
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(), 
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=9), panel.grid = element_blank(), 
        legend.title = element_blank())
old.ed

mid.ed <- ggplot(data=subset(prev.sum.sex.mid, bat_species=="Eidolon dupreanum")) +
  geom_rect(data=demo.dat.mid.ed, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values=colz1) + scale_fill_manual(values=filz) +
  facet_grid(bat_sex~bat_species+region) + 
  theme_bw()+
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(), 
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=9), panel.grid = element_blank(), 
        legend.title = element_blank())
mid.ed


##Now do rousettus
old.rm <- ggplot(data=subset(prev.sum.sex.old, bat_species=="Rousettus madagascariensis")) +
  geom_rect(data=demo.dat.old.rm, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values=colz1) + scale_fill_manual(values=filz) +
  #geom_line(aes(x=mid_date, y=prev, color=virus), linetype=2) +
  facet_grid(bat_sex~bat_species+region) + 
  theme_bw()+
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(), 
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=9), panel.grid = element_blank(), 
        legend.title = element_blank())
old.rm 

mid.rm <- ggplot(data=subset(prev.sum.sex.mid, bat_species=="Rousettus madagascariensis")) +
  geom_rect(data=demo.dat.mid.rm, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values=colz1) + scale_fill_manual(values=filz) +
  #geom_line(aes(x=mid_date, y=prev, color=virus), linetype=2) +
  facet_grid(bat_sex~bat_species+region) + 
  theme_bw()+
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(), 
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=9), panel.grid = element_blank(), 
        legend.title = element_blank())
mid.rm 

##Now do pteropus
old.pr <- ggplot(data=subset(prev.sum.sex.old, bat_species=="Pteropus rufus")) +
  geom_rect(data=demo.dat.old.pr, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values=colz1) + scale_fill_manual(values=filz) +
  #geom_line(aes(x=mid_date, y=prev, color=virus), linetype=2) +
  facet_grid(bat_sex~bat_species+region) + 
  theme_bw()+
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(), 
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=9), panel.grid = element_blank(), 
        legend.title = element_blank())
old.pr 

mid.pr <- ggplot(data=subset(prev.sum.sex.mid, bat_species=="Pteropus rufus")) +
  geom_rect(data=demo.dat.mid.pr, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values=colz1) + scale_fill_manual(values=filz) +
  #geom_line(aes(x=mid_date, y=prev, color=virus), linetype=2) +
  facet_grid(bat_sex~bat_species+region) + 
  theme_bw()+
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(), 
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=9), panel.grid = element_blank(), 
        legend.title = element_blank())
mid.pr 


#all plots together

#add category for the time period
prev.sum.sex.old$dataset <- "2014-2016"
prev.sum.sex.mid$dataset <- "2018-2020"
prev.sum.old$dataset <- "2014-2016"
prev.sum.mid$dataset <- "2018-2020"

demo.dat.old$dataset<-"2014-2016"
demo.dat.mid$dataset<-"2018-2020"

#bind the time periods together for plotting
prev.sum.sex.plot <- rbind(prev.sum.sex.old, prev.sum.sex.mid)
prev.sum.plot <- rbind(prev.sum.old, prev.sum.mid)

#bind together the demographic events
demo.dat.plot <- rbind(demo.dat.old, demo.dat.mid)

mid <- ggplot(data=subset(prev.sum.plot)) +
  geom_rect(data=demo.dat.plot, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values=colz1) + scale_fill_manual(values=filz) +
  geom_line(aes(x=mid_date, y=prev, color=virus), linetype=2) +
  facet_grid(bat_species~dataset, scales="free") + 
  #scale_x_date(limits = as.Date(c("2018-02-18","2019-09-20")))+
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
