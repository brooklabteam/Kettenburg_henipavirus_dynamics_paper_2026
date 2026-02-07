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

# # remove bats from the dataset with no age class listed
molecular.dat = subset(merge.dat.pos, !is.na(bat_age_class))
molecular.dat = subset(molecular.dat, !is.na(bat_sex))
molecular.dat = subset(molecular.dat, !is.na(bat_species))

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
#month
lact.dat.month <- cbind.data.frame(bat_species=rep(c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"), each=2), 
                                   start = c(11.5,1, 11.4,1, 8.1,1), 
                                   stop=c(12,2.05,12,2.7,12.9,1))
lact.dat.month$type = "lactation"
lact.dat.month$bat_sex="female"
gest.dat.month <- cbind.data.frame(bat_species=c("Eidolon dupreanum","Rousettus madagascariensis", "Pteropus rufus"), 
                                   start = c(8.15, 9.45, 7.2),
                                   stop=c(11.5, 11.39, 8.1))
gest.dat.month$type = "gestation"
gest.dat.month$bat_sex = "female"
#and winter
cal.dat.month <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"), 
                                   start = c(5.5,5.5,5.5),
                                   stop=c(10.3,10.3,10.3))
cal.dat.month$type <- "winter"
cal.dat.month$bat_sex="male"

library(reshape2)

repro.dat.month <- rbind(lact.dat.month, gest.dat.month)
demo.dat.month <- rbind(gest.dat.month, lact.dat.month, cal.dat.month)

#for plotting
merge.plot <- ddply(molecular.dat,.(bat_species, year, month, virus), summarise, N=length(pos), pos=sum(pos))
merge.plot$molecular <-   merge.plot$pos/merge.plot$N

merge.plot.noyear <- ddply(molecular.dat,.(bat_species, month, virus), summarise, N=length(pos), pos=sum(pos))
merge.plot.noyear$molecular <-   merge.plot.noyear$pos/merge.plot.noyear$N

merge.plot.sex <- ddply(molecular.dat,.(bat_species, bat_sex, year, month, virus), summarise, N=length(pos), pos=sum(pos))
merge.plot.sex$molecular <-   merge.plot.sex$pos/merge.plot.sex$N

merge.plot.sex.noyear <- ddply(molecular.dat,.(bat_species, bat_sex, month, virus), summarise, N=length(pos), pos=sum(pos))
merge.plot.sex.noyear$molecular <-   merge.plot.sex.noyear$pos/merge.plot.sex.noyear$N

shapez=c("2013"=16,"2014"=17,"2015"=15, "2016"=3, "2018"=7, "2019"=8, "2020"=18, "2022"=9, "2023"=10)

###plotting by day of year
#all datasets
intra <- 
  ggplot(merge.plot) + 
  #ggplot(merge.plot.noyear) + 
  geom_rect(data=demo.dat.month, aes(xmin=start, xmax=stop, ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(data=merge.plot, aes(x=month, y=molecular, size=N, color=virus, shape=year)) +
  #geom_point(data=merge.plot.noyear, aes(x=month, y=molecular, size=N, color=virus, shape=year)) +
  facet_grid(bat_species~.) + 
  theme_bw() +
  scale_color_manual(values=colz1) +
  scale_fill_manual(values=filz) +
  scale_shape_manual(values=shapez) +
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(), 
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=9), panel.grid = element_blank(), legend.title = element_blank()) +
  guides(fill="none")+
  coord_cartesian(expand = F) + ylab("Prevalence") + xlab("Month")
intra

intra.sex <- 
  ggplot(merge.plot.sex) + 
  #ggplot(merge.plot.sex.noyear) + 
  geom_rect(data=demo.dat.month, aes(xmin=start, xmax=stop, ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(data=merge.plot.sex, aes(x=month, y=molecular, size=N, color=virus, shape=year)) +
  #geom_point(data=merge.plot.sex.noyear, aes(x=month, y=molecular, size=N, color=virus, shape=year)) +
  facet_grid(bat_species~bat_sex) + 
  theme_bw() +
  scale_color_manual(values=colz1) +
  scale_fill_manual(values=filz) +
  scale_shape_manual(values=shapez) +
  theme(axis.text = element_text(size=9), #axis.title.x = element_blank(), 
        axis.title = element_text(size=9), strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=9), panel.grid = element_blank(), legend.title = element_blank()) +
  guides(fill="none")+
  coord_cartesian(expand = F) + ylab("Prevalence") + xlab("Month")
intra.sex



