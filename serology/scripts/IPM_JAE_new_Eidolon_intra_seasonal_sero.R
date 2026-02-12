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

homewd = "/Users/gwenddolenkettenburg/Desktop/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
batwd = "/Users/gwenddolenkettenburg/Desktop/Developer/Madagascar-Bat-Data"

# plot of seasonality of antigens for analysis
# also GAMs of predictors of seropositivity and/or seroprevalence
################################################################################
################################################################################
## Section 1. Prep the data

# First, load the luminex data after determining who is seropositive/
# seronegative in the other scripts
lum.dat <- read.csv(file = paste0(homewd, "/serology/data/mix_sub_data_75_90_96_long.csv"))
head(lum.dat)
lum.dat = subset(lum.dat, bat_species=="Eidolon dupreanum"| bat_species=="Rousettus madagascariensis" & antigen=="GhV")
  
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
#names(merge.dat)[names(merge.dat)=="mid_date"] <- "mid_date_old"
#names(merge.dat)[names(merge.dat)=="processing_date"] <- "mid_date"

# #curate date once more
 sampling.date <- ddply(merge.dat, .(bat_species, roost_site, mid_date), summarise)
 sampling.date# examine
# #write over those that are super close
 merge.dat$mid_date <- as.character(merge.dat$mid_date)
 merge.dat$mid_date[merge.dat$bat_species=="Eidolon dupreanum" & merge.dat$roost_site=="Moramanga" & merge.dat$mid_date=="2018-02-15"] <- "2018-02-14"
 merge.dat$mid_date[merge.dat$bat_species=="Rousettus madagascariensis" & merge.dat$roost_site=="Moramanga" & merge.dat$mid_date=="2019-01-09"] <- "2019-01-05"
 merge.dat$mid_date <- as.Date(merge.dat$mid_date)
 
 
####################################################################################################
####################################################################################################
## SECTION 1: Seasonality of seroprevalence 

# now look at seroprevalence by sampling session by species and antigen
# plot data first
seroprev.sum <- ddply(merge.dat, .(bat_species, roost_site, antigen, mid_date), summarise, N=length(seropos), seropos=sum(seropos))
#seroprev.sum <- ddply(merge.dat, .(bat_species, roost_site, antigen, quarter_year), summarise, N=length(seropos), seropos=sum(seropos))
#names(seroprev.sum)[names(seroprev.sum)=="quarter_year"] <- "mid_date"
#seroprev.sum <- ddply(merge.dat, .(bat_species, roost_site, antigen, mid_date), summarise, N=length(seropos), seropos=sum(seropos))
seroprev.sum$seroneg <- seroprev.sum$N - seroprev.sum$seropos
seroprev.sum$seroprev <- seroprev.sum$seropos/seroprev.sum$N
seroprev.sum$seroprev_lci <-binom.exact(x=seroprev.sum$seropos, n=seroprev.sum$N)$lower
seroprev.sum$seroprev_uci <-binom.exact(x=seroprev.sum$seropos, n=seroprev.sum$N)$upper
head(seroprev.sum)

#and plot

# it's pretty evident there are correlations in seasonality 
# between these antigens

colorz = c('AngV' = "#F8766D", 'HeV' ="#A58AFF", 
           'NiV' = "#DB72FB", 'GhV' = "#00C094", 
           'CedV'= "#53B400", "MojV" = "#619CFF")

p1 <- ggplot(data=seroprev.sum) + geom_point(aes(x=mid_date, y=seroprev, size=N, color=antigen)) +
  geom_line(aes(x=mid_date, y=seroprev, color=antigen)) + scale_color_manual(values = colorz) +
  facet_grid(bat_species~roost_site) + theme_bw() 
p1

# time series is super short for Ankarana and does not span a whole year,
# so let's remove it from the data for seasonal analysis
# we can come back to it for age and day of year

# add in the timing of the reproductive season in Moramanga for both species
# from angelo and santino's paper
# eidolon, gestation is from: Aug-3 - Nov 15; lactation Nov 16 - Feb 2
# rousettus, gestation is from: Sept 11 - Dec 11; lactation Dec 12 - Feb 19
# dry season = May 15 - October 10

gest.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis", "Rousettus madagascariensis"), start = c(as.Date("2018-08-03"), as.Date("2019-08-03"), as.Date("2018-09-11"), as.Date("2019-09-11")), stop=c(as.Date("2018-11-15"), as.Date("2019-11-15"), as.Date("2018-12-11"), as.Date("2019-12-11")))
gest.dat$season = "gestation"
lact.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis"), start = c(as.Date("2018-01-01"), as.Date("2018-11-16"), as.Date("2019-11-16"), as.Date("2018-01-01"), as.Date("2018-12-12"), as.Date("2019-12-12")), stop=c(as.Date("2018-02-02"), as.Date("2019-02-02"), as.Date("2020-02-02"), as.Date("2018-02-19"), as.Date("2019-02-19"), as.Date("2020-02-19")))
lact.dat$season = "lactation"

#and winter
cal.dat <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis", "Rousettus madagascariensis"), start = c(as.Date("2018-05-15"), as.Date("2019-05-15"), as.Date("2018-05-15"), as.Date("2019-05-15")), stop=c(as.Date("2018-10-10"), as.Date("2019-10-10"), as.Date("2018-10-10"), as.Date("2019-10-10")))
cal.dat$season <- "winter"

cal.dat <- rbind(cal.dat, gest.dat, lact.dat)
cal.dat.save=cal.dat

#now add to plot

p2 <- ggplot(data=subset(seroprev.sum, roost_site=="Moramanga")) + 
  geom_rect(data=cal.dat, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=mid_date, y=seroprev, size=N, color=antigen)) +
  geom_line(aes(x=mid_date, y=seroprev, color=antigen)) + scale_color_manual(values = colorz) +
  facet_grid(bat_species~.) + theme_bw()

p2

#what if you did males and females separately?
cal.dat$bat_sex <- "female"
cal.dat$bat_sex[cal.dat$season=="winter"] <- "male"

seroprev.sum.sex <- ddply(merge.dat, .(bat_species, roost_site, bat_sex, antigen,  mid_date), summarise, N=length(seropos), seropos=sum(seropos))
#seroprev.sum.sex <- ddply(merge.dat, .(bat_species, roost_site, bat_sex, antigen,  quarter_year), summarise, N=length(seropos), seropos=sum(seropos))
#names(seroprev.sum.sex)[names(seroprev.sum.sex)=="quarter_year"] <- "mid_date"
#seroprev.sum.sex <- ddply(merge.dat, .(bat_species, roost_site, bat_sex, antigen,  mid_date), summarise, N=length(seropos), seropos=sum(seropos))
seroprev.sum.sex$seroneg <- seroprev.sum.sex$N - seroprev.sum.sex$seropos
seroprev.sum.sex$seroprev <- seroprev.sum.sex$seropos/seroprev.sum.sex$N
seroprev.sum.sex$seroprev_lci <-binom.exact(x=seroprev.sum.sex$seropos, n=seroprev.sum.sex$N)$lower
seroprev.sum.sex$seroprev_uci <-binom.exact(x=seroprev.sum.sex$seropos, n=seroprev.sum.sex$N)$upper

head(seroprev.sum.sex)
seroprev.sum.sex = subset(seroprev.sum.sex, bat_sex=="female" | bat_sex=="male")

p3 <- ggplot(data=subset(seroprev.sum.sex, roost_site=="Moramanga")) + 
  geom_rect(data=cal.dat, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=mid_date, y=seroprev, size=N, color=antigen)) + scale_color_manual(values = colorz) +
  geom_line(aes(x=mid_date, y=seroprev, color=antigen)) +
  facet_grid(bat_species~bat_sex) + theme_bw()

p3
# interesting differences here

#check out Ankarana
p3b <- ggplot(data=subset(seroprev.sum.sex, roost_site=="Ankarana")) + 
  geom_rect(data=cal.dat, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=mid_date, y=seroprev, size=N, color=antigen)) +
  geom_line(aes(x=mid_date, y=seroprev, color=antigen)) + scale_color_manual(values = colorz) +
  facet_grid(bat_species~bat_sex) + theme_bw()

p3b

# time series is too short here, so we will allow the model to zero out data 
# from the Ankarana site if it does not find it useful. 

# we will need to make a bonferonni correction on p-values
# since we are doing multiple hypothesis tests on these data
# 4 correlated antigens for the eidolon
.001/length(unique(seroprev.sum.sex$antigen[seroprev.sum.sex$bat_species=="Eidolon dupreanum"])) #0.00025 threshold for significance for Eidolon
# 3 for the rousettus
.001/length(unique(seroprev.sum.sex$antigen[seroprev.sum.sex$bat_species=="Rousettus madagascariensis"]))#0.00025 threshold for significance for Rousettus


# you want to test whether including sex or antigen matters in driving seroprevalence
seroprev.sum.sex$antigen <- as.factor(seroprev.sum.sex$antigen)
seroprev.sum.sex$bat_sex <- as.factor(seroprev.sum.sex$bat_sex)
seroprev.sum.sex$roost_site <- as.factor(seroprev.sum.sex$roost_site)

# prior work suggests important sex-based differences, so we will 
# model sex interactions using the "by" term and as a fixed effect
# we will also model all antigens together but allow for antigen-specific 
# effects of sex and site

# research question to be answered: is there sex-specific seasonality in
# the serological data for each antigen?

# for each species and antigen, we first test: 
# is there seasonality in seroprevalence by antigen without accounting for sex?
# we can include both sites (Ankarana and Moramanga) but allow the model to drop 
# terms internally if one dataset (e.g. Ankarana) is not robust enough.

seroprev.sum$year <- year(seroprev.sum$mid_date)
seroprev.sum$doy <- yday(seroprev.sum$mid_date)
seroprev.sum$antigen <- as.factor(seroprev.sum$antigen)
seroprev.sum$roost_site <- as.factor(seroprev.sum$roost_site)
seroprev.sum$year_f <- as.factor(seroprev.sum$year)


#GhV
Eid.sum.GhV <- subset(seroprev.sum, bat_species=="Eidolon dupreanum" & antigen=="GhV")
Eid.sum.GhV$month <- month(Eid.sum.GhV$mid_date)
Eid.sum.GhV$roost_site <- droplevels(Eid.sum.GhV$roost_site)
Eid.sum.GhV$antigen <- droplevels(Eid.sum.GhV$antigen)
Eid.sum.GhV$year_f <- droplevels(Eid.sum.GhV$year_f)


gamEidNoSexGhV <- gam(cbind(seropos, seroneg) ~ 
                         roost_site + 
                         s(as.numeric(doy), bs = "cc", k =10, by = roost_site) +
                         s(as.numeric(mid_date),  bs="tp", k=6, by= roost_site),
                       data= Eid.sum.GhV, 
                       family = "binomial", 
                       method = "REML", 
                       gamma = 1,
                       select = TRUE) #



summary(gamEidNoSexGhV) 
#plot(gamEidNoSexGhV)
# sig seasonality in Moramanga inter and intra 

#now try with sex too
seroprev.sum.sex$doy <- yday(seroprev.sum.sex$mid_date)
seroprev.sum.sex$year <- year(seroprev.sum.sex$mid_date)

Eid.sex.GhV = subset(seroprev.sum.sex, bat_species=="Eidolon dupreanum" & antigen=="GhV")
Eid.sex.GhV$year_f <- as.factor(Eid.sex.GhV$year)
Eid.sex.GhV$antigen <- droplevels(Eid.sex.GhV$antigen)
Eid.sex.GhV$roost_site <- droplevels(Eid.sex.GhV$roost_site)
Eid.sex.GhV$bat_sex <- droplevels(Eid.sex.GhV$bat_sex)

group <- interaction( Eid.sex.GhV$bat_sex, Eid.sex.GhV$roost_site, drop = TRUE)

Eid.sex.GhV$group <- group
Eid.sex.GhV$year <- as.numeric(as.character(Eid.sex.GhV$year))

gamEidSexGhV <- gam(cbind(seropos, seroneg)~ 
                   roost_site + bat_sex +
                   s(as.numeric(doy), bs = c("cc"), k = 10, by = group) +
                   s(as.numeric(mid_date), bs="tp", k=6, by=group),
                 data=Eid.sex.GhV, 
                 family = "binomial", 
                 method = "REML", 
                 gamma =1,
                 select = TRUE) #this term allows non-impactful terms to essentially be removed from model - check edf in summary to see if this happened to any
summary(gamEidSexGhV)
# significant trends female Moramanga interannual, male moramanga intra borderline

AIC(gamEidNoSexGhV, gamEidSexGhV) # a lot better with no sex

#now try for AngV as well
Eid.sum.AngV <- subset(seroprev.sum, bat_species=="Eidolon dupreanum" & antigen=="AngV")
# 
Eid.sum.AngV$year <- year(Eid.sum.AngV$mid_date)
Eid.sum.AngV$year_f <- as.factor(Eid.sum.AngV$year)
Eid.sum.AngV$month <- month(Eid.sum.AngV$mid_date)
Eid.sum.AngV$roost_site <- droplevels(Eid.sum.AngV$roost_site)
Eid.sum.AngV$antigen <- droplevels(Eid.sum.AngV$antigen)
Eid.sum.AngV$year_f <- droplevels(Eid.sum.AngV$year_f)

gamEidNoSexAngV <- gam(cbind(seropos, seroneg) ~ 
                        roost_site + 
                        s(as.numeric(doy), bs = "cc", k =10, by = roost_site) +
                         s(as.numeric(mid_date), bs="tp", k=6, by= roost_site),
                      data= Eid.sum.AngV, 
                      family = "binomial", 
                      method = "REML", 
                      gamma = 1,
                      select = TRUE) #
summary(gamEidNoSexAngV) #Moramanga inter boderline


gamEidNoSexAngVPlot <- gam(cbind(seropos, seroneg) ~ 
                         roost_site + 
                         s(as.numeric(doy), bs = "cc", k =10, by = roost_site) +
                         s(as.numeric(mid_date), bs="tp", k=6, by= roost_site),
                       data= Eid.sum.AngV, 
                       family = "binomial", 
                       method = "REML", 
                       gamma = .48,
                       select = TRUE) #

summary(gamEidNoSexAngV) #Moramanga interannual boderline
#plot(gamEidNoSexAngV)

#now try with sex too

seroprev.sum.sex$doy <- yday(seroprev.sum.sex$mid_date)
seroprev.sum.sex$year <- year(seroprev.sum.sex$mid_date)

Eid.sex.AngV = subset(seroprev.sum.sex, bat_species=="Eidolon dupreanum" & antigen=="AngV")
Eid.sex.AngV$antigen <- droplevels(Eid.sex.AngV$antigen)
Eid.sex.AngV$roost_site <- droplevels(Eid.sex.AngV$roost_site)
Eid.sex.AngV$bat_sex <- droplevels(Eid.sex.AngV$bat_sex)

group <- interaction( Eid.sex.AngV$bat_sex, Eid.sex.AngV$roost_site, drop = TRUE)

Eid.sex.AngV$group <- group
Eid.sex.AngV$year_f <- as.factor(Eid.sex.AngV$year)
Eid.sex.AngV$year <- as.numeric(as.character(Eid.sex.AngV$year))

gamEidSexAngV <- gam(cbind(seropos, seroneg)~ 
                      roost_site + bat_sex +
                      s(as.numeric(doy), bs = c("cc"), k = 10, by = group) +
                      s(as.numeric(mid_date), bs="tp", k=6, by=group),
                    data=Eid.sex.AngV, 
                    family = "binomial", 
                    method = "REML", 
                    gamma =1,
                    select = TRUE) #this term allows non-impactful terms to essentially be removed from model - check edf in summary to see if this happened to any
summary(gamEidSexAngV) #Ankarana males inter boderline

AIC(gamEidNoSexAngV, gamEidSexAngV) # better with no sex


#and HeV
Eid.sum.HeV <- subset(seroprev.sum, bat_species=="Eidolon dupreanum" & antigen=="HeV")
# 
Eid.sum.HeV$year <- year(Eid.sum.HeV$mid_date)
Eid.sum.HeV$year_f <- as.factor(Eid.sum.HeV$year)
Eid.sum.HeV$month <- month(Eid.sum.HeV$mid_date)
Eid.sum.HeV$roost_site <- droplevels(Eid.sum.HeV$roost_site)
Eid.sum.HeV$antigen <- droplevels(Eid.sum.HeV$antigen)
Eid.sum.HeV$year_f <- droplevels(Eid.sum.HeV$year_f)

gamEidNoSexHeV <- gam(cbind(seropos, seroneg) ~ 
                        roost_site + 
                        s(as.numeric(doy), bs = "cc", k =10, by = roost_site) +
                        s(as.numeric(mid_date),  bs="tp", k=6, by= roost_site),
                      data= Eid.sum.HeV, 
                      family = "binomial", 
                      method = "REML", 
                      gamma = 1,
                      select = TRUE) #


gamEidNoSexHeVPlot <- gam(cbind(seropos, seroneg) ~ 
                            roost_site + 
                            s(as.numeric(doy), bs = "cc", k =10, by = roost_site) +
                            s(as.numeric(mid_date),  bs="tp", k=6, by= roost_site),
                          data= Eid.sum.HeV, 
                          family = "binomial", 
                          method = "REML", 
                          gamma = .22,
                          select = TRUE) #

summary(gamEidNoSexHeV) #nothing sig
#plot(gamEidNoSexHeV)


#now try with sex too

seroprev.sum.sex$doy <- yday(seroprev.sum.sex$mid_date)
seroprev.sum.sex$year <- year(seroprev.sum.sex$mid_date)

Eid.sex.HeV = subset(seroprev.sum.sex, bat_species=="Eidolon dupreanum" & antigen=="HeV")
Eid.sex.HeV$antigen <- droplevels(Eid.sex.HeV$antigen)
Eid.sex.HeV$roost_site <- droplevels(Eid.sex.HeV$roost_site)
Eid.sex.HeV$bat_sex <- droplevels(Eid.sex.HeV$bat_sex)

group <- interaction( Eid.sex.HeV$bat_sex, Eid.sex.HeV$roost_site, drop = TRUE)

Eid.sex.HeV$group <- group
Eid.sex.HeV$year_f <- as.factor(Eid.sex.HeV$year)
Eid.sex.HeV$year <- as.numeric(as.character(Eid.sex.HeV$year))

gamEidSexHeV <- gam(cbind(seropos, seroneg)~ 
                      roost_site + bat_sex +
                      s(as.numeric(doy), bs = c("cc"), k = 10, by = group) +
                      s(as.numeric(mid_date), bs="tp", k=6, by=group),
                    data=Eid.sex.HeV, 
                    family = "binomial", 
                    method = "REML", 
                    gamma =1,
                    select = TRUE) #this term allows non-impactful terms to essentially be removed from model - check edf in summary to see if this happened to any
summary(gamEidSexHeV)
# significant interannual trends for female Moramanga 
# most Ankarana  data reduced to edf = 0, especial interannual


AIC(gamEidNoSexHeV, gamEidSexHeV) # better with no sex


#and NiV
#now try for NiV as well
Eid.sum.NiV <- subset(seroprev.sum, bat_species=="Eidolon dupreanum" & antigen=="NiV")
# 
Eid.sum.NiV$year <- year(Eid.sum.NiV$mid_date)
Eid.sum.NiV$year_f <- as.factor(Eid.sum.NiV$year)
Eid.sum.NiV$month <- month(Eid.sum.NiV$mid_date)
Eid.sum.NiV$roost_site <- droplevels(Eid.sum.NiV$roost_site)
Eid.sum.NiV$antigen <- droplevels(Eid.sum.NiV$antigen)
Eid.sum.NiV$year_f <- droplevels(Eid.sum.NiV$year_f)

gamEidNoSexNiV <- gam(cbind(seropos, seroneg) ~ 
                        roost_site + 
                        s(as.numeric(doy), bs = "cc", k =10, by = roost_site) +
                        s(as.numeric(mid_date),  bs="tp", k=6, by= roost_site),
                      data= Eid.sum.NiV, 
                      family = "binomial", 
                      method = "REML", 
                      gamma = 1,
                      select = TRUE) #


gamEidNoSexNiVPlot <- gam(cbind(seropos, seroneg) ~ 
                            roost_site + 
                            s(as.numeric(doy), bs = "cc", k =10, by = roost_site) +
                            s(as.numeric(mid_date),  bs="tp", k=6, by= roost_site),
                          data= Eid.sum.NiV, 
                          family = "binomial", 
                          method = "REML", 
                          gamma = .97,
                          select = TRUE) #

summary(gamEidNoSexNiV) #nothing sig
#plot(gamEidNoSexNiV)
# Ankarana edf = 0


#now try with sex too

seroprev.sum.sex$doy <- yday(seroprev.sum.sex$mid_date)
seroprev.sum.sex$year <- year(seroprev.sum.sex$mid_date)

Eid.sex.NiV = subset(seroprev.sum.sex, bat_species=="Eidolon dupreanum" & antigen=="NiV")
Eid.sex.NiV$antigen <- droplevels(Eid.sex.NiV$antigen)
Eid.sex.NiV$roost_site <- droplevels(Eid.sex.NiV$roost_site)
Eid.sex.NiV$bat_sex <- droplevels(Eid.sex.NiV$bat_sex)

group <- interaction( Eid.sex.NiV$bat_sex, Eid.sex.NiV$roost_site, drop = TRUE)

Eid.sex.NiV$group <- group
Eid.sex.NiV$year_f <- as.factor(Eid.sex.NiV$year)
Eid.sex.NiV$year <- as.numeric(as.character(Eid.sex.NiV$year))

gamEidSexNiV <- gam(cbind(seropos, seroneg)~ 
                      roost_site + bat_sex +
                      s(as.numeric(doy), bs = c("cc"), k = 10, by = group) +
                      s(as.numeric(mid_date), bs="tp", k=6, by=group),
                    data=Eid.sex.NiV, 
                    family = "binomial", 
                    method = "REML", 
                    gamma =1,
                    select = TRUE) #this term allows non-impactful terms to essentially be removed from model - check edf in summary to see if this happened to any
summary(gamEidSexNiV)
#nothing sig
# most Ankarana  data reduced to edf = 0, especial interannual


AIC(gamEidNoSexNiV, gamEidSexNiV) # better with no sex





#and predict just within a year
pred.yr <- cbind.data.frame(doy= 1:365, mid_date = seq(from = as.Date("2019-01-01"), to= as.Date("2019-12-31"), by=1))
head(pred.yr)
pred.yr$antigen <- "GhV"
pred.yr$year_f <- "2019"
pred.yr$roost_site <- "Moramanga"
pred.yr2 <- pred.yr
pred.yr3 <- pred.yr
pred.yr4 <- pred.yr
pred.yr2$antigen <- "AngV"
pred.yr3$antigen <- "HeV"
pred.yr4$antigen <- "NiV"
pred.yr <- rbind(pred.yr, pred.yr2, pred.yr3, pred.yr4)


pred.yr$seroprevalence <-NA
pred.yr$seroprev_uci <- NA
pred.yr$seroprev_lci <- NA

pred.yr$seroprevalence[pred.yr$antigen=="GhV"] <- predict(gamEidNoSexGhV, type = "response", newdata = pred.yr[pred.yr$antigen=="GhV",], exclude = "s(mid_date)")
pred.yr$seroprev_uci[pred.yr$antigen=="GhV"] <- predict(gamEidNoSexGhV, type = "response", newdata = pred.yr[pred.yr$antigen=="GhV",], exclude = "s(mid_date)") + 1.96*predict(gamEidNoSexGhV, type = "response", newdata = pred.yr[pred.yr$antigen=="GhV",], exclude = "s(mid_date)", se.fit = T)$se
pred.yr$seroprev_lci[pred.yr$antigen=="GhV"] <- predict(gamEidNoSexGhV, type = "response", newdata = pred.yr[pred.yr$antigen=="GhV",], exclude = "s(mid_date)") - 1.96*predict(gamEidNoSexGhV, type = "response", newdata = pred.yr[pred.yr$antigen=="GhV",], exclude = "s(mid_date)", se.fit = T)$se
pred.yr$bat_species <- "Eidolon dupreanum"


pred.yr$seroprevalence[pred.yr$antigen=="AngV"] <- predict(gamEidNoSexAngVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="AngV",], exclude = "s(mid_date)")
pred.yr$seroprev_uci[pred.yr$antigen=="AngV"] <- predict(gamEidNoSexAngVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="AngV",], exclude = "s(mid_date)") + 1.96*predict(gamEidNoSexAngVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="AngV",], exclude = "s(mid_date)", se.fit = T)$se
pred.yr$seroprev_lci[pred.yr$antigen=="AngV"] <- predict(gamEidNoSexAngVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="AngV",], exclude = "s(mid_date)") - 1.96*predict(gamEidNoSexAngVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="AngV",], exclude = "s(mid_date)", se.fit = T)$se



pred.yr$seroprevalence[pred.yr$antigen=="HeV"] <- predict(gamEidNoSexHeVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="HeV",], exclude = "s(mid_date)")
pred.yr$seroprev_uci[pred.yr$antigen=="HeV"] <- predict(gamEidNoSexHeVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="HeV",], exclude = "s(mid_date)") + 1.96*predict(gamEidNoSexHeVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="HeV",], exclude = "s(mid_date)", se.fit = T)$se
pred.yr$seroprev_lci[pred.yr$antigen=="HeV"] <- predict(gamEidNoSexHeVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="HeV",], exclude = "s(mid_date)") - 1.96*predict(gamEidNoSexHeVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="HeV",], exclude = "s(mid_date)", se.fit = T)$se



pred.yr$seroprevalence[pred.yr$antigen=="NiV"] <- predict(gamEidNoSexNiVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="NiV",], exclude = "s(mid_date)")
pred.yr$seroprev_uci[pred.yr$antigen=="NiV"] <- predict(gamEidNoSexNiVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="NiV",], exclude = "s(mid_date)") + 1.96*predict(gamEidNoSexNiVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="NiV",], exclude = "s(mid_date)", se.fit = T)$se
pred.yr$seroprev_lci[pred.yr$antigen=="NiV"] <- predict(gamEidNoSexNiVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="NiV",], exclude = "s(mid_date)") - 1.96*predict(gamEidNoSexNiVPlot, type = "response", newdata = pred.yr[pred.yr$antigen=="NiV",], exclude = "s(mid_date)", se.fit = T)$se



pred.yr$seroprev_lci[pred.yr$seroprev_lci<0] <- 0
pred.yr$seroprev_lci[pred.yr$seroprev_uci>1] <- 1

cal.dat$year <- year(cal.dat$stop)
cal.sub = subset(cal.dat, bat_species =="Eidolon dupreanum" & year==2019)
cal.sub$doy_start <- yday(cal.sub$start)
cal.sub$doy_stop <- yday(cal.sub$stop)
cal.sub$year <- year(cal.sub$stop)
cal.sub$doy_stop[cal.sub$season=="lactation"] <- 365
cal.sub <- rbind(cal.sub, cal.sub[nrow(cal.sub),])
cal.sub$doy_start[4] <- 1
cal.sub$doy_stop[4] <- yday("2019-02-02")



filz <- c(colorz, 'gestation' = "yellow", 'lactation' = "greenyellow", 'winter' = "skyblue")


filz.season <- filz[(length(filz)-2):length(filz)]
filz.antigen <- filz[1:(length(filz)-2)]


pIntra <- ggplot(data=subset(seroprev.sum, roost_site=="Moramanga" & bat_species =="Eidolon dupreanum")) +
  geom_rect(data=cal.sub, aes(xmin=as.Date(doy_start), xmax=as.Date(doy_stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=doy, y=seroprev, size=N, color=antigen, shape =year_f)) + 
  scale_color_manual(values=colorz) + 
  scale_fill_manual(values=filz.season) +
  ggnewscale::new_scale_fill() +
  geom_linerange(aes(x=doy, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen), linetype=2, alpha=.3, show.legend = F) +
  geom_ribbon(data=subset(pred.yr), aes(x=doy, ymin= seroprev_lci, ymax=seroprev_uci, fill=antigen), alpha=.3, show.legend = F) +
  scale_fill_manual(values=filz.antigen) +
  geom_line(data=subset(pred.yr), aes(x=doy, y= seroprevalence, color=antigen)) +
  scale_shape_discrete(name="year") +
  facet_grid(~bat_species) + theme_bw() +
  theme(axis.text.x = element_text(size=14), 
        axis.title.x = element_text(size=16),
        axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),
        strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  coord_cartesian(expand = F)  + xlab("day of year") + guides(alpha="none")
pIntra





# And do the same for Rousettus
# "GhV" only

Rou.sum.GhV <- subset(seroprev.sum, bat_species=="Rousettus madagascariensis" & antigen=="GhV")
Rou.sum.GhV$antigen <- droplevels(Rou.sum.GhV$antigen)
Rou.sum.GhV$roost_site <- droplevels(Rou.sum.GhV$roost_site)
Rou.sum.GhV$year_f <- as.factor(Rou.sum.GhV$year)



gamRouNoSexGhV <- gam(cbind(seropos, seroneg) ~ roost_site +
                     s(as.numeric(doy), k=10, bs="cc", by = roost_site) +
                     s(as.numeric(mid_date), k=6, bs="fs", by = roost_site), 
                     data=Rou.sum.GhV, 
                     family = "binomial", 
                     method = "REML", 
                     gamma = 1,
                     select = TRUE) #this term allows non-impactful terms to essentially be removed from model - check edf in summary to see if this happened to any
summary(gamRouNoSexGhV)  
# no sig

#and with sex
Rou.sex.GhV = subset(seroprev.sum.sex, bat_species=="Rousettus madagascariensis" & antigen=="GhV" )
Rou.sex.GhV$antigen <- droplevels(Rou.sex.GhV$antigen)
Rou.sex.GhV$roost_site <- droplevels(Rou.sex.GhV$roost_site)
Rou.sex.GhV$bat_sex <- droplevels(Rou.sex.GhV$bat_sex)
Rou.sex.GhV$year_f <- as.factor(Rou.sex.GhV$year)

group4 <- interaction(Rou.sex.GhV$roost_site, Rou.sex.GhV$bat_sex)

gamRouSexGhV <- gam(cbind(seropos, seroneg)~ roost_site + bat_sex +
                   s(as.numeric(doy),  k=10, bs="cc", by = group4) +
                   s(as.numeric(mid_date), k=6, bs="fs", by = group4), 
                   data=Rou.sex.GhV, 
                   family = "binomial", 
                   method = "REML", 
                   gamma =1,
                   select = TRUE) #this term allows non-impactful terms to essentially be removed from model - check edf in summary to see if this happened to any
summary(gamRouSexGhV) 
# no sig

# compare model with and without sex
AIC(gamRouNoSexGhV, gamRouSexGhV) # simpler model without sex is still more parsimonious


#use above to predict and overlay on the data
#and predict just within a year
pred.yr.rou <- cbind.data.frame(doy= 1:365, mid_date = seq(from = as.Date("2019-01-01"), to= as.Date("2019-12-31"), by=1))
head(pred.yr.rou)
pred.yr.rou$antigen <- "GhV"
pred.yr.rou$year_f <- "2019"
pred.yr.rou$roost_site <- "Moramanga"

pred.yr.rou$seroprevalence <-NA
pred.yr.rou$seroprev_uci <- NA
pred.yr.rou$seroprev_lci <- NA

pred.yr.rou$seroprevalence[pred.yr.rou$antigen=="GhV"] <- predict(gamRouNoSexGhV, type = "response", newdata = pred.yr.rou[pred.yr.rou$antigen=="GhV",], exclude = "s(mid_date)")
pred.yr.rou$seroprev_uci[pred.yr.rou$antigen=="GhV"] <- predict(gamRouNoSexGhV, type = "response", newdata = pred.yr.rou[pred.yr.rou$antigen=="GhV",], exclude = "s(mid_date)") + 1.96*predict(gamRouNoSexGhV, type = "response", newdata = pred.yr.rou[pred.yr.rou$antigen=="GhV",], exclude = "s(mid_date)", se.fit = T)$se
pred.yr.rou$seroprev_lci[pred.yr.rou$antigen=="GhV"] <- predict(gamRouNoSexGhV, type = "response", newdata = pred.yr.rou[pred.yr.rou$antigen=="GhV",], exclude = "s(mid_date)") - 1.96*predict(gamRouNoSexGhV, type = "response", newdata = pred.yr.rou[pred.yr.rou$antigen=="GhV",], exclude = "s(mid_date)", se.fit = T)$se
pred.yr.rou$bat_species <- "Rousettus madagascariensis"

rou.sub = subset(cal.dat, bat_species=="Rousettus madagascariensis" & year==2019)
rou.sub$doy_start <- yday(rou.sub$start)
rou.sub$doy_stop <- yday(rou.sub$stop)

p7 <- ggplot(data=subset(seroprev.sum, roost_site=="Moramanga" & bat_species =="Rousettus madagascariensis")) +
  geom_rect(data=rou.sub, aes(xmin=as.Date(doy_start), xmax=as.Date(doy_stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=doy, y=seroprev, size=N, color=antigen, shape =year_f)) + scale_color_manual(values=colorz) + scale_fill_manual(values=filz) +
  geom_linerange(aes(x=doy, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen), linetype=2, alpha=.3, show.legend = F) +
  geom_ribbon(data=subset(pred.yr.rou), aes(x=doy, ymin= seroprev_lci, ymax=seroprev_uci, fill=antigen), alpha=.3) +
  geom_line(data=subset(pred.yr.rou), aes(x=doy, y= seroprevalence, color=antigen)) +
  facet_grid(~bat_species) + theme_bw()
p7






 ################################################################################
 ################################################################################
 ## Section 3. Add the historical IPM time series for comparison
 
 # load the IPM data and the JAE data
 
IPM.dat <- read.csv(file = paste0(batwd,"/Serological-Data/IPM-data-2006-2007.csv"),header = T, stringsAsFactors = F)
head(IPM.dat) 

# and do the same male/female specific serology for both for the antigens of 
# interest: HeV, NiV 

unique(IPM.dat$prev)
IPM.dat = subset(IPM.dat, prev=="0" | prev=="1")
unique(IPM.dat$species)
unique(IPM.dat$type)
IPM.dat = subset(IPM.dat, type=="ELISA_HeV" | type =="ELISA_NiV")

sort(unique(IPM.dat$date)) #each of these is a discrete sampling event
unique(IPM.dat$site)#"AngavoKely" "AngavoBe"  
subset(IPM.dat, date=="2007-04-05" ) #angavokely
subset(IPM.dat, date=="2007-04-11" ) #angavobe

#summarise to seropositivity
IPM.dat$prev <- as.numeric(IPM.dat$prev)
IPM.sero <- ddply(IPM.dat,.(type, date, sex),summarise, N_seropos=sum(prev), N_tot=length(prev))
IPM.sero.nosex <- ddply(IPM.dat,.(type, date),summarise, N_seropos=sum(prev), N_tot=length(prev))
IPM.sero$N_seroneg = IPM.sero$N_tot-IPM.sero$N_seropos
IPM.sero$seroprev <- IPM.sero$N_seropos/IPM.sero$N_tot
IPM.sero$seroprev_lci <-binom.exact(x=IPM.sero$N_seropos, n=IPM.sero$N_tot)$lower
IPM.sero$seroprev_uci <-binom.exact(x=IPM.sero$N_seropos, n=IPM.sero$N_tot)$upper
IPM.sero$date <- as.Date(IPM.sero$date)
head(IPM.sero)
unique(IPM.sero$type)
IPM.sero$type[IPM.sero$type=="ELISA_HeV"] <- "HeV"
IPM.sero$type[IPM.sero$type=="ELISA_NiV"] <- "NiV"
names(IPM.sero)[names(IPM.sero)=="type"] <- "antigen"
IPM.sero$antigen <- as.factor(IPM.sero$antigen)
unique(IPM.sero$sex)

IPM.sero$sex[IPM.sero$sex=="F"] <- "female"
IPM.sero$sex[IPM.sero$sex=="M"] <- "male"

IPM.sero$sex <- as.factor(IPM.sero$sex)


IPM.sero.nosex$N_seroneg = IPM.sero.nosex$N_tot-IPM.sero.nosex$N_seropos
IPM.sero.nosex$seroprev <- IPM.sero.nosex$N_seropos/IPM.sero.nosex$N_tot
IPM.sero.nosex$seroprev_lci <-binom.exact(x=IPM.sero.nosex$N_seropos, n=IPM.sero.nosex$N_tot)$lower
IPM.sero.nosex$seroprev_uci <-binom.exact(x=IPM.sero.nosex$N_seropos, n=IPM.sero.nosex$N_tot)$upper
IPM.sero.nosex$date <- as.Date(IPM.sero.nosex$date)
head(IPM.sero.nosex)
unique(IPM.sero.nosex$type)
IPM.sero.nosex$type[IPM.sero.nosex$type=="ELISA_HeV"] <- "HeV"
IPM.sero.nosex$type[IPM.sero.nosex$type=="ELISA_NiV"] <- "NiV"
names(IPM.sero.nosex)[names(IPM.sero.nosex)=="type"] <- "antigen"
IPM.sero.nosex$antigen <- as.factor(IPM.sero.nosex$antigen)
IPM.sero.nosex$year <- year(IPM.sero.nosex$date)

IPM.sero.nosex$year_f <-as.factor(IPM.sero.nosex$year)

IPM.sero.nosex$doy <- yday(IPM.sero.nosex$date)
IPM.sero$doy <- yday(IPM.sero$date)
IPM.sero$year <- year(IPM.sero$date)
IPM.sero$year_f <- as.factor(IPM.sero$year)

IPM.sero = subset(IPM.sero, antigen=="NiV")
IPM.sero.nosex = subset(IPM.sero.nosex, antigen=="NiV")

# first, is there seasoyear()# first, is there seasonality without accounting for sex?
gamIPMNoSex <- gam(cbind(N_seropos, N_seroneg)~ 
                    s(as.numeric(doy), k=10, bs="cc") +
                     s(as.numeric(date), k=6, bs="tp"),
                     data=IPM.sero.nosex, family = "binomial", select = TRUE, 
                    method = "REML")
summary(gamIPMNoSex) # yes, very significant for both antigens and both timescales


# now, with sex
gamIPMSex <- gam(cbind(N_seropos, N_seroneg)~sex +
                   s(as.numeric(doy), k=10, bs="cc", by=sex) +
                   s(as.numeric(date), k=6, bs="tp", by = sex),
                   data=IPM.sero, family = "binomial", select = TRUE,
                  method = "REML")
summary(gamIPMSex) # yes, very significant for both antigens and both timescales
#plot(gamIPMSex)
#intraanual for both sexes
# interannual for males and borderline for females

AIC(gamIPMNoSex, gamIPMSex) # simpler, no sex model prevails


#and predict just within a year
pred.yr.IPM <- cbind.data.frame(doy= 1:365, date=seq(from = as.Date("2006-01-01"), to=as.Date("2006-12-31"), by=1))
head(pred.yr.IPM)
pred.yr.IPM$antigen <- "NiV"
pred.yr.IPM$year_f <- "2006"
pred.yr.IPM$roost_site <- "Moramanga"
# pred.yr.IPM2 <- pred.yr.IPM
# pred.yr.IPM2$antigen <- "HeV"
# pred.yr.IPM <- rbind(pred.yr.IPM, pred.yr.IPM2)

pred.yr.IPM$seroprevalence <-NA
pred.yr.IPM$seroprev_uci <- NA
pred.yr.IPM$seroprev_lci <- NA

pred.yr.IPM$seroprevalence[pred.yr.IPM$antigen=="NiV"] <- predict(gamIPMNoSex, type = "response", newdata = pred.yr.IPM[pred.yr.IPM$antigen=="NiV",], exclude = "s(date)")
pred.yr.IPM$seroprev_uci[pred.yr.IPM$antigen=="NiV"] <- predict(gamIPMNoSex, type = "response", newdata = pred.yr.IPM[pred.yr.IPM$antigen=="NiV",], exclude = "s(date)") + 1.96*predict(gamIPMNoSex, type = "response", newdata = pred.yr.IPM[pred.yr.IPM$antigen=="NiV",], exclude = "s(date)", se.fit = T)$se
pred.yr.IPM$seroprev_lci[pred.yr.IPM$antigen=="NiV"] <- predict(gamIPMNoSex, type = "response", newdata = pred.yr.IPM[pred.yr.IPM$antigen=="NiV",], exclude = "s(date)") - 1.96*predict(gamIPMNoSex, type = "response", newdata = pred.yr.IPM[pred.yr.IPM$antigen=="NiV",], exclude = "s(date)", se.fit = T)$se
pred.yr.IPM$bat_species <- "Eidolon dupreanum"


# pred.yr.IPM$seroprevalence[pred.yr.IPM$antigen=="HeV"] <- predict(gamIPMNoSex, type = "response", newdata = pred.yr.IPM[pred.yr.IPM$antigen=="HeV",], exclude = "s(year_corrected_f)")
# pred.yr.IPM$seroprev_uci[pred.yr.IPM$antigen=="HeV"] <- predict(gamIPMNoSex, type = "response", newdata = pred.yr.IPM[pred.yr.IPM$antigen=="HeV",], exclude = "s(year_corrected_f)") + 1.96*predict(gamIPMNoSex, type = "response", newdata = pred.yr.IPM[pred.yr.IPM$antigen=="HeV",], exclude = "s(year_corrected_f)", se.fit = T)$se
# pred.yr.IPM$seroprev_lci[pred.yr.IPM$antigen=="HeV"] <- predict(gamIPMNoSex, type = "response", newdata = pred.yr.IPM[pred.yr.IPM$antigen=="HeV",], exclude = "s(year_corrected_f)") - 1.96*predict(gamIPMNoSex, type = "response", newdata = pred.yr.IPM[pred.yr.IPM$antigen=="HeV",], exclude = "s(year_corrected_f)", se.fit = T)$se


IPM.sero$bat_species=="Eidolon dupreanum"
IPM.sero$sig <- 0
IPM.sero$sig[IPM.sero$antigen=="NiV"] <- 1
IPM.sero$sig[IPM.sero$antigen=="HeV"] <- 0
IPM.sero$sig <- as.factor(IPM.sero$sig)

alphaz = c('0' = 0, '1' = 1, '2' = 0.5)

filz1 <- filz[1:(length(filz)-3)]
filz2 <- filz[(length(filz)-2):length(filz)]

IPM.sero.nosex$sig <- 1
IPM.sero.nosex$sig[IPM.sero.nosex$antigen=="HeV"] <- 2
IPM.sero.nosex$sig <- as.factor(IPM.sero.nosex$sig)

#with legend
pIPMLeg <- ggplot(data=IPM.sero.nosex) +
  scale_color_manual(values=colorz) +
  scale_fill_manual(values=filz2) +
  scale_size_continuous(name="N bats", breaks = c(10,20,30)) +
  geom_rect(data=cal.sub, aes(xmin=as.Date(doy_start), xmax=as.Date(doy_stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=doy, y=seroprev, size=N_tot, color=antigen, alpha = sig, shape=year_f)) +
  #geom_line(aes(x=date, y=seroprev, color=antigen, alpha = sig), linetype=2) +
  geom_linerange(aes(x=doy, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen, alpha=sig, group=year_f), linetype=2, show.legend = F) +
  scale_alpha_manual(values=alphaz) +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values=filz1) +
  geom_ribbon(data=subset(pred.yr.IPM, antigen=="HeV"), aes(x=doy, ymin= seroprev_lci, ymax=seroprev_uci, fill=antigen), alpha=.03, show.legend = F) +
  geom_line(data=subset(pred.yr.IPM, antigen=="HeV"), aes(x=doy, y= seroprevalence, color=antigen), alpha=.3) +
  geom_ribbon(data=subset(pred.yr.IPM, antigen=="NiV"), aes(x=doy, ymin= seroprev_lci, ymax=seroprev_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.yr.IPM, antigen=="NiV"), aes(x=doy, y= seroprevalence, color=antigen)) +
  facet_grid(~bat_species) + theme_bw()+
  theme(axis.text = element_text(size=14), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16), strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines")) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + guides(alpha="none"); pIPMLeg



#strip
pIPM <- ggplot(data=IPM.sero.nosex) +
  scale_color_manual(values=colorz) +
  scale_fill_manual(values=filz2) +
  scale_size_continuous(name="N bats", breaks = c(10,20,30)) +
  geom_rect(data=cal.sub, aes(xmin=as.Date(doy_start), xmax=as.Date(doy_stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=doy, y=seroprev, size=N_tot, color=antigen, alpha = sig, shape=year_f)) +
  #geom_line(aes(x=date, y=seroprev, color=antigen, alpha = sig), linetype=2) +
  geom_linerange(aes(x=doy, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen, alpha=sig, group=year_f), linetype=2, show.legend = F) +
  scale_alpha_manual(values=alphaz) +
  guides(fill="none") +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values=filz1) +
  geom_ribbon(data=subset(pred.yr.IPM, antigen=="HeV"), aes(x=doy, ymin= seroprev_lci, ymax=seroprev_uci, fill=antigen), alpha=.03, show.legend = F) +
  geom_line(data=subset(pred.yr.IPM, antigen=="HeV"), aes(x=doy, y= seroprevalence, color=antigen), alpha=.3) +
  geom_ribbon(data=subset(pred.yr.IPM, antigen=="NiV"), aes(x=doy, ymin= seroprev_lci, ymax=seroprev_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.yr.IPM, antigen=="NiV"), aes(x=doy, y= seroprevalence, color=antigen)) +
  facet_grid(~bat_species) + theme_bw()+
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=16), strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines")) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + xlab("day of year") +
  guides(alpha="none", color="none", shape="none", size="none", fill="none"); pIPM


################################################################################
################################################################################
## Section 4. Add the historical JAE time series for comparison

JAE.dat <-  read.csv(file = paste0(homewd,"/serology/data/mix_2019_lum_dat.csv"),header = T, stringsAsFactors = F)

# new.dat <-  read.csv(file = paste0(homewd,"/data/mix_2019_lum_dat.csv"),header = T, stringsAsFactors = F)
head(JAE.dat)


head(catch.dat)
catch.merge <- dplyr::select(catch.dat,sampleid, processing_date, bat_sex, roost_site, sampling_session)

JAE.long <- merge(JAE.dat, catch.merge, by="sampleid", all.x = T)
head(JAE.long)
JAE.long$processing_date <- as.Date(JAE.long$processing_date)
unique(JAE.long$antigen)
unique(JAE.long$roost_site)
JAE.long$roost_site[is.na(JAE.long$roost_site)] <- "Moramanga"
JAE.long$roost_site[JAE.long$roost_site=="Ankarana_Canyon_Edup"] <- "Ankarana"
JAE.long$roost_site[JAE.long$roost_site=="Mahabo_Edup"] <- "Mahabo"
JAE.long$roost_site[JAE.long$roost_site=="Ankarana_Canyon_Rmad"] <- "Ankarana"
JAE.long$roost_site[JAE.long$roost_site=="Makira_Rmad"] <- "Makira"
JAE.long$roost_site[JAE.long$roost_site!="Makira" &JAE.long$roost_site!="Ankarana" & JAE.long$roost_site!="Mahabo" ] <- "Moramanga"

#get sampling session
sample.df <- ddply(JAE.long,.(sampling_session), summarise, min_date = min(processing_date), max_date =max(processing_date))
sample.df <- sample.df[complete.cases(sample.df),]
sample.df$mid_date <- ((sample.df$max_date-sample.df$min_date)/2) + sample.df$min_date

JAE.long <- merge(JAE.long, sample.df, by="sampling_session")

#now, sub to just NiVin Eidolon
JAE.long$plot <- 0
JAE.long$plot[JAE.long$bat_species=="Eidolon dupreanum" & JAE.long$antigen == "NiV"]<- 1


#and do sex- and antigen-specific GAMS
JAE.merge = subset(JAE.long, plot==1)
head(JAE.merge)

###here better compile data

#if in the same month and year, make the same datapoint
JAE.merge$month <- as.character(month(JAE.merge$mid_date))
JAE.merge$month[JAE.merge$month!="10" &JAE.merge$month!="11" & JAE.merge$month!="12"] <- paste0("0", JAE.merge$month[JAE.merge$month!="10" &JAE.merge$month!="11" & JAE.merge$month!="12"])
JAE.merge$year <- year(JAE.merge$mid_date)
JAE.merge$plot_date <- as.Date(as.character(paste0(JAE.merge$year, "-", JAE.merge$month, "-15")))

#unique(JAE.merge$site)
unique(JAE.merge$roost_site)
# JAE.merge$roost_site[JAE.merge$roost_site=="Angavokely_Edup"| JAE.merge$roost_site=="Angavobe_Edup"|JAE.merge$roost_site=="Maromizaha_Rmad"] <- "Moramanga"
# JAE.merge = subset(JAE.merge, roost_site=="Moramanga" )
JAE.sum <- ddply(JAE.merge,.(bat_species, roost_site, bat_sex, antigen, plot_date), summarise, N_seropos=sum(seropos), N_tot=length(seropos))
head(JAE.sum)
names(JAE.sum)[names(JAE.sum)=="plot_date"] <- "mid_date"


JAE.sum.all <- ddply(JAE.merge,.(bat_species, roost_site, antigen, plot_date), summarise, N_seropos=sum(seropos), N_tot=length(seropos))
head(JAE.sum.all)
names(JAE.sum.all)[names(JAE.sum.all)=="plot_date"] <-"mid_date"

JAE.sum.all$N_seroneg = JAE.sum.all$N_tot-JAE.sum.all$N_seropos

JAE.sum$N_seroneg = JAE.sum$N_tot-JAE.sum$N_seropos
JAE.sum$seroprev <- JAE.sum$N_seropos/JAE.sum$N_tot
JAE.sum$seroprev_lci <- binom.exact(JAE.sum$N_seropos, JAE.sum$N_tot)$lower
JAE.sum$seroprev_uci <- binom.exact(JAE.sum$N_seropos, JAE.sum$N_tot)$upper
JAE.sum$bat_sex[JAE.sum$bat_sex=="M"]<- "male"
JAE.sum$bat_sex[JAE.sum$bat_sex=="F"]<- "female"
JAE.sum$bat_sex <- as.factor(JAE.sum$bat_sex)
JAE.sum$mid_date <- as.Date(JAE.sum$mid_date)
JAE.sum$antigen <- as.factor(JAE.sum$antigen)
JAE.sum$year <- year(JAE.sum$mid_date)
JAE.sum$doy <- yday(JAE.sum$mid_date)
JAE.sum$year_f <- as.factor(JAE.sum$year)

#is there seasonality without sex by antigen?
#largest possible k that runs (k=-1 is actually 10)
JAE.sum.all$antigen <- as.factor(JAE.sum.all$antigen)
JAE.sum.all$roost_site <- as.factor(JAE.sum.all$roost_site)
JAE.sum.all$doy <- yday(JAE.sum.all$mid_date)
JAE.sum.all$year <- year(JAE.sum.all$mid_date)
JAE.sum.all$year_f <- as.factor(JAE.sum.all$year)

# all antigens
JAE.Eid <- subset(JAE.sum.all, bat_species=="Eidolon dupreanum")
JAE.Eid$antigen <- droplevels(JAE.Eid$antigen)
JAE.Eid$roost_site <- droplevels(JAE.Eid$roost_site)

#below is very slow to run
# gamJEidNoSex <- gam(cbind(N_seropos, N_seroneg)~ antigen + roost_site + 
#                       s(as.numeric(doy), year_f, k=7, bs="fs", by = interaction(antigen, roost_site)) +
#                       s(as.numeric(year), k=3, bs="tp",  by = interaction(antigen, roost_site)), 
#                       #s(as.numeric(mid_date), k=7, bs="cr"), 
#                       data=JAE.Eid, 
#                       family = "binomial", 
#                       select = TRUE,
#                       method = "REML")

gamJEidNoSex <- gam(cbind(N_seropos, N_seroneg)~  roost_site + 
                      s(as.numeric(doy), k=9, bs="cc",  by = roost_site) +
                      s(as.numeric(mid_date), k=6, bs="tp",  by =  roost_site),
                    data=JAE.Eid, 
                    family = "binomial", 
                    select = TRUE,
                    method = "REML")


gamJEidNoSexPlot <- gam(cbind(N_seropos, N_seroneg)~  roost_site + 
                      s(as.numeric(doy), k=9, bs="cc",  by = roost_site) +
                      s(as.numeric(mid_date), k=6, bs="tp",  by =  roost_site),
                    data=JAE.Eid, 
                    gamma = .25,
                    family = "binomial", 
                    select = TRUE,
                    method = "REML")



summary(gamJEidNoSex)
#plot(gamJEidNoSex)

# and with sex
# all antigens
JAE.Eid.sex <- subset(JAE.sum, bat_species=="Eidolon dupreanum")
JAE.Eid.sex$antigen <- droplevels(JAE.Eid.sex$antigen)
JAE.Eid.sex$roost_site <- as.factor(JAE.Eid.sex$roost_site)
JAE.Eid.sex$roost_site <- droplevels(JAE.Eid.sex$roost_site)
JAE.Eid.sex$bat_sex <- droplevels(JAE.Eid.sex$bat_sex)
# 


gamJEidSex <- gam(cbind(N_seropos, N_seroneg)~  roost_site + bat_sex +
                    s(as.numeric(doy),  k=9, bs="cc", by = interaction(roost_site, bat_sex)) +
                    s(as.numeric(mid_date), k=6, bs="tp",  by = interaction(antigen, roost_site, bat_sex)), 
                  #s(as.numeric(mid_date), k=7, bs="cr"), 
                  data=JAE.Eid.sex, 
                  family = "binomial", 
                  select = TRUE,
                  method = "REML")


summary(gamJEidSex)

AIC(gamJEidNoSex, gamJEidSex) #no sex is better


#and plot the predictions
pred.yr.JAE <- cbind.data.frame(doy= 1:365, mid_date=seq(from=as.Date("2015-01-01"), to=as.Date("2015-12-31"), by=1))
head(pred.yr.JAE)
pred.yr.JAE$antigen <- "NiV"
pred.yr.JAE$year_f <- "2015"
pred.yr.JAE$roost_site <- "Moramanga"
# pred.yr.JAE2 <- pred.yr.JAE
# pred.yr.JAE2$antigen <- "HeV"
# pred.yr.JAE <- rbind(pred.yr.JAE, pred.yr.JAE2)

pred.yr.JAE$seroprevalence <-NA
pred.yr.JAE$seroprev_uci <- NA
pred.yr.JAE$seroprev_lci <- NA

pred.yr.JAE$seroprevalence[pred.yr.JAE$antigen=="NiV"] <- predict(gamJEidNoSexPlot, type = "response", newdata = pred.yr.JAE[pred.yr.JAE$antigen=="NiV",], exclude = "s(mid_date)")
pred.yr.JAE$seroprev_uci[pred.yr.JAE$antigen=="NiV"] <- predict(gamJEidNoSexPlot, type = "response", newdata = pred.yr.JAE[pred.yr.JAE$antigen=="NiV",], exclude = "s(mid_date)") + 1.96*predict(gamJEidNoSexPlot, type = "response", newdata = pred.yr.JAE[pred.yr.JAE$antigen=="NiV",], exclude = "s(mid_date)", se.fit = T)$se
pred.yr.JAE$seroprev_lci[pred.yr.JAE$antigen=="NiV"] <- predict(gamJEidNoSexPlot, type = "response", newdata = pred.yr.JAE[pred.yr.JAE$antigen=="NiV",], exclude = "s(mid_date)") - 1.96*predict(gamJEidNoSexPlot, type = "response", newdata = pred.yr.JAE[pred.yr.JAE$antigen=="NiV",], exclude = "s(mid_date)", se.fit = T)$se
pred.yr.JAE$bat_species <- "Eidolon dupreanum"


# pred.yr.JAE$seroprevalence[pred.yr.JAE$antigen=="HeV"] <- predict(gamJEidNoSex, type = "response", newdata = pred.yr.JAE[pred.yr.JAE$antigen=="HeV",], exclude = "s(year_corrected_f)")
# pred.yr.JAE$seroprev_uci[pred.yr.JAE$antigen=="HeV"] <- predict(gamJEidNoSex, type = "response", newdata = pred.yr.JAE[pred.yr.JAE$antigen=="HeV",], exclude = "s(year_corrected_f)") + 1.96*predict(gamJEidNoSex, type = "response", newdata = pred.yr.JAE[pred.yr.JAE$antigen=="HeV",], exclude = "s(year_corrected_f)", se.fit = T)$se
# pred.yr.JAE$seroprev_lci[pred.yr.JAE$antigen=="HeV"] <- predict(gamJEidNoSex, type = "response", newdata = pred.yr.JAE[pred.yr.JAE$antigen=="HeV",], exclude = "s(year_corrected_f)") - 1.96*predict(gamJEidNoSex, type = "response", newdata = pred.yr.JAE[pred.yr.JAE$antigen=="HeV",], exclude = "s(year_corrected_f)", se.fit = T)$se


pred.yr.JAE$seroprev_lci[pred.yr.JAE$seroprev_lci<0]<-0
pred.yr.JAE$seroprev_uci[pred.yr.JAE$seroprev_uci>1]<-1
pred.yr.JAE$bat_species="Eidolon dupreanum"


JAE.sum.all$sig <- 2
JAE.sum.all$sig[JAE.sum.all$antigen=="NiV"] <- 1

JAE.sum.all$sig<- as.factor(JAE.sum.all$sig)

JAE.sum.all$seroprevalence <- JAE.sum.all$N_seropos/JAE.sum.all$N_tot
JAE.sum.all$seroprev_uci <- binom.exact(JAE.sum.all$N_seropos,JAE.sum.all$N_tot)$upper
JAE.sum.all$seroprev_lci <- binom.exact(JAE.sum.all$N_seropos,JAE.sum.all$N_tot)$lower
#only plot the Eidolon data since the rousettus did not test the same antigens
#and came out negative for any seasonal signature

#strip
pJAE_new <- ggplot(data=subset(JAE.sum.all, bat_species=="Eidolon dupreanum")) +
  scale_color_manual(values=colorz) +
  scale_fill_manual(values=filz2) +
  scale_size_continuous(name="N bats", breaks = c(10,20,30)) +
  geom_rect(data=cal.sub, aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3, show.legend = F) +
  geom_point(aes(x=doy, y=seroprevalence, size=N_tot, color=antigen, alpha=sig), show.legend = F) +
  #geom_line(aes(x=doy, y=seroprev, color=antigen, alpha=sig), linetype=2, show.legend = F) +
  geom_linerange(aes(x=doy, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen), linetype=2, show.legend = F) +
  scale_alpha_manual(values=alphaz) +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values=filz1) +
  geom_ribbon(data=subset(pred.yr.JAE, bat_species=="Eidolon dupreanum" & antigen=="HeV") , aes(x=doy, ymin= seroprev_lci, ymax=seroprev_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.yr.JAE, bat_species=="Eidolon dupreanum"& antigen=="HeV"), aes(x=doy, y= seroprevalence, color=antigen), alpha=.3, show.legend = F) +
  geom_ribbon(data=subset(pred.yr.JAE, bat_species=="Eidolon dupreanum" & antigen=="NiV") , aes(x=doy, ymin= seroprev_lci, ymax=seroprev_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.yr.JAE, bat_species=="Eidolon dupreanum"& antigen=="NiV"), aes(x=doy, y= seroprevalence, color=antigen), show.legend = F) +
  facet_grid(~bat_species, space = "free", scales = "free") + theme_bw() +
  theme(axis.text.x = element_text(size=14), 
        axis.title.x = element_text(size=16), 
        strip.background.x = element_rect(fill="white"), 
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_blank(), strip.background.y = element_blank(), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,1,.5,1), "lines")) +
  coord_cartesian(expand = F) + xlab("day of year") + ylab("seroprevalence") + guides(alpha="none"); pJAE_new




# 
# #now put it all together
# Fig2 <- cowplot::plot_grid(pIPM, pJAE, pInterAnnual, ncol = 3, nrow = 1, rel_widths = c(.5,.51,1.21), labels = c("A", "B", "C"), label_size = 22, label_x = c(0,-.04,-.02));Fig2

Fig2 <- cowplot::plot_grid(pIPM, pJAE_new, pIntra, ncol = 3, nrow = 1, rel_widths = c(1,.9,1.1), labels = c("A", "B", "C"), label_size = 22, label_x = c(0,-.04,-.02));Fig2

