rm(list=ls())

#this script uses recapture data and aged data to put some restraints on
#the duration of both adult and maternal Ab-mediated immunity.
#this precedes the fitting of our mechanistic models

library(plyr)
library(dplyr)   
library(ggplot2)
library(stats)
library(mgcv)
library(reshape2) 
library(ggh4x)
library(minpack.lm)
library(devtools)
library(rstan)
library(ggplot2)

#fit declining maternal immunity function to data


################################################################################
################################################################################
homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

################################################################################
################################################################################


# First, load the luminex data after determining who is seropositive/
# seronegative in the other scripts
lum.dat <- read.csv(file = paste0(homewd, "/master_data_files/mix_sub_data_75_90_96.csv"))
head(lum.dat)

lum.dat = subset(lum.dat, bat_species=="Eidolon dupreanum" & antigen=="GhV"  |bat_species=="Eidolon dupreanum" & antigen=="AngV" | bat_species=="Rousettus madagascariensis" & antigen=="GhV")

lum.dat <- dplyr::select(lum.dat, sampleid, antigen, MFI, log10MFI, seropos, seropos_lci, seropos_uci)

##merge with catching data


# link this to the catching data to plot by sampling session
catch.dat <- read.csv(file = paste0(homewd, "/master_data_files/catching_data.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))


head(catch.dat)
names(catch.dat)

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, pup_or_mom, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)

#and merge all three datasets
merge.dat <- merge(lum.dat, catch.dat, by="sampleid", all.x = T)
head(merge.dat)


#change roost site to just Moramanga region or Ankarana
unique(merge.dat$roost_site)
merge.dat$roost_site[merge.dat$roost_site=="Angavobe_Edup" |merge.dat$roost_site=="Angavokely_Edup"| merge.dat$roost_site=="Maromizaha_Rmad"] <- "Moramanga"
merge.dat$roost_site[merge.dat$roost_site!="Moramanga"] <- "Ankarana"

head(merge.dat)

#subset to mom-pup pairs
mom.df <- subset(merge.dat, pupmom_pair=="yes")
unique(mom.df$pup_or_mom) #"mom" "pup"

#and build a paired dataset for both
pup.df = subset(mom.df, pup_or_mom=="pup")
names(pup.df)[names(pup.df)=="sampleid"] <- "pup_sampleid"
names(pup.df)[names(pup.df)=="MFI"] <- "pup_MFI"
names(pup.df)[names(pup.df)=="log10MFI"] <- "pup_log10MFI"
names(pup.df)[names(pup.df)=="seropos"] <- "pup_seropos"


pup.df <- dplyr::select(pup.df, pup_sampleid, antigen, pup_MFI, pup_log10MFI, pup_seropos)

#and merge
mom.df = subset(mom.df, pup_or_mom=="mom")
names(mom.df)[names(mom.df)=="sampleid_pupmom_pair"] <- "pup_sampleid"

mom.df <- merge(mom.df, pup.df, by=c("pup_sampleid", "antigen"), all.x = T)
head(mom.df)
mom.df = subset(mom.df, !is.na(log10MFI) & !is.na(pup_log10MFI))
length(unique(mom.df$sampleid[mom.df$antigen=="GhV"])) #29
length(mom.df$sampleid[mom.df$antigen=="GhV"]) #29
length(mom.df$sampleid[mom.df$bat_species=="Rousettus madagascariensis" & mom.df$antigen=="GhV"]) #6. not enough

mom.df = subset(mom.df, bat_species =="Eidolon dupreanum")
mom.df$seropos[mom.df$seropos==0] <- "seronegative"
mom.df$seropos[mom.df$seropos==1] <- "seropositive"

mom.df$pup_seropos[mom.df$pup_seropos==0] <- "seronegative"
mom.df$pup_seropos[mom.df$pup_seropos==1] <- "seropositive"


mom.df$seropos <- as.factor(mom.df$seropos)
mom.df$pup_seropos <- as.factor(mom.df$pup_seropos)

#and plot by species and antigen

p1 <- ggplot(mom.df) + geom_point(aes(x=log10MFI, 
                       y=pup_log10MFI, color=seropos, 
                       shape = pup_seropos), size=2) + facet_grid(~antigen) + theme_bw() +
                       scale_color_manual("mom serostatus", values = c("seronegative" = "cyan3", "seropositive"="tomato"))+
                       scale_shape_discrete("pup serostatus") + ylab("log10 MFI (pup)") + xlab("log10 MFI (mom)") +
                      theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
                      legend.text=element_text(size=12),legend.title=element_text(size=12),
                      strip.text = element_text(size=14), axis.title = element_text(size=14), axis.text = element_text(size=12)) 
p1  


#and do the regression on the data for which at least one is seropos
pos.df <- subset(mom.df, antigen=="GhV" & seropos=="seropositive" |antigen=="GhV" & pup_seropos=="seropositive") #19 individuals
m1  <- lm(pup_log10MFI~log10MFI, data = pos.df)
summary(m1)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.65674    0.27544   6.015 1.39e-05 ***
#   log10MFI     0.65816    0.08064   8.162 2.78e-07 ***
#R-sq = 0.78

pred.df <- cbind.data.frame(log10MFI=seq(2.5,4.5, by=.1))
pred.df$pup_log10MFI <- NA
pred.df$pup_log10MFI <- predict(m1, newdata = pred.df, type = "response")
pred.df$lci <- predict(m1, newdata = pred.df, type = "response") - 1.96*predict(m1, newdata = pred.df, type = "response", se.fit = T)$se
pred.df$uci <- predict(m1, newdata = pred.df, type = "response") + 1.96*predict(m1, newdata = pred.df, type = "response", se.fit = T)$se
pred.df$antigen = "GhV"

p2 <- ggplot(mom.df) + geom_point(aes(x=log10MFI, 
                                      y=pup_log10MFI, color=seropos, 
                                      shape = pup_seropos), size=2) + facet_grid(~antigen) + theme_bw() +
  geom_ribbon(data=pred.df, aes(x=log10MFI, ymin=lci, ymax=uci), alpha=.3) +
  geom_line(data=pred.df, aes(x=log10MFI, y=pup_log10MFI)) +
  scale_color_manual("mom serostatus", values = c("seronegative" = "cyan3", "seropositive"="tomato"))+
  scale_shape_discrete("pup serostatus") + ylab("log10 MFI (pup)") + xlab("log10 MFI (mom)") +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
        legend.text=element_text(size=12),legend.title=element_text(size=12),
        strip.text = element_text(size=14), axis.title = element_text(size=14), axis.text = element_text(size=12)) 
p2

text.df <- cbind.data.frame(x=0.7, y=4.6, antigen="GhV", label ="italic(R)^2 == 0.78")

p3 <- ggplot(mom.df) + geom_point(aes(x=log10MFI, 
                                      y=pup_log10MFI, color=seropos, 
                                      shape = pup_seropos), size=2) + facet_grid(~antigen) + theme_bw() +
  geom_ribbon(data=pred.df, aes(x=log10MFI, ymin=lci, ymax=uci), alpha=.3) +
  geom_line(data=pred.df, aes(x=log10MFI, y=pup_log10MFI)) +
  geom_label(data=text.df, aes(x=x, y=y, label=label), parse=TRUE, label.size = 0) +
  scale_color_manual("mom serostatus", values = c("seronegative" = "cyan3", "seropositive"="tomato"))+
  scale_shape_discrete("pup serostatus") + ylab("log10 MFI (pup)") + xlab("log10 MFI (mom)") +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
        legend.text=element_text(size=12),legend.title=element_text(size=12),
        strip.text = element_text(size=14), axis.title = element_text(size=14), axis.text = element_text(size=12)) 
p3
