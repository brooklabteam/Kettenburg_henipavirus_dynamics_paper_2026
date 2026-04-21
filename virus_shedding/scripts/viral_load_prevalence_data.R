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


homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/"
batwd = "/Users/carabrook/Developer/Madagascar-Bat-Data"

# First, load the shedding data 
shed.dat <- read.csv(file = paste0(homewd, "/virus_shedding/data/henipavirus_shedding_data.csv"))
head(shed.dat)
names(shed.dat)


shed.dat <- dplyr::select(shed.dat, sampleid, GhV_shedding_ur, AngV_shedding_ur, GhV_shedding_ct_ur, AngV_shedding_ct_ur)

names(shed.dat)[names(shed.dat) == "GhV_shedding_ur"] <- "GhV_pos"
names(shed.dat)[names(shed.dat) == "AngV_shedding_ur"] <- "AngV_pos"

names(shed.dat)[names(shed.dat) == "GhV_shedding_ct_ur"] <- "GhV_ct"
names(shed.dat)[names(shed.dat) == "AngV_shedding_ct_ur"] <- "AngV_ct"


# link this to the catching data to plot by sampling session
catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))

head(catch.dat)
names(catch.dat)
catch.dat$net_date[!is.na(catch.dat$net_date)] <- as.Date(catch.dat$net_date[!is.na(catch.dat$net_date)], format = "%m/%d/%Y")
catch.dat$net_date[is.na(catch.dat$net_date)] <- catch.dat$processing_date[is.na(catch.dat$net_date)]

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, net_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)
names(catch.dat)[names(catch.dat)=="net_date"] <- "processing_date"

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

#subset to eidolon
merge.eid <- subset(merge.dat, bat_species=="Eidolon dupreanum")

#now get prevalence for every day of the dataset
shed.prev <- ddply(merge.eid,.(processing_date), summarise, NGhVpos = sum(GhV_pos), NAngVpos = sum(AngV_pos), Ntot=length(GhV_pos))
head(shed.prev)
shed.prev$GhVprev <- shed.prev$NGhVpos/shed.prev$Ntot
shed.prev$AngVprev <- shed.prev$NAngVpos/shed.prev$Ntot
shed.prev <- dplyr::select(shed.prev,-(NGhVpos), -(NAngVpos))

merge.eid <- merge(merge.eid, shed.prev, by="processing_date")
head(merge.eid)
length(unique(merge.eid$barcode_id[merge.eid$GhV_ct>0]))#68
length(merge.eid$barcode_id[merge.eid$GhV_ct>0])#71. there are three repeated
merge.eid$barcode_id[duplicated(merge.eid$barcode_id) & merge.eid$GhV_ct>0]
#"921", "1086_969"  ,"956000005671203" - these are three bats sampled at multiple timepoints
# not so many - can rerun without them as well 

m1 <- lm(GhV_ct ~ GhVprev:roost_site, data=subset(merge.eid, GhV_ct>0))
summary(m1) #moramanga and ankarana
m1b <- lm(GhV_ct ~ GhVprev, data=subset(merge.eid, GhV_ct>0 & roost_site=="Moramanga"))
summary(m1b) #improves

#or the pearson correlation
with(subset(merge.eid, GhV_ct>0 & roost_site=="Moramanga"), cor.test(GhVprev, GhV_ct, method = "pearson"))
#t = -1.8863, df = 38, p-value = 0.06691
with(subset(merge.eid, GhV_ct>0 & roost_site=="Ankarana"), cor.test(GhVprev, GhV_ct, method = "pearson"))
#t = -0.92411, df = 27, p-value = 0.3636

#spearman
with(subset(merge.eid, GhV_ct>0 & roost_site=="Moramanga"), cor.test(GhVprev, GhV_ct, method = "spearman"))
#S = 15069, p-value = 0.007975

with(subset(merge.eid, GhV_ct>0 & roost_site=="Ankarana"), cor.test(GhVprev, GhV_ct, method = "spearman"))
#S = 5258.3, p-value = 0.1201

ggplot(data= subset(merge.eid, GhV_ct>0)) + geom_point(aes(x=GhV_ct, y=GhVprev, color=roost_site, size=Ntot)) + theme_bw()

library(sjPlot)
plot_model(m1, type="pred")
plot_model(m1b, type="pred")

m2 <- lm(AngV_ct ~ AngVprev:roost_site, data=subset(merge.eid, AngV_ct>0))
summary(m2)#Ankarana is sig
m2b <- lm(AngV_ct ~ AngVprev, data=subset(merge.eid, AngV_ct>0 & roost_site=="Moramanga"))
summary(m2b) #not sig. only 12 data points
plot_model(m2, type="pred")
plot_model(m2b, type="pred")

#or the pearson correlation
with(subset(merge.eid, AngV_ct>0 & roost_site=="Moramanga"), cor.test(AngVprev, AngV_ct, method = "pearson"))
#t = 0.50217, df = 10, p-value = 0.6264
with(subset(merge.eid, AngV_ct>0 & roost_site=="Ankarana"), cor.test(AngVprev, AngV_ct, method = "pearson"))
# t = -1.0181, df = 8, p-value = 0.3384

#spearman
with(subset(merge.eid, AngV_ct>0& roost_site=="Moramanga"), cor.test(AngVprev, AngV_ct, method = "spearman"))
#S = 244.34, p-value = 0.6515
with(subset(merge.eid, AngV_ct>0& roost_site=="Ankarana"), cor.test(AngVprev, AngV_ct, method = "spearman"))
#S = 232.23, p-value = 0.2425

ggplot(data= subset(merge.eid, AngV_ct>0)) + geom_point(aes(x=AngV_ct, y=AngVprev, color=roost_site, size=Ntot)) + theme_bw()

#try a permutation test to shuffle the Ct and prevalence relationships

perm_cor_stratified <- function(x, y, site, n_perm = 10000, method = "spearman", seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # observed correlation (pooled, but structure respected later)
  obs_cor <- cor(x, y, method = method, use = "complete.obs")
  
  # permutation distribution
  perm_stats <- replicate(n_perm, {
    
    y_perm <- y
    
    # shuffle within each site
    for (s in unique(site)) {
      idx <- which(site == s)
      y_perm[idx] <- sample(y[idx])
    }
    
    cor(x, y_perm, method = method, use = "complete.obs")
  })
  
  p_val <- mean(abs(perm_stats) >= abs(obs_cor))
  
  return(list(
    observed_correlation = obs_cor,
    p_value = p_val
  ))
}
perm_cor_test <- function(x, y, n_perm = 10000, method = "spearman", seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # observed correlation
  obs_cor <- cor(x, y, method = method, use = "complete.obs")
  
  # permutation distribution
  perm_stats <- replicate(n_perm, {
    y_perm <- sample(y)  # shuffle y
    cor(x, y_perm, method = method, use = "complete.obs")
  })
  
  # two-sided p-value
  p_val <- mean(abs(perm_stats) >= abs(obs_cor))
  
  return(list(
    observed_correlation = obs_cor,
    p_value = p_val,
    permuted_correlations = perm_stats
  ))
}#no site stratification

#first with GhV
out.GhV <- perm_cor_stratified(x= subset(merge.eid, GhV_ct>0)$GhVprev, y = subset(merge.eid, GhV_ct>0)$GhV_ct, site=subset(merge.eid, GhV_ct>0)$roost_site)
#$observed_correlation: -0.3497536; p_value: 0.003
#prevalence and Ct are negatively correlated, after accounting for site

#if we remove Mahabo:
perm_cor_stratified(x= subset(merge.eid, GhV_ct>0 & roost_site!="Mahabo")$GhVprev, y = subset(merge.eid, GhV_ct>0& roost_site!="Mahabo")$GhV_ct, site=subset(merge.eid, GhV_ct>0& roost_site!="Mahabo")$roost_site)
#$observed_correlation: -0.3633688; p_value: 9e-04
#no change to results - just more sig

#and with AngV
out.AngV <- perm_cor_stratified(x= subset(merge.eid, AngV_ct>0)$AngVprev, y = subset(merge.eid, AngV_ct>0)$AngV_ct, site=subset(merge.eid, AngV_ct>0)$roost_site)
#$observed_correlation -0.4112351; $p_value 0.3451
#not sig for angv


#and is the correlation improved by a lag?
#summarize by sampling session and site
 
#merge.eid.mora = subset(merge.eid, roost_site=="Moramanga")
#head(merge.eid.mora)

#summarize by mean/median ct and prevalence per sampling session
#first get middate of sampling session
merge.eid$processing_date <- as.Date(merge.eid$processing_date)

sampling.session <- ddply(merge.eid,.(sampling_session), summarise, min_date =min(processing_date), max_date=max(processing_date))

sampling.session$mid_date <- sampling.session$min_date + (sampling.session$max_date-sampling.session$min_date)/2

sampling.session$duration <- sampling.session$max_date-sampling.session$min_date
sampling.session <- dplyr::select(sampling.session, -(min_date), -(max_date))

merge.eid <-  merge(merge.eid, sampling.session, by="sampling_session")

head(merge.eid)
merge.eid$GhV_ct[merge.eid$GhV_ct==0] <- NA
merge.eid$AngV_ct[merge.eid$AngV_ct==0] <- NA
mora.sum <- ddply(merge.eid, .(sampling_session, mid_date, roost_site), summarise, GhV_pos = sum(GhV_pos), AngV_pos=sum(AngV_pos), N_tot=length(sampleid), meanGhVCt = mean(GhV_ct, na.rm=T), meanAngVCt=mean(AngV_ct, na.rm=T))
head(mora.sum)
tail(mora.sum)
mora.sum <- arrange(mora.sum, mid_date)
mora.sum <- subset(mora.sum, mid_date <"2021-01-01" & roost_site!="Mahabo") #only one mahabo dataset so no use in including it

mora.sum$prevGhV <- mora.sum$GhV_pos/mora.sum$N_tot
mora.sum$prevAngV <- mora.sum$AngV_pos/mora.sum$N_tot

ggplot(mora.sum) + geom_point(aes(x = meanGhVCt, y=prevGhV, color=mid_date)) + facet_grid(~roost_site)
ggplot(mora.sum) + geom_point(aes(x = meanAngVCt, y=prevGhV, color=mid_date)) + facet_grid(~roost_site)

#order by date and site --
mora.sum <- arrange(mora.sum, roost_site, mid_date)

#now explore lags between sampling sessions
