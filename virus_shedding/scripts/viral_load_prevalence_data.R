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

#subset to eidolon
merge.eid <- subset(merge.dat, bat_species=="Eidolon dupreanum")

#now get prevalence for every sampling session
#first get midpoint date for sampling session

sampling.session <- ddply(merge.dat,.(sampling_session), summarise, min_date =min(processing_date), max_date=max(processing_date))
sampling.session$mid_date <- sampling.session$min_date + ((sampling.session$max_date-sampling.session$min_date)/2)

merge.eid <- merge(merge.eid, sampling.session, by="sampling_session")
merge.eid$GhV_ct[merge.eid$GhV_ct==0] <- NA
merge.eid$AngV_ct[merge.eid$AngV_ct==0] <- NA

shed.prev <- ddply(merge.eid,.(sampling_session, roost_site, mid_date), summarise, NGhVpos = sum(GhV_pos), NAngVpos = sum(AngV_pos), Ntot=length(GhV_pos), meanGhVCt=mean(GhV_ct, na.rm=T), meanAngVCt=mean(AngV_ct, na.rm=T), sdGhVCt=sd(GhV_ct, na.rm=T), sdAngVCt=sd(AngV_ct, na.rm=T))
head(shed.prev)
shed.prev$GhVprev <- shed.prev$NGhVpos/shed.prev$Ntot
shed.prev$AngVprev <- shed.prev$NAngVpos/shed.prev$Ntot

shed.prev$NGhVneg <- shed.prev$Ntot - shed.prev$NGhVpos
shed.prev$NAngVneg <- shed.prev$Ntot - shed.prev$NAngVpos
  
shed.prev$GhVprevlci <- shed.prev$GhVprevuci <- shed.prev$AngVprevlci <- shed.prev$AngVprevuci <-  NA

for(i in 1:nrow(shed.prev)){
  tmp <- binom.test(x=c(shed.prev$NGhVpos[i], shed.prev$NGhVneg[i]))$conf.int  
  tmp2 <- binom.test(x=c(shed.prev$NAngVpos[i], shed.prev$NAngVneg[i]))$conf.int  
  shed.prev$GhVprevlci[i] <- tmp[1]
  shed.prev$GhVprevuci[i] <- tmp[2]
  shed.prev$AngVprevlci[i] <- tmp2[1]
  shed.prev$AngVprevuci[i] <- tmp2[2]
}

head(shed.prev)

m1 <- glm(cbind(NGhVpos, NGhVneg) ~ meanGhVCt:roost_site,
           family = binomial,
           data = shed.prev)

#and as a sensitivity analysis, try a version that weights those 
#sessions with more individuals contributing to Ct score
m1a <- glm(cbind(NGhVpos, NGhVneg) ~ meanGhVCt:roost_site,
    family = binomial,
    weights = NGhVpos,
    data = shed.prev)

summary(m1) #moramanga is sig
summary(m1a) #all three sites are significant


#and AngV
m2 <- glm(cbind(NAngVpos, NAngVneg) ~ meanAngVCt:roost_site,
          family = binomial,
          data = shed.prev)

m2a <- glm(cbind(NAngVpos, NAngVneg) ~ meanAngVCt:roost_site,
           family = binomial,
           weights = NAngVpos,
           data = shed.prev)

summary(m2) #not sig
summary(m2a) #still not sig


#now try a permutation test - shuffle Cts across sampling sessions within a site
#first drop mahabo because it had only one visit


weighted_cor <- function(x, y, w) {
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]
} #spearman
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
    
  }
  
  p_val <- mean(abs(perm_stats) >= abs(obs))
  
  list(observed = obs, p_value = p_val)
}
perm_cor_unweighted <- function(df, virus, n_perm = 10000, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  if(virus=="GhV"){
    
    
    df=subset(df, !is.na(meanGhVCt))
    #w <- df$NGhVpos  # or n_bats
    obs <- cor(df$meanGhVCt, df$GhVprev, method = "spearman")
    
    perm_stats <- replicate(n_perm, {
      
      df_perm <- df
      
      for (s in unique(df$roost_site)) {
        idx <- which(df$roost_site == s)
        df_perm$meanGhVCt[idx] <- sample(df$meanGhVCt[idx])
      }
      
      cor(df_perm$meanGhVCt, df_perm$GhVprev,method = "spearman")
    })
    
  }else if (virus=="AngV"){
    
    df=subset(df, !is.na(meanAngVCt))
    #w <- df$NAngVpos  # or n_bats
    obs <- cor(df$meanAngVCt, df$AngVprev, method = "spearman")
    
    perm_stats <- replicate(n_perm, {
      
      df_perm <- df
      
      for (s in unique(df$roost_site)) {
        idx <- which(df$roost_site == s)
        df_perm$meanAngVCt[idx] <- sample(df$meanAngVCt[idx])
      }
      
      cor(df_perm$meanAngVCt, df_perm$AngVprev, method = "spearman")
    })
    
  }
  
  p_val <- mean(abs(perm_stats) >= abs(obs))
  
  list(observed = obs, p_value = p_val)
}
perm_cor_weighted(df = subset(shed.prev, roost_site!="Mahabo"), virus="GhV",n_perm = 10000)
#$observed: -0.3563853; $p_value= 0.1066. almost significant
perm_cor_weighted(df = subset(shed.prev, roost_site!="Mahabo"), virus="AngV",n_perm = 10000)
#$observed: -0.3293304; $p_value= 0.6732. not significant

#and unweighted
perm_cor_unweighted(df = subset(shed.prev, roost_site!="Mahabo"), virus="GhV",n_perm = 10000)
#$observed: -0.3280877  $p_value= 0.1056. nearly sig
perm_cor_unweighted(df = subset(shed.prev, roost_site!="Mahabo"), virus="AngV",n_perm = 10000)
#$observed: -0.3414698; $p_value= 0.6138. not significant

# now make a plot for the actual paper (either supp or one component of main):
shed.prev$GhVCt_uci <- shed.prev$meanGhVCt+1.96*shed.prev$sdGhVCt
shed.prev$GhVCt_lci <- shed.prev$meanGhVCt-1.96*shed.prev$sdGhVCt

shed.prev$AngVCt_uci <- shed.prev$meanAngVCt+1.96*shed.prev$sdAngVCt
shed.prev$AngVCt_lci <- shed.prev$meanAngVCt-1.96*shed.prev$sdAngVCt

head(shed.prev)

pred.df <- cbind.data.frame(meanGhVCt = rep(25:39, 2), roost_site=rep(c("Moramanga", "Ankarana"), each=15))
pred.df$prevalence <- predict(m1, newdata = pred.df, type = "response")
pred.df$prevalence_lci <- predict(m1, newdata = pred.df, type = "response") -1.96*predict(m1, newdata = pred.df, type = "response", se.fit = T)$se
pred.df$prevalence_lci <- predict(m1, newdata = pred.df, type = "response") +1.96*predict(m1, newdata = pred.df, type = "response", se.fit = T)$se

#I think I like it better without the model, so just plotting data here
p1 <- ggplot(data= subset(shed.prev, !is.na(meanGhVCt))) + 
      geom_point(aes(x=meanGhVCt, y=GhVprev, color=roost_site, size=Ntot)) + theme_bw() +
      geom_linerange(aes(x=meanGhVCt, ymin=GhVprevlci, ymax=GhVprevuci, color=roost_site), linetype=3, alpha=.5) +
      geom_linerange(aes(xmin=GhVCt_lci, xmax=GhVCt_uci, y=GhVprev, color=roost_site), linetype=3, alpha=.5) +
      #geom_line(data=pred.df, aes(x=meanGhVCt, y=prevalence, color = roost_site)) +
      xlab("Ct (GhV)") + ylab("prevalence (GhV)") +
      scale_size_continuous(name="N bats") +
      scale_color_discrete(name="roost site") +
      theme(panel.grid = element_blank(), axis.text = element_text(size=14),
        axis.title = element_text(size=18), legend.position = "inside",
        legend.position.inside = c(.15,.73), legend.text = element_text(size=12),
        legend.title = element_text(size=14));p1

ggsave(file = paste0(homewd, "final_figures_tables/tmp/FigX_prev_by_Ct.png"),
       plot= p1,
       units="mm",  
       width=55, 
       height=45, 
       scale=3, 
       dpi=300)

#save this one as supplementary and pair with below lag plot
p2 <- ggplot(data= subset(shed.prev, !is.na(meanAngVCt))) + 
  geom_point(aes(x=meanAngVCt, y=AngVprev, color=roost_site, size=Ntot)) + theme_bw() +
  geom_linerange(aes(x=meanAngVCt, ymin=AngVprevlci, ymax=AngVprevuci, color=roost_site), linetype=3, alpha=.5) +
  geom_linerange(aes(xmin=AngVCt_lci, xmax=AngVCt_uci, y=AngVprev, color=roost_site), linetype=3, alpha=.5) +
  xlab("Ct (AngV)") + ylab("\nprevalence (AngV)") +
  scale_size_continuous(name="N bats") +
  scale_color_manual(name="roost site", values=c("Ankarana" = "#F8766D", "Moramanga"="#619CFF")) +
  theme(panel.grid = element_blank(), axis.text = element_text(size=14),
        axis.title = element_text(size=18), legend.position = "inside",
        legend.position.inside = c(.15,.66), legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        plot.margin = unit(c(.5,.5,1,1.1), "lines") );p2


 

#and is the correlation improved by a lag?
#summarize by mean/median ct and prevalence per sampling session
#first get middate of sampling session
head(shed.prev)
tail(shed.prev)
shed.prev <- subset(shed.prev, mid_date <"2021-01-01" & roost_site!="Mahabo") #only one singleton in mahabo dataset so no use in including it


ggplot(shed.prev) + geom_point(aes(x = meanGhVCt, y=GhVprev, color=mid_date)) + facet_grid(~roost_site)
ggplot(shed.prev) + geom_point(aes(x = meanAngVCt, y=AngVprev, color=mid_date)) + facet_grid(~roost_site)

#order by date and site --
shed.prev <- arrange(shed.prev, roost_site, mid_date)

#now explore lags between sampling sessions
#data gets much thinner, so let's only look at GhV for now
lagged_pairs <- function(df) {
  df = df[!is.na(df$meanGhVCt),]
  results <- data.frame()
  
  for (i in 1:(nrow(df)-1)) {
    for (j in (i+1):nrow(df)) {
      
      dt <- as.numeric(df$mid_date[j] - df$mid_date[i])  # in days
      
      results <- rbind(results, data.frame(
        lag_days = dt,
        meanGhVCt = df$meanGhVCt[i],
        #meanAngVCt = df$meanAngVCt[i],
        prevGhV = df$GhVprev[j],
        #prevAngV = df$prevAngV[j],
        roost_site = df$roost_site[i]
      ))
    }
  }
  
  return(results)
}


#and permute - use "lag_limit" to only look at lags within a biologically relevant window
perm_test_short_lag <- function(df, lag_limit, n_perm = 10000) {
  df = df[!is.na(df$meanGhVCt),]
  
  
  site.split <- dlply(df,.(roost_site))
  
  pairs <- data.table::rbindlist(lapply(site.split, lagged_pairs))
  #sort(unique(lag.sites$lag_days))
  
  # observed
  pairs <- subset(pairs, lag_days > 0 & lag_days <= lag_limit)
  
  obs <- cor(pairs$meanGhVCt, pairs$prevGhV, method = "spearman")  
  
  perm_stats <- replicate(n_perm, {
    
    df_perm <- df
    
    # shuffle Ct within site
    for (s in unique(df$roost_site)) {
      idx <- which(df$roost_site == s)
      df_perm$meanGhVCt[idx] <- sample(df$meanGhVCt[idx])
    }
    
    pairs_perm <- lagged_pairs(df_perm)
    pairs_perm <- subset(pairs_perm, lag_days > 0 & lag_days <= lag_limit & !is.na(meanGhVCt))
    
    cor(pairs_perm$meanGhVCt, pairs_perm$prevGhV, method = "spearman")
  })
  
  perm_stats <- perm_stats[!is.na(perm_stats)]
  mean(abs(perm_stats) >= abs(obs))
}

perm_test_short_lag(df=shed.prev, lag_limit = 40, n_perm = 10000)
#p-val = 0.5236
#no support for lag of 40 days

perm_test_short_lag(df=shed.prev, lag_limit = 30, n_perm = 10000)
#p-val = 0.6003
#no support for lag of 30 days 

perm_test_short_lag(df=shed.prev, lag_limit = 15, n_perm = 10000)
# p-val = 0.088
# this is slightly significant. it is actually a lag of 2-3 days because
# that is all we have in the dataset that is <15 (so reducing to lower lags is
# not helpful)


## CONCLUSION
# no support for lagged data as more predictive... but we don't have the fine-scale
# sampling that is needed to parse this. biologically, if there was a lag, it would only be
# in the ~5-7-10 day window reflecting the lag in peak viral load and population-level transmission

#let's plot the lag duration vs. correlation
bin_lags <- function(pairs) {
  
  pairs <- arrange(pairs, lag_days)
  pairs$lag_bin <- cut(
    pairs$lag_days,
    breaks = c(0, 15, 30, 50, 70, 90),
    labels = c("0–15d", "15–30d", "30–50d", "50–70d", "70–90d"),
    right = TRUE
  )
  
  pairs[!is.na(pairs$lag_bin), ]
}

bin.df <- bin_lags(lagged_pairs(shed.prev))

#get correlation
summarise_bins <- function(pairs) {
  
  aggregate(cbind(meanGhVCt, prevGhV) ~ lag_bin, data = pairs, function(x) NA)  # placeholder
  
  library(dplyr)
  
  pairs %>%
    group_by(lag_bin) %>%
    summarise(
      r = cor(meanGhVCt, prevGhV, method = "spearman"),#, use = "complete.obs"),
      n = n()
    )
}

#just do Moramanga data because not enough pairs in Ankarana
bin.sum <- summarise_bins(bin.df[bin.df$roost_site=="Moramanga",])

plot_lag_correlation <- function(summary_df) {
  
  ggplot(data=summary_df)+
    geom_point(aes(x=lag_bin, y=r, size=n)) +
    scale_size_continuous(name="N bats") +
    theme_bw() + xlab("Lag window") +
    ylab("Spearman correlation for GhV\n(Ct vs prevalence)") +
    geom_hline(aes(yintercept=0), linetype=2) +
    theme(panel.grid = element_blank(), axis.text = element_text(size=14),
          axis.title = element_text(size=18), legend.position = "inside",
          legend.position.inside = c(.1,.75), legend.text = element_text(size=12),
          legend.title = element_text(size=14),
          plot.margin = unit(c(1.5,.5,.5,1), "lines")) 
    
}

p3 <- plot_lag_correlation(bin.sum)
p3
# we can see that the negative correlation is strongest at lower lags, indicating
# that this is a temporally tight relationship - but again, we lack the data to really
# find the lag specifically

plot.supp <- cowplot::plot_grid(p2,p3, ncol=1, nrow=2, labels=c("A", "B"), label_size = 22, label_x = -.02, label_y = 1.02); plot.supp

ggsave(file = paste0(homewd, "final_figures_tables/tmp/FigSX_prev_by_Ct_AngV_plus_lag.png"),
       plot= plot.supp,
       units="mm",  
       width=55, 
       height=75, 
       scale=3, 
       dpi=300)


#now repeat all for rousettus
rm(list=ls())


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

#subset to rousettus
merge.eid <- subset(merge.dat, bat_species=="Rousettus madagascariensis")

#now get prevalence for every sampling session
#first get midpoint date for sampling session

sampling.session <- ddply(merge.dat,.(sampling_session), summarise, min_date =min(processing_date), max_date=max(processing_date))
sampling.session$mid_date <- sampling.session$min_date + ((sampling.session$max_date-sampling.session$min_date)/2)

merge.eid <- merge(merge.eid, sampling.session, by="sampling_session")
merge.eid$GhV_ct[merge.eid$GhV_ct==0] <- NA
merge.eid$AngV_ct[merge.eid$AngV_ct==0] <- NA

shed.prev <- ddply(merge.eid,.(sampling_session, roost_site, mid_date), summarise, NGhVpos = sum(GhV_pos), NAngVpos = sum(AngV_pos), Ntot=length(GhV_pos), meanGhVCt=mean(GhV_ct, na.rm=T), meanAngVCt=mean(AngV_ct, na.rm=T), sdGhVCt=sd(GhV_ct, na.rm=T), sdAngVCt=sd(AngV_ct, na.rm=T))
head(shed.prev)
shed.prev$GhVprev <- shed.prev$NGhVpos/shed.prev$Ntot
shed.prev$AngVprev <- shed.prev$NAngVpos/shed.prev$Ntot

shed.prev$NGhVneg <- shed.prev$Ntot - shed.prev$NGhVpos
shed.prev$NAngVneg <- shed.prev$Ntot - shed.prev$NAngVpos

shed.prev$GhVprevlci <- shed.prev$GhVprevuci <- shed.prev$AngVprevlci <- shed.prev$AngVprevuci <-  NA

for(i in 1:nrow(shed.prev)){
  tmp <- binom.test(x=c(shed.prev$NGhVpos[i], shed.prev$NGhVneg[i]))$conf.int  
  tmp2 <- binom.test(x=c(shed.prev$NAngVpos[i], shed.prev$NAngVneg[i]))$conf.int  
  shed.prev$GhVprevlci[i] <- tmp[1]
  shed.prev$GhVprevuci[i] <- tmp[2]
  shed.prev$AngVprevlci[i] <- tmp2[1]
  shed.prev$AngVprevuci[i] <- tmp2[2]
}

head(shed.prev)

m1 <- glm(cbind(NGhVpos, NGhVneg) ~ meanGhVCt:roost_site,
          family = binomial,
          data = shed.prev)

#and as a sensitivity analysis, try a version that weights those 
#sessions with more individuals contributing to Ct score
m1a <- glm(cbind(NGhVpos, NGhVneg) ~ meanGhVCt:roost_site,
           family = binomial,
           weights = NGhVpos,
           data = shed.prev)

summary(m1) #no sig
summary(m1a) #Moramanga is significant


#and AngV
m2 <- glm(cbind(NAngVpos, NAngVneg) ~ meanAngVCt:roost_site,
          family = binomial,
          data = shed.prev)

m2a <- glm(cbind(NAngVpos, NAngVneg) ~ meanAngVCt:roost_site,
           family = binomial,
           weights = NAngVpos,
           data = shed.prev)

summary(m2) #not sig
summary(m2a) #still not sig


#now try a permutation test - shuffle Cts across sampling sessions within a site
#first drop mahabo because it had only one visit


weighted_cor <- function(x, y, w) {
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]
} #spearman
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
    
  }
  
  p_val <- mean(abs(perm_stats) >= abs(obs))
  
  list(observed = obs, p_value = p_val)
}
perm_cor_unweighted <- function(df, virus, n_perm = 10000, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  if(virus=="GhV"){
    
    
    df=subset(df, !is.na(meanGhVCt))
    #w <- df$NGhVpos  # or n_bats
    obs <- cor(df$meanGhVCt, df$GhVprev, method = "spearman")
    
    perm_stats <- replicate(n_perm, {
      
      df_perm <- df
      
      for (s in unique(df$roost_site)) {
        idx <- which(df$roost_site == s)
        df_perm$meanGhVCt[idx] <- sample(df$meanGhVCt[idx])
      }
      
      cor(df_perm$meanGhVCt, df_perm$GhVprev,method = "spearman")
    })
    
  }else if (virus=="AngV"){
    
    df=subset(df, !is.na(meanAngVCt))
    #w <- df$NAngVpos  # or n_bats
    obs <- cor(df$meanAngVCt, df$AngVprev, method = "spearman")
    
    perm_stats <- replicate(n_perm, {
      
      df_perm <- df
      
      for (s in unique(df$roost_site)) {
        idx <- which(df$roost_site == s)
        df_perm$meanAngVCt[idx] <- sample(df$meanAngVCt[idx])
      }
      
      cor(df_perm$meanAngVCt, df_perm$AngVprev, method = "spearman")
    })
    
  }
  
  p_val <- mean(abs(perm_stats) >= abs(obs))
  
  list(observed = obs, p_value = p_val)
}
perm_cor_weighted(df = subset(shed.prev, roost_site!="Makira"), virus="GhV",n_perm = 10000)
#$observed: 0.5241467; $p_value= 0.0461. significant... but backwards
perm_cor_weighted(df = subset(shed.prev, roost_site!="Makira"), virus="AngV",n_perm = 10000)
#$observed: -0.608751; $p_value= 0.1964. not significant but close

#and unweighted
perm_cor_unweighted(df = subset(shed.prev, roost_site!="Makira"), virus="GhV",n_perm = 10000)
#$observed: 0.5222545  $p_value= 0.0902. sig
perm_cor_unweighted(df = subset(shed.prev, roost_site!="Mahabo"), virus="AngV",n_perm = 10000)
#$observed: -0.7807201; $p_value= 0.0564. significant here

# now make a plot for the actual paper (either supp or one component of main):
shed.prev$GhVCt_uci <- shed.prev$meanGhVCt+1.96*shed.prev$sdGhVCt
shed.prev$GhVCt_lci <- shed.prev$meanGhVCt-1.96*shed.prev$sdGhVCt

shed.prev$AngVCt_uci <- shed.prev$meanAngVCt+1.96*shed.prev$sdAngVCt
shed.prev$AngVCt_lci <- shed.prev$meanAngVCt-1.96*shed.prev$sdAngVCt

head(shed.prev)

pred.df <- cbind.data.frame(meanGhVCt = rep(25:39, 2), roost_site=rep(c("Moramanga", "Ankarana"), each=15))
pred.df$prevalence <- predict(m1, newdata = pred.df, type = "response")
pred.df$prevalence_lci <- predict(m1, newdata = pred.df, type = "response") -1.96*predict(m1, newdata = pred.df, type = "response", se.fit = T)$se
pred.df$prevalence_lci <- predict(m1, newdata = pred.df, type = "response") +1.96*predict(m1, newdata = pred.df, type = "response", se.fit = T)$se

#I think I like it better without the model, so just plotting data here
p1 <- ggplot(data= subset(shed.prev, !is.na(meanGhVCt))) + 
  geom_point(aes(x=meanGhVCt, y=GhVprev, color=roost_site, size=Ntot)) + theme_bw() +
  geom_linerange(aes(x=meanGhVCt, ymin=GhVprevlci, ymax=GhVprevuci, color=roost_site), linetype=3, alpha=.5) +
  geom_linerange(aes(xmin=GhVCt_lci, xmax=GhVCt_uci, y=GhVprev, color=roost_site), linetype=3, alpha=.5) +
  #geom_line(data=pred.df, aes(x=meanGhVCt, y=prevalence, color = roost_site)) +
  xlab("Ct (GhV)") + ylab("prevalence (GhV)") +
  scale_size_continuous(name="N bats") +
  scale_color_discrete(name="roost site") +
  theme(panel.grid = element_blank(), axis.text = element_text(size=14),
        axis.title = element_text(size=18), legend.position = "inside",
        legend.position.inside = c(.15,.73), legend.text = element_text(size=12),
        legend.title = element_text(size=14));p1

# ggsave(file = paste0(homewd, "final_figures_tables/tmp/FigX_prev_by_Ct_Rousettus.png"),
#        plot= p1,
#        units="mm",  
#        width=55, 
#        height=45, 
#        scale=3, 
#        dpi=300)

#save this one as supplementary and pair with below lag plot
p2 <- ggplot(data= subset(shed.prev, !is.na(meanAngVCt))) + 
  geom_point(aes(x=meanAngVCt, y=AngVprev, color=roost_site, size=Ntot), show.legend = F) + theme_bw() +
  geom_linerange(aes(x=meanAngVCt, ymin=AngVprevlci, ymax=AngVprevuci, color=roost_site), linetype=3, alpha=.5, show.legend = F) +
  geom_linerange(aes(xmin=AngVCt_lci, xmax=AngVCt_uci, y=AngVprev, color=roost_site), linetype=3, alpha=.5, show.legend = F) +
  xlab("Ct (AngV)") + ylab("\nprevalence (AngV)") +
  scale_size_continuous(name="N bats",breaks = c(5,10,15), labels = c(5,10,15)) +
  scale_color_manual(name="roost site", values=c("Ankarana" = "#F8766D", "Moramanga"="#619CFF")) +
  theme(panel.grid = element_blank(), axis.text = element_text(size=14),
        axis.title = element_text(size=18), legend.position = "inside",
        legend.position.inside = c(.15,.66), legend.text = element_text(size=12),
        legend.title = element_text(size=14)) + coord_cartesian(xlim=c(25,37.5), ylim = c(0,1));p2




#and is the correlation improved by a lag?
#summarize by mean/median ct and prevalence per sampling session
#first get middate of sampling session
head(shed.prev)
tail(shed.prev)
shed.prev <- subset(shed.prev, mid_date <"2021-01-01" & roost_site!="Makira") #only one singleton in mahabo dataset so no use in including it


ggplot(shed.prev) + geom_point(aes(x = meanGhVCt, y=GhVprev, color=mid_date)) + facet_grid(~roost_site)
ggplot(shed.prev) + geom_point(aes(x = meanAngVCt, y=AngVprev, color=mid_date)) + facet_grid(~roost_site)

#order by date and site --
shed.prev <- arrange(shed.prev, roost_site, mid_date)

#now explore lags between sampling sessions
#data gets much thinner, so let's only look at GhV for now
lagged_pairs <- function(df) {
  df = df[!is.na(df$meanGhVCt),]
  results <- data.frame()
  
  for (i in 1:(nrow(df)-1)) {
    for (j in (i+1):nrow(df)) {
      
      dt <- as.numeric(df$mid_date[j] - df$mid_date[i])  # in days
      
      results <- rbind(results, data.frame(
        lag_days = dt,
        meanGhVCt = df$meanGhVCt[i],
        #meanAngVCt = df$meanAngVCt[i],
        prevGhV = df$GhVprev[j],
        #prevAngV = df$prevAngV[j],
        roost_site = df$roost_site[i]
      ))
    }
  }
  
  return(results)
}


#and permute - use "lag_limit" to only look at lags within a biologically relevant window
perm_test_short_lag <- function(df, lag_limit, n_perm = 10000) {
  df = df[!is.na(df$meanGhVCt),]
  
  
  site.split <- dlply(df,.(roost_site))
  
  pairs <- data.table::rbindlist(lapply(site.split, lagged_pairs))
  #sort(unique(lag.sites$lag_days))
  
  # observed
  pairs <- subset(pairs, lag_days > 0 & lag_days <= lag_limit)
  
  obs <- cor(pairs$meanGhVCt, pairs$prevGhV, method = "spearman")  
  
  perm_stats <- replicate(n_perm, {
    
    df_perm <- df
    
    # shuffle Ct within site
    for (s in unique(df$roost_site)) {
      idx <- which(df$roost_site == s)
      df_perm$meanGhVCt[idx] <- sample(df$meanGhVCt[idx])
    }
    
    pairs_perm <- lagged_pairs(df_perm)
    pairs_perm <- subset(pairs_perm, lag_days > 0 & lag_days <= lag_limit & !is.na(meanGhVCt))
    
    cor(pairs_perm$meanGhVCt, pairs_perm$prevGhV, method = "spearman")
  })
  
  perm_stats <- perm_stats[!is.na(perm_stats)]
  mean(abs(perm_stats) >= abs(obs))
}

#none of these run because we have no pairs of GhV shedding data for 
#rousettus with < 40 day lag
perm_test_short_lag(df=shed.prev, lag_limit = 40, n_perm = 10000)
#NA

#so let's just make the base plots as supplemental for Rousettus

plot.out <- cowplot::plot_grid(p1,p2, nrow = 1, ncol=2, labels = c("A", "B"), label_size = 22); plot.out
# now add the title

ggsave(file = paste0(homewd, "final_figures_tables/tmp/FigX_prev_by_Ct_Rousettus.png"),
        plot= plot.out,
        units="mm",  
        width=95, 
        height=45, 
        scale=3, 
        dpi=300)
