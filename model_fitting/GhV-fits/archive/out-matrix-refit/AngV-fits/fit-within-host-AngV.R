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



################################################################################
################################################################################
homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
batwd = "/Users/carabrook/Developer/Madagascar-Bat-Data"
################################################################################
################################################################################

#fit declining maternal immunity function to data

dat <- read.csv(paste0(homewd,"/model_fitting/AngV-fits/fit-1-AngV/age-seroprevalence-data.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)[names(dat)=="species"] <- "bat_species"
dat$log10MFI <- log10(dat$MFI)

#try AngV
dat.mat = subset(dat, antigen=="AngV" & age <=1 & bat_species=="Eidolon dupreanum")
sort(dat.mat$log10MFI)
sort(dat.mat$prev) #our, only 1 seropos by this threshold
sort(dat.mat$prev_lci) #still just one
subset(dat, sampleid=="WAY045") #also seropositive to GhV

# it is clearly possible to inherit matABs to AngV but not necessarily a sure thing
# so we will fit one rate of waning immunity and not guarantee any classes that 
# are protected. 



################################################################################
################################################################################

#now  fit adult Ab trajectory to data

# first, collect all the adults that have a 3-prong peak


#load the recaptures
catch.dat <- read.csv(paste0(batwd,"/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)
head(catch.dat)
recap = subset(catch.dat, recap_new=="yes")
head(recap)
names(recap)
#slim to columns of interest
recap <- dplyr::select(recap, barcode_id, sampleid, roost_site, processing_date, bat_species, bat_age_class, bat_sex, mass_forearm_residual)
head(recap)


#just take rousettus and eidolon
recap=subset(recap, bat_species!="Pteropus rufus")
recap$processing_date <- as.Date(recap$processing_date)
min(recap$processing_date)#"2013-11-24"
max(recap$processing_date)#"2026-01-24"



#and now load the new luminex data 
lum.dat <- read.csv(file = paste0(homewd, "/serology/data/mix_sub_data_75_90_96_long.csv"))
head(lum.dat)

#slim
lum.dat <- dplyr::select(lum.dat, sampleid,bat_species, antigen, MFI, seropos, seropos_lci, seropos_uci)


head(lum.dat)
length(unique(lum.dat$sampleid)) #1055
head(recap)

lum.dat = subset(lum.dat, bat_species=="Eidolon dupreanum" & antigen=="AngV")
recap = subset(recap, bat_species=="Eidolon dupreanum")

recap <- merge(recap, lum.dat, by=c("bat_species", "sampleid"), all.x = T)
head(recap)
min(recap$processing_date)#"2013-11-24"
max(recap$processing_date)#"2026-01-22"


recap = subset(recap, antigen=="AngV")

# split by antigen and barcode and make sure there are multiple timepoints. 
# color as bats that stay neg, seroconvert, serorevert or do both
id.recap <- function(df){
  if(nrow(df)>1){
    df <- arrange(df, processing_date)
    #df.new <- dplyr::select(df, -(keep))
    df.new <- df
    df.new$class <- NA
    if(sum(as.numeric(as.character(df$seropos)))==0){
      df.new$class <- "always seronegative" # all seroneg
    }else if(sum(as.numeric(as.character(df$seropos)))==length(df$seropos)){
      df.new$class <- "always seropositive" #all seropos
    }else{
      vec.sero <- unique(diff(as.numeric(as.character(df$seropos))))
      if(sum(unique(vec.sero))==0){
        df.new$class <- "seroconversion-seroreversion" #both
      }else if((sum(unique(match(vec.sero,-1)),na.rm=T)==0) & (sum(unique(match(vec.sero,1)),na.rm=T)==1)){
        df.new$class <- "seroconversion" #both
      }else if((sum(unique(match(vec.sero,1)),na.rm=T)==0) & (sum(unique(match(vec.sero,-1)),na.rm=T)==1)){
        df.new$class <- "seroreversion" #both
      }
    }
    df.new$N_cap<- nrow(df.new)
    return(df.new)
  }
  
}
recap.list <- dlply(recap, .(antigen, barcode_id))

recap.sum <- data.table::rbindlist(lapply(recap.list, id.recap))

recap.sum$seropos <- as.factor(recap.sum$seropos)

#and plot
recap.sum$is_seroneg <- 1
recap.sum$is_seroneg[recap.sum$class!="always seronegative"] <- 0
recap.sum$is_seroneg <- as.factor(recap.sum$is_seroneg)

unique(recap.sum$class)
recap.sum$seropos <- as.character(recap.sum$seropos)
recap.sum$seropos[recap.sum$seropos==0] <- "seronegative"
recap.sum$seropos[recap.sum$seropos==1] <- "seropositive"
recap.sum$seropos <- factor(recap.sum$seropos, levels=c("seronegative", "seropositive"))

filz = c('seronegative' = "navy",'seropositive'= "violetred")

colorz1 = c('always seronegative' = "navy", 
            "always seropositive" = "violetred4",
            "seroconversion-seroreversion"="goldenrod2", 
            "seroreversion" = "lightseagreen",
            "seroconversion" = "orangered")

shapez = c('female' = 21, 'male' = 24)
alphaz <- c('0' = 1, '1'=.1)


recap.sum$class <- factor(recap.sum$class, levels = c("always seronegative", "always seropositive", "seroconversion","seroreversion", "seroconversion-seroreversion" ))

Eid.recap <- ggplot(subset(recap.sum, bat_species=="Eidolon dupreanum")) + theme_bw()+ 
  geom_line(aes(x=processing_date, y=barcode_id, group=barcode_id, color=class, alpha=is_seroneg), size=1) + 
  geom_point(aes(x=processing_date, y=barcode_id, group=barcode_id, shape=bat_sex, fill=seropos, alpha=is_seroneg), size=3, color="gray65") + 
  scale_shape_manual(name="bat sex", values = shapez) + 
  scale_fill_manual(name = "serostatus", values = filz) +
  scale_color_manual(name= "serological\ntrajectory", values = colorz1) + ylab("individual barcode id\n\n") +
  facet_nested(bat_sex~antigen, scales = "free", space = "free") + scale_alpha_manual(values = alphaz) +
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        legend.box = "vertical",
        legend.position = "left",
        panel.spacing.y = unit(.1, "lines"),
        panel.spacing.x = unit(1, "lines"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        plot.margin = unit(c(.3,.3,0,2),"lines"),
        axis.text.x = element_text(size=14),  strip.text.x = element_text(size=16),
        strip.text.y = element_text(size=16),#, face="italic"),
        panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
  scale_x_continuous(breaks = c(as.Date("2014-07-01"),as.Date("2015-07-01"), as.Date("2016-07-01"), as.Date("2017-07-01"), as.Date("2018-07-01"), as.Date("2019-07-01")),labels = c(2014,2015,2016,2017,2018,2019)) +
  guides(alpha="none",  shape=guide_legend(order=1, override.aes = list(size=3, color="black")), fill=guide_legend(order=2, override.aes = list(shape=21,size=3)), color=guide_legend(order=3, ncol=1));Eid.recap




head(recap.sum)

#get those with < 1 year bewteen sampling dates
recap.list <- dlply(recap.sum, .(barcode_id))

count_spacing <- function(df, time_days){
  df <- arrange(df, processing_date)
  diff= c(diff(df$processing_date))
  
  df$days_diff <- c(NA, diff)
  
  df$time_series <- cumsum(c(1,diff))
  if(sum(df$days_diff[!is.na(df$days_diff)]<time_days)){
    return(df)
  }
  
}

recap.short <-  data.table::rbindlist(lapply(recap.list, count_spacing, time_days=6*31))

#recap.short <- subset(recap.short, N_cap>2)


recap.short$simple_class <- recap.short$class
recap.short$simple_class <- as.character(recap.short$simple_class)
recap.short$simple_class[recap.short$simple_class!="always seronegative"] <- "dynamic"
recap.short = subset(recap.short, simple_class=="dynamic")
alphaz2 = c("always seronegative" = .3, "dynamic" = 1 )

recap.short$era <- "2014-2016"
recap.short$era[recap.short$processing_date>as.Date("2018-01-01")] <- "2018-2020"

# add <- recap.short[1,]
# add$era <- "2014-2016"
# add$processing_date <- as.Date("2017-01-01")
# add$MFI <- NA
# add$antigen <- "NiV"
# add2 <- add
# add2$antigen <- "HeV"
# add <- rbind(add, add2)
# add2 <- add
# add2$era <- "2018-2020"
# add2$processing_date <- as.Date("2017-01-02")
# 
# recap.short <- rbind(recap.short, add, add2)

shapez = c("female" = 21, "male" =24)

#and include plot of the cutoff line by species/antigen/era
cutoff.df <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 6), rep("Rousettus madagascariensis",3)), antigen=c("NiV", "HeV", "NiV", "HeV", "GhV", "AngV", "GhV", "BdbV", "RavV"), era=c(rep("2014-2016",2), rep("2018-2020",7)))
cutoff.df$MFI = c(757,215,936,2851,1141,1763,2642,5647,6617)
cutoff.df$antigen <- factor(cutoff.df$antigen, levels=c("NiV", "HeV", "GhV", "AngV", "BdbV", "RavV"))
cutoff.df = subset(cutoff.df, antigen=="GhV")

EidMFIall <- ggplot(data=subset(recap.short, bat_species=="Eidolon dupreanum")) + 
  facet_nested(bat_sex~., scales = "free", space = "free") + theme_bw() +
  scale_shape_manual(name="bat sex", values = shapez) + 
  scale_alpha_manual(values = alphaz2) + ylab("mean fluorescence intensity (MFI)")+
  theme(axis.title.y = element_text(size=16), axis.title.x = element_blank(),
        legend.box = "horizontal", strip.text.y = element_text(size=16),#, face = "italic"),
        strip.background.y = element_rect(fill="white"),
        strip.background.x = element_blank(),
        axis.text = element_text(size=14), 
        panel.spacing.y = unit(.1, "lines"),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,.3,.3,.3),"lines"),
        panel.spacing.x = unit(c(0,1,0,1,1),"lines"),
        panel.grid = element_blank()) +
  geom_hline(data = subset(cutoff.df, bat_species=="Eidolon dupreanum"), aes(yintercept = MFI), color="red", linetype=2) + 
  geom_line(aes(x=processing_date, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1)+
  geom_point(aes(x=processing_date, y=MFI, fill=barcode_id, shape = bat_sex, alpha = simple_class), show.legend = F, size=3) +
  scale_x_continuous(breaks = c(as.Date("2014-07-01"),as.Date("2015-07-01"), as.Date("2016-07-01"), as.Date("2017-07-01"), as.Date("2018-07-01"), as.Date("2019-07-01")),labels = c(2014,2015,2016,2017,2018,2019)) +
  guides(alpha="none",  shape=guide_legend(order=1));EidMFIall


#now try restricting to only triptics that go "low-high-low". center on the peak
#also need the first and second point to be within 6 months of each other
recap.list <- dlply(recap.short, .(barcode_id))
find.tryptic <- function(df){
  if(nrow(df)>=3){
    
    if(df$MFI[1]<df$MFI[2]){
      df.out = df[1:3,]
      
      
      df.out$new_ts <- NA
      df.out$new_ts[df.out$MFI==max(df.out$MFI)] <- 0
      
      df.out$new_ts[is.na(df.out$new_ts)] <- df.out$time_series[is.na(df.out$new_ts)] + -df.out$time_series[df.out$MFI==max(df.out$MFI)]
      
      return(df.out)
    }else if(df$MFI[2]<df$MFI[3] & nrow(df)>3){
      df.out = df[2:nrow(df),]
      
      #remove anything separated by >6 months and recompute diff
      df.out =subset(df.out, days_diff<(6*31))
      df.out$days_diff <- c(NA, diff(df.out$processing_date))
      
      if(sum(df.out$days_diff[!is.na(df.out$days_diff)]<(6*31))>1){ #otherwise you lose it
        
        
        
        df.out$new_ts <- NA
        df.out$new_ts[df.out$MFI==max(df.out$MFI)] <- 0
        
        df.out$new_ts[is.na(df.out$new_ts)] <- df.out$time_series[is.na(df.out$new_ts)] + -df.out$time_series[df.out$MFI==max(df.out$MFI)]
        
        return(df.out)
      }     }
    
    
  }
  
  
}


triptic_df <- data.table::rbindlist(lapply(recap.list, find.tryptic))

triptic_df$seropos <- as.character(triptic_df$seropos)

tryptic.plot <- ggplot(data=triptic_df) + 
  facet_nested(bat_sex~., scales = "free", space = "free") + theme_bw() +
  scale_shape_manual(name="bat sex", values = shapez) + 
  scale_alpha_manual(values = alphaz2) + ylab("mean fluorescence intensity (MFI)")+
  theme(axis.title.y = element_text(size=16), axis.title.x = element_blank(),
        legend.box = "horizontal", strip.text.y = element_text(size=16),#, face = "italic"),
        strip.background.y = element_rect(fill="white"),
        strip.background.x = element_blank(),
        axis.text = element_text(size=14), 
        panel.spacing.y = unit(.1, "lines"),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,.3,.3,.3),"lines"),
        panel.spacing.x = unit(c(0,1,0,1,1),"lines"),
        panel.grid = element_blank()) +
  geom_hline(data = subset(cutoff.df, bat_species=="Eidolon dupreanum"), aes(yintercept = MFI), color="red", linetype=2) + 
  geom_line(aes(x=processing_date, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1)+
  geom_point(aes(x=processing_date, y=MFI, fill=barcode_id, shape = bat_sex, alpha = simple_class), show.legend = F, size=3) +
  scale_x_continuous(breaks = c(as.Date("2014-07-01"),as.Date("2015-07-01"), as.Date("2016-07-01"), as.Date("2017-07-01"), as.Date("2018-07-01"), as.Date("2019-07-01")),labels = c(2014,2015,2016,2017,2018,2019)) +
  guides(alpha="none",  shape=guide_legend(order=1));tryptic.plot

#and plot with centered peak
tryptic.plot <- ggplot(data=triptic_df) + 
  facet_nested(bat_sex~., scales = "free", space = "free") + theme_bw() +
  scale_shape_manual(name="bat sex", values = shapez) + 
  xlab("Time relative to observed peak (t_rel)") +
  ylab("MFI") +
  scale_alpha_manual(values = alphaz2) + ylab("mean fluorescence intensity (MFI)")+
  theme(axis.title.y = element_text(size=16), 
        legend.box = "horizontal", strip.text.y = element_text(size=16),#, face = "italic"),
        strip.background.y = element_rect(fill="white"),
        strip.background.x = element_blank(),
        axis.text = element_text(size=14), 
        panel.spacing.y = unit(.1, "lines"),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,.3,.3,.3),"lines"),
        panel.spacing.x = unit(c(0,1,0,1,1),"lines"),
        panel.grid = element_blank()) +
  geom_hline(data = subset(cutoff.df, bat_species=="Eidolon dupreanum"), aes(yintercept = MFI), color="red", linetype=2) + 
  geom_line(aes(x=new_ts, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1)+
  geom_point(aes(x=new_ts, y=MFI, fill=barcode_id, shape = bat_sex, alpha = simple_class), show.legend = F, size=3) +
  guides(alpha="none",  shape=guide_legend(order=1));tryptic.plot

#interpolate that curve for the top male peak

plot1 <- ggplot(data=subset(triptic_df, barcode_id=="956000010868919")) + 
  facet_nested(bat_sex~., scales = "free", space = "free") + theme_bw() +
  scale_shape_manual(name="bat sex", values = shapez) + 
  xlab("Time relative to observed peak (t_rel)") +
  ylab("MFI") +
  scale_alpha_manual(values = alphaz2) + ylab("mean fluorescence intensity (MFI)")+
  theme(axis.title.y = element_text(size=16), 
        legend.box = "horizontal", strip.text.y = element_text(size=16),#, face = "italic"),
        strip.background.y = element_rect(fill="white"),
        strip.background.x = element_blank(),
        axis.text = element_text(size=14), 
        panel.spacing.y = unit(.1, "lines"),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,.3,.3,.3),"lines"),
        panel.spacing.x = unit(c(0,1,0,1,1),"lines"),
        panel.grid = element_blank()) +
  geom_hline(data = subset(cutoff.df, bat_species=="Eidolon dupreanum"), aes(yintercept = MFI), color="red", linetype=2) + 
  geom_line(aes(x=new_ts, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1)+
  geom_point(aes(x=new_ts, y=MFI, fill=barcode_id, shape = bat_sex, alpha = simple_class), show.legend = F, size=3) +
  guides(alpha="none",  shape=guide_legend(order=1));plot1

pre.MFI <- lm(MFI~new_ts, data = subset(triptic_df, barcode_id=="956000010868919" & new_ts<1))


#MFI = 1141 is the cutoff for seropositivity
add.pseudo <- triptic_df[triptic_df$barcode_id=="956000010868919",]
add.pseudo <- add.pseudo[1,]
add.pseudo$MFI <- c(predict(pre.MFI, newdata = data.frame(new_ts = -110)))
add.pseudo$new_ts <- c(-110)
add.pseudo$seropos <- "seronegative"

triptic_df <- rbind(triptic_df, add.pseudo)


# and plot with centered peak again 
tryptic.plot2 <- ggplot(data=triptic_df) + 
  facet_nested(bat_sex~., scales = "free", space = "free") + theme_bw() +
  scale_shape_manual(name="bat sex", values = shapez) + 
  xlab("Time relative to observed peak (t_rel)") +
  ylab("MFI") +
  scale_alpha_manual(values = alphaz2) + ylab("mean fluorescence intensity (MFI)")+
  theme(axis.title.y = element_text(size=16), 
        legend.box = "horizontal", strip.text.y = element_text(size=16),#, face = "italic"),
        strip.background.y = element_rect(fill="white"),
        strip.background.x = element_blank(),
        axis.text = element_text(size=14), 
        panel.spacing.y = unit(.1, "lines"),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,.3,.3,.3),"lines"),
        panel.spacing.x = unit(c(0,1,0,1,1),"lines"),
        panel.grid = element_blank()) +
  geom_hline(data = subset(cutoff.df, bat_species=="Eidolon dupreanum"), aes(yintercept = MFI), color="red", linetype=2) + 
  geom_line(aes(x=new_ts, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1)+
  geom_point(aes(x=new_ts, y=MFI, fill=barcode_id, shape = bat_sex, alpha = simple_class), show.legend = F, size=3) +
  guides(alpha="none",  shape=guide_legend(order=1)) +
  geom_vline(aes(xintercept = 28))+ #duration 2 seropositive timesteps
  geom_vline(aes(xintercept = 42))+ #duration 3 seropositive timesteps
  geom_vline(aes(xintercept = 135)) #duration 9 seropositive timesteps - this is the mean
tryptic.plot2

#save for stan and include marker of the last seroneg and first seropos for each bat

names(triptic_df)
triptic.list <- dlply(triptic_df,.(barcode_id))

prep.data.stan <- function(df){
  df$t_obs <- df$new_ts + abs(min(df$new_ts))+1
  df$last_neg <- max(df$t_obs[df$new_ts<0 & df$seropos=="seronegative"])
  df$first_pos <- min(df$t_obs[df$seropos=="seropositive"])
  
  return(df)
}

triptic_df <- data.table::rbindlist(lapply(triptic.list, prep.data.stan))

#now run stan
# Map barcode_id to consecutive integers for Stan
bat_levels <- unique(triptic_df$barcode_id)
triptic_df <- triptic_df %>%
  mutate(bat_idx = as.numeric(factor(barcode_id, levels = bat_levels)))

# Prepare Stan data
stan_data <- list(
  N = nrow(triptic_df),
  N_bats = length(bat_levels),
  bat_id = triptic_df$bat_idx,
  t_obs = triptic_df$t_obs,
  MFI = triptic_df$MFI,
  t_exp_min = triptic_df$last_neg,
  t_exp_max = triptic_df$first_pos
)

# Fit Stan model
fit <- stan(file = paste0(homewd,"/model_fitting/GhV-fits/latent_exposure_model.stan"),
            data = stan_data,
            iter = 10000, chains = 4, cores = 4, warmup = 2000,
            control = list(adapt_delta = 0.99))
print(fit, pars = c("A", "B", "sigma", "t_exp"), probs = c(0.05, 0.5, 0.95))

posterior <- as.data.frame(fit)

#get the hyperparameters (population-level)
hyperpars <- posterior %>%
  dplyr::select(mu_A, sigma_A, mu_B, sigma_B, sigma)

# quick summary
summary(hyperpars)

#get the bat-specific parameters (A and B for each bat)
bat_params <- posterior %>%
  dplyr::select(matches("^A\\[|^B\\["))  # matches A[1], A[2], ... B[1], B[2], ...

# latent exposure times
t_exp_post <- posterior %>%
  dplyr::select(matches("^t_exp\\["))     # t_exp[1], t_exp[2], ...

bat_summary <- bat_params %>%
  dplyr::summarise(across(everything(), list(
    mean = ~mean(.),
    median = ~median(.),
    q05 = ~quantile(., 0.05),
    q95 = ~quantile(., 0.95)
  )))

head(bat_summary)

# to plot the individual bat fits:

sim_fun <- function(A, B, t_exp_val, sigma, t_seq) {
  tau <- pmax(t_seq - t_exp_val, 0)
  A * tau * exp(-B * tau) + rnorm(length(tau), 0, sigma)
}

# simulate trajectory for bat 1
sim_traj <- sim_fun(A=mean(posterior$`A[1]`), B=mean(posterior$`B[1]`), t_exp_val=mean(posterior$`t_exp[1]`), sigma=mean(posterior$sigma), t_seq=0:150)

plot(0:150, sim_traj) #trajectory for a single bat
t_seq = 0:150

sim_list <- lapply(1:4, function(i) {
  
  A_hat <- mean(posterior[[paste0("A[", i, "]")]])
  B_hat <- mean(posterior[[paste0("B[", i, "]")]])
  sigma_hat <- mean(posterior[["sigma"]])
  t_exp_val <- mean(posterior[[paste0("t_exp[", i, "]")]])
  data.frame(
    t = t_seq,
    MFI_sim = sim_fun(A_hat, B_hat, t_exp_val, sigma_hat, t_seq),
    bat = i
  )
})

sim_df <- do.call(rbind, sim_list)

stan_df <- as.data.frame(stan_data)

p.ind.bat <- ggplot(data=sim_df)+
  geom_line(aes(x=t, y=MFI_sim,group=bat, color=as.factor(bat) )) + theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "inside", legend.position.inside = c(.85,.85)) +
  geom_point(data=stan_df, (aes(x=t_obs, y=MFI, color=as.factor(bat_id)))); p.ind.bat

#now simulate many new bats from the population-level hyperparameters

hyperpars <- posterior %>%
  dplyr::select(mu_A, sigma_A, mu_B, sigma_B, sigma)



simulate_new_bat_stochastic <- function(mu_A, sigma_A, mu_B, sigma_B, sigma, t_seq, t_exp = 0, MFI_min = 0, MFI_max = 50000, shrink_factor) {
  
  # draw bat-level parameters from population
  #shrink_factor = 0 is fully shrunk to mean
  #shrink_factor = 1 is full sigma variance
  A <- rnorm(1, mean = mu_A, sd = sigma_A*shrink_factor)
  B <- rnorm(1, mean = mu_B, sd = sigma_B*shrink_factor)
  
  tau <- pmax(t_seq - t_exp, 0)
  
  # simulate raw MFI
  MFI <- A * tau * exp(-B * tau) + rnorm(length(tau), 0, sigma)
  
  # constrain MFI to [MFI_min, MFI_max]
  MFI <- pmin(pmax(MFI, MFI_min), MFI_max)
  
  return(MFI)
}


n_new_bats <- 300
t_seq <- seq(0, 150, by = 1)   # time in days

# sample 300 posterior rows for variation
posterior_rows <- hyperpars[sample(nrow(hyperpars), n_new_bats, replace = TRUE), ]

sim_list <- lapply(1:n_new_bats, function(i) {
  row <- posterior_rows[i, ]
  data.frame(
    t = t_seq,
    MFI_sim = simulate_new_bat_stochastic(row$mu_A, row$sigma_A, row$mu_B, row$sigma_B, row$sigma, t_seq, shrink_factor = 0),
    bat = i
  )
})


sim_df <- do.call(rbind, sim_list)

#and the deterministic
simulate_new_bat_deterministic <- function(A, B, t_seq, t_exp = 0) {
  
  A <- mean(posterior$mu_A)
  B <- mean(posterior$mu_B)
  t_exp <- 0   # or mean(posterior$t_exp[bat_id]) if you want average exposure
  
  
  tau <- pmax(t_seq - t_exp, 0)
  
  # simulate raw MFI
  MFI <- A * tau * exp(-B * tau) 
  
  
  
  return(MFI)
}

det.bat <- cbind.data.frame(t=t_seq,
                            MFI=simulate_new_bat_deterministic(A=row$mu_A, B=row$mu_B, t_seq))


ggplot() +
  geom_line(data= sim_df, aes(x = t, y = MFI_sim, group = bat),alpha = 0.2, color = "steelblue") +
  labs(title = "Simulated antibody trajectories for 300 new bats",
       x = "Days since exposure",
       y = "MFI") +
  geom_line(data=det.bat, aes(x=t, y=MFI),size=1) +
  theme(panel.grid = element_blank()) +
  theme_bw() + coord_cartesian(ylim=c(0,40000)) + 
  geom_hline(aes(yintercept =cutoff.df$MFI[cutoff.df$bat_species=="Eidolon dupreanum"]), linetype=2)

min(det.bat$t[det.bat$MFI<cutoff.df$MFI[cutoff.df$bat_species=="Eidolon dupreanum"] & det.bat$t>25])
# 119 days (compared with 135 estimated for GhV) duration Ab immunity for AngV
# but there is a lot of variation in the stochastic runs
# that's 3.8 months - but we start seeing many dropoff at ~60 days (2 months)

# just like with GhV, we'll keep two protected classes for AngV R class (duration of 28
# days of total protection) and then estimate the rate of waning immunity out of the third.