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


dat <- read.csv(paste0(homewd,"/model_fitting/GhV-fits/matrix-fit-1a/age-seroprevalence-data.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)[names(dat)=="species"] <- "bat_species"
dat$log10MFI <- log10(dat$MFI)

# catch.dat <- read.csv(file = paste0(homewd, "/master_data_files/catching_data.csv"), header = T, stringsAsFactors = F)
# lum.dat <- read.csv(file = paste0(homewd, "/serology/data/mix_sub_data_75_90_96_long.csv"))
# 
# #now  fit adult Ab trajectory to data
# #this functions does it all and saves the data for plotting
# 
#   
#   #load the recaptures
#   recap = subset(catch.dat, recap_new=="yes")
#   # head(recap)
#   # names(recap)
#   #slim to columns of interest
#   recap <- dplyr::select(recap, barcode_id, sampleid, roost_site, processing_date, bat_species, bat_age_class, bat_sex, mass_forearm_residual)
#   # head(recap)
#   
#   
#   #just take  eidolon
#   recap=subset(recap, bat_species=="Eidolon dupreanum")
#   recap$processing_date <- as.Date(recap$processing_date)
#   # min(recap$processing_date)#"2013-11-24"
#   # max(recap$processing_date)#"2026-01-24"
#   
#   
#   
#   #and now load the new luminex data 
#   
#   #slim
#   lum.dat <- dplyr::select(lum.dat, sampleid,bat_species, antigen, MFI, seropos, seropos_lci, seropos_uci)
#   
#   #head(lum.dat)
#   #length(unique(lum.dat$sampleid)) #1055
#   #head(recap)
#   
#   lum.dat = subset(lum.dat, bat_species=="Eidolon dupreanum" & antigen=="GhV")
#   
#   
#   recap <- merge(recap, lum.dat, by=c("bat_species", "sampleid"), all.x = T)
#   #head(recap)
#   # min(recap$processing_date)#"2013-11-24"
#   # max(recap$processing_date)#"2026-01-22"
#   # 
#   
#   recap = subset(recap, antigen=="GhV")
#   
#   # split by antigen and barcode and make sure there are multiple timepoints. 
#   # color as bats that stay neg, seroconvert, serorevert or do both
#   id.recap <- function(df){
#     if(nrow(df)>1){
#       df <- arrange(df, processing_date)
#       #df.new <- dplyr::select(df, -(keep))
#       df.new <- df
#       df.new$class <- NA
#       if(sum(as.numeric(as.character(df$seropos)))==0){
#         df.new$class <- "always seronegative" # all seroneg
#       }else if(sum(as.numeric(as.character(df$seropos)))==length(df$seropos)){
#         df.new$class <- "always seropositive" #all seropos
#       }else{
#         vec.sero <- unique(diff(as.numeric(as.character(df$seropos))))
#         if(sum(unique(vec.sero))==0){
#           df.new$class <- "seroconversion-seroreversion" #both
#         }else if((sum(unique(match(vec.sero,-1)),na.rm=T)==0) & (sum(unique(match(vec.sero,1)),na.rm=T)==1)){
#           df.new$class <- "seroconversion" #both
#         }else if((sum(unique(match(vec.sero,1)),na.rm=T)==0) & (sum(unique(match(vec.sero,-1)),na.rm=T)==1)){
#           df.new$class <- "seroreversion" #both
#         }
#       }
#       df.new$N_cap<- nrow(df.new)
#       return(df.new)
#     }
#     
#   }
#   recap.list <- dlply(recap, .(antigen, barcode_id))
#   
#   recap.sum <- data.table::rbindlist(lapply(recap.list, id.recap))
#   
#   recap.sum$seropos <- as.factor(recap.sum$seropos)
#   
#   #and plot
#   recap.sum$is_seroneg <- 1
#   recap.sum$is_seroneg[recap.sum$class!="always seronegative"] <- 0
#   recap.sum$is_seroneg <- as.factor(recap.sum$is_seroneg)
#   
#   #unique(recap.sum$class)
#   recap.sum$seropos <- as.character(recap.sum$seropos)
#   recap.sum$seropos[recap.sum$seropos==0] <- "seronegative"
#   recap.sum$seropos[recap.sum$seropos==1] <- "seropositive"
#   recap.sum$seropos <- factor(recap.sum$seropos, levels=c("seronegative", "seropositive"))
#   
#   
#   #get those with < 1 year bewteen sampling dates
#   recap.list <- dlply(recap.sum, .(barcode_id))
#   
#   count_spacing <- function(df, time_days){
#     df <- arrange(df, processing_date)
#     diff= c(diff(df$processing_date))
#     
#     df$days_diff <- c(NA, diff)
#     
#     df$time_series <- cumsum(c(1,diff))
#     if(sum(df$days_diff[!is.na(df$days_diff)]<time_days)){
#       return(df)
#     }
#     
#   }
#   
#   recap.short <-  data.table::rbindlist(lapply(recap.list, count_spacing, time_days=6*31))
#   
#   #recap.short <- subset(recap.short, N_cap>2)
#   
#   
#   recap.short$simple_class <- recap.short$class
#   recap.short$simple_class <- as.character(recap.short$simple_class)
#   recap.short$simple_class[recap.short$simple_class!="always seronegative"] <- "dynamic"
#   recap.short = subset(recap.short, simple_class=="dynamic")
#   #alphaz2 = c("always seronegative" = .3, "dynamic" = 1 )
#   
#   recap.short$era <- "2014-2016"
#   recap.short$era[recap.short$processing_date>as.Date("2018-01-01")] <- "2018-2020"
#   
#   
#   
#   #now try restricting to only triptics that go "low-high-low". center on the peak
#   #also need the first and second point to be within 6 months of each other
#   recap.list <- dlply(recap.short, .(barcode_id))
#   find.tryptic <- function(df){
#     if(nrow(df)>=3){
#       
#       if(df$MFI[1]<df$MFI[2]){
#         df.out = df[1:3,]
#         
#         
#         df.out$new_ts <- NA
#         df.out$new_ts[df.out$MFI==max(df.out$MFI)] <- 0
#         
#         df.out$new_ts[is.na(df.out$new_ts)] <- df.out$time_series[is.na(df.out$new_ts)] + -df.out$time_series[df.out$MFI==max(df.out$MFI)]
#         
#         return(df.out)
#       }else if(df$MFI[2]<df$MFI[3] & nrow(df)>3){
#         df.out = df[2:nrow(df),]
#         
#         #remove anything separated by >6 months and recompute diff
#         df.out =subset(df.out, days_diff<(6*31))
#         df.out$days_diff <- c(NA, diff(df.out$processing_date))
#         
#         if(sum(df.out$days_diff[!is.na(df.out$days_diff)]<(6*31))>1){ #otherwise you lose it
#           
#           
#           
#           df.out$new_ts <- NA
#           df.out$new_ts[df.out$MFI==max(df.out$MFI)] <- 0
#           
#           df.out$new_ts[is.na(df.out$new_ts)] <- df.out$time_series[is.na(df.out$new_ts)] + -df.out$time_series[df.out$MFI==max(df.out$MFI)]
#           
#           return(df.out)
#         }     }
#       
#       
#     }
#     
#     
#   }
#   
#   
#   triptic_df <- data.table::rbindlist(lapply(recap.list, find.tryptic))
#   
#   triptic_df$seropos <- as.character(triptic_df$seropos)
#   
#   
#   pre.MFI <- lm(MFI~new_ts, data = subset(triptic_df, barcode_id=="956000006105427" & new_ts<1))
#   post.MFI <- lm(MFI~new_ts, data = subset(triptic_df, barcode_id=="956000006105427" & new_ts>-1))
#   # 
#   #MFI = 1141 is the cutoff for seropositivity
#   add.pseudo <- triptic_df[triptic_df$barcode_id=="956000006105427",]
#   add.pseudo <- add.pseudo[1:2,]
#   add.pseudo$MFI <- c(predict(pre.MFI, newdata = data.frame(new_ts = -90)), predict(post.MFI, newdata = data.frame(new_ts = 100)))
#   add.pseudo$new_ts <- c(-90, 100)
#   add.pseudo$seropos <- "seronegative"
#   
#   triptic_df <- rbind(triptic_df, add.pseudo)
#   
#   
#   #names(triptic_df)
#   triptic.list <- dlply(triptic_df,.(barcode_id))
#   
#   prep.data.stan <- function(df){
#     df$t_obs <- df$new_ts + abs(min(df$new_ts))+1
#     df$last_neg <- max(df$t_obs[df$new_ts<0 & df$seropos=="seronegative"])
#     df$first_pos <- min(df$t_obs[df$seropos=="seropositive"])
#     
#     return(df)
#   }
#   
#   triptic_df <- data.table::rbindlist(lapply(triptic.list, prep.data.stan))
#   
#   #now run stan
#   # Map barcode_id to consecutive integers for Stan
#   bat_levels <- unique(triptic_df$barcode_id)
#   triptic_df <- triptic_df %>%
#     mutate(bat_idx = as.numeric(factor(barcode_id, levels = bat_levels)))
#   
#   
#below is all exploratory and produces some temporary plots, but above is all that is needed

# first, collect all the adults that have a 3-prong peak


#load the recaptures
catch.dat <- read.csv(paste0(homewd,"/master_data_files/catching_data.csv"), header = T, stringsAsFactors = F)
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

lum.dat = subset(lum.dat, bat_species=="Eidolon dupreanum" & antigen=="GhV")
recap = subset(recap, bat_species=="Eidolon dupreanum")

recap <- merge(recap, lum.dat, by=c("bat_species", "sampleid"), all.x = T)
head(recap)
min(recap$processing_date)#"2013-11-24"
max(recap$processing_date)#"2026-01-22"


recap = subset(recap, antigen=="GhV")

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

#interpolate that curve for the top peak

plot1 <- ggplot(data=subset(triptic_df, barcode_id=="956000006105427")) + 
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

pre.MFI <- lm(MFI~new_ts, data = subset(triptic_df, barcode_id=="956000006105427" & new_ts<1))
post.MFI <- lm(MFI~new_ts, data = subset(triptic_df, barcode_id=="956000006105427" & new_ts>-1))

#MFI = 1141 is the cutoff for seropositivity
add.pseudo <- triptic_df[triptic_df$barcode_id=="956000006105427",]
add.pseudo <- add.pseudo[1:2,]
add.pseudo$MFI <- c(predict(pre.MFI, newdata = data.frame(new_ts = -90)), predict(post.MFI, newdata = data.frame(new_ts = 100)))
add.pseudo$new_ts <- c(-90, 100)
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
  guides(alpha="none",  shape=guide_legend(order=1)) 
  # geom_vline(aes(xintercept = 28))+ #duration 2 seropositive timesteps
  # geom_vline(aes(xintercept = 42))+ #duration 3 seropositive timesteps
  # geom_vline(aes(xintercept = 135)) #duration 9 seropositive timesteps - this is the mean
tryptic.plot2


# and plot with centered peak again 
#and add in a way to distinguish the extrapolated points from the observed
triptic_df$bat_idx <- as.numeric(as.factor(triptic_df$barcode_id))
triptic_df$status <- "observed"
triptic_df$status[triptic_df$bat_idx==2 & triptic_df$new_ts < -80] <- "inferred"
triptic_df$status[triptic_df$bat_idx==2 & triptic_df$new_ts > 80] <- "inferred"

shapez <- c('1'=16, '2'=17, '3'=15, '4' = 3)
tryptic.plot3 <- ggplot(data=triptic_df) + 
  #facet_nested(bat_sex~., scales = "free", space = "free") + 
  theme_bw() +
  scale_shape_manual( values = shapez) + 
  xlab("Time relative to observed peak (days)") +
  ylab("MFI") +
  scale_alpha_manual(values = alphaz2) + 
  ylab("mean fluorescence intensity (MFI)")+
  theme(axis.title = element_text(size=18), 
        legend.box = "horizontal", strip.text.y = element_text(size=16),#, face = "italic"),
        strip.background.y = element_rect(fill="white"),
        strip.background.x = element_blank(),
        axis.text = element_text(size=14), 
        panel.spacing.y = unit(.1, "lines"),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,.3,.3,.3),"lines"),
        panel.spacing.x = unit(c(0,1,0,1,1),"lines"),
        panel.grid = element_blank()) +
  geom_hline(data = subset(cutoff.df, bat_species=="Eidolon dupreanum"), aes(yintercept = MFI), linetype=2) + 
  geom_line(aes(x=new_ts, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1)+
  geom_point(data = subset(triptic_df, status=="inferred"), aes(x=new_ts, y=MFI, shape = as.factor(bat_idx)), show.legend = F, size=4, color="black") +
  geom_point(aes(x=new_ts, y=MFI, color=barcode_id, shape = as.factor(as.numeric(bat_idx)), alpha = simple_class), show.legend = F, size=3) +
  guides(alpha="none",  shape=guide_legend(order=1)) #+
tryptic.plot3

#
ggsave(file = paste0(homewd, "/final_figures_tables/supplemental/SuppFig7.png"),
       plot = tryptic.plot3,
       units="mm",  
       width=65, 
       height=50, 
       scale=2.5, 
       dpi=300)

