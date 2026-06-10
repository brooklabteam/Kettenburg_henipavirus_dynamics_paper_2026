rm(list=ls())
#set wd
#.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4")

library(epitools)
library(plyr)
library(dplyr)
library(cowplot)
library(deSolve)
library(mgcv)
library(lubridate)
library(matrixcalc)
library(Matrix)
library(ggplot2)
library(mvtnorm)

#load data
homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

#load the completed models with age-seroprevalence and time series
age.sero.df <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/final_age_sero_df_GhV.csv"), header = T, stringsAsFactors = F)
biwk.ts.df <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/final_biweekly_ts_GhV.csv"), header = T, stringsAsFactors = F)

#and load the data
dat <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/fit-1/age-seroprevalence-data.csv"), stringsAsFactors = F, header = T)
head(dat)
sero.dat <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/fit-1/biweekly-serology-data.csv"), stringsAsFactors = F, header = T)
shed.dat <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/fit-1/biweekly-shedding-data.csv"), stringsAsFactors = F, header = T)
head(dat)
head(sero.dat)
head(shed.dat)

#plotting scripts
add.age <- function(df){
  df <- arrange(df, age)
  bins = c(c(0,.5,1,2), c(3:20))
  
  df$age_bin <- NA
  for (i in 1:length(df$age)){
    df$age_bin[i] <- max(bins[df$age[i]>=bins], na.rm = T)
  }
  
  
  
  #you want the midpoint between either end of the age class now
  vect.age.yr = sort(unique(df$age_bin))
  vect.age.plot =   rep(NA, (length(vect.age.yr)-1))
  for (i in 1:length(vect.age.plot)){
    vect.age.plot[i] <- ((vect.age.yr[i+1] - vect.age.yr[i])/2) + vect.age.yr[i]
  }
  
  for (i in 1:length(df$age_bin)){
    df$age_bin[i] <- min(vect.age.plot[df$age_bin[i]<= vect.age.plot], na.rm = T)
  }
  
  
  df$age_bin[df$age_bin == Inf] =  (ceiling(max(df$age)) - vect.age.yr[length(vect.age.yr)])/2  + vect.age.yr[length(vect.age.yr)]
  
  return(df)
  
}
get.seroprev.dat = function(data, vis_split, cutoff){
  
  visbin = seq(3,floor(max(data$age)), vis_split)
  #but breakdown the early years into more
  visbin = c(c(0,.5, 1, 2),   visbin)
  data$age_year <- NA
  
  for (i in 1:length(data$age_year)){
    tmp = visbin[data$age[i] > visbin]
    data$age_year[i] <- tmp[length(tmp)] 
  }
  
  if(cutoff=="mean"){
    dat.sum2 <- ddply(data, .(age_year), summarise,  prevalence=sum(prev)/length(prev), count=length(prev))  
  }else if (cutoff=="uci"){
    dat.sum2 <- ddply(data, .(age_year), summarise,  prevalence=sum(prev_uci)/length(prev_uci), count=length(prev_uci)) 
  }else if (cutoff=="lci"){
    dat.sum2 <- ddply(data, .(age_year), summarise,  prevalence=sum(prev_lci)/length(prev_lci), count=length(prev_lci)) 
  }
  
  #you want the midpoint between either end of the age class now
  vect.age.yr = sort(unique(data$age_year))
  dat.sum2$age_plot = NA
  
  for (i in 1:(length(dat.sum2$age_year)-1)){
    dat.sum2$age_plot [i] = ((dat.sum2$age_year[i + 1] - dat.sum2$age_year[i])/2) + dat.sum2$age_year[i]
  }
  dat.sum2$age_plot[length(dat.sum2$age_plot)] =  (ceiling(max(data$age)) - dat.sum2$age_year[length(dat.sum2$age_year)])/2  + dat.sum2$age_year[length(dat.sum2$age_year)]
  
  names(dat.sum2)  <- c("real_age", "prevalence", "count", "age_year") #<- names(dat.sum2.tmp)
  
  # dat.sum2 <-  rbind(dat.sum2.tmp, dat.sum2)
  
  dat.sum3 = dat.sum2
  dat.sum3$class = "seropositive"
  rownames(dat.sum3) <- c()
  
  return(dat.sum3)
  
}
plot.models.data.CI.squiggle <- function(CI.age.fits, CI.ts.fits, dat1, dat.sero1, dat.prev1, species1, antigen1, cutoff1, vis_split,  do.plot, do.save, filename1, filename2, filename3){
  
  mod.all = subset(CI.age.fits, antigen==antigen1 & cutoff==cutoff1 & species==species1)
  CI.ts.fits1 = subset(CI.ts.fits, antigen==antigen1 & cutoff==cutoff1 & species==species1)
  
  mod.all$model <- factor(mod.all$model, levels=c("MSIR", "MSIR-mu-sick", "MSIR-mu-rec", "MSIR-mu-sick-mu-rec",
                                                  "MSIRS", "MSIRS-mu-sick", "MSIRS-mu-rec", "MSIRS-mu-sick-mu-rec",
                                                  "MSIRN", "MSIRN-mu-sick", "MSIRN-mu-rec", "MSIRN-mu-sick-mu-rec",
                                                  "MSIRNI", "MSIRNI-mu-sick", "MSIRNI-mu-rec", "MSIRNI-mu-sick-mu-rec",
                                                  "MSIRNS", "MSIRNS-mu-sick", "MSIRNS-mu-rec", "MSIRNS-mu-sick-mu-rec",
                                                  "MSIRNNI", "MSIRNNI-mu-sick", "MSIRNNI-mu-rec", "MSIRNNI-mu-sick-mu-rec",
                                                  "MSIRNNS", "MSIRNNS-mu-sick", "MSIRNNS-mu-rec", "MSIRNNS-mu-sick-mu-rec",
                                                  "MSIRNIR", "MSIRNIR-mu-sick", "MSIRNIR-mu-rec", "MSIRNIR-mu-sick-mu-rec",
                                                  "MSIRNIRN", "MSIRNIRN-mu-sick", "MSIRNIRN-mu-sick-mu-rec","MSIRNIRN-mu-rec", 
                                                  "MSIRNIRN-LT", "MSIRNIRN-LT-mu-sick", "MSIRNIRN-LT-mu-sick-mu-rec","MSIRNIRN-LT-mu-rec", 
                                                  "MSIRNSIRN", "MSIRNSIRN-mu-sick", "MSIRNSIRN-mu-rec", "MSIRNSIRN-mu-sick-mu-rec",
                                                  "MSIRNSIRN-LT", "MSIRNSIRN-LT-mu-sick", "MSIRNSIRN-LT-mu-rec", "MSIRNSIRN-LT-mu-sick-mu-rec",
                                                  "MSIRSIR", "MSIRSIR-mu-sick", "MSIRSIR-mu-rec", "MSIRSIR-mu-sick-mu-rec",
                                                  "MSIRNSIRNS", "MSIRNSIRNS-mu-sick", "MSIRNSIRNS-mu-rec", "MSIRNSIRNS-mu-sick-mu-rec"))
  
  CI.ts.fits1$model <- factor(CI.ts.fits1$model, levels=c("MSIR", "MSIR-mu-sick", "MSIR-mu-rec", "MSIR-mu-sick-mu-rec",
                                                  "MSIRS", "MSIRS-mu-sick", "MSIRS-mu-rec", "MSIRS-mu-sick-mu-rec",
                                                  "MSIRN", "MSIRN-mu-sick", "MSIRN-mu-rec", "MSIRN-mu-sick-mu-rec",
                                                  "MSIRNI", "MSIRNI-mu-sick", "MSIRNI-mu-rec", "MSIRNI-mu-sick-mu-rec",
                                                  "MSIRNS", "MSIRNS-mu-sick", "MSIRNS-mu-rec", "MSIRNS-mu-sick-mu-rec",
                                                  "MSIRNNI", "MSIRNNI-mu-sick", "MSIRNNI-mu-rec", "MSIRNNI-mu-sick-mu-rec",
                                                  "MSIRNNS", "MSIRNNS-mu-sick", "MSIRNNS-mu-rec", "MSIRNNS-mu-sick-mu-rec",
                                                  "MSIRNIR", "MSIRNIR-mu-sick", "MSIRNIR-mu-rec", "MSIRNIR-mu-sick-mu-rec",
                                                  "MSIRNIRN", "MSIRNIRN-mu-sick", "MSIRNIRN-mu-sick-mu-rec","MSIRNIRN-mu-rec", 
                                                  "MSIRNIRN-LT", "MSIRNIRN-LT-mu-sick", "MSIRNIRN-LT-mu-sick-mu-rec","MSIRNIRN-LT-mu-rec", 
                                                  "MSIRNSIRN", "MSIRNSIRN-mu-sick", "MSIRNSIRN-mu-rec", "MSIRNSIRN-mu-sick-mu-rec",
                                                  "MSIRNSIRN-LT", "MSIRNSIRN-LT-mu-sick", "MSIRNSIRN-LT-mu-rec", "MSIRNSIRN-LT-mu-sick-mu-rec",
                                                  "MSIRSIR", "MSIRSIR-mu-sick", "MSIRSIR-mu-rec", "MSIRSIR-mu-sick-mu-rec",
                                                  "MSIRNSIRNS", "MSIRNSIRNS-mu-sick", "MSIRNSIRNS-mu-rec", "MSIRNSIRNS-mu-sick-mu-rec"))
  
  
  #first pull out your data...
  dat.tmp <- subset(dat1, species==species1 & antigen==antigen1)
  dat.tmp <- subset(dat.tmp, !is.na(age))
  
  #get data in age-seroprev form
  dat.sum = get.seroprev.dat(data=dat.tmp[!is.na(dat.tmp$age),], vis_split = vis_split, cutoff = cutoff1)
  dat.sum$seropos = dat.sum$prevalence*dat.sum$count
  
  dat.sum$seroprev_lci <-dat.sum$seroprev_uci <-NA
  for(i in 1:length(dat.sum$real_age)){
    dat.sum$seroprev_lci[i] <- binom.test(x=dat.sum$seropos[i], n=dat.sum$count[i])$conf.int[1]  
    dat.sum$seroprev_uci[i] <-  binom.test(x=dat.sum$seropos[i], n=dat.sum$count[i])$conf.int[2]
  }
  
  dat.sero = subset(dat.sero1, bat_species==species1 & antigen==antigen1)
  dat.sero$year<- year(dat.sero$processing_date)
  dat.sero.sum <- ddply(dat.sero,.(year, biwk), summarize, Nseropos = sum(seropos), Ntot=length(seropos))
  dat.sero.sum$seroprevalence <- dat.sero.sum$Nseropos/dat.sero.sum$Ntot
  dat.sero.sum$seroprev_lci <-dat.sero.sum$seroprev_uci <-NA
  for(i in 1:length(dat.sero.sum$year)){
    dat.sero.sum$seroprev_lci[i] <- binom.test(x=dat.sero.sum$Nseropos[i], n=dat.sero.sum$Ntot[i])$conf.int[1]  
    dat.sero.sum$seroprev_uci[i] <-  binom.test(x=dat.sero.sum$Nseropos[i], n=dat.sero.sum$Ntot[i])$conf.int[2]
  }
  
  dat.prev = subset(dat.prev1, bat_species==species1 & antigen==antigen1)
  dat.prev$year<- year(dat.prev$processing_date)
  dat.prev.sum <- ddply(dat.prev,.(year, biwk), summarize, Npos = sum(shed_pos), Ntot=length(shed_pos))
  dat.prev.sum$prevalence <- dat.prev.sum$Npos/dat.prev.sum$Ntot
  dat.prev.sum$prev_lci <-dat.prev.sum$prev_uci <-NA
  for(i in 1:length(dat.prev.sum$year)){
    dat.prev.sum$prev_lci[i] <- binom.test(x=dat.prev.sum$Npos[i], n=dat.prev.sum$Ntot[i])$conf.int[1]  
    dat.prev.sum$prev_uci[i] <-  binom.test(x=dat.prev.sum$Npos[i], n=dat.prev.sum$Ntot[i])$conf.int[2]
  }
  
  #plot with data
  max.a = ceiling(max(dat.sum$age_year))
  p1 <- ggplot() +  theme_bw() +
    geom_line(data=dat.sum, aes(x=age_year, y=prevalence), color ="black", linetype=2) +
    geom_linerange(data=dat.sum, aes(x=age_year, ymin = seroprev_lci, ymax=seroprev_uci), linetype=3, alpha=.5)+ 
    geom_point(data=dat.sum, aes(x=age_year, y=prevalence, size=count), fill="white", color ="black", pch=21) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + 
    scale_size_continuous(name="N bats") +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), 
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=mod.all, aes(x=age, y=seroprevalence), color="cornflowerblue") + 
    geom_ribbon(data=mod.all, aes(x=age, ymin=seroprev_lci, ymax = seroprev_uci), fill="cornflowerblue", alpha=.1) + 
    geom_ribbon(data=mod.all, aes(x=age, ymin=seroprev_lci_25, ymax = seroprev_uci_75), fill="cornflowerblue", alpha=.3) + 
    facet_wrap(model~., ncol=4) 
  
  
  
  print(p1)
  
  
  
  #and plot the others
  
  #sero 
  p2 <- ggplot() +  theme_bw() +
    geom_linerange(data=dat.sero.sum, aes(x=biwk, ymin=seroprev_lci, ymax=seroprev_uci, color = as.factor(year)), linetype=3, alpha=.5)+ 
    geom_point(data=dat.sero.sum, aes(x=biwk, y=seroprevalence, size=Ntot, shape = year, color = as.factor(year)), fill="white",  pch=21) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, 26)) + ylab("seroprevalence") + xlab("biweek of year")  + 
    scale_size_continuous(name="N bats") +
    scale_color_discrete(name="year") +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), 
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=CI.ts.fits1, aes(x=biwk, y=seroprevalence), color="cornflowerblue") + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=seroprev_lci, ymax = seroprev_uci), fill="cornflowerblue", alpha=.1) + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=seroprev_lci_25, ymax = seroprev_uci_75), fill="cornflowerblue", alpha=.3) + 
    facet_wrap(model~., ncol=4) 
  
  #prev
  p3 <- ggplot() +  theme_bw() +
    geom_linerange(data=dat.prev.sum, aes(x=biwk, ymin=prev_lci, ymax=prev_uci, color = as.factor(year)), linetype=3, alpha=.5)+ 
    geom_point(data=dat.prev.sum, aes(x=biwk, y=prevalence, size=Ntot, shape = year, color = as.factor(year)), fill="white",  pch=21) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, 26)) + ylab("prevalence") + xlab("biweek of year")  + 
    scale_size_continuous(name="N bats") +
    scale_color_discrete(name="year") +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), 
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=CI.ts.fits1, aes(x=biwk, y=prevalence), color="red") + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=prev_lci, ymax = seroprev_uci), fill="red", alpha=.1) + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=prev_lci_25, ymax = seroprev_uci_75), fill="red", alpha=.3) + 
    facet_wrap(model~., ncol=4) 
  
  if(do.save==TRUE){
    ggsave(file = filename1,
           plot = p1,
           units="mm",  
           width=100, 
           height=130, 
           scale=3, 
           dpi=300)
    
    ggsave(file = filename2,
           plot = p2,
           units="mm",  
           width=100, 
           height=130, 
           scale=3, 
           dpi=300)
    
    ggsave(file = filename3,
           plot = p3,
           units="mm",  
           width=100, 
           height=130, 
           scale=3, 
           dpi=300)
  }
  
  
  
  
}
#this version does not include the squiggles
plot.models.data.CI <- function(CI.age.fits, CI.ts.fits, dat1, dat.sero1, dat.prev1, species1, antigen1, cutoff1, vis_split,  do.plot, do.save, filename1, filename2, filename3){
  
  mod.all = subset(CI.age.fits, antigen==antigen1 & cutoff==cutoff1 & species==species1)
  
  CI.ts.fits1 = subset(CI.ts.fits, antigen==antigen1 & cutoff==cutoff1 & species==species1)
  
  mod.all$model <- factor(mod.all$model, levels=c("MSIR", "MSIR-mu-sick", "MSIR-mu-rec", "MSIR-mu-sick-mu-rec",
                                                  "MSIRS", "MSIRS-mu-sick", "MSIRS-mu-rec", "MSIRS-mu-sick-mu-rec",
                                                  "MSIRN", "MSIRN-mu-sick", "MSIRN-mu-rec", "MSIRN-mu-sick-mu-rec",
                                                  "MSIRNI", "MSIRNI-mu-sick", "MSIRNI-mu-rec", "MSIRNI-mu-sick-mu-rec",
                                                  "MSIRNS", "MSIRNS-mu-sick", "MSIRNS-mu-rec", "MSIRNS-mu-sick-mu-rec",
                                                  "MSIRNNI", "MSIRNNI-mu-sick", "MSIRNNI-mu-rec", "MSIRNNI-mu-sick-mu-rec",
                                                  "MSIRNNS", "MSIRNNS-mu-sick", "MSIRNNS-mu-rec", "MSIRNNS-mu-sick-mu-rec",
                                                  "MSIRNIR", "MSIRNIR-mu-sick", "MSIRNIR-mu-rec", "MSIRNIR-mu-sick-mu-rec",
                                                  "MSIRNIRN", "MSIRNIRN-mu-sick", "MSIRNIRN-mu-sick-mu-rec","MSIRNIRN-mu-rec", 
                                                  "MSIRNIRN-LT", "MSIRNIRN-LT-mu-sick", "MSIRNIRN-LT-mu-sick-mu-rec","MSIRNIRN-LT-mu-rec", 
                                                  "MSIRNSIRN", "MSIRNSIRN-mu-sick", "MSIRNSIRN-mu-rec", "MSIRNSIRN-mu-sick-mu-rec",
                                                  "MSIRNSIRN-LT", "MSIRNSIRN-LT-mu-sick", "MSIRNSIRN-LT-mu-rec", "MSIRNSIRN-LT-mu-sick-mu-rec",
                                                  "MSIRSIR", "MSIRSIR-mu-sick", "MSIRSIR-mu-rec", "MSIRSIR-mu-sick-mu-rec",
                                                  "MSIRNSIRNS", "MSIRNSIRNS-mu-sick", "MSIRNSIRNS-mu-rec", "MSIRNSIRNS-mu-sick-mu-rec"))
  
  CI.ts.fits1$model <- factor(CI.ts.fits1$model, levels=c("MSIR", "MSIR-mu-sick", "MSIR-mu-rec", "MSIR-mu-sick-mu-rec",
                                                          "MSIRS", "MSIRS-mu-sick", "MSIRS-mu-rec", "MSIRS-mu-sick-mu-rec",
                                                          "MSIRN", "MSIRN-mu-sick", "MSIRN-mu-rec", "MSIRN-mu-sick-mu-rec",
                                                          "MSIRNI", "MSIRNI-mu-sick", "MSIRNI-mu-rec", "MSIRNI-mu-sick-mu-rec",
                                                          "MSIRNS", "MSIRNS-mu-sick", "MSIRNS-mu-rec", "MSIRNS-mu-sick-mu-rec",
                                                          "MSIRNNI", "MSIRNNI-mu-sick", "MSIRNNI-mu-rec", "MSIRNNI-mu-sick-mu-rec",
                                                          "MSIRNNS", "MSIRNNS-mu-sick", "MSIRNNS-mu-rec", "MSIRNNS-mu-sick-mu-rec",
                                                          "MSIRNIR", "MSIRNIR-mu-sick", "MSIRNIR-mu-rec", "MSIRNIR-mu-sick-mu-rec",
                                                          "MSIRNIRN", "MSIRNIRN-mu-sick", "MSIRNIRN-mu-sick-mu-rec","MSIRNIRN-mu-rec", 
                                                          "MSIRNIRN-LT", "MSIRNIRN-LT-mu-sick", "MSIRNIRN-LT-mu-sick-mu-rec","MSIRNIRN-LT-mu-rec", 
                                                          "MSIRNSIRN", "MSIRNSIRN-mu-sick", "MSIRNSIRN-mu-rec", "MSIRNSIRN-mu-sick-mu-rec",
                                                          "MSIRNSIRN-LT", "MSIRNSIRN-LT-mu-sick", "MSIRNSIRN-LT-mu-rec", "MSIRNSIRN-LT-mu-sick-mu-rec",
                                                          "MSIRSIR", "MSIRSIR-mu-sick", "MSIRSIR-mu-rec", "MSIRSIR-mu-sick-mu-rec",
                                                          "MSIRNSIRNS", "MSIRNSIRNS-mu-sick", "MSIRNSIRNS-mu-rec", "MSIRNSIRNS-mu-sick-mu-rec"))
  
  
  #now summarize across age bins
  # #split model and add ages
  mod.split <- dlply(mod.all, .(model))
  # 
  # 
  # 
  mod.all <- data.table::rbindlist(lapply(mod.split, add.age))
  # 
  mod.all$seropos_lci <- mod.all$seroprev_lci*mod.all$N
  mod.all$seropos_uci <- mod.all$seroprev_uci*mod.all$N
  mod.all$seropos_lci_25 <- mod.all$seroprev_lci_25*mod.all$N
  mod.all$seropos_uci_75 <- mod.all$seroprev_uci_75*mod.all$N
  # 
  #  
  model.sum = ddply(mod.all,.(model, age_bin), summarise, N_pos = sum(seropos, na.rm = F), N_pos_lci = sum(seropos_lci, na.rm = F), N_pos_uci = sum(seropos_uci, na.rm = F),  N_pos_lci_25 = sum(seropos_lci_25, na.rm = F), N_pos_uci_75 = sum(seropos_uci_75, na.rm = F), Ntot=sum(N, na.rm = F))
  model.sum$seroprevalence <- model.sum$N_pos/model.sum$Ntot
  model.sum$seroprev_lci <- model.sum$N_pos_lci/model.sum$Ntot
  model.sum$seroprev_uci <- model.sum$N_pos_uci/model.sum$Ntot
  model.sum$seroprev_lci_25 <- model.sum$N_pos_lci_25/model.sum$Ntot
  model.sum$seroprev_uci_75 <- model.sum$N_pos_uci_75/model.sum$Ntot
  # 
  #first pull out your data...
  dat.tmp <- subset(dat1, species==species1 & antigen==antigen1)
  dat.tmp <- subset(dat.tmp, !is.na(age))
  
  #get data in age-seroprev form
  dat.sum = get.seroprev.dat(data=dat.tmp[!is.na(dat.tmp$age),], vis_split = vis_split, cutoff = cutoff1)
  dat.sum$seropos = dat.sum$prevalence*dat.sum$count
  
  dat.sum$seroprev_lci <-dat.sum$seroprev_uci <-NA
  for(i in 1:length(dat.sum$real_age)){
    dat.sum$seroprev_lci[i] <- binom.test(x=dat.sum$seropos[i], n=dat.sum$count[i])$conf.int[1]  
    dat.sum$seroprev_uci[i] <-  binom.test(x=dat.sum$seropos[i], n=dat.sum$count[i])$conf.int[2]
  }
  
  dat.sero = subset(dat.sero1, bat_species==species1 & antigen==antigen1)
  dat.sero$year<- year(dat.sero$processing_date)
  dat.sero.sum <- ddply(dat.sero,.(year, biwk), summarise, Nseropos = sum(seropos), Ntot=length(seropos))
  dat.sero.sum$seroprevalence <- dat.sero.sum$Nseropos/dat.sero.sum$Ntot
  dat.sero.sum$seroprev_lci <-dat.sero.sum$seroprev_uci <-NA
  for(i in 1:length(dat.sero.sum$year)){
    dat.sero.sum$seroprev_lci[i] <- binom.test(x=dat.sero.sum$Nseropos[i], n=dat.sero.sum$Ntot[i])$conf.int[1]  
    dat.sero.sum$seroprev_uci[i] <-  binom.test(x=dat.sero.sum$Nseropos[i], n=dat.sero.sum$Ntot[i])$conf.int[2]
  }
  
  dat.prev = subset(dat.prev1, bat_species==species1 & antigen==antigen1)
  dat.prev$year<- year(dat.prev$processing_date)
  dat.prev.sum <- ddply(dat.prev,.(year, biwk), summarise, Npos = sum(shed_pos), Ntot=length(shed_pos))
  dat.prev.sum$prevalence <- dat.prev.sum$Npos/dat.prev.sum$Ntot
  dat.prev.sum$prev_lci <-dat.prev.sum$prev_uci <-NA
  for(i in 1:length(dat.prev.sum$year)){
    dat.prev.sum$prev_lci[i] <- binom.test(x=dat.prev.sum$Npos[i], n=dat.prev.sum$Ntot[i])$conf.int[1]  
    dat.prev.sum$prev_uci[i] <-  binom.test(x=dat.prev.sum$Npos[i], n=dat.prev.sum$Ntot[i])$conf.int[2]
  }
  
  #plot with data
  max.a = ceiling(max(dat.sum$age_year))
  p1 <- ggplot() +  theme_bw() +
    geom_line(data=dat.sum, aes(x=age_year, y=prevalence), color ="black", linetype=2) +
    geom_linerange(data=dat.sum, aes(x=age_year, ymin = seroprev_lci, ymax=seroprev_uci), linetype=3, alpha=.5)+ 
    geom_point(data=dat.sum, aes(x=age_year, y=prevalence, size=count), fill="white", color ="black", pch=21) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + 
    scale_size_continuous(name="N bats") +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), 
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=model.sum, aes(x=age_bin, y=seroprevalence), color="cornflowerblue") + 
    geom_ribbon(data=model.sum, aes(x=age_bin, ymin=seroprev_lci, ymax = seroprev_uci), fill="cornflowerblue", alpha=.1) + 
    geom_ribbon(data=model.sum, aes(x=age_bin, ymin=seroprev_lci_25, ymax = seroprev_uci_75), fill="cornflowerblue", alpha=.3) + 
    facet_wrap(model~., ncol=4) 
  
  
  
  print(p1)
  
  
  
  #and plot the others
  
  #sero 
  p2 <- ggplot() +  theme_bw() +
    geom_linerange(data=dat.sero.sum, aes(x=biwk, ymin=seroprev_lci, ymax=seroprev_uci, color = as.factor(year)), linetype=3, alpha=.5)+ 
    geom_point(data=dat.sero.sum, aes(x=biwk, y=seroprevalence, size=Ntot, color = as.factor(year)), fill="white", pch=21) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, 26)) + ylab("seroprevalence") + xlab("biweek of year")  + 
    scale_size_continuous(name="N bats") +
    scale_color_discrete(name="year") +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), 
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=CI.ts.fits1, aes(x=biwk, y=seroprevalence), color="cornflowerblue") + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=seroprev_lci, ymax = seroprev_uci), fill="cornflowerblue", alpha=.1) + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=seroprev_lci_25, ymax = seroprev_uci_75), fill="cornflowerblue", alpha=.3) + 
    facet_wrap(model~., ncol=4) 
  
  #prev
  p3 <- ggplot() +  theme_bw() +
    geom_linerange(data=dat.prev.sum, aes(x=biwk, ymin=prev_lci, ymax=prev_uci, color = as.factor(year)), linetype=3, alpha=.5)+ 
    geom_point(data=dat.prev.sum, aes(x=biwk, y=prevalence, size=Ntot, color = as.factor(year)), fill="white", pch=21) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, 26)) + ylab("prevalence") + xlab("biweek of year")  + 
    scale_size_continuous(name="N bats") +
    scale_color_discrete(name="year") +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), 
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=CI.ts.fits1, aes(x=biwk, y=prevalence), color="red") + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=prev_lci, ymax = seroprev_uci), fill="red", alpha=.1) + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=prev_lci_25, ymax = seroprev_uci_75), fill="red", alpha=.3) + 
    facet_wrap(model~., ncol=4) 
  
  if(do.save==TRUE){
    ggsave(file = filename1,
           plot = p1,
           units="mm",  
           width=100, 
           height=130, 
           scale=3, 
           dpi=300)
    
    ggsave(file = filename2,
           plot = p2,
           units="mm",  
           width=100, 
           height=130, 
           scale=3, 
           dpi=300)
    
    ggsave(file = filename3,
           plot = p3,
           units="mm",  
           width=100, 
           height=130, 
           scale=3, 
           dpi=300)
  }
  
  
  
  
}
plot.models.data.CI.squiggle.one <- function(CI.age.fits, CI.ts.fits, dat1, dat.sero1, dat.prev1, species1, antigen1, cutoff1, vis_split,  do.plot, do.save, filename1, filename2, filename3){
  
  mod.all = subset(CI.age.fits, antigen==antigen1 & cutoff==cutoff1 & species==species1)
  CI.ts.fits1 = subset(CI.ts.fits, antigen==antigen1 & cutoff==cutoff1 & species==species1)
  
  mod.all$model <- factor(mod.all$model, levels=c("MSIR", "MSIR-mu-sick", "MSIR-mu-rec", "MSIR-mu-sick-mu-rec",
                                                  "MSIRS", "MSIRS-mu-sick", "MSIRS-mu-rec", "MSIRS-mu-sick-mu-rec",
                                                  "MSIRN", "MSIRN-mu-sick", "MSIRN-mu-rec", "MSIRN-mu-sick-mu-rec",
                                                  "MSIRNI", "MSIRNI-mu-sick", "MSIRNI-mu-rec", "MSIRNI-mu-sick-mu-rec",
                                                  "MSIRNS", "MSIRNS-mu-sick", "MSIRNS-mu-rec", "MSIRNS-mu-sick-mu-rec",
                                                  "MSIRNNI", "MSIRNNI-mu-sick", "MSIRNNI-mu-rec", "MSIRNNI-mu-sick-mu-rec",
                                                  "MSIRNNS", "MSIRNNS-mu-sick", "MSIRNNS-mu-rec", "MSIRNNS-mu-sick-mu-rec",
                                                  "MSIRNIR", "MSIRNIR-mu-sick", "MSIRNIR-mu-rec", "MSIRNIR-mu-sick-mu-rec",
                                                  "MSIRNIRN", "MSIRNIRN-mu-sick", "MSIRNIRN-mu-sick-mu-rec","MSIRNIRN-mu-rec", 
                                                  "MSIRNIRN-LT", "MSIRNIRN-LT-mu-sick", "MSIRNIRN-LT-mu-sick-mu-rec","MSIRNIRN-LT-mu-rec", 
                                                  "MSIRNSIRN", "MSIRNSIRN-mu-sick", "MSIRNSIRN-mu-rec", "MSIRNSIRN-mu-sick-mu-rec",
                                                  "MSIRNSIRN-LT", "MSIRNSIRN-LT-mu-sick", "MSIRNSIRN-LT-mu-rec", "MSIRNSIRN-LT-mu-sick-mu-rec",
                                                  "MSIRSIR", "MSIRSIR-mu-sick", "MSIRSIR-mu-rec", "MSIRSIR-mu-sick-mu-rec",
                                                  "MSIRNSIRNS", "MSIRNSIRNS-mu-sick", "MSIRNSIRNS-mu-rec", "MSIRNSIRNS-mu-sick-mu-rec"))
  
  CI.ts.fits1$model <- factor(CI.ts.fits1$model, levels=c("MSIR", "MSIR-mu-sick", "MSIR-mu-rec", "MSIR-mu-sick-mu-rec",
                                                          "MSIRS", "MSIRS-mu-sick", "MSIRS-mu-rec", "MSIRS-mu-sick-mu-rec",
                                                          "MSIRN", "MSIRN-mu-sick", "MSIRN-mu-rec", "MSIRN-mu-sick-mu-rec",
                                                          "MSIRNI", "MSIRNI-mu-sick", "MSIRNI-mu-rec", "MSIRNI-mu-sick-mu-rec",
                                                          "MSIRNS", "MSIRNS-mu-sick", "MSIRNS-mu-rec", "MSIRNS-mu-sick-mu-rec",
                                                          "MSIRNNI", "MSIRNNI-mu-sick", "MSIRNNI-mu-rec", "MSIRNNI-mu-sick-mu-rec",
                                                          "MSIRNNS", "MSIRNNS-mu-sick", "MSIRNNS-mu-rec", "MSIRNNS-mu-sick-mu-rec",
                                                          "MSIRNIR", "MSIRNIR-mu-sick", "MSIRNIR-mu-rec", "MSIRNIR-mu-sick-mu-rec",
                                                          "MSIRNIRN", "MSIRNIRN-mu-sick", "MSIRNIRN-mu-sick-mu-rec","MSIRNIRN-mu-rec", 
                                                          "MSIRNIRN-LT", "MSIRNIRN-LT-mu-sick", "MSIRNIRN-LT-mu-sick-mu-rec","MSIRNIRN-LT-mu-rec", 
                                                          "MSIRNSIRN", "MSIRNSIRN-mu-sick", "MSIRNSIRN-mu-rec", "MSIRNSIRN-mu-sick-mu-rec",
                                                          "MSIRNSIRN-LT", "MSIRNSIRN-LT-mu-sick", "MSIRNSIRN-LT-mu-rec", "MSIRNSIRN-LT-mu-sick-mu-rec",
                                                          "MSIRSIR", "MSIRSIR-mu-sick", "MSIRSIR-mu-rec", "MSIRSIR-mu-sick-mu-rec",
                                                          "MSIRNSIRNS", "MSIRNSIRNS-mu-sick", "MSIRNSIRNS-mu-rec", "MSIRNSIRNS-mu-sick-mu-rec"))
  
  
  #first pull out your data...
  dat.tmp <- subset(dat1, species==species1 & antigen==antigen1)
  dat.tmp <- subset(dat.tmp, !is.na(age))
  
  #get data in age-seroprev form
  dat.sum = get.seroprev.dat(data=dat.tmp[!is.na(dat.tmp$age),], vis_split = vis_split, cutoff = cutoff1)
  dat.sum$seropos = dat.sum$prevalence*dat.sum$count
  
  dat.sum$seroprev_lci <-dat.sum$seroprev_uci <-NA
  for(i in 1:length(dat.sum$real_age)){
    dat.sum$seroprev_lci[i] <- binom.test(x=dat.sum$seropos[i], n=dat.sum$count[i])$conf.int[1]  
    dat.sum$seroprev_uci[i] <-  binom.test(x=dat.sum$seropos[i], n=dat.sum$count[i])$conf.int[2]
  }
  
  dat.sero = subset(dat.sero1, bat_species==species1 & antigen==antigen1)
  dat.sero$year<- year(dat.sero$processing_date)
  dat.sero.sum <- ddply(dat.sero,.(year, biwk), summarize, Nseropos = sum(seropos), Ntot=length(seropos))
  dat.sero.sum$seroprevalence <- dat.sero.sum$Nseropos/dat.sero.sum$Ntot
  dat.sero.sum$seroprev_lci <-dat.sero.sum$seroprev_uci <-NA
  for(i in 1:length(dat.sero.sum$year)){
    dat.sero.sum$seroprev_lci[i] <- binom.test(x=dat.sero.sum$Nseropos[i], n=dat.sero.sum$Ntot[i])$conf.int[1]  
    dat.sero.sum$seroprev_uci[i] <-  binom.test(x=dat.sero.sum$Nseropos[i], n=dat.sero.sum$Ntot[i])$conf.int[2]
  }
  
  dat.prev = subset(dat.prev1, bat_species==species1 & antigen==antigen1)
  dat.prev$year<- year(dat.prev$processing_date)
  dat.prev.sum <- ddply(dat.prev,.(year, biwk), summarize, Npos = sum(shed_pos), Ntot=length(shed_pos))
  dat.prev.sum$prevalence <- dat.prev.sum$Npos/dat.prev.sum$Ntot
  dat.prev.sum$prev_lci <-dat.prev.sum$prev_uci <-NA
  for(i in 1:length(dat.prev.sum$year)){
    dat.prev.sum$prev_lci[i] <- binom.test(x=dat.prev.sum$Npos[i], n=dat.prev.sum$Ntot[i])$conf.int[1]  
    dat.prev.sum$prev_uci[i] <-  binom.test(x=dat.prev.sum$Npos[i], n=dat.prev.sum$Ntot[i])$conf.int[2]
  }
  
  #plot with data
  max.a = ceiling(max(dat.sum$age_year))
  p1 <- ggplot() +  theme_bw() +
    geom_line(data=dat.sum, aes(x=age_year, y=prevalence), color ="black", linetype=2) +
    geom_linerange(data=dat.sum, aes(x=age_year, ymin = seroprev_lci, ymax=seroprev_uci), linetype=3, alpha=.5)+ 
    geom_point(data=dat.sum, aes(x=age_year, y=prevalence, size=count), fill="white", color ="black", pch=21) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + 
    scale_size_continuous(name="N bats") +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), 
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=mod.all, aes(x=age, y=seroprevalence), color="cornflowerblue") + 
    geom_ribbon(data=mod.all, aes(x=age, ymin=seroprev_lci, ymax = seroprev_uci), fill="cornflowerblue", alpha=.1) + 
    geom_ribbon(data=mod.all, aes(x=age, ymin=seroprev_lci_25, ymax = seroprev_uci_75), fill="cornflowerblue", alpha=.3) + 
    facet_wrap(model~., ncol=4) 
  
  
  
  print(p1)
  
  
  
  #and plot the others
  
  #sero 
  p2 <- ggplot() +  theme_bw() +
    geom_linerange(data=dat.sero.sum, aes(x=biwk, ymin=seroprev_lci, ymax=seroprev_uci, color = as.factor(year)), linetype=3, alpha=.5)+ 
    geom_point(data=dat.sero.sum, aes(x=biwk, y=seroprevalence, size=Ntot, shape = year, color = as.factor(year)), fill="white",  pch=21) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, 26)) + ylab("seroprevalence") + xlab("biweek of year")  + 
    scale_size_continuous(name="N bats") +
    scale_color_discrete(name="year") +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), 
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=CI.ts.fits1, aes(x=biwk, y=seroprevalence), color="cornflowerblue") + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=seroprev_lci, ymax = seroprev_uci), fill="cornflowerblue", alpha=.1) + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=seroprev_lci_25, ymax = seroprev_uci_75), fill="cornflowerblue", alpha=.3) + 
    facet_wrap(model~., ncol=4) 
  
  #prev
  p3 <- ggplot() +  theme_bw() +
    geom_linerange(data=dat.prev.sum, aes(x=biwk, ymin=prev_lci, ymax=prev_uci, color = as.factor(year)), linetype=3, alpha=.5)+ 
    geom_point(data=dat.prev.sum, aes(x=biwk, y=prevalence, size=Ntot, shape = year, color = as.factor(year)), fill="white",  pch=21) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, 26)) + ylab("prevalence") + xlab("biweek of year")  + 
    scale_size_continuous(name="N bats") +
    scale_color_discrete(name="year") +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), 
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=CI.ts.fits1, aes(x=biwk, y=prevalence), color="red") + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=prev_lci, ymax = seroprev_uci), fill="red", alpha=.1) + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=prev_lci_25, ymax = seroprev_uci_75), fill="red", alpha=.3) + 
    facet_wrap(model~., ncol=4) 
  
  if(do.save==TRUE){
    ggsave(file = filename1,
           plot = p1,
           units="mm",  
           width=55, 
           height=45,
           scale=3, 
           dpi=300)
    
    ggsave(file = filename2,
           plot = p2,
           units="mm",  
           width=55, 
           height=45, 
           scale=3, 
           dpi=300)
    
    ggsave(file = filename3,
           plot = p3,
           units="mm",  
           width=55, 
           height=45, 
           scale=3, 
           dpi=300)
  }
  
  
  
  
}
#this version does not include the squiggles
plot.models.data.CI.one <- function(CI.age.fits, CI.ts.fits, dat1, dat.sero1, dat.prev1, species1, antigen1, cutoff1, vis_split,  do.plot, do.save, filename1, filename2, filename3){
  
  mod.all = subset(CI.age.fits, antigen==antigen1 & cutoff==cutoff1 & species==species1)
  
  CI.ts.fits1 = subset(CI.ts.fits, antigen==antigen1 & cutoff==cutoff1 & species==species1)
  
  mod.all$model <- factor(mod.all$model, levels=c("MSIR", "MSIR-mu-sick", "MSIR-mu-rec", "MSIR-mu-sick-mu-rec",
                                                  "MSIRS", "MSIRS-mu-sick", "MSIRS-mu-rec", "MSIRS-mu-sick-mu-rec",
                                                  "MSIRN", "MSIRN-mu-sick", "MSIRN-mu-rec", "MSIRN-mu-sick-mu-rec",
                                                  "MSIRNI", "MSIRNI-mu-sick", "MSIRNI-mu-rec", "MSIRNI-mu-sick-mu-rec",
                                                  "MSIRNS", "MSIRNS-mu-sick", "MSIRNS-mu-rec", "MSIRNS-mu-sick-mu-rec",
                                                  "MSIRNNI", "MSIRNNI-mu-sick", "MSIRNNI-mu-rec", "MSIRNNI-mu-sick-mu-rec",
                                                  "MSIRNNS", "MSIRNNS-mu-sick", "MSIRNNS-mu-rec", "MSIRNNS-mu-sick-mu-rec",
                                                  "MSIRNIR", "MSIRNIR-mu-sick", "MSIRNIR-mu-rec", "MSIRNIR-mu-sick-mu-rec",
                                                  "MSIRNIRN", "MSIRNIRN-mu-sick", "MSIRNIRN-mu-sick-mu-rec","MSIRNIRN-mu-rec", 
                                                  "MSIRNIRN-LT", "MSIRNIRN-LT-mu-sick", "MSIRNIRN-LT-mu-sick-mu-rec","MSIRNIRN-LT-mu-rec", 
                                                  "MSIRNSIRN", "MSIRNSIRN-mu-sick", "MSIRNSIRN-mu-rec", "MSIRNSIRN-mu-sick-mu-rec",
                                                  "MSIRNSIRN-LT", "MSIRNSIRN-LT-mu-sick", "MSIRNSIRN-LT-mu-rec", "MSIRNSIRN-LT-mu-sick-mu-rec",
                                                  "MSIRSIR", "MSIRSIR-mu-sick", "MSIRSIR-mu-rec", "MSIRSIR-mu-sick-mu-rec",
                                                  "MSIRNSIRNS", "MSIRNSIRNS-mu-sick", "MSIRNSIRNS-mu-rec", "MSIRNSIRNS-mu-sick-mu-rec"))
  
  CI.ts.fits1$model <- factor(CI.ts.fits1$model, levels=c("MSIR", "MSIR-mu-sick", "MSIR-mu-rec", "MSIR-mu-sick-mu-rec",
                                                          "MSIRS", "MSIRS-mu-sick", "MSIRS-mu-rec", "MSIRS-mu-sick-mu-rec",
                                                          "MSIRN", "MSIRN-mu-sick", "MSIRN-mu-rec", "MSIRN-mu-sick-mu-rec",
                                                          "MSIRNI", "MSIRNI-mu-sick", "MSIRNI-mu-rec", "MSIRNI-mu-sick-mu-rec",
                                                          "MSIRNS", "MSIRNS-mu-sick", "MSIRNS-mu-rec", "MSIRNS-mu-sick-mu-rec",
                                                          "MSIRNNI", "MSIRNNI-mu-sick", "MSIRNNI-mu-rec", "MSIRNNI-mu-sick-mu-rec",
                                                          "MSIRNNS", "MSIRNNS-mu-sick", "MSIRNNS-mu-rec", "MSIRNNS-mu-sick-mu-rec",
                                                          "MSIRNIR", "MSIRNIR-mu-sick", "MSIRNIR-mu-rec", "MSIRNIR-mu-sick-mu-rec",
                                                          "MSIRNIRN", "MSIRNIRN-mu-sick", "MSIRNIRN-mu-sick-mu-rec","MSIRNIRN-mu-rec", 
                                                          "MSIRNIRN-LT", "MSIRNIRN-LT-mu-sick", "MSIRNIRN-LT-mu-sick-mu-rec","MSIRNIRN-LT-mu-rec", 
                                                          "MSIRNSIRN", "MSIRNSIRN-mu-sick", "MSIRNSIRN-mu-rec", "MSIRNSIRN-mu-sick-mu-rec",
                                                          "MSIRNSIRN-LT", "MSIRNSIRN-LT-mu-sick", "MSIRNSIRN-LT-mu-rec", "MSIRNSIRN-LT-mu-sick-mu-rec",
                                                          "MSIRSIR", "MSIRSIR-mu-sick", "MSIRSIR-mu-rec", "MSIRSIR-mu-sick-mu-rec",
                                                          "MSIRNSIRNS", "MSIRNSIRNS-mu-sick", "MSIRNSIRNS-mu-rec", "MSIRNSIRNS-mu-sick-mu-rec"))
  
  
  #now summarize across age bins
  # #split model and add ages
  mod.split <- dlply(mod.all, .(model))
  # 
  # 
  # 
  mod.all <- data.table::rbindlist(lapply(mod.split, add.age))
  # 
  mod.all$seropos_lci <- mod.all$seroprev_lci*mod.all$N
  mod.all$seropos_uci <- mod.all$seroprev_uci*mod.all$N
  mod.all$seropos_lci_25 <- mod.all$seroprev_lci_25*mod.all$N
  mod.all$seropos_uci_75 <- mod.all$seroprev_uci_75*mod.all$N
  # 
  #  
  model.sum = ddply(mod.all,.(model, age_bin), summarise, N_pos = sum(seropos, na.rm = F), N_pos_lci = sum(seropos_lci, na.rm = F), N_pos_uci = sum(seropos_uci, na.rm = F),  N_pos_lci_25 = sum(seropos_lci_25, na.rm = F), N_pos_uci_75 = sum(seropos_uci_75, na.rm = F), Ntot=sum(N, na.rm = F))
  model.sum$seroprevalence <- model.sum$N_pos/model.sum$Ntot
  model.sum$seroprev_lci <- model.sum$N_pos_lci/model.sum$Ntot
  model.sum$seroprev_uci <- model.sum$N_pos_uci/model.sum$Ntot
  model.sum$seroprev_lci_25 <- model.sum$N_pos_lci_25/model.sum$Ntot
  model.sum$seroprev_uci_75 <- model.sum$N_pos_uci_75/model.sum$Ntot
  # 
  #first pull out your data...
  dat.tmp <- subset(dat1, species==species1 & antigen==antigen1)
  dat.tmp <- subset(dat.tmp, !is.na(age))
  
  #get data in age-seroprev form
  dat.sum = get.seroprev.dat(data=dat.tmp[!is.na(dat.tmp$age),], vis_split = vis_split, cutoff = cutoff1)
  dat.sum$seropos = dat.sum$prevalence*dat.sum$count
  
  dat.sum$seroprev_lci <-dat.sum$seroprev_uci <-NA
  for(i in 1:length(dat.sum$real_age)){
    dat.sum$seroprev_lci[i] <- binom.test(x=dat.sum$seropos[i], n=dat.sum$count[i])$conf.int[1]  
    dat.sum$seroprev_uci[i] <-  binom.test(x=dat.sum$seropos[i], n=dat.sum$count[i])$conf.int[2]
  }
  
  dat.sero = subset(dat.sero1, bat_species==species1 & antigen==antigen1)
  dat.sero$year<- year(dat.sero$processing_date)
  dat.sero.sum <- ddply(dat.sero,.(year, biwk), summarize, Nseropos = sum(seropos), Ntot=length(seropos))
  dat.sero.sum$seroprevalence <- dat.sero.sum$Nseropos/dat.sero.sum$Ntot
  dat.sero.sum$seroprev_lci <-dat.sero.sum$seroprev_uci <-NA
  for(i in 1:length(dat.sero.sum$year)){
    dat.sero.sum$seroprev_lci[i] <- binom.test(x=dat.sero.sum$Nseropos[i], n=dat.sero.sum$Ntot[i])$conf.int[1]  
    dat.sero.sum$seroprev_uci[i] <-  binom.test(x=dat.sero.sum$Nseropos[i], n=dat.sero.sum$Ntot[i])$conf.int[2]
  }
  
  dat.prev = subset(dat.prev1, bat_species==species1 & antigen==antigen1)
  dat.prev$year<- year(dat.prev$processing_date)
  dat.prev.sum <- ddply(dat.prev,.(year, biwk), summarize, Npos = sum(shed_pos), Ntot=length(shed_pos))
  dat.prev.sum$prevalence <- dat.prev.sum$Npos/dat.prev.sum$Ntot
  dat.prev.sum$prev_lci <-dat.prev.sum$prev_uci <-NA
  for(i in 1:length(dat.prev.sum$year)){
    dat.prev.sum$prev_lci[i] <- binom.test(x=dat.prev.sum$Npos[i], n=dat.prev.sum$Ntot[i])$conf.int[1]  
    dat.prev.sum$prev_uci[i] <-  binom.test(x=dat.prev.sum$Npos[i], n=dat.prev.sum$Ntot[i])$conf.int[2]
  }
  
  #plot with data
  max.a = ceiling(max(dat.sum$age_year))
  p1 <- ggplot() +  theme_bw() +
    geom_line(data=dat.sum, aes(x=age_year, y=prevalence), color ="black", linetype=2) +
    geom_linerange(data=dat.sum, aes(x=age_year, ymin = seroprev_lci, ymax=seroprev_uci), linetype=3, alpha=.5)+ 
    geom_point(data=dat.sum, aes(x=age_year, y=prevalence, size=count), fill="white", color ="black", pch=21) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + 
    scale_size_continuous(name="N bats") +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), 
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=model.sum, aes(x=age, y=seroprevalence), color="cornflowerblue") + 
    geom_ribbon(data=model.sum, aes(x=age, ymin=seroprev_lci, ymax = seroprev_uci), fill="cornflowerblue", alpha=.1) + 
    geom_ribbon(data=model.sum, aes(x=age, ymin=seroprev_lci_25, ymax = seroprev_uci_75), fill="cornflowerblue", alpha=.3) + 
    facet_wrap(model~., ncol=4) 
  
  
  
  print(p1)
  
  
  
  #and plot the others
  
  #sero 
  p2 <- ggplot() +  theme_bw() +
    geom_linerange(data=dat.sero.sum, aes(x=biwk, ymin=seroprev_lci, ymax=seroprev_uci, color = as.factor(year)), linetype=3, alpha=.5)+ 
    geom_point(data=dat.sero.sum, aes(x=biwk, y=seroprevalence, size=Ntot, color = as.factor(year)), fill="white", pch=21) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, 26)) + ylab("seroprevalence") + xlab("biweek of year")  + 
    scale_size_continuous(name="N bats") +
    scale_color_discrete(name="year") +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), 
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=CI.ts.fits1, aes(x=biwk, y=seroprevalence), color="cornflowerblue") + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=seroprev_lci, ymax = seroprev_uci), fill="cornflowerblue", alpha=.1) + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=seroprev_lci_25, ymax = seroprev_uci_75), fill="cornflowerblue", alpha=.3) + 
    facet_wrap(model~., ncol=4) 
  
  #prev
  p3 <- ggplot() +  theme_bw() +
    geom_linerange(data=dat.prev.sum, aes(x=biwk, ymin=prev_lci, ymax=prev_uci, color = as.factor(year)), linetype=3, alpha=.5)+ 
    geom_point(data=dat.prev.sum, aes(x=biwk, y=prevalence, size=Ntot, color = as.factor(year)), fill="white", pch=21) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, 26)) + ylab("prevalence") + xlab("biweek of year")  + 
    scale_size_continuous(name="N bats") +
    scale_color_discrete(name="year") +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), 
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=CI.ts.fits1, aes(x=biwk, y=prevalence), color="red") + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=prev_lci, ymax = seroprev_uci), fill="red", alpha=.1) + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=prev_lci_25, ymax = seroprev_uci_75), fill="red", alpha=.3) + 
    facet_wrap(model~., ncol=4) 
  
  if(do.save==TRUE){
    ggsave(file = filename1,
           plot = p1,
           units="mm",  
           width=55, 
           height=45, 
           scale=3, 
           dpi=300)
    
    ggsave(file = filename2,
           plot = p2,
           units="mm",  
           width=55, 
           height=45, 
           scale=3, 
           dpi=300)
    
    ggsave(file = filename3,
           plot = p3,
           units="mm",  
           width=55, 
           height=45, 
           scale=3, 
           dpi=300)
  }
  
  
  
  
}


plot.models.data.CI.squiggle(CI.age.fits = age.sero.df,
                             CI.ts.fits = biwk.ts.df,
                             dat1=dat, 
                             dat.sero1 = sero.dat,
                             dat.prev1 = shed.dat,
                             antigen1 = "GhV", 
                             cutoff1="mean", 
                             species1 = "Eidolon dupreanum",  
                             vis_split = 3, 
                             do.plot=TRUE, 
                             do.save=TRUE, 
                             filename1=paste0(homewd,"/final_figures_tables/tmp/GhV-fits-age-seroprev.png"),
                             filename2=paste0(homewd,"/final_figures_tables/tmp/GhV-fits-biweekly-serology.png"),
                             filename3=paste0(homewd,"/final_figures_tables/tmp/GhV-fits-biweekly-prevalence.png"))



#and if you want no squiggles
plot.models.data.CI(CI.age.fits = age.sero.df,
                             CI.ts.fits = biwk.ts.df,
                             dat1=dat, 
                             dat.sero1 = sero.dat,
                             dat.prev1 = shed.dat,
                             antigen1 = "GhV", 
                             cutoff1="mean", 
                             species1 = "Eidolon dupreanum",  
                             vis_split = 3, 
                             do.plot=TRUE, 
                             do.save=TRUE, 
                             filename1=paste0(homewd,"/final_figures_tables/tmp/GhV-fits-age-seroprev-no-squiggle.png"),
                             filename2=paste0(homewd,"/final_figures_tables/tmp/GhV-fits-biweekly-serology-no-squiggle.png"),
                             filename3=paste0(homewd,"/final_figures_tables/tmp/GhV-fits-biweekly-prevalence-no-squiggle.png"))



#and just the best fit

plot.models.data.CI.squiggle.one(CI.age.fits = subset(age.sero.df, model=="MSIRNS-mu-rec"),
                             CI.ts.fits = subset(biwk.ts.df, model=="MSIRNS-mu-rec"),
                             dat1=dat, 
                             dat.sero1 = sero.dat,
                             dat.prev1 = shed.dat,
                             antigen1 = "GhV", 
                             cutoff1="mean", 
                             species1 = "Eidolon dupreanum",  
                             vis_split = 3, 
                             do.plot=TRUE, 
                             do.save=TRUE, 
                             filename1=paste0(homewd,"/final_figures_tables/tmp/GhV-fits-age-seroprev-best-fit.png"),
                             filename2=paste0(homewd,"/final_figures_tables/tmp/GhV-fits-biweekly-serology-best-fit.png"),
                             filename3=paste0(homewd,"/final_figures_tables/tmp/GhV-fits-biweekly-prevalence-best-fit.png"))



#and if you want no squiggles
plot.models.data.CI.one(CI.age.fits = subset(age.sero.df, model=="MSIRNS-mu-rec"),
                    CI.ts.fits = subset(biwk.ts.df, model=="MSIRNS-mu-rec"),
                    dat1=dat, 
                    dat.sero1 = sero.dat,
                    dat.prev1 = shed.dat,
                    antigen1 = "GhV", 
                    cutoff1="mean", 
                    species1 = "Eidolon dupreanum",  
                    vis_split = 3, 
                    do.plot=TRUE, 
                    do.save=TRUE, 
                    filename1=paste0(homewd,"/final_figures_tables/tmp/GhV-fits-age-seroprev-no-squiggle-best-fit.png"),
                    filename2=paste0(homewd,"/final_figures_tables/tmp/GhV-fits-biweekly-serology-no-squiggle-best-fit.png"),
                    filename3=paste0(homewd,"/final_figures_tables/tmp/GhV-fits-biweekly-prevalence-no-squiggle-best-fit.png"))
