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
library(magick)
library(grid)
library(pdftools)
#load data
homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
setwd(homewd)

#load the completed models with age-seroprevalence and time series
age.sero.df <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/final_age_sero_df_GhV.csv"), header = T, stringsAsFactors = F)
biwk.ts.df <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/final_biweekly_ts_GhV.csv"), header = T, stringsAsFactors = F)

#and load the data
dat <- read.csv(file = paste0(homewd,"/master_data_files/age-seroprevalence-data.csv"), stringsAsFactors = F, header = T)
head(dat)
sero.dat <- read.csv(file = paste0(homewd,"/master_data_files/biweekly-serology-data.csv"), stringsAsFactors = F, header = T)
shed.dat <- read.csv(file = paste0(homewd,"/master_data_files/biweekly-shedding-data.csv"), stringsAsFactors = F, header = T)
head(dat)
head(sero.dat)
head(shed.dat)
#shed.dat = subset(shed.dat, processing_date <"2022-01-01")

#plotting scripts
add.age <- function(df){
  df <- arrange(df, age)
  bins = c(c(0,.5,1,2), c(3:20))
  #bins = c(c(0,.25,.5,.75,1,2), c(3:20))
  #bins = c(c(0,.6,.9,1,2), c(3:20))
  
  df$age_bin <- NA
  for (i in 1:length(df$age)){
    df$age_bin[i] <- max(bins[df$age[i]>bins], na.rm = T)
  }
  
  df$age_plot<- NA
  
  #you want the midpoint between either end of the age class now
  vect.age.yr = sort(unique(df$age_bin))
  vect.age.plot =   rep(NA, (length(vect.age.yr)-1))
  for (i in 1:length(vect.age.plot)){
    vect.age.plot[i] <- ((vect.age.yr[i+1] - vect.age.yr[i])/2) + vect.age.yr[i]
  }
  
  vect.age.plot <- c(  vect.age.plot, 20)
  
  #now add the plotting bin to the actual age bins
  
  for (i in 1:length(df$age_bin)){
    df$age_plot[i] <- min(vect.age.plot[df$age_bin[i]<= vect.age.plot], na.rm = T)
  }
  
  
  #df$age_bin[df$age_bin == Inf] =  (ceiling(max(df$age)) - vect.age.yr[length(vect.age.yr)])/2  + vect.age.yr[length(vect.age.yr)]
  
  return(df)
  
}
get.seroprev.dat = function(data, vis_split, cutoff){
  
  visbin = seq(3,floor(max(data$age)), vis_split)
  #but breakdown the early years into more
  visbin = c(c(0,.5, 1,2),   visbin)
  #visbin = c(c(0,.25,.5,.75, 1,2),   visbin)
  #visbin = c(c(0,.6,.9, 1,2),   visbin)
  data$age_year <- NA
  
  for (i in 1:length(data$age_year)){
    tmp = visbin[data$age[i] > visbin]
    data$age_year[i] <- tmp[length(tmp)] 
  }
  
  if(cutoff=="mean"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev)/length(prev), count=length(prev))  
  }else if (cutoff=="uci"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev_uci)/length(prev_uci), count=length(prev_uci)) 
  }else if (cutoff=="lci"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev_lci)/length(prev_lci), count=length(prev_lci)) 
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
return.models.data.CI.one <- function(CI.age.fits, CI.ts.fits, dat1, dat.sero1, dat.prev1, species1, antigen1, cutoff1, vis_split){
  
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
  model.sum = ddply(mod.all,.(model, age_plot), summarise, N_pos = sum(seropos, na.rm = F), N_pos_lci = sum(seropos_lci, na.rm = F), N_pos_uci = sum(seropos_uci, na.rm = F),  N_pos_lci_25 = sum(seropos_lci_25, na.rm = F), N_pos_uci_75 = sum(seropos_uci_75, na.rm = F), Ntot=sum(N, na.rm = F))
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
    geom_line(data=dat.sum, aes(x=age_year, y=prevalence), color ="black", linetype=2, show.legend = F) +
    geom_linerange(data=dat.sum, aes(x=age_year, ymin = seroprev_lci, ymax=seroprev_uci), linetype=3, alpha=.5, show.legend = F)+ 
    geom_point(data=dat.sum, aes(x=age_year, y=prevalence, size=count), fill="white", color ="black", pch=21, show.legend = F) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + 
    scale_size_continuous(name="N bats", breaks = c(0,10,20,30,40,50)) +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), strip.text = element_blank(),
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=model.sum, aes(x=age_plot, y=seroprevalence), color="cornflowerblue") + 
    geom_ribbon(data=model.sum, aes(x=age_plot, ymin=seroprev_lci, ymax = seroprev_uci), fill="cornflowerblue", alpha=.1) + 
    geom_ribbon(data=model.sum, aes(x=age_plot, ymin=seroprev_lci_25, ymax = seroprev_uci_75), fill="cornflowerblue", alpha=.3) + 
    facet_wrap(model~., ncol=4) ;p1
  
  
  
  yearz = c('2013'="#F8766D", '2014' ="#C49A00", '2015'= "#53B400",'2016'= "#00C094", '2018'= "#00B6EB", '2019'= "#A58AFF", '2020'= "#FB61D7")
  #and plot the others
  
  #sero 
  p2 <- ggplot() +  theme_bw() +
    geom_linerange(data=dat.sero.sum, aes(x=biwk, ymin=seroprev_lci, ymax=seroprev_uci, color = as.factor(year)), linetype=3, alpha=.5, show.legend = F)+ 
    geom_point(data=dat.sero.sum, aes(x=biwk, y=seroprevalence, size=Ntot, color = as.factor(year)), fill="white", pch=21, show.legend = F) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, 26)) + ylab("seroprevalence") + xlab("biweek of year")  + 
    scale_size_continuous(name="N bats", breaks = c(0,10,20,30,40,50)) +
    scale_color_manual(name="year", values =   yearz) +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), strip.text = element_blank(),
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=CI.ts.fits1, aes(x=biwk, y=seroprevalence), color="cornflowerblue") + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=seroprev_lci, ymax = seroprev_uci), fill="cornflowerblue", alpha=.1) + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=seroprev_lci_25, ymax = seroprev_uci_75), fill="cornflowerblue", alpha=.3) + 
    facet_wrap(model~., ncol=4) ;p2
  
  #prev
  p3 <- ggplot() +  theme_bw() +
    geom_linerange(data=dat.prev.sum, aes(x=biwk, ymin=prev_lci, ymax=prev_uci, color = as.factor(year)), linetype=3, alpha=.5)+ 
    geom_point(data=dat.prev.sum, aes(x=biwk, y=prevalence, size=Ntot, color = as.factor(year)), fill="white", pch=21) +
    theme(legend.position = "right",  legend.spacing = unit(.03,"cm"), 
          panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, 26)) + ylab("prevalence") + xlab("biweek of year")  + 
    scale_size_continuous(name="N bats", breaks = c(0,10,20,30,40,50)) +
    scale_color_manual(name="year", values =   yearz) +
    guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=18), legend.direction = "vertical", legend.box = "vertical",strip.text = element_blank(),
          axis.text = element_text(size=14), panel.grid = element_blank()) +
    geom_line(data=CI.ts.fits1, aes(x=biwk, y=prevalence), color="red") + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=prev_lci, ymax = seroprev_uci), fill="red", alpha=.1) + 
    geom_ribbon(data=CI.ts.fits1, aes(x=biwk, ymin=prev_lci_25, ymax = seroprev_uci_75), fill="red", alpha=.3) + 
    facet_wrap(model~., ncol=4) ; p3
  
  
  pbottom <- cowplot::plot_grid(p1, p2, p3, ncol=3, nrow=1, rel_widths = c(1,1,1.2), labels = c("D", "E", "F"), label_size = 22, align = "h");pbottom
 
  
  return(pbottom)
}

plot.matAb <- function(dat){
  
  names(dat)[names(dat)=="species"] <- "bat_species"
  dat$log10MFI <- log10(dat$MFI)
  
  #focus just on GhV for now
  dat.mat = subset(dat, antigen=="GhV" & age <=1 & bat_species=="Eidolon dupreanum")
  #sort(dat.mat$log10MFI)

  
  #now get your titer trajectory for maternal immunity
  #want to remove those that might have been born seronegative -
  #remove any with age <1 week that are seroneg
  dat.mat[dat.mat$age<= (1/52) & dat.mat$prev==0,] # just 1
  dat.mat = dat.mat[dat.mat$sampleid!="KEL402",]
  
  #some cutoffs were established by mixture modeling
  cutoff.df <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 3)), antigen=rep(c( "GhV"),3))
  cutoff.df$MFI = c(1141,1763, 1619.25)#, 2037.75, 905.25, 1667.25)
  cutoff.df$type = rep(c("mean", "uci", "lci"),each=1)
  cutoff.df$log10MFI<- log10(cutoff.df$MFI)
  
  
  #the viral kinetics model is embedded here in the fitted script
  # alpha = the initial concentration of maternal antibodies that can wane.
  # lambda = the decay rate (kinetics). It dictates how fast those maternal antibodies disappear.
  # eta = the baseline or background antibody level. In maternal immunity models, this usually represents either the infant's own early endogenous antibody production or the technical detection limit (background noise) of the assay used to measure the antibod
  
  fit <- nls(log10MFI ~ alpha * exp(-lambda * age) + eta,
             data = dat.mat,
             start = list(alpha = max(dat.mat$log10MFI, na.rm=TRUE),
                          lambda = 0.5,
                          eta = min(dat.mat$log10MFI, na.rm=TRUE)))
  
  
  
  # Compute deterministic mean trajectory (black line)
  x <- seq(0, max(dat.mat$log10MFI), length.out = 100)
  T <- 365   # number of time steps (days)
  mu <- coef(fit)[["eta"]] + exp(-coef(fit)[["lambda"]] * ((0:(T-1))/365)) * (max(x) - coef(fit)[["eta"]])
  
  df.mean <- cbind.data.frame(age=(0:(T-1))/365, log10MFI=mu)
  
  #and the upper and lower lci
  CI.dat <- confint(fit)
  mu_lci <- CI.dat['eta','2.5%'] + exp(-CI.dat['lambda','97.5%']* ((0:(T-1))/365)) * (max(x) - CI.dat['eta','2.5%'])
  mu_uci <- CI.dat['eta','97.5%'] + exp(-CI.dat['lambda','2.5%'] * ((0:(T-1))/365)) * (max(x) - CI.dat['eta','97.5%'])
  
  
  df.mean$log10MFI_lci <-  mu_lci
  df.mean$log10MFI_uci <-  mu_uci
  df.mean$log10MFI_lci[df.mean$log10MFI_lci<0] <- 0
  
  #replot with just the deterministic model
  p1 <- ggplot(dat.mat) + #ggtitle("Simulated maternal antibody trajectories for 100 bats") +
    geom_point(data = dat.mat, aes(x=age,y=log10MFI), size=2) + 
    geom_hline(aes(yintercept =cutoff.df$log10MFI[cutoff.df$type=="mean"]), linetype=2, color="red") +
    theme_bw() + theme(panel.grid = element_blank(),
                       axis.title = element_text(size=18),
                       plot.title = element_text(size=18),
                       axis.text = element_text(size=14)) +
    ylab(bquote(log[10]~MFI))+
    geom_ribbon(data=df.mean, aes(x=age, ymin=log10MFI_lci, ymax=log10MFI_uci), fill="cornflowerblue", alpha=.2) +
    geom_line(data=df.mean, aes(x=age, y=log10MFI), color="cornflowerblue", size=1)
  
  return(p1)  
}
plot.adultAB <- function(det_df, stan_df){
  #and include plot of the cutoff line by species/antigen/era
  cutoff.df <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 6), rep("Rousettus madagascariensis",3)), antigen=c("NiV", "HeV", "NiV", "HeV", "GhV", "AngV", "GhV", "BdbV", "RavV"), era=c(rep("2014-2016",2), rep("2018-2020",7)))
  cutoff.df$MFI = c(757,215,936,2851,1141,1763,2642,5647,6617)
  cutoff.df$antigen <- factor(cutoff.df$antigen, levels=c("NiV", "HeV", "GhV", "AngV", "BdbV", "RavV"))
  cutoff.df = subset(cutoff.df, antigen=="GhV")
  # 
  

  p4 <- ggplot() +
    labs(
      x = "Days since exposure",
      y = "MFI") +
    geom_ribbon(data=det_df, aes(x=t, ymin=MFI_lci, ymax=MFI_uci),alpha=.2, fill="cornflowerblue") +
    geom_line(data=det_df, aes(x=t, y=MFI),size=1, color="cornflowerblue") +
    theme_bw() + coord_cartesian(ylim=c(0,30000), xlim=c(0,150)) + 
    theme(panel.grid = element_blank(),
          plot.margin = margin(t = 1,r = 70,b = 1,l = 1),
          axis.title = element_text(size=18),
          plot.title = element_text(size=18),
          axis.text = element_text(size=14)) +
    geom_hline(aes(yintercept =cutoff.df$MFI[cutoff.df$bat_species=="Eidolon dupreanum" ]), linetype=2, color="red") +
    geom_point(data=stan_df, (aes(x=new_ts, y=MFI, shape=as.factor(bat_id))), size=2, show.legend = F); p4
  #geom_line(data=stan_df, aes(x=t_obs, y=MFI, color=as.factor(bat_id)), linetype=2); p4
  
  return(p4)
  
}

make.Fig5 <- function(panel1, CI.age.fits, CI.ts.fits, dat1, dat.sero1, dat.prev1, species1, antigen1, cutoff1, vis_split,det.df, stan_df, filename){
  
  img <- image_read_pdf(panel1, density = 300)
  g <- rasterGrob(as.raster(img), interpolate = TRUE)
  p1 <- ggplot() + annotation_custom(g,xmin = 0,xmax = 1,ymin = 0,ymax = 1) + coord_cartesian(xlim = c(-0.03, 1.03),ylim = c(-0.03, 1.03),clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0.01,0.01,0.01,0.01), plot.background = element_rect(fill="white"))
  p2 <- plot.matAb(dat = dat1)
  p3 <- plot.adultAB(det_df = det.df, stan_df = stan_df)
  p4 <- return.models.data.CI.one(CI.age.fits = CI.age.fits,  CI.ts.fits = CI.ts.fits, dat1=dat1,  dat.sero1 = dat.sero1,dat.prev1 =dat.prev1, antigen1 =antigen1, cutoff1=cutoff1,species1=species1, vis_split=vis_split)
  
  
  
  ptop <- cowplot::plot_grid(p1, p2,p3, ncol=3, nrow=1, rel_widths = c(1.1,1,1.3), labels = c("A", "B", "C"), label_size = 22, label_x = c(0,0,-.02))
  
  Fig5 <- cowplot:: plot_grid(ptop, p4, ncol=1, nrow=2)
  
  ggsave(file =filename,
         plot = Fig5,
         units="mm",  
         width=125, 
         height=75, 
         scale=3, 
         dpi=300)
  
}


make.Fig5(panel1 = paste0(homewd, "/final_figures_tables/tmp/MSIRNS.pdf"),
          CI.age.fits = subset(age.sero.df, model=="MSIRNNS-mu-rec"),
          CI.ts.fits = subset(biwk.ts.df, model=="MSIRNNS-mu-rec"),
          dat1=dat, 
          dat.sero1 = sero.dat,
          dat.prev1 = shed.dat,
          antigen1 = "GhV", 
          cutoff1="mean", 
          species1 = "Eidolon dupreanum",  
          vis_split = 3,
          det.df=read.csv(file =paste0(homewd, "/serology/data/stan_fit_deterministic_df_adultAB.csv"), header = T, stringsAsFactors = F), 
          stan_df=read.csv(file = paste0(homewd, "/serology/data/stan_fitted_bats.csv"), header = T, stringsAsFactors = F),
          filename =  paste0(homewd, "/final_figures_tables/tmp/Fig5.png"))


