rm(list=ls())
#set wd
#.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4")

library(plyr)
library(dplyr)   
library(epitools)
library(cowplot)
library(deSolve)
library(mgcv)
library(lubridate)
library(matrixcalc)
library(Matrix)
library(ggplot2)
library(mvtnorm)
library(corpcor)

#load data
#homewd = "/Users/carabrook/Developer/luminex-serology-madabat/scripts/age-sero-fitting/make-CIs/mean-GhV"
#setwd(homewd)
#dat <- read.csv(file =  paste0(homewd,"/data/age-seroprevalence-data.csv"), header = T, stringsAsFactors = F)
dat <- read.csv(file =  "age-seroprevalence-data.csv", header = T, stringsAsFactors = F)
head(dat)

#helper functions

transform.vect <- function(vec, s, c){
  vec2 <- t(commutation.matrix(r=s,c=c))%*%vec
  return(vec2)
}
find.biweek = function(t, times){
  
  biwks <- sort(unique(round(revtrunc(times),4)))
  this.wk = round(revtrunc(times[t]), 4)
  this.biwk <- which(biwks==this.wk)
  
  
  return(this.biwk)
  
  
}
get.age.struct = function(pop, s){ #function that collapses population vector in disease state form to age structure
  age.mat <- mat_split(M=matrix(pop, ncol=1), r=s, c=1)
  age.mat.list <- c()
  for (i in 1:dim(age.mat)[3]){
    age.mat.list[[i]] <- age.mat[,,i]
  }
  age.mat <- Reduce('+', age.mat.list)
  return(age.mat)
}
get.age.struct.M = function(pop, c){ #function that collapses population vector in disease state form to age structure
  age.mat <- mat_split(M=matrix(pop, ncol=1), r=c, c=1)
  age.mat.list <- c()
  for (i in 1:dim(age.mat)[3]){
    age.mat.list[[i]] <- age.mat[,,i]
  }
  #zage.mat <- Reduce('+', age.mat.list)
  age.mat <- lapply(age.mat.list, sum)
  age.mat.dat <- c( c(1:length(age.mat)), unlist(age.mat))
  names(age.mat.dat) <- c("age", "pop")
  return(age.mat)
}
prev.by.age <- function(dat){
  dat$n_age = sum(dat$count)
  dat$prev = dat$count/dat$n_age
  #and seroprev
  dat$seropos = sum(dat$count[dat$class=="M"], dat$count[dat$class=="R"])
  dat$seroprev = dat$seropos/dat$n_age
  return(dat)
}
revtrunc = function(x){
  newx = x - floor(x)
  return(newx)
}
spline.fit = function(data){
  #replace any Inf vals only the okay values
  #data <- data[complete.cases(data),]
  data$seroprevalence[data$seroprevalence==Inf] <- 1
  data$seroprevalence[data$seroprevalence<0] <- 0
  # if(length(data$age)>1){
  spline1 = with(data, smooth.spline(x=age, y=seroprevalence)) 
  #  return(spline1)
  # }else{
  #  return(NA)
  #  }
}
pred.fit <- function(spline1, ages.new){
  prev.comp = predict(spline1, x=ages.new)
  return(prev.comp$y)
}
get.seas.seroprev = function(dat){
  #get doy in bat calendar
  #convert to biweek
  #calc seroprev by biweek
  
  dat$doy <- yday(dat$date)
  #now correct for birthday of bat in questiob
  if(unique(dat$species=="Pteropus rufus")){
    bday <- yday("2015-10-01") #320
  } else if(unique(dat$species=="Eidolon dupreanum")){
    bday = yday("2015-11-01")
  }
  #write over for doy
  dat$new_doy = NA
  for (i in 1:length(dat$doy)){
    if(dat$doy[i] > bday){
      dat$new_doy[i] <- dat$doy[i] - bday
    } else if (dat$doy[i] <= bday){
      dat$new_doy[i] <- dat$doy[i] + (365-bday)
    }
  }
  
  dat$doy <- dat$new_doy
  dat <- dplyr::select(dat, -(new_doy))
  
  #now convert to biweek
  dat$biwk <- NA
  brk <- seq(0,365, by=14)
  for (i in 1:length(dat$biwk)){
    tmp <- brk[dat$doy[i] >=  brk]
    tmp2 <- tmp[length(tmp)] 
    dat$biwk[i] <- which(brk==tmp2)
  }
  
  #now calc seroprev by biweek
  biwk.sero <- ddply(dat, .(biwk), summarize, seropos = sum(prev), sero_lci = sum(prev_lci), sero_uci = sum(prev_uci), n=length(prev))  
  biwk.sero$seroprev = biwk.sero$seropos/biwk.sero$n
  biwk.sero$seroprev_lci = biwk.sero$sero_lci/biwk.sero$n
  biwk.sero$seroprev_uci = biwk.sero$sero_uci/biwk.sero$n
  biwk.sero$biwk = as.numeric(biwk.sero$biwk)
  
  return(biwk.sero)
}
get.mod.seas = function(mod.out){
  #convert time to doy
  mod.out$doy = mod.out$time*365
  #then to biwk
  brk <- seq(0,365, by=14) #doy breaks by biweek
  brk <- brk[-length(brk)]
  #brk.seq = 1:length(brk)
  mod.out$biwk <- NA
  for (i in 1:length(mod.out$biwk)){
    tmp <- brk[mod.out$doy[i] >=  brk]
    tmp2 <- tmp[length(tmp)] 
    mod.out$biwk[i] <- which(brk==tmp2)
  }
  #then return
  return(mod.out)
}
age.sero.plot = function(dat){
  with(dat, plot(age, seroprevalence, type = "b", ylim=c(0,1)))
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
build.pop.mat = function(surv, surv_juv, s, adult_fec){
  pop.mat = matrix(0,  nrow=(s-1), ncol = (s-1))
  diag(pop.mat) = surv
  diag(pop.mat)[1] = surv_juv
  col_s = c(rep(0, s-2), surv)
  pop.mat = cbind(pop.mat, col_s)
  row1 = c(0, rep((adult_fec*surv), (s-1))) #bats reproduce for the first time at the end of the second year of life. good for E. dup and P. ruf
  pop.mat = rbind(row1,pop.mat)
  return(pop.mat)
  
}#for stable age distribution
mat_split <- function(M, r, c){
  nr <- ceiling(nrow(M)/r)
  nc <- ceiling(ncol(M)/c)
  newM <- matrix(NA, nr*r, nc*c)
  newM[1:nrow(M), 1:ncol(M)] <- M
  
  div_k <- kronecker(matrix(seq_len(nr*nc), nr, byrow = TRUE), matrix(1, r, c))
  matlist <- split(newM, div_k)
  N <- length(matlist)
  mats <- unlist(matlist)
  dim(mats)<-c(r, c, N)
  return(mats)
}#for matrix slicing. give it the number of rows and columns you want in the resulting matrices
stack.age = function(dat, s){
  dat.new= list()
  for (i in 1:s){
    dat.new[[i]] = dat[i,]
  }
  dat.new = c(unlist(dat.new))
  return(dat.new)
}
stack.class = function(mat.dat, c){
  
  new.dat = list()
  for (i in 1:c){
    new.dat[[i]] = mat.dat[i,]
  }
  new.dat = c(unlist(new.dat))
  return(new.dat)
}
get.seroprev.biwk.bin.model = function(data, bin){
  
  #add in every second or third biwk
  
  biwk.df = cbind.data.frame(biwk=1:26)
  tmp=sort(seq(1,26, by=bin))
  biwk.df$biwk_bin <- sort(rep(tmp, length.out=length(biwk.df$biwk)))
  
  data <- merge(data, biwk.df, by="biwk")
  
  #if(cutoff=="mean"){
  dat.sum2 <- ddply(data, .(biwk_bin, age), summarize,  seroprevalence=sum(seropos)/length(seropos), count=length(prev))  
  #}else if (cutoff=="uci"){
  #  dat.sum2 <- ddply(data, .(biwk_bin, age), summarize,  seroprevalence=sum(seropos)/length(prev_uci), count=length(prev_uci)) 
  #}else if (cutoff=="lci"){
  # dat.sum2 <- ddply(data, .(biwk_bin, age), summarize,  seroprevalence=sum(seropos)/length(prev_lci), count=length(prev_lci)) 
  #}
  # 
  # #you want the midpoint between either end of the age class now
  # vect.age.yr = sort(unique(data$age_year))
  # dat.sum2$age_plot = NA
  # 
  # for (i in 1:(length(dat.sum2$age_year)-1)){
  #   dat.sum2$age_plot [i] = ((dat.sum2$age_year[i + 1] - dat.sum2$age_year[i])/2) + dat.sum2$age_year[i]
  # }
  # dat.sum2$age_plot[length(dat.sum2$age_plot)] =  (ceiling(max(data$age)) - dat.sum2$age_year[length(dat.sum2$age_year)])/2  + dat.sum2$age_year[length(dat.sum2$age_year)]
  # 
  # names(dat.sum2)  <- c("biwk_bin", "real_age", "prevalence", "count", "age_year") #<- names(dat.sum2.tmp)
  # 
  # # dat.sum2 <-  rbind(dat.sum2.tmp, dat.sum2)
  # 
  dat.sum3 = dat.sum2
  dat.sum3$class = "seropositive"
  rownames(dat.sum3) <- c()
  
  return(dat.sum3)
  
}

run.model <- function(par.df, s, biwk_bin, sim_pop, yrs, ntyr, N_stat, age.brk, slope.wane, burnin, do.plot, do.save, filename){
  
  
  #convert input pars to variables
  model = sapply(strsplit(par.df$model, "-mu"),'[',1)
  recov = par.df$recov
  sero = par.df$sero
  wane = par.df$wane
  sigma=par.df$sigma
  epsilon=par.df$epsilon
  sigma2=par.df$sigma2
  beta2=par.df$beta2
  mort = par.df$mort
  mort_juv = par.df$mort_juv
  recru = par.df$recru
  recru2 = par.df$recru2
  adult_fec = par.df$adult_fec
  mu.sick = par.df$mu_sick
  mu.rec = par.df$mu_rec
  add.inf.mort = (if(par.df$mu_sick>1){TRUE}else{FALSE})
  add.rec.mort = (if(par.df$mu_rec>1){TRUE}else{FALSE})
  beta = par.df$beta
  
  
  
  # choose and run model:
  
  if (model=="MSIR"){
    
    
    #run it to get the mean
    model.out  = sim.met.MSIR.age(burnin=burnin,
                                  sim_pop=sim_pop,
                                  yrs=yrs,
                                  ntyr=ntyr,
                                  s=s,
                                  beta=beta,
                                  age.brk = age.brk,
                                  sero=sero,
                                  recov=recov,
                                  mort=mort, 
                                  mort_juv=mort_juv, 
                                  adult_fec=adult_fec, 
                                  wane=wane, 
                                  slope.wane=slope.wane,
                                  sigma=sigma,
                                  model=model,
                                  mu.sick = mu.sick,
                                  mu.rec = mu.rec,
                                  add.inf.mort = add.inf.mort,
                                  add.rec.mort = add.rec.mort)  
    
  }else if (model=="MSIRS"){
    
    
    
    model.out  = sim.met.MSIR.age(burnin=burnin,
                                  sim_pop=sim_pop,
                                  yrs=yrs,
                                  ntyr=ntyr,
                                  s=s,
                                  beta=beta,
                                  age.brk = age.brk,
                                  sero=sero,
                                  recov=recov,
                                  mort=mort, 
                                  mort_juv=mort_juv, 
                                  adult_fec=adult_fec, 
                                  wane=wane, 
                                  slope.wane=slope.wane,
                                  sigma=sigma,
                                  model=model,
                                  mu.sick = mu.sick,
                                  mu.rec = mu.rec,
                                  add.inf.mort = add.inf.mort,
                                  add.rec.mort = add.rec.mort) 
    
    
  }else if (model =="MSIRN" | model =="MSIRNI"){
    
    model.out  = sim.met.MSIRN.age(burnin=burnin,
                                   sim_pop=sim_pop,
                                   yrs=yrs,
                                   ntyr=ntyr,
                                   s=s,
                                   beta=beta,
                                   age.brk = age.brk,
                                   recov=recov,
                                   mort=mort, 
                                   mort_juv=mort_juv,
                                   adult_fec=adult_fec, 
                                   wane=wane,
                                   recru=recru,
                                   slope.wane=slope.wane,
                                   sigma=sigma,
                                   model=model,
                                   mu.sick = mu.sick,
                                   mu.rec = mu.rec,
                                   add.inf.mort = add.inf.mort,
                                   N_stat=N_stat,
                                   add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRNS"){
    
    
    model.out  = sim.met.MSIRNS.age(burnin=burnin,
                                    sim_pop=sim_pop,
                                    yrs=yrs,
                                    ntyr=ntyr,
                                    s=s,
                                    beta=beta,
                                    age.brk = age.brk,
                                    recov=recov,
                                    mort=mort, 
                                    mort_juv=mort_juv,
                                    adult_fec=adult_fec, 
                                    wane=wane,
                                    recru=recru,
                                    slope.wane=slope.wane,
                                    sigma=sigma,
                                    model=model,
                                    mu.sick = mu.sick,
                                    mu.rec = mu.rec,
                                    add.inf.mort = add.inf.mort,
                                    N_stat=N_stat,
                                    add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRNNS"){
    
    
    model.out  = sim.met.MSIRNNS.age(burnin=burnin,
                                     sim_pop=sim_pop,
                                     yrs=yrs,
                                     ntyr=ntyr,
                                     s=s,
                                     beta=beta,
                                     age.brk = age.brk,
                                     recov=recov,
                                     mort=mort, 
                                     mort_juv=mort_juv,
                                     adult_fec=adult_fec, 
                                     wane=wane,
                                     recru=recru,
                                     epsilon=epsilon,
                                     slope.wane=slope.wane,
                                     sigma=sigma,
                                     model=model,
                                     mu.sick = mu.sick,
                                     mu.rec = mu.rec,
                                     add.inf.mort = add.inf.mort,
                                     N_stat=N_stat,
                                     add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRNNI"){
    
    
    model.out  = sim.met.MSIRNNI.age(burnin=burnin,
                                     sim_pop=sim_pop,
                                     yrs=yrs,
                                     ntyr=ntyr,
                                     s=s,
                                     beta=beta,
                                     age.brk = age.brk,
                                     recov=recov,
                                     mort=mort, 
                                     mort_juv=mort_juv,
                                     adult_fec=adult_fec, 
                                     wane=wane,
                                     recru=recru,
                                     epsilon=epsilon,
                                     slope.wane=slope.wane,
                                     sigma=sigma,
                                     model=model,
                                     mu.sick = mu.sick,
                                     mu.rec = mu.rec,
                                     add.inf.mort = add.inf.mort,
                                     N_stat=N_stat,
                                     add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRNIRN"){
    
    
    model.out  = sim.met.MSIRNIRN.age(burnin=burnin,
                                      sim_pop=sim_pop,
                                      yrs=yrs,
                                      ntyr=ntyr,
                                      s=s,
                                      beta=beta,
                                      age.brk = age.brk,
                                      recov=recov,
                                      mort=mort, 
                                      mort_juv=mort_juv,
                                      adult_fec=adult_fec, 
                                      wane=wane,
                                      recru=recru,
                                      slope.wane=slope.wane,
                                      sigma=sigma,
                                      sigma2=sigma2,
                                      model=model,
                                      mu.sick = mu.sick,
                                      mu.rec = mu.rec,
                                      add.inf.mort = add.inf.mort,
                                      N_stat=N_stat,
                                      add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRNSIRN"){
    
    
    model.out  = sim.met.MSIRNSIRN.age(burnin=burnin,
                                       sim_pop=sim_pop,
                                       yrs=yrs,
                                       ntyr=ntyr,
                                       s=s,
                                       beta=beta,
                                       age.brk = age.brk,
                                       recov=recov,
                                       mort=mort, 
                                       mort_juv=mort_juv,
                                       adult_fec=adult_fec, 
                                       wane=wane,
                                       recru=recru,
                                       slope.wane=slope.wane,
                                       sigma=sigma,
                                       sigma2=sigma2,
                                       beta2=beta2,
                                       model=model,
                                       mu.sick = mu.sick,
                                       mu.rec = mu.rec,
                                       add.inf.mort = add.inf.mort,
                                       N_stat=N_stat,
                                       add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRSIR"){
    
    
    model.out  = sim.met.MSIRSIR.age(burnin=burnin,
                                     sim_pop=sim_pop,
                                     yrs=yrs,
                                     ntyr=ntyr,
                                     s=s,
                                     beta=beta,
                                     age.brk = age.brk,
                                     recov=recov,
                                     mort=mort, 
                                     mort_juv=mort_juv,
                                     adult_fec=adult_fec, 
                                     wane=wane,
                                     slope.wane=slope.wane,
                                     sigma=sigma,
                                     beta2=beta2,
                                     model=model,
                                     mu.sick = mu.sick,
                                     mu.rec = mu.rec,
                                     add.inf.mort = add.inf.mort,
                                     add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRNIRN-LT"){
    
    
    model.out  = sim.met.MSIRNIRN.LT.age(burnin=burnin,
                                         sim_pop=sim_pop,
                                         yrs=yrs,
                                         ntyr=ntyr,
                                         s=s,
                                         beta=beta,
                                         age.brk = age.brk,
                                         recov=recov,
                                         mort=mort, 
                                         mort_juv=mort_juv,
                                         adult_fec=adult_fec, 
                                         wane=wane,
                                         recru=recru,
                                         slope.wane=slope.wane,
                                         sigma=sigma,
                                         sigma2=sigma2,
                                         epsilon=epsilon,
                                         model=model,
                                         mu.sick = mu.sick,
                                         mu.rec = mu.rec,
                                         add.inf.mort = add.inf.mort,
                                         N_stat=N_stat,
                                         add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRNSIRN-LT"){
    
    
    model.out  = sim.met.MSIRNSIRN.LT.age(burnin=burnin,
                                          sim_pop=sim_pop,
                                          yrs=yrs,
                                          ntyr=ntyr,
                                          s=s,
                                          beta=beta,
                                          age.brk = age.brk,
                                          recov=recov,
                                          mort=mort, 
                                          mort_juv=mort_juv,
                                          adult_fec=adult_fec, 
                                          wane=wane,
                                          recru=recru,
                                          slope.wane=slope.wane,
                                          sigma=sigma,
                                          sigma2=sigma2,
                                          beta2=beta2,
                                          epsilon=epsilon,
                                          model=model,
                                          mu.sick = mu.sick,
                                          mu.rec = mu.rec,
                                          add.inf.mort = add.inf.mort,
                                          N_stat=N_stat,
                                          add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRNIR"){
    
    model.out  = sim.met.MSIRNIR.age(burnin=burnin,
                                     sim_pop=sim_pop,
                                     yrs=yrs,
                                     ntyr=ntyr,
                                     s=s,
                                     beta=beta,
                                     age.brk = age.brk,
                                     recov=recov,
                                     mort=mort, 
                                     mort_juv=mort_juv,
                                     adult_fec=adult_fec, 
                                     wane=wane,
                                     recru=recru,
                                     slope.wane=slope.wane,
                                     sigma=sigma,
                                     sigma2=sigma2,
                                     model=model,
                                     mu.sick = mu.sick,
                                     mu.rec = mu.rec,
                                     add.inf.mort = add.inf.mort,
                                     N_stat=N_stat,
                                     add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRNSIRNS"){
    
    
    
    model.out  = sim.met.MSIRNSIRNS.age(burnin=burnin,
                                        sim_pop=sim_pop,
                                        yrs=yrs,
                                        ntyr=ntyr,
                                        s=s,
                                        beta=beta,
                                        age.brk = age.brk,
                                        recov=recov,
                                        mort=mort, 
                                        mort_juv=mort_juv,
                                        adult_fec=adult_fec, 
                                        wane=wane,
                                        recru=recru,
                                        slope.wane=slope.wane,
                                        sigma=sigma,
                                        sigma2=sigma2,
                                        beta2=beta2,
                                        recru2=recru2,
                                        model=model,
                                        mu.sick = mu.sick,
                                        mu.rec = mu.rec,
                                        add.inf.mort = add.inf.mort,
                                        N_stat=N_stat,
                                        add.rec.mort = add.rec.mort) 
    
  } 
  
  
  
  
  #and plot
  if(do.plot==TRUE){
    
    p1 <- ggplot() + theme_bw() +
      geom_line(data=model.out, aes(x=age, y=seroprevalence, color=biwk)) + #facet_wrap(model~.) +
      theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
      coord_cartesian(ylim=c(0,1))
    
    print(p1)
    
    
    if(do.save==TRUE){
      ggsave(file = filename,
             units="mm",  
             width=80, 
             height=60, 
             scale=3, 
             dpi=300)
    }
    
    
  }
  
  #return your model stats, including estimated pars, the lci/uci and the output seroprev with CIs
  
  if(add.inf.mort==TRUE & add.rec.mort==FALSE){
    model = paste0(model,"-mu-sick")
  }else if (add.inf.mort==FALSE & add.rec.mort==TRUE){
    model = paste0(model,"-mu-rec")
  }else if (add.inf.mort==TRUE &  add.rec.mort==TRUE){
    model = paste0(model,"-mu-sick-mu-rec")
  }
  
  model.out$model = model
  
  
  
  
  return(model.out)
}


#build functions
buildTMat_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	beta, recov, sigma, sero, wane, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  mort.vect.sick <- c((1-surv.juv.biwk.sick), rep((1-surv.biwk.sick), (nage-1))) 
  mort.vect.rec <- c((1-surv.juv.biwk.rec), rep((1-surv.biwk.rec), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #psi = 1-exp(-sero*sum(I_m)) #force of seroconversion, also depends on density of I class
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi))) #vector
  psi = 1-exp(-sero*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,4,4) # MSIR
  
  Tmat <- matrix(0,4*nage,4*nage) #MSIR for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]-psi
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,2] <- psi
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[2,4] <- sigma[j]
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],4);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- 1-mort.vect.sick[j]
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too
      surv[4] <- 1-mort.vect.rec[j]
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*4+1):(j*4+4),((j-1)*4+1):(j*4)] <- mat1*surv*age.rate[j]
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up
      Tmat[((j-1)*4+1):(j*4),((j-1)*4+1):(j*4)] <- mat2*surv*(1-age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*4+1):(j*4),((j-1)*4+1):(j*4)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildTMat_MSIRN_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	 beta, recru, recov, sigma, wane, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  mort.vect.sick <- c((1-surv.juv.biwk.sick), rep((1-surv.biwk.sick), (nage-1))) 
  mort.vect.rec <- c((1-surv.juv.biwk.rec), rep((1-surv.biwk.rec), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(recru)==1) recru <- rep(recru,nage)
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  #Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,5,5) # MSIRN
  
  Tmat <- matrix(0,5*nage,5*nage) #MSIRN for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1- recru[j]
    mat1[3,5] <- recru[j] #recrudescence. moving from N back to I
    
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],5);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- 1-mort.vect.sick[j]
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too
      surv[4] <- 1-mort.vect.rec[j]
      surv[5] <- 1-mort.vect.rec[j]
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*5+1):(j*5+5),((j-1)*5+1):(j*5)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildTMat_MSIRNS_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	 beta, recru, recov, sigma, wane, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  mort.vect.sick <- c((1-surv.juv.biwk.sick), rep((1-surv.biwk.sick), (nage-1))) 
  mort.vect.rec <- c((1-surv.juv.biwk.rec), rep((1-surv.biwk.rec), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(recru)==1) recru <- rep(recru,nage)
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  #Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,5,5) # MSIRN
  
  Tmat <- matrix(0,5*nage,5*nage) #MSIRN for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1- recru[j]
    mat1[2,5] <- recru[j] #recrudescence. moving from N back to I
    
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],5);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- 1-mort.vect.sick[j]
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too
      surv[4] <- 1-mort.vect.rec[j]
      surv[5] <- 1-mort.vect.rec[j]
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*5+1):(j*5+5),((j-1)*5+1):(j*5)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildTMat_MSIRNNI_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	 beta, recru, recov, sigma, wane, epsilon, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  mort.vect.sick <- c((1-surv.juv.biwk.sick), rep((1-surv.biwk.sick), (nage-1))) 
  mort.vect.rec <- c((1-surv.juv.biwk.rec), rep((1-surv.biwk.rec), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(recru)==1) recru <- rep(recru,nage)
  if (length(epsilon)==1) epsilon <- rep(epsilon,nage)
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  #Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,6,6) # MSIRNN
  
  Tmat <- matrix(0,6*nage,6*nage) #MSIRNN for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1- recru[j] - epsilon[j]
    mat1[3,5] <- recru[j] #recrudescence. moving from N1 back to I
    mat1[6,5] <- epsilon[j] #moving from N1 to N2. 2 years
    mat1[6,6] <- 1 #once here, you stay
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],6);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- 1-mort.vect.sick[j]
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too, including N1 and N2
      surv[4] <- 1-mort.vect.rec[j]
      surv[5] <- 1-mort.vect.rec[j]
      surv[6] <- 1-mort.vect.rec[j]
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*6+1):(j*6+6),((j-1)*6+1):(j*6)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*6+1):(j*6),((j-1)*6+1):(j*6)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*6+1):(j*6),((j-1)*6+1):(j*6)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildTMat_MSIRNNS_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	 beta, recru, recov, sigma, wane, epsilon, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  mort.vect.sick <- c((1-surv.juv.biwk.sick), rep((1-surv.biwk.sick), (nage-1))) 
  mort.vect.rec <- c((1-surv.juv.biwk.rec), rep((1-surv.biwk.rec), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(recru)==1) recru <- rep(recru,nage)
  if (length(epsilon)==1) epsilon <- rep(epsilon,nage)
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  #Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,6,6) # MSIRNN
  
  Tmat <- matrix(0,6*nage,6*nage) #MSIRNN for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1- recru[j] - epsilon[j]
    mat1[2,5] <- recru[j] #recrudescence. moving from N1 back to S
    mat1[6,5] <- epsilon[j] #moving from N1 to N2. 2 years
    mat1[6,6] <- 1 #once here, you stay
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],6);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- 1-mort.vect.sick[j]
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too, including N1 and N2
      surv[4] <- 1-mort.vect.rec[j]
      surv[5] <- 1-mort.vect.rec[j]
      surv[6] <- 1-mort.vect.rec[j]
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*6+1):(j*6+6),((j-1)*6+1):(j*6)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*6+1):(j*6),((j-1)*6+1):(j*6)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*6+1):(j*6),((j-1)*6+1):(j*6)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildTMat_MSIRNIRN_LT_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk,surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	 beta, recru, recov, sigma, wane, sigma2, epsilon, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  mort.vect.sick <- c((1-surv.juv.biwk.sick), rep((1-surv.biwk.sick), (nage-1))) 
  mort.vect.rec <- c((1-surv.juv.biwk.rec), rep((1-surv.biwk.rec), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(epsilon)==1) epsilon <- rep(epsilon,nage)
  if (length(recru)==1) recru <- rep(recru,nage)
  if (length(sigma2)==1) sigma2 <- rep(sigma2,nage)
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  I_m2 = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,6] #I2 is #6 
  foi = 1-exp(-beta*((sum(I_m) + sum(I_m2))/sum(Npop_epi))) #includes both types of infecteds
  #Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,8,8) # MSIRNIRN
  
  Tmat <- matrix(0,8*nage,8*nage) #MSIRNIRNfor every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1- recru[j] - epsilon[j]
    mat1[6,5] <- recru[j] #recrudescence. moving from N1 to I2
    mat1[6,6] <- 1-recov #those staying in I2
    mat1[7,6] <-  recov #moving from I2 to R2
    mat1[7,7] <-  1-sigma2[j] #staying in R2
    mat1[8,7] <-  sigma2[j] # loss of ABs, from R2 to N2
    mat1[8,8] <- 1 #once in N2, you stay
    mat1[8,5] <- epsilon[j] #advance to LT immunity (N2)
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],8);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- 1-mort.vect.sick[j]
      surv[6] <- 1-mort.vect.sick[j]
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too, including N1 and N2
      surv[4] <- 1-mort.vect.rec[j]
      surv[5] <- 1-mort.vect.rec[j]
      surv[6] <- 1-mort.vect.rec[j]
      surv[7] <- 1-mort.vect.rec[j]
      surv[8] <- 1-mort.vect.rec[j]
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*8+1):(j*8+8),((j-1)*8+1):(j*8)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*8+1):(j*8),((j-1)*8+1):(j*8)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*8+1):(j*8),((j-1)*8+1):(j*8)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}# this version allows for advance to LT immunity
buildTMat_MSIRNIRN_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk,surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	 beta, recru, recov, sigma, wane, sigma2, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  mort.vect.sick <- c((1-surv.juv.biwk.sick), rep((1-surv.biwk.sick), (nage-1))) 
  mort.vect.rec <- c((1-surv.juv.biwk.rec), rep((1-surv.biwk.rec), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(recru)==1) recru <- rep(recru,nage)
  if (length(sigma2)==1) sigma2 <- rep(sigma2,nage)
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I4 comes #3
  I_m2 = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,6] #I2 is #6 
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  #foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  foi = 1-exp(-beta*((sum(I_m) + sum(I_m2))/sum(Npop_epi))) #includes both types of infecteds
  #Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,8,8) # MSIRNIRN
  
  Tmat <- matrix(0,8*nage,8*nage) #MSIRNIRN for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1- recru[j] 
    mat1[6,5] <- recru[j] #recrudescence. moving from N1 to I2
    mat1[6,6] <- 1-recov #those staying in I2
    mat1[7,6] <-  recov #moving from I2 to R2
    mat1[7,7] <-  1-sigma2[j] #staying in R2
    mat1[8,7] <-  sigma2[j] # loss of ABs, from R2 to N2
    mat1[8,8] <- 1 #once in N2, you stay
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],8);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- 1-mort.vect.sick[j]
      surv[6] <- 1-mort.vect.sick[j]
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too, including N1 and N2
      surv[4] <- 1-mort.vect.rec[j]
      surv[5] <- 1-mort.vect.rec[j]
      surv[6] <- 1-mort.vect.rec[j]
      surv[7] <- 1-mort.vect.rec[j]
      surv[8] <- 1-mort.vect.rec[j]
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*8+1):(j*8+8),((j-1)*8+1):(j*8)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*8+1):(j*8),((j-1)*8+1):(j*8)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*8+1):(j*8),((j-1)*8+1):(j*8)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildTMat_MSIRNIR_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk,surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	 beta, recru, recov, sigma, wane, sigma2, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  mort.vect.sick <- c((1-surv.juv.biwk.sick), rep((1-surv.biwk.sick), (nage-1))) 
  mort.vect.rec <- c((1-surv.juv.biwk.rec), rep((1-surv.biwk.rec), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(recru)==1) recru <- rep(recru,nage)
  if (length(sigma2)==1) sigma2 <- rep(sigma2,nage)
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I4 comes #3
  I_m2 = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,6] #I2 is #6 
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  #foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  foi = 1-exp(-beta*((sum(I_m) + sum(I_m2))/sum(Npop_epi))) #includes both types of infecteds
  #Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,7,7) # MSIRNIR
  
  Tmat <- matrix(0,7*nage,7*nage) #MSIRNIR for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1- recru[j] 
    mat1[6,5] <- recru[j] #recrudescence. moving from N to I2
    mat1[6,6] <- 1-recov #those staying in I2
    mat1[7,6] <-  recov #moving from I2 to R2
    mat1[7,7] <-  1-sigma2[j] #staying in R2. those that do not wane back to N
    mat1[5,7] <-  sigma2[j] # loss of ABs, from R2 back to N
    
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],7);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- 1-mort.vect.sick[j]
      surv[6] <- 1-mort.vect.sick[j]
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too, including N1 and N2
      surv[4] <- 1-mort.vect.rec[j]
      surv[5] <- 1-mort.vect.rec[j]
      surv[6] <- 1-mort.vect.rec[j]
      surv[7] <- 1-mort.vect.rec[j]
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*7+1):(j*7+7),((j-1)*7+1):(j*7)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*7+1):(j*7),((j-1)*7+1):(j*7)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*7+1):(j*7),((j-1)*7+1):(j*7)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildTMat_MSIRNSIRNS_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	 beta, recru, recru2, recov, sigma, wane, sigma2, beta2, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  mort.vect.sick <- c((1-surv.juv.biwk.sick), rep((1-surv.biwk.sick), (nage-1))) 
  mort.vect.rec <- c((1-surv.juv.biwk.rec), rep((1-surv.biwk.rec), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  
  #and beta2
  if (length(beta2)==1) beta2 <- rep(beta2,nage) # 
  if (length(beta2)==s) beta2 <- beta2
  if (length(beta2) > 1 & length(beta2)<s){
    beta2.list <- list()
    for (i in 1:length(beta2)){
      beta2.list[[i]] = rep(beta2[i], age.brk[i])
    }
    beta2 = c(unlist(beta2.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(recru)==1) recru <- rep(recru,nage)
  if (length(recru2)==1) recru2 <- rep(recru2,nage)
  if (length(sigma2)==1) sigma2 <- rep(sigma2,nage)
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I1 comes #3
  I_m2 = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,7] #always. for all model forms, I2 comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  foi = 1-exp(-beta*((sum(I_m) + sum(I_m2))/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  foi2 = 1-exp(-beta2*((sum(I_m) + sum(I_m2))/sum(Npop_epi))) #uses different beta because susceptibles in the second class might be less suscetible
  #Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,9,9) # MSIRNSIRNS
  
  Tmat <- matrix(0,9*nage,9*nage) #MSIRNSIRNS for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1- recru[j] 
    mat1[6,5] <- recru[j] #recrudescence. moving from N1 to S2
    mat1[6,6] <- 1-foi2[j] #those staying in S2 who are not infected
    mat1[6,9] <- recru2[j] #recrudescence. moving from N2 to S2
    mat1[7,6] <-  foi2[j] #moving from S2 to I2
    mat1[7,7] <-  1-recov #staying in I2
    mat1[8,7] <-  recov #moving I2 to R2
    mat1[8,8] <-  1-sigma2[j] #staying in R2
    mat1[9,8] <-  sigma2[j] # loss of ABs, from R2 to N2
    mat1[9,9] <- 1-recru2[j] #those remaining in N2 who don't recrudesce back to S2
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],9);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- 1-mort.vect.sick[j]
      surv[7] <- 1-mort.vect.sick[j]
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too, including N1 and N2
      surv[4] <- 1-mort.vect.rec[j]
      surv[5] <- 1-mort.vect.rec[j]
      surv[6] <- 1-mort.vect.rec[j]
      surv[7] <- 1-mort.vect.rec[j]
      surv[8] <- 1-mort.vect.rec[j]
      surv[9] <- 1-mort.vect.rec[j]
    } 
    
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*9+1):(j*9+9),((j-1)*9+1):(j*9)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*9+1):(j*9),((j-1)*9+1):(j*9)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*9+1):(j*9),((j-1)*9+1):(j*9)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildTMat_MSIRNSIRN_LT_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	 beta, recru, recov, sigma, wane, sigma2, beta2, epsilon, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  mort.vect.sick <- c((1-surv.juv.biwk.sick), rep((1-surv.biwk.sick), (nage-1))) 
  mort.vect.rec <- c((1-surv.juv.biwk.rec), rep((1-surv.biwk.rec), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  
  #and beta2
  if (length(beta2)==1) beta2 <- rep(beta2,nage) # 
  if (length(beta2)==s) beta2 <- beta2
  if (length(beta2) > 1 & length(beta2)<s){
    beta2.list <- list()
    for (i in 1:length(beta2)){
      beta2.list[[i]] = rep(beta2[i], age.brk[i])
    }
    beta2 = c(unlist(beta2.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(epsilon)==1) epsilon <- rep(epsilon,nage)
  if (length(recru)==1) recru <- rep(recru,nage)
  if (length(sigma2)==1) sigma2 <- rep(sigma2,nage)
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I1 comes #3
  I_m2 = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,7] #always. for all model forms, I2 comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  foi = 1-exp(-beta*((sum(I_m) + sum(I_m2))/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  foi2 = 1-exp(-beta2*((sum(I_m) + sum(I_m2))/sum(Npop_epi))) #uses different beta because susceptibles in the second class might be less suscetible
  #Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,9,9) # MSIRNSIRN
  
  Tmat <- matrix(0,9*nage,9*nage) #MSIRNSIRN for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1- recru[j] -epsilon[j]
    mat1[6,5] <- recru[j] #recrudescence. moving from N1 to S2
    mat1[6,6] <- 1-foi2[j] #those staying in S2 who are not infected
    mat1[7,6] <-  foi2[j] #moving from S2 to I2
    mat1[7,7] <-  1-recov #staying in I2
    mat1[8,7] <-  recov #moving I2 to R2
    mat1[8,8] <-  1-sigma2[j] #staying in R2
    mat1[9,8] <-  sigma2[j] # loss of ABs, from R2 to N2
    mat1[9,9] <- 1 #once in N2, you stay
    mat1[9,5] <- epsilon[j] #advance directly from N1 to N2
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],9);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- 1-mort.vect.sick[j]
      surv[7] <- 1-mort.vect.sick[j]
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too, including N1 and N2
      surv[4] <- 1-mort.vect.rec[j]
      surv[5] <- 1-mort.vect.rec[j]
      surv[6] <- 1-mort.vect.rec[j]
      surv[7] <- 1-mort.vect.rec[j]
      surv[8] <- 1-mort.vect.rec[j]
      surv[9] <- 1-mort.vect.rec[j]
    } 
    
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*9+1):(j*9+9),((j-1)*9+1):(j*9)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*9+1):(j*9),((j-1)*9+1):(j*9)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*9+1):(j*9),((j-1)*9+1):(j*9)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildTMat_MSIRNSIRN_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	 beta, recru, recov, sigma, wane, sigma2, beta2, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  mort.vect.sick <- c((1-surv.juv.biwk.sick), rep((1-surv.biwk.sick), (nage-1))) 
  mort.vect.rec <- c((1-surv.juv.biwk.rec), rep((1-surv.biwk.rec), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  
  #and beta2
  if (length(beta2)==1) beta2 <- rep(beta2,nage) # 
  if (length(beta2)==s) beta2 <- beta2
  if (length(beta2) > 1 & length(beta2)<s){
    beta2.list <- list()
    for (i in 1:length(beta2)){
      beta2.list[[i]] = rep(beta2[i], age.brk[i])
    }
    beta2 = c(unlist(beta2.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(recru)==1) recru <- rep(recru,nage)
  if (length(sigma2)==1) sigma2 <- rep(sigma2,nage)
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I1 comes #3
  I_m2 = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,7] #always. for all model forms, I2 comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  foi = 1-exp(-beta*((sum(I_m) + sum(I_m2))/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  foi2 = 1-exp(-beta2*((sum(I_m) + sum(I_m2))/sum(Npop_epi))) #uses different beta because susceptibles in the second class might be less suscetible
  #Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,9,9) # MSIRNSIRN
  
  Tmat <- matrix(0,9*nage,9*nage) #MSIRNSIRN for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1- recru[j] 
    mat1[6,5] <- recru[j] #recrudescence. moving from N1 to S2
    mat1[6,6] <- 1-foi2[j] #those staying in S2 who are not infected
    mat1[7,6] <-  foi2[j] #moving from S2 to I2
    mat1[7,7] <-  1-recov #staying in I2
    mat1[8,7] <-  recov #moving I2 to R2
    mat1[8,8] <-  1-sigma2[j] #staying in R2
    mat1[9,8] <-  sigma2[j] # loss of ABs, from R2 to N2
    mat1[9,9] <- 1 #once in N2, you stay
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],9);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- 1-mort.vect.sick[j]
      surv[7] <- 1-mort.vect.sick[j]
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too, including N1 and N2
      surv[4] <- 1-mort.vect.rec[j]
      surv[5] <- 1-mort.vect.rec[j]
      surv[6] <- 1-mort.vect.rec[j]
      surv[7] <- 1-mort.vect.rec[j]
      surv[8] <- 1-mort.vect.rec[j]
      surv[9] <- 1-mort.vect.rec[j]
    } 
    
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*9+1):(j*9+9),((j-1)*9+1):(j*9)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*9+1):(j*9),((j-1)*9+1):(j*9)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*9+1):(j*9),((j-1)*9+1):(j*9)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildTMat_MSIRSIR_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	beta, recov, sigma, wane, beta2, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  mort.vect.sick <- c((1-surv.juv.biwk.sick), rep((1-surv.biwk.sick), (nage-1))) 
  mort.vect.rec <- c((1-surv.juv.biwk.rec), rep((1-surv.biwk.rec), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  
  #and beta2
  if (length(beta2)==1) beta2 <- rep(beta2,nage) # 
  if (length(beta2)==s) beta2 <- beta2
  if (length(beta2) > 1 & length(beta2)<s){
    beta2.list <- list()
    for (i in 1:length(beta2)){
      beta2.list[[i]] = rep(beta2[i], age.brk[i])
    }
    beta2 = c(unlist(beta2.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  #if (length(recru)==1) recru <- rep(recru,nage)
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I1 comes #3
  I_m2 = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,6] #I2 at #6
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  foi = 1-exp(-beta*((sum(I_m) + sum(I_m2))/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  foi2 = 1-exp(-beta2*((sum(I_m) + sum(I_m2))/sum(Npop_epi))) #uses different beta because susceptibles in the second class might be less suscetible
  #Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,7,7) # MSIRSIR
  
  Tmat <- matrix(0,7*nage,7*nage) #MSIRSIR for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #those the remain Ab+ in R1
    mat1[5,4] <- sigma[j] #Ab wane into S2 class
    mat1[5,5] <- 1-foi2[j] #those staying in S2 who are not infected
    mat1[6,5] <-  foi2[j] #moving from S2 to I2
    mat1[6,6] <-  1-recov #staying in I2
    mat1[7,6] <-  recov #moving I2 to R2
    mat1[7,7] <- 1 #once in R2, you stay
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],7);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    
    
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- 1-mort.vect.sick[j]
      surv[6] <- 1-mort.vect.sick[j]
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too, including N1 and N2
      surv[4] <- 1-mort.vect.rec[j]
      surv[5] <- 1-mort.vect.rec[j]
      surv[6] <- 1-mort.vect.rec[j]
      surv[7] <- 1-mort.vect.rec[j]
    } 
    
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*7+1):(j*7+7),((j-1)*7+1):(j*7)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*7+1):(j*7),((j-1)*7+1):(j*7)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*7+1):(j*7),((j-1)*7+1):(j*7)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}

buildFMatrix <- function(age.classes, adult_fec, surv.biwk, biwk){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = adult_fec*surv.biwk*.2
  }else if (biwk==3|biwk==25){
    new.fec = adult_fec*surv.biwk*.1
  }else if (biwk==4|biwk ==24){
    new.fec = adult_fec*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,4*nage,4*nage)
  
  for (j in 1:nage) {
    Fmat[1,((j-1)*4+1):(j*4)] <- c(0,0,fert.biwk[j],fert.biwk[j]) #the mat immune
    Fmat[2,((j-1)*4+1):(j*4)]<- c(fert.biwk[j],fert.biwk[j],0,0) #the susceptible (mom didn't get sick)
  }
  
  return(Fmat)
}
buildFMatrix_MSIRN <- function(age.classes, adult_fec, surv.biwk, biwk, N_stat){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = adult_fec*surv.biwk*.25
  }else if (biwk==3 | biwk==25){
    new.fec = adult_fec*surv.biwk*.1
    #}else if (biwk==4|biwk ==24){
    #new.fec = adult_fec*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,5*nage,5*nage)
  
  for (j in 1:nage) { #no fertility in first age class (mat immune)
    if(N_stat=="matAB"){
      #try it assuming N moms produce maternally immune pups
      Fmat[1,((j-1)*5+1):(j*5)] <- c(0,0,fert.biwk[j],fert.biwk[j], fert.biwk[j]) #the mat immune offspring (so mom = I, R, and N)
      Fmat[2,((j-1)*5+1):(j*5)]<-  c(fert.biwk[j],fert.biwk[j],0,0, 0) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S)
      
    }else if (N_stat=="matSus"){
      Fmat[1,((j-1)*5+1):(j*5)] <- c(0,0,fert.biwk[j],fert.biwk[j], 0) #the mat immune offpring (so mom = I and R but not N)
      Fmat[2,((j-1)*5+1):(j*5)]<- c(fert.biwk[j],fert.biwk[j],0,0, fert.biwk[j]) #the susceptible offspring (so mom= M, S, and N)
    }
    
    
    
  }
  
  return(Fmat)
}
buildFMatrix_MSIRN <- function(age.classes, adult_fec, surv.biwk, biwk, N_stat){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = adult_fec*surv.biwk*.25
  }else if (biwk==3 | biwk==25){
    new.fec = adult_fec*surv.biwk*.1
    #}else if (biwk==4|biwk ==24){
    #new.fec = adult_fec*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,5*nage,5*nage)
  
  for (j in 1:nage) { #no fertility in first age class (mat immune)
    if(N_stat=="matAB"){
      #try it assuming N moms produce maternally immune pups
      Fmat[1,((j-1)*5+1):(j*5)] <- c(0,0,fert.biwk[j],fert.biwk[j], fert.biwk[j]) #the mat immune offspring (so mom = I, R, and N)
      Fmat[2,((j-1)*5+1):(j*5)]<-  c(fert.biwk[j],fert.biwk[j],0,0, 0) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S)
      
    }else if (N_stat=="matSus"){
      Fmat[1,((j-1)*5+1):(j*5)] <- c(0,0,fert.biwk[j],fert.biwk[j], 0) #the mat immune offpring (so mom = I and R but not N)
      Fmat[2,((j-1)*5+1):(j*5)]<- c(fert.biwk[j],fert.biwk[j],0,0, fert.biwk[j]) #the susceptible offspring (so mom= M, S, and N)
    }
    
    
    
  }
  
  return(Fmat)
}
buildFMatrix_MSIRNN <- function(age.classes, adult_fec, surv.biwk, biwk, N_stat){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = adult_fec*surv.biwk*.25
  }else if (biwk==3 | biwk==25){
    new.fec = adult_fec*surv.biwk*.1
    #}else if (biwk==4|biwk ==24){
    #new.fec = adult_fec*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,6*nage,6*nage)
  
  for (j in 1:nage) { #no fertility in first age class (mat immune)
    if(N_stat=="matAB"){
      #try it assuming N moms produce maternally immune pups
      Fmat[1,((j-1)*6+1):(j*6)] <- c(0,0,fert.biwk[j],fert.biwk[j], fert.biwk[j], fert.biwk[j]) #the mat immune offspring (so mom = I, R, N1, N2)
      Fmat[2,((j-1)*6+1):(j*6)]<-  c(fert.biwk[j],fert.biwk[j],0,0, 0,0) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S)
      
    }else if (N_stat=="matSus"){
      Fmat[1,((j-1)*6+1):(j*6)] <- c(0,0,fert.biwk[j],fert.biwk[j], 0,0) #the mat immune offpring (so mom = I and R but not N1 or N2)
      Fmat[2,((j-1)*6+1):(j*6)]<- c(fert.biwk[j],fert.biwk[j],0,0, fert.biwk[j], fert.biwk[j]) #the susceptible offspring (so mom= M, S, and N1, N2)
    }
    
    
    
  }
  
  return(Fmat)
}
buildFMatrix_MSIRNIRN <- function(age.classes, adult_fec, surv.biwk, biwk, N_stat){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = adult_fec*surv.biwk*.25
  }else if (biwk==3 | biwk==25){
    new.fec = adult_fec*surv.biwk*.1
    #}else if (biwk==4|biwk ==24){
    #new.fec = adult_fec*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,8*nage,8*nage)
  
  for (j in 1:nage) { #no fertility in first age class (mat immune)
    if(N_stat=="matAB"){
      #try it assuming any mom who has previously seen infection will produce maternally immune pups
      Fmat[1,((j-1)*8+1):(j*8)] <- c(0,0,fert.biwk[j],fert.biwk[j], fert.biwk[j], fert.biwk[j], fert.biwk[j], fert.biwk[j]) #the mat immune offspring (so mom = I1, R1, N1,I2,R2,N2)
      Fmat[2,((j-1)*8+1):(j*8)]<-  c(fert.biwk[j],fert.biwk[j],0,0, 0,0, 0,0) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S)
      
    }else if (N_stat=="matSus"){ #only moms currently Ab+ produce maternally immune pups
      Fmat[1,((j-1)*8+1):(j*8)] <- c(0,0,fert.biwk[j],fert.biwk[j], 0,fert.biwk[j],fert.biwk[j], 0) #the mat immune offpring (so mom = I1 and R1 and I2 and R2 but not N1 or N2)
      Fmat[2,((j-1)*8+1):(j*8)]<- c(fert.biwk[j],fert.biwk[j],0,0, fert.biwk[j],0,0, fert.biwk[j]) #the susceptible offspring (so mom= M, S, and N1, N2)
    }
    
    
    
  }
  
  return(Fmat)
}
buildFMatrix_MSIRNIR <- function(age.classes, adult_fec, surv.biwk, biwk, N_stat){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = adult_fec*surv.biwk*.25
  }else if (biwk==3 | biwk==25){
    new.fec = adult_fec*surv.biwk*.1
    #}else if (biwk==4|biwk ==24){
    #new.fec = adult_fec*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,7*nage,7*nage)
  
  for (j in 1:nage) { #no fertility in first age class (mat immune)
    if(N_stat=="matAB"){
      #try it assuming any mom who has previously seen infection will produce maternally immune pups
      Fmat[1,((j-1)*7+1):(j*7)] <- c(0,0,fert.biwk[j],fert.biwk[j], fert.biwk[j], fert.biwk[j], fert.biwk[j]) #the mat immune offspring (so mom = I1, R1, N1,I2,R2,N2)
      Fmat[2,((j-1)*7+1):(j*7)]<-  c(fert.biwk[j],fert.biwk[j],0,0, 0,0, 0) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S)
      
    }else if (N_stat=="matSus"){ #only moms currently Ab+ produce maternally immune pups
      Fmat[1,((j-1)*7+1):(j*7)] <- c(0,0,fert.biwk[j],fert.biwk[j], 0,fert.biwk[j],fert.biwk[j]) #the mat immune offpring (so mom = I1 and R1 and I2 and R2 but not N1 or N2)
      Fmat[2,((j-1)*7+1):(j*7)]<- c(fert.biwk[j],fert.biwk[j],0,0, fert.biwk[j],0,0) #the susceptible offspring (so mom= M, S, and N1, N2)
    }
    
    
    
  }
  
  return(Fmat)
}
buildFMatrix_MSIRNSIRN <- function(age.classes, adult_fec, surv.biwk, biwk, N_stat){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = adult_fec*surv.biwk*.25
  }else if (biwk==3 | biwk==25){
    new.fec = adult_fec*surv.biwk*.1
    #}else if (biwk==4|biwk ==24){
    #new.fec = adult_fec*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,9*nage,9*nage)
  
  for (j in 1:nage) { #no fertility in first age class (mat immune)
    if(N_stat=="matAB"){
      #try it assuming any mom who has previously seen infection will produce maternally immune pups
      Fmat[1,((j-1)*9+1):(j*9)] <- c(0,0,fert.biwk[j],fert.biwk[j], fert.biwk[j],0, fert.biwk[j], fert.biwk[j], fert.biwk[j]) #the mat immune offspring (so mom = I1, R1, N1,I2,R2,N2)
      Fmat[2,((j-1)*9+1):(j*9)]<-  c(fert.biwk[j],fert.biwk[j],0,0,0, fert.biwk[j], 0, 0,0) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S1 or S2)
      
    }else if (N_stat=="matSus"){ #only moms currently Ab+ produce maternally immune pups
      Fmat[1,((j-1)*9+1):(j*9)] <- c(0,0,fert.biwk[j],fert.biwk[j], 0,0, fert.biwk[j],fert.biwk[j], 0) #the mat immune offpring (so mom = I1 and R1 and I2 and R2 but not N1 or N2)
      Fmat[2,((j-1)*9+1):(j*9)]<- c(fert.biwk[j],fert.biwk[j],0,0, fert.biwk[j], fert.biwk[j],0,0, fert.biwk[j]) #the susceptible offspring (so mom= M, S1,S2, and N1, N2)
    }
    
    
    
  }
  
  return(Fmat)
}
buildFMatrix_MSIRSIR <- function(age.classes, adult_fec, surv.biwk, biwk){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = adult_fec*surv.biwk*.25
  }else if (biwk==3 | biwk==25){
    new.fec = adult_fec*surv.biwk*.1
    #}else if (biwk==4|biwk ==24){
    #new.fec = adult_fec*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,7*nage,7*nage)
  
  for (j in 1:nage) { #no fertility in first age class (mat immune)
    
    #try it assuming any mom who has previously seen infection will produce maternally immune pups
    Fmat[1,((j-1)*7+1):(j*7)] <- c(0,0,fert.biwk[j],fert.biwk[j],0, fert.biwk[j], fert.biwk[j]) #the mat immune offspring (so mom = I1, R1, I2,R2)
    Fmat[2,((j-1)*7+1):(j*7)]<-  c(fert.biwk[j],fert.biwk[j],0,0, fert.biwk[j], 0, 0) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S1 or S2)
    
    
    
  }
  
  return(Fmat)
}

#and sim models
sim.met.MSIR.age <- function(burnin, sim_pop, yrs, ntyr, age.brk, s, beta, sero, recov, mort, mort_juv, adult_fec, wane, slope.wane, sigma, mu.sick, mu.rec, add.inf.mort,add.rec.mort, model){
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  
  #first pull out your data...
  c=4
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  # plot(stab.struct, xlab="Age", ylab="Proportion")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion")
  
  
  #print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #print(i)
    Tmat <- buildTMat_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr), surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr), surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, sigma=(1-exp(-(sigma*(1/ntyr)))), recov=recov, sero=sero, wane=(1-exp(-(wane*(1/ntyr)))), add.rec.mort=add.rec.mort,  add.inf.mort=add.inf.mort, age.rate=(1/ntyr), slope.wane=slope.wane)
    #print(i)
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #print(biwk1)
    
    #feed into fertility matrix:
    Fmat <- buildFMatrix(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=20, c=4)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S", "I", "R")
  dat.tot$N = rowSums((dat.tot[,2:ncol(dat.tot)]))
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$N
  prop.tot$S = prop.tot$S/prop.tot$N
  prop.tot$I = prop.tot$I/prop.tot$N
  prop.tot$R = prop.tot$R/prop.tot$N
  
  
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5])))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S", "I", "R"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S", "I", "R"))
  
  # colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  dat.tot$proportion = dat.tot$count/sim_pop
  
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S", "I", "R"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$times <- round(eq.dat.sero$times,6)
  eq.dat.sero$age[eq.dat.sero$times==50] <- eq.dat.sero$age[eq.dat.sero$times==50]+1
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #remove age at 0
  #eq.sero.tot = subset(eq.sero.tot, age >0)
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  
  list.by.T = dlply(eq.sero.tot, .(times))
  #convert them into a matrix by biweek to keep track for fitting
  out.mat = do.call("rbind", list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  
  return(out.mat)#, seas.prev))
}
sim.met.MSIRN.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, age.brk, recov, recru, mort, mort_juv, adult_fec, wane, slope.wane, sigma, mu.sick, mu.rec, add.inf.mort, add.rec.mort, model, N_stat){
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  
  #first pull out your data...
  c=5
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRN_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr), surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr), surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, sigma=(1-exp(-(sigma*(1/ntyr)))), recru=(1-exp(-(recru*(1/ntyr)))), recov=recov, wane=(1-exp(-(wane*(1/ntyr)))), add.inf.mort=add.inf.mort, add.rec.mort=add.rec.mort, age.rate=(1/ntyr), slope.wane=slope.wane)
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRN(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
    
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S", "I", "R", "N")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S = prop.tot$S/prop.tot$pop
  prop.tot$I = prop.tot$I/prop.tot$pop
  prop.tot$R = prop.tot$R/prop.tot$pop
  prop.tot$N = prop.tot$N/prop.tot$pop
  
 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5]), cbind(dat.tot[,1], dat.tot[,6])))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S", "I", "R", "N"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S", "I", "R", "N"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot[,1],prop.tot[,2]), cbind(prop.tot[,1], prop.tot[,3]), cbind(prop.tot[,1], prop.tot[,4]), cbind(prop.tot[,1], prop.tot[,5]), cbind(prop.tot[,1], prop.tot[,6])))
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S", "I", "R", "N"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S", "I", "R", "N"))
  
  
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S", "I", "R", "N"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  #remove age at 0 - got rid of this bug
  #eq.sero.tot = subset(eq.sero.tot, age >0)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))
  #convert them into a matrix by biweek to keep track for fitting
  out.mat = data.table::rbindlist(list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}
sim.met.MSIRNS.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, age.brk, recov, recru, mort, mort_juv, adult_fec, wane, slope.wane, sigma, mu.sick, mu.rec, add.inf.mort, add.rec.mort, model, N_stat){
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  
  #first pull out your data...
  c=5
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRNS_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr), surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr), surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, sigma=(1-exp(-(sigma*(1/ntyr)))), recru=(1-exp(-(recru*(1/ntyr)))), recov=recov, wane=(1-exp(-(wane*(1/ntyr)))), add.inf.mort=add.inf.mort, add.rec.mort=add.rec.mort, age.rate=(1/ntyr), slope.wane=slope.wane)
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRN(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
    
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S", "I", "R", "N")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  # par(mfrow = c(1,1))
  # with(dat.tot, plot(time, pop, type="l")) #ylim =c(0,1.2*max(N))))
  # #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  # 
  # par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="navy", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # with(dat.tot,
  #      plot(times, pop, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(pop))))
  # 
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S = prop.tot$S/prop.tot$pop
  prop.tot$I = prop.tot$I/prop.tot$pop
  prop.tot$R = prop.tot$R/prop.tot$pop
  prop.tot$N = prop.tot$N/prop.tot$pop
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #     plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # with(prop.tot, lines(times, N, col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #     plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], N[(length(times)-27*2): (length(times))], col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5]), cbind(dat.tot[,1], dat.tot[,6])))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S", "I", "R", "N"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S", "I", "R", "N"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot[,1],prop.tot[,2]), cbind(prop.tot[,1], prop.tot[,3]), cbind(prop.tot[,1], prop.tot[,4]), cbind(prop.tot[,1], prop.tot[,5]), cbind(prop.tot[,1], prop.tot[,6])))
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S", "I", "R", "N"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S", "I", "R", "N"))
  
  
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S", "I", "R", "N"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  #remove age at 0 - got rid of this bug
  #eq.sero.tot = subset(eq.sero.tot, age >0)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))

  #convert them into a matrix by biweek to keep track for fitting
  out.mat = data.table::rbindlist(list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  # some months get doubles biwks
  #  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~month, ncol=3) 
  # # # 
  # #  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  # #    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~biwk, ncol=4) 
  # # 
  # #names(eq.sero.tot)[1] ="biwk"
  #plot one of them to get a sense
  #plot.age.seroprev = list.by.T[[5]]
  
  #age.sero.plot(plot.age.seroprev)
  
  #if you want to, plot all of them and save
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}
sim.met.MSIRNNI.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, age.brk, recov, recru, mort, mort_juv, adult_fec, wane, slope.wane, sigma, epsilon, mu.sick, mu.rec, add.inf.mort, add.rec.mort, model, N_stat){
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  #first pull out your data...
  c=6
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  N2_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init, N2_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRNNI_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr), surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr), surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, sigma=(1-exp(-(sigma*(1/ntyr)))), recru=(1-exp(-(recru*(1/ntyr)))), recov=recov, wane=(1-exp(-(wane*(1/ntyr)))), add.inf.mort=add.inf.mort, add.rec.mort=add.rec.mort, age.rate=(1/ntyr), slope.wane=slope.wane, epsilon = (1-exp(-(epsilon*(1/ntyr)))))
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRNN(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
    
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S", "I", "R", "N1", "N2")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  # par(mfrow = c(1,1))
  # with(dat.tot, plot(time, pop, type="l")) #ylim =c(0,1.2*max(N))))
  # #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  # 
  # par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="navy", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # with(dat.tot,
  #      plot(times, pop, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(pop))))
  # 
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S = prop.tot$S/prop.tot$pop
  prop.tot$I = prop.tot$I/prop.tot$pop
  prop.tot$R = prop.tot$R/prop.tot$pop
  prop.tot$N1 = prop.tot$N1/prop.tot$pop
  prop.tot$N2 = prop.tot$N2/prop.tot$pop
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #     plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # with(prop.tot, lines(times, N, col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #     plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], N[(length(times)-27*2): (length(times))], col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5]), cbind(dat.tot[,1], dat.tot[,6]), cbind(dat.tot[,1], dat.tot[,7])))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S", "I", "R", "N1", "N2"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S", "I", "R", "N1", "N2"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N1' = "blue", 'N2'="navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot[,1],prop.tot[,2]), cbind(prop.tot[,1], prop.tot[,3]), cbind(prop.tot[,1], prop.tot[,4]), cbind(prop.tot[,1], prop.tot[,5]), cbind(prop.tot[,1], prop.tot[,6]), cbind(prop.tot[,1], prop.tot[,7])))
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S", "I", "R", "N1", "N2"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S", "I", "R", "N1", "N2"))
  
  
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S", "I", "R", "N1", "N2"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  #remove age at 0 - got rid of this bug
  #eq.sero.tot = subset(eq.sero.tot, age >0)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))

  #convert them into a matrix by biweek to keep track for fitting
  out.mat = data.table::rbindlist(list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  # some months get doubles biwks
  #  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~month, ncol=3) 
  # # # # 
  #   ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #     geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~biwk, ncol=4) 
  # # # 
  # # #names(eq.sero.tot)[1] ="biwk"
  #plot one of them to get a sense
  #plot.age.seroprev = list.by.T[[5]]
  
  #age.sero.plot(plot.age.seroprev)
  
  #if you want to, plot all of them and save
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}
sim.met.MSIRNNS.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, age.brk, recov, recru, mort, mort_juv, adult_fec, wane, slope.wane, sigma, epsilon, mu.sick, mu.rec, add.inf.mort, add.rec.mort, model, N_stat){
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  #first pull out your data...
  c=6
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  N2_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init, N2_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRNNS_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr), surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr), surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, sigma=(1-exp(-(sigma*(1/ntyr)))), recru=(1-exp(-(recru*(1/ntyr)))), recov=recov, wane=(1-exp(-(wane*(1/ntyr)))), add.inf.mort=add.inf.mort, add.rec.mort=add.rec.mort, age.rate=(1/ntyr), slope.wane=slope.wane, epsilon = (1-exp(-(epsilon*(1/ntyr)))))
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRNN(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
    
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S", "I", "R", "N1", "N2")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  # par(mfrow = c(1,1))
  # with(dat.tot, plot(time, pop, type="l")) #ylim =c(0,1.2*max(N))))
  # #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  # 
  # par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="navy", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # with(dat.tot,
  #      plot(times, pop, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(pop))))
  # 
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S = prop.tot$S/prop.tot$pop
  prop.tot$I = prop.tot$I/prop.tot$pop
  prop.tot$R = prop.tot$R/prop.tot$pop
  prop.tot$N1 = prop.tot$N1/prop.tot$pop
  prop.tot$N2 = prop.tot$N2/prop.tot$pop
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #     plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # with(prop.tot, lines(times, N, col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #     plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], N[(length(times)-27*2): (length(times))], col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5]), cbind(dat.tot[,1], dat.tot[,6]), cbind(dat.tot[,1], dat.tot[,7])))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S", "I", "R", "N1", "N2"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S", "I", "R", "N1", "N2"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N1' = "blue", 'N2'="navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot[,1],prop.tot[,2]), cbind(prop.tot[,1], prop.tot[,3]), cbind(prop.tot[,1], prop.tot[,4]), cbind(prop.tot[,1], prop.tot[,5]), cbind(prop.tot[,1], prop.tot[,6]), cbind(prop.tot[,1], prop.tot[,7])))
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S", "I", "R", "N1", "N2"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S", "I", "R", "N1", "N2"))
  
  
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S", "I", "R", "N1", "N2"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  #remove age at 0 - got rid of this bug
  #eq.sero.tot = subset(eq.sero.tot, age >0)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))

  #convert them into a matrix by biweek to keep track for fitting
  out.mat = data.table::rbindlist(list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  # some months get doubles biwks
  #  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~month, ncol=3) 
  # # # # 
  #   ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #     geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~biwk, ncol=4) 
  # # # 
  # # #names(eq.sero.tot)[1] ="biwk"
  #plot one of them to get a sense
  #plot.age.seroprev = list.by.T[[5]]
  
  #age.sero.plot(plot.age.seroprev)
  
  #if you want to, plot all of them and save
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}
sim.met.MSIRSIR.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, beta2, age.brk, recov,  mort, mort_juv, adult_fec, wane, slope.wane, sigma, mu.sick, mu.rec, add.inf.mort, add.rec.mort, model){
  
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  #first pull out your data...
  c=7 #state variables
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  S2_init = rep(0, s)
  I2_init = rep(0, s)
  R2_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, S2_init, I2_init, R2_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRSIR_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr), surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr), surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, beta2=beta2, sigma=(1-exp(-(sigma*(1/ntyr)))), recov=recov, wane=(1-exp(-(wane*(1/ntyr)))), add.inf.mort=add.inf.mort, add.rec.mort=add.rec.mort, age.rate=(1/ntyr), slope.wane=slope.wane)
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRSIR(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
    
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S1", "I1", "R1", "S2",  "I2", "R2")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  # par(mfrow = c(1,1))
  # with(dat.tot, plot(time, pop, type="l")) #ylim =c(0,1.2*max(N))))
  # #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  # 
  # par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="navy", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # with(dat.tot,
  #      plot(times, pop, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(pop))))
  # 
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S1 = prop.tot$S1/prop.tot$pop
  prop.tot$I1 = prop.tot$I1/prop.tot$pop
  prop.tot$R1 = prop.tot$R1/prop.tot$pop
  prop.tot$S2 = prop.tot$S2/prop.tot$pop
  prop.tot$I2 = prop.tot$I2/prop.tot$pop
  prop.tot$R2 = prop.tot$R2/prop.tot$pop
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #     plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # with(prop.tot, lines(times, N, col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #     plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], N[(length(times)-27*2): (length(times))], col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot$time,dat.tot$M), cbind(dat.tot$time,dat.tot$S1), cbind(dat.tot$time,dat.tot$I1), cbind(dat.tot$time,dat.tot$R1), cbind(dat.tot$time,dat.tot$S2), cbind(dat.tot$time,dat.tot$I2), cbind(dat.tot$time,dat.tot$R2)))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S1", "I1", "R1", "S2", "I2", "R2"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S1", "I1", "R1","S2", "I2", "R2"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N1' = "blue", 'N2'="navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot$time,prop.tot$M), cbind(prop.tot$time,prop.tot$S1), cbind(prop.tot$time,prop.tot$I1), cbind(prop.tot$time,prop.tot$R1), cbind(prop.tot$time,prop.tot$S2), cbind(prop.tot$time,prop.tot$I2), cbind(prop.tot$time,prop.tot$R2)))
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S1", "I1", "R1", "S2","I2", "R2"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S1", "I1", "R1", "S2","I2", "R2"))
  
  
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S1", "I1", "R1", "S2","I2", "R2"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R2" |eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R2" |eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  #remove age at 0 - got rid of this bug
  #eq.sero.tot = subset(eq.sero.tot, age >0)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))
  
  #convert them into a matrix by biweek to keep track for fitting
  out.mat = data.table::rbindlist(list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  # some months get doubles biwks
  #  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~month, ncol=3) 
  # # # # # # 
  #     ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #       geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~biwk, ncol=4) 
  # # # # # 
  # # # # #names(eq.sero.tot)[1] ="biwk"
  # #plot one of them to get a sense
  # #plot.age.seroprev = list.by.T[[5]]
  
  #age.sero.plot(plot.age.seroprev)
  
  #if you want to, plot all of them and save
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}
sim.met.MSIRNIRN.LT.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, age.brk, recov, recru, mort, mort_juv, adult_fec, wane, slope.wane, sigma, sigma2, mu.sick, mu.rec, epsilon, add.inf.mort, add.rec.mort, model, N_stat){
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  #first pull out your data...
  c=8 #state variables
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  N2_init = rep(0, s)
  I2_init = rep(0, s)
  R2_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init, I2_init, R2_init, N2_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRNIRN_LT_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr), surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr), surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, sigma=(1-exp(-(sigma*(1/ntyr)))), recru=(1-exp(-(recru*(1/ntyr)))), recov=recov, wane=(1-exp(-(wane*(1/ntyr)))), epsilon=(1-exp(-(epsilon*(1/ntyr)))),  add.inf.mort=add.inf.mort, add.rec.mort=add.rec.mort, age.rate=(1/ntyr), slope.wane=slope.wane, sigma2 = (1-exp(-(sigma2*(1/ntyr)))))
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRNIRN(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
    
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S1", "I1", "R1", "N1", "I2", "R2", "N2")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  # par(mfrow = c(1,1))
  # with(dat.tot, plot(time, pop, type="l")) #ylim =c(0,1.2*max(N))))
  # #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  # 
  # par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="navy", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # with(dat.tot,
  #      plot(times, pop, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(pop))))
  # 
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S1 = prop.tot$S1/prop.tot$pop
  prop.tot$I1 = prop.tot$I1/prop.tot$pop
  prop.tot$R1 = prop.tot$R1/prop.tot$pop
  prop.tot$N1 = prop.tot$N1/prop.tot$pop
  prop.tot$I2 = prop.tot$I2/prop.tot$pop
  prop.tot$R2 = prop.tot$R2/prop.tot$pop
  prop.tot$N2 = prop.tot$N2/prop.tot$pop
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #     plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # with(prop.tot, lines(times, N, col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #     plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], N[(length(times)-27*2): (length(times))], col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot$time,dat.tot$M), cbind(dat.tot$time,dat.tot$S1), cbind(dat.tot$time,dat.tot$I1), cbind(dat.tot$time,dat.tot$R1), cbind(dat.tot$time,dat.tot$N1), cbind(dat.tot$time,dat.tot$I2), cbind(dat.tot$time,dat.tot$R2), cbind(dat.tot$time,dat.tot$N2)))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S1", "I1", "R1", "N1","I2", "R2", "N2"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S1", "I1", "R1", "N1","I2", "R2", "N2"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N1' = "blue", 'N2'="navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot$time,prop.tot$M), cbind(prop.tot$time,prop.tot$S1), cbind(prop.tot$time,prop.tot$I1), cbind(prop.tot$time,prop.tot$R1), cbind(prop.tot$time,prop.tot$N1), cbind(prop.tot$time,prop.tot$I2), cbind(prop.tot$time,prop.tot$R2), cbind(prop.tot$time,prop.tot$N2)))
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S1", "I1", "R1", "N1","I2", "R2", "N2"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S1", "I1", "R1", "N1","I2", "R2", "N2"))
  
  
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S1", "I1", "R1", "N1","I2", "R2", "N2"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R2" | eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R2" | eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  #remove age at 0 - got rid of this bug
  #eq.sero.tot = subset(eq.sero.tot, age >0)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))
  
  #convert them into a matrix by biweek to keep track for fitting
  out.mat = data.table::rbindlist(list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  # some months get doubles biwks
  #  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~month, ncol=3) 
  # # # # # # 
  #     ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #       geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~biwk, ncol=4) 
  # # # # # 
  # # # #names(eq.sero.tot)[1] ="biwk"
  # #plot one of them to get a sense
  # #plot.age.seroprev = list.by.T[[5]]
  
  #age.sero.plot(plot.age.seroprev)
  
  #if you want to, plot all of them and save
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}
sim.met.MSIRNIR.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, age.brk, recov, recru, mort, mort_juv, adult_fec, wane, slope.wane, sigma, sigma2, mu.sick, mu.rec, add.inf.mort, add.rec.mort, model, N_stat){
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  #first pull out your data...
  c=7 #state variables
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  I2_init = rep(0, s)
  R2_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init, I2_init, R2_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRNIR_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr), surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr), surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, sigma=(1-exp(-(sigma*(1/ntyr)))), recru=(1-exp(-(recru*(1/ntyr)))), recov=recov, wane=(1-exp(-(wane*(1/ntyr)))), add.inf.mort=add.inf.mort, add.rec.mort=add.rec.mort, age.rate=(1/ntyr), slope.wane=slope.wane, sigma2 = (1-exp(-(sigma2*(1/ntyr)))))
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRNIR(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
    
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S1", "I1", "R1", "N1", "I2", "R2")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  # par(mfrow = c(1,1))
  # with(dat.tot, plot(time, pop, type="l")) #ylim =c(0,1.2*max(N))))
  # #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  # 
  # par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="navy", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # with(dat.tot,
  #      plot(times, pop, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(pop))))
  # 
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S1 = prop.tot$S1/prop.tot$pop
  prop.tot$I1 = prop.tot$I1/prop.tot$pop
  prop.tot$R1 = prop.tot$R1/prop.tot$pop
  prop.tot$N1 = prop.tot$N1/prop.tot$pop
  prop.tot$I2 = prop.tot$I2/prop.tot$pop
  prop.tot$R2 = prop.tot$R2/prop.tot$pop
  
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #     plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # with(prop.tot, lines(times, N, col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #     plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], N[(length(times)-27*2): (length(times))], col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot$time,dat.tot$M), cbind(dat.tot$time,dat.tot$S1), cbind(dat.tot$time,dat.tot$I1), cbind(dat.tot$time,dat.tot$R1), cbind(dat.tot$time,dat.tot$N1), cbind(dat.tot$time,dat.tot$I2), cbind(dat.tot$time,dat.tot$R2)))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S1", "I1", "R1", "N1","I2", "R2"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S1", "I1", "R1", "N1","I2", "R2"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N1' = "blue", 'N2'="navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot$time,prop.tot$M), cbind(prop.tot$time,prop.tot$S1), cbind(prop.tot$time,prop.tot$I1), cbind(prop.tot$time,prop.tot$R1), cbind(prop.tot$time,prop.tot$N1), cbind(prop.tot$time,prop.tot$I2), cbind(prop.tot$time,prop.tot$R2)))
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S1", "I1", "R1", "N1","I2", "R2"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S1", "I1", "R1", "N1","I2", "R2"))
  
  
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S1", "I1", "R1", "N1","I2", "R2"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R2" | eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R2" | eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  #remove age at 0 - got rid of this bug
  #eq.sero.tot = subset(eq.sero.tot, age >0)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))

  #convert them into a matrix by biweek to keep track for fitting
  out.mat = data.table::rbindlist(list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  # some months get doubles biwks
  #  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~month, ncol=3) 
  # # # # # # 
  #     ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #       geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~biwk, ncol=4) 
  # # # # # 
  # # # #names(eq.sero.tot)[1] ="biwk"
  # #plot one of them to get a sense
  # #plot.age.seroprev = list.by.T[[5]]
  
  #age.sero.plot(plot.age.seroprev)
  
  #if you want to, plot all of them and save
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}
sim.met.MSIRNSIRN.LT.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, beta2, age.brk, recov, recru, mort, mort_juv, adult_fec, wane, slope.wane, sigma, sigma2, mu.sick, mu.rec, epsilon, add.inf.mort, add.rec.mort, model, N_stat){
  
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  
  #first pull out your data...
  c=9 #state variables
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  N2_init = rep(0, s)
  S2_init = rep(0, s)
  I2_init = rep(0, s)
  R2_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init, S2_init, I2_init, R2_init, N2_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRNSIRN_LT_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr), surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr), surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, beta2=beta2, sigma=(1-exp(-(sigma*(1/ntyr)))), recru=(1-exp(-(recru*(1/ntyr)))), recov=recov, wane=(1-exp(-(wane*(1/ntyr)))), epsilon=(1-exp(-(epsilon*(1/ntyr)))), add.inf.mort=add.inf.mort, add.rec.mort=add.rec.mort, age.rate=(1/ntyr), slope.wane=slope.wane, sigma2 = (1-exp(-(sigma2*(1/ntyr)))))
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRNSIRN(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
    
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S1", "I1", "R1", "N1","S2",  "I2", "R2", "N2")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  # par(mfrow = c(1,1))
  # with(dat.tot, plot(time, pop, type="l")) #ylim =c(0,1.2*max(N))))
  # #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  # 
  # par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="navy", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # with(dat.tot,
  #      plot(times, pop, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(pop))))
  # 
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S1 = prop.tot$S1/prop.tot$pop
  prop.tot$I1 = prop.tot$I1/prop.tot$pop
  prop.tot$R1 = prop.tot$R1/prop.tot$pop
  prop.tot$N1 = prop.tot$N1/prop.tot$pop
  prop.tot$S2 = prop.tot$S2/prop.tot$pop
  prop.tot$I2 = prop.tot$I2/prop.tot$pop
  prop.tot$R2 = prop.tot$R2/prop.tot$pop
  prop.tot$N2 = prop.tot$N2/prop.tot$pop
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #     plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # with(prop.tot, lines(times, N, col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #     plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], N[(length(times)-27*2): (length(times))], col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot$time,dat.tot$M), cbind(dat.tot$time,dat.tot$S1), cbind(dat.tot$time,dat.tot$I1), cbind(dat.tot$time,dat.tot$R1), cbind(dat.tot$time,dat.tot$N1), cbind(dat.tot$time,dat.tot$S2), cbind(dat.tot$time,dat.tot$I2), cbind(dat.tot$time,dat.tot$R2), cbind(dat.tot$time,dat.tot$N2)))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S1", "I1", "R1", "N1","S2", "I2", "R2", "N2"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S1", "I1", "R1", "N1","S2", "I2", "R2", "N2"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N1' = "blue", 'N2'="navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot$time,prop.tot$M), cbind(prop.tot$time,prop.tot$S1), cbind(prop.tot$time,prop.tot$I1), cbind(prop.tot$time,prop.tot$R1), cbind(prop.tot$time,prop.tot$N1), cbind(prop.tot$time,prop.tot$S2), cbind(prop.tot$time,prop.tot$I2), cbind(prop.tot$time,prop.tot$R2), cbind(prop.tot$time,prop.tot$N2)))
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S1", "I1", "R1", "N1","S2","I2", "R2", "N2"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S1", "I1", "R1", "N1","S2","I2", "R2", "N2"))
  
  
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S1", "I1", "R1", "N1","S2","I2", "R2", "N2"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R2" |eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R2" |eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  #remove age at 0 - got rid of this bug
  #eq.sero.tot = subset(eq.sero.tot, age >0)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))

  #convert them into a matrix by biweek to keep track for fitting
  out.mat = data.table::rbindlist(list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  # some months get doubles biwks
  #  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~month, ncol=3) 
  # # # # # # 
  #     ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #       geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~biwk, ncol=4) 
  # # # # # 
  # # # # #names(eq.sero.tot)[1] ="biwk"
  # #plot one of them to get a sense
  # #plot.age.seroprev = list.by.T[[5]]
  
  #age.sero.plot(plot.age.seroprev)
  
  #if you want to, plot all of them and save
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}
sim.met.MSIRNIRN.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, age.brk, recov, recru, mort, mort_juv, adult_fec, wane, slope.wane, sigma, sigma2, mu.sick, mu.rec, add.inf.mort, add.rec.mort, model, N_stat){
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  #first pull out your data...
  c=8 #state variables
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  N2_init = rep(0, s)
  I2_init = rep(0, s)
  R2_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init, I2_init, R2_init, N2_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRNIRN_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr), surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr), surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, sigma=(1-exp(-(sigma*(1/ntyr)))), recru=(1-exp(-(recru*(1/ntyr)))), recov=recov, wane=(1-exp(-(wane*(1/ntyr)))), add.inf.mort=add.inf.mort, add.rec.mort=add.rec.mort, age.rate=(1/ntyr), slope.wane=slope.wane, sigma2 = (1-exp(-(sigma2*(1/ntyr)))))
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRNIRN(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
    
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S1", "I1", "R1", "N1", "I2", "R2", "N2")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  # par(mfrow = c(1,1))
  # with(dat.tot, plot(time, pop, type="l")) #ylim =c(0,1.2*max(N))))
  # #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  # 
  # par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="navy", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # with(dat.tot,
  #      plot(times, pop, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(pop))))
  # 
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S1 = prop.tot$S1/prop.tot$pop
  prop.tot$I1 = prop.tot$I1/prop.tot$pop
  prop.tot$R1 = prop.tot$R1/prop.tot$pop
  prop.tot$N1 = prop.tot$N1/prop.tot$pop
  prop.tot$I2 = prop.tot$I2/prop.tot$pop
  prop.tot$R2 = prop.tot$R2/prop.tot$pop
  prop.tot$N2 = prop.tot$N2/prop.tot$pop
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #     plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # with(prop.tot, lines(times, N, col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #     plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], N[(length(times)-27*2): (length(times))], col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot$time,dat.tot$M), cbind(dat.tot$time,dat.tot$S1), cbind(dat.tot$time,dat.tot$I1), cbind(dat.tot$time,dat.tot$R1), cbind(dat.tot$time,dat.tot$N1), cbind(dat.tot$time,dat.tot$I2), cbind(dat.tot$time,dat.tot$R2), cbind(dat.tot$time,dat.tot$N2)))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S1", "I1", "R1", "N1","I2", "R2", "N2"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S1", "I1", "R1", "N1","I2", "R2", "N2"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N1' = "blue", 'N2'="navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot$time,prop.tot$M), cbind(prop.tot$time,prop.tot$S1), cbind(prop.tot$time,prop.tot$I1), cbind(prop.tot$time,prop.tot$R1), cbind(prop.tot$time,prop.tot$N1), cbind(prop.tot$time,prop.tot$I2), cbind(prop.tot$time,prop.tot$R2), cbind(prop.tot$time,prop.tot$N2)))
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S1", "I1", "R1", "N1","I2", "R2", "N2"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S1", "I1", "R1", "N1","I2", "R2", "N2"))
  
  
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S1", "I1", "R1", "N1","I2", "R2", "N2"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R2" | eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R2" | eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  #remove age at 0 - got rid of this bug
  #eq.sero.tot = subset(eq.sero.tot, age >0)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))

  #convert them into a matrix by biweek to keep track for fitting
  out.mat = data.table::rbindlist(list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  # some months get doubles biwks
  #  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~month, ncol=3) 
  # # # # # # 
  #     ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #       geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~biwk, ncol=4) 
  # # # # # 
  # # # #names(eq.sero.tot)[1] ="biwk"
  # #plot one of them to get a sense
  # #plot.age.seroprev = list.by.T[[5]]
  
  #age.sero.plot(plot.age.seroprev)
  
  #if you want to, plot all of them and save
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}
sim.met.MSIRNSIRN.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, beta2, age.brk, recov, recru, mort, mort_juv, adult_fec, wane, slope.wane, sigma, sigma2, mu.sick, mu.rec, add.inf.mort, add.rec.mort, model, N_stat){
  
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  
  #first pull out your data...
  c=9 #state variables
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  N2_init = rep(0, s)
  S2_init = rep(0, s)
  I2_init = rep(0, s)
  R2_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init, S2_init, I2_init, R2_init, N2_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRNSIRN_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr), surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr), surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, beta2=beta2, sigma=(1-exp(-(sigma*(1/ntyr)))), recru=(1-exp(-(recru*(1/ntyr)))), recov=recov, wane=(1-exp(-(wane*(1/ntyr)))), add.inf.mort=add.inf.mort, add.rec.mort=add.rec.mort, age.rate=(1/ntyr), slope.wane=slope.wane, sigma2 = (1-exp(-(sigma2*(1/ntyr)))))
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRNSIRN(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
    
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S1", "I1", "R1", "N1","S2",  "I2", "R2", "N2")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  # par(mfrow = c(1,1))
  # with(dat.tot, plot(time, pop, type="l")) #ylim =c(0,1.2*max(N))))
  # #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  # 
  # par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="navy", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # with(dat.tot,
  #      plot(times, pop, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(pop))))
  # 
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S1 = prop.tot$S1/prop.tot$pop
  prop.tot$I1 = prop.tot$I1/prop.tot$pop
  prop.tot$R1 = prop.tot$R1/prop.tot$pop
  prop.tot$N1 = prop.tot$N1/prop.tot$pop
  prop.tot$S2 = prop.tot$S2/prop.tot$pop
  prop.tot$I2 = prop.tot$I2/prop.tot$pop
  prop.tot$R2 = prop.tot$R2/prop.tot$pop
  prop.tot$N2 = prop.tot$N2/prop.tot$pop
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #     plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # with(prop.tot, lines(times, N, col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #     plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], N[(length(times)-27*2): (length(times))], col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot$time,dat.tot$M), cbind(dat.tot$time,dat.tot$S1), cbind(dat.tot$time,dat.tot$I1), cbind(dat.tot$time,dat.tot$R1), cbind(dat.tot$time,dat.tot$N1), cbind(dat.tot$time,dat.tot$S2), cbind(dat.tot$time,dat.tot$I2), cbind(dat.tot$time,dat.tot$R2), cbind(dat.tot$time,dat.tot$N2)))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S1", "I1", "R1", "N1","S2", "I2", "R2", "N2"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S1", "I1", "R1", "N1","S2", "I2", "R2", "N2"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N1' = "blue", 'N2'="navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot$time,prop.tot$M), cbind(prop.tot$time,prop.tot$S1), cbind(prop.tot$time,prop.tot$I1), cbind(prop.tot$time,prop.tot$R1), cbind(prop.tot$time,prop.tot$N1), cbind(prop.tot$time,prop.tot$S2), cbind(prop.tot$time,prop.tot$I2), cbind(prop.tot$time,prop.tot$R2), cbind(prop.tot$time,prop.tot$N2)))
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S1", "I1", "R1", "N1","S2","I2", "R2", "N2"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S1", "I1", "R1", "N1","S2","I2", "R2", "N2"))
  
  
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S1", "I1", "R1", "N1","S2","I2", "R2", "N2"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R2" |eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R2" |eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  #remove age at 0 - got rid of this bug
  #eq.sero.tot = subset(eq.sero.tot, age >0)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))
  
  #convert them into a matrix by biweek to keep track for fitting
  out.mat = data.table::rbindlist(list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  # some months get doubles biwks
  #  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~month, ncol=3) 
  # # # # # # 
  #     ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #       geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~biwk, ncol=4) 
  # # # # # 
  # # # # #names(eq.sero.tot)[1] ="biwk"
  # #plot one of them to get a sense
  # #plot.age.seroprev = list.by.T[[5]]
  
  #age.sero.plot(plot.age.seroprev)
  
  #if you want to, plot all of them and save
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}
sim.met.MSIRNSIRNS.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, beta2, age.brk, recov, recru, recru2, mort, mort_juv, adult_fec, wane, slope.wane, sigma, sigma2, mu.sick, mu.rec, add.inf.mort, add.rec.mort, model, N_stat){
  
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  
  #first pull out your data...
  c=9 #state variables
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  N2_init = rep(0, s)
  S2_init = rep(0, s)
  I2_init = rep(0, s)
  R2_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init, S2_init, I2_init, R2_init, N2_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRNSIRNS_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr), surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr), surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, beta2=beta2, sigma=(1-exp(-(sigma*(1/ntyr)))), recru=(1-exp(-(recru*(1/ntyr)))), recru2=(1-exp(-(recru2*(1/ntyr)))), recov=recov, wane=(1-exp(-(wane*(1/ntyr)))), add.inf.mort=add.inf.mort, add.rec.mort=add.rec.mort, age.rate=(1/ntyr), slope.wane=slope.wane, sigma2 = (1-exp(-(sigma2*(1/ntyr)))))
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRNSIRN(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
    
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S1", "I1", "R1", "N1","S2",  "I2", "R2", "N2")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  # par(mfrow = c(1,1))
  # with(dat.tot, plot(time, pop, type="l")) #ylim =c(0,1.2*max(N))))
  # #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  # 
  # par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="navy", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # with(dat.tot,
  #      plot(times, pop, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(pop))))
  # 
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S1 = prop.tot$S1/prop.tot$pop
  prop.tot$I1 = prop.tot$I1/prop.tot$pop
  prop.tot$R1 = prop.tot$R1/prop.tot$pop
  prop.tot$N1 = prop.tot$N1/prop.tot$pop
  prop.tot$S2 = prop.tot$S2/prop.tot$pop
  prop.tot$I2 = prop.tot$I2/prop.tot$pop
  prop.tot$R2 = prop.tot$R2/prop.tot$pop
  prop.tot$N2 = prop.tot$N2/prop.tot$pop
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #     plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # with(prop.tot, lines(times, N, col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #     plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], N[(length(times)-27*2): (length(times))], col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot$time,dat.tot$M), cbind(dat.tot$time,dat.tot$S1), cbind(dat.tot$time,dat.tot$I1), cbind(dat.tot$time,dat.tot$R1), cbind(dat.tot$time,dat.tot$N1), cbind(dat.tot$time,dat.tot$S2), cbind(dat.tot$time,dat.tot$I2), cbind(dat.tot$time,dat.tot$R2), cbind(dat.tot$time,dat.tot$N2)))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S1", "I1", "R1", "N1","S2", "I2", "R2", "N2"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S1", "I1", "R1", "N1","S2", "I2", "R2", "N2"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N1' = "blue", 'N2'="navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot$time,prop.tot$M), cbind(prop.tot$time,prop.tot$S1), cbind(prop.tot$time,prop.tot$I1), cbind(prop.tot$time,prop.tot$R1), cbind(prop.tot$time,prop.tot$N1), cbind(prop.tot$time,prop.tot$S2), cbind(prop.tot$time,prop.tot$I2), cbind(prop.tot$time,prop.tot$R2), cbind(prop.tot$time,prop.tot$N2)))
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S1", "I1", "R1", "N1","S2","I2", "R2", "N2"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S1", "I1", "R1", "N1","S2","I2", "R2", "N2"))
  
  
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S1", "I1", "R1", "N1","S2","I2", "R2", "N2"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R2" |eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R2" |eq.dat.sero$class=="R1" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  #remove age at 0 - got rid of this bug
  #eq.sero.tot = subset(eq.sero.tot, age >0)
  #ggplot(eq.sero.tot) + geom_line(aes(x=age,y=seroprevalence)) + facet_wrap(~times, ncol=4)
  #length(unique(eq.dat.sero$times)) #27
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))
 
  #convert them into a matrix by biweek to keep track for fitting
  out.mat = data.table::rbindlist(list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  # some months get doubles biwks
  #  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~month, ncol=3) 
  # # # # # # 
  #     ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
  #       geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~biwk, ncol=4) 
  # # # # # 
  # # # # #names(eq.sero.tot)[1] ="biwk"
  # #plot one of them to get a sense
  # #plot.age.seroprev = list.by.T[[5]]
  
  #age.sero.plot(plot.age.seroprev)
  
  #if you want to, plot all of them and save
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}



#likelihood function
#this for a biweekly observation process
#this for a biweekly observation process
log.lik.age <- function(par, data, cutoff, burnin, sim_pop, yrs, ntyr, age.brk, s, par.df1, model, N_stat){
  
  if (par.df1$mu_sick>1 & par.df1$mu_rec==1){
    add.inf.mort <- TRUE
    add.rec.mort <- FALSE
    mu.sick = exp(par[length(par)])
    par <- par[1:(length(par)-1)]
    mu.rec = 1
    
    #mu sick can't be <1 or much >4 because you'd have 0 survivorship
    if(mu.sick<1 | mu.sick > 4){
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
    }
    
  }else if (par.df1$mu_sick ==1 & par.df1$mu_rec>1){
    add.rec.mort <- TRUE
    add.inf.mort <- FALSE
    mu.rec = exp(par[length(par)])
    par <- par[1:(length(par)-1)]
    mu.sick = 1
    
    #mu rec can't be <1 and can't be >4 because you'd have 0 survivorship
    if(mu.rec<1 | mu.rec>4){
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
    }
    
  }else if (par.df1$mu_sick > 1 & par.df1$mu_rec>1){
    add.inf.mort <- TRUE
    add.rec.mort <- TRUE
    
    mu.rec = exp(par[length(par)])# last
    par <- par[1:(length(par)-1)]
    mu.sick = exp(par[length(par)]) #second to last
    par <- par[1:(length(par)-1)]
    
    #mu rec and mu sick can't be <1
    if(mu.rec<1 | mu.sick <1 | mu.rec>4 | mu.sick>4){
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
    }
    
  }else{
    
    add.inf.mort <- FALSE
    add.rec.mort <- FALSE
    mu.rec = 1
    mu.sick = 1
  }
  
  # and a few additional constraints within each model formula:
  # rates of Ab waning (sigma), maternal immune waning (wane), 
  # and recrudescence (recru) can't be faster than a duration
  # of one biweek of the year (1/.0385) = 25.97 = 26
  # no constraints on beta
  
  
  
  #run your model of choice
  if(par.df1$model=="MSIR"){
    
    #  par = c(beta)  
    
    out = sim.met.MSIR.age(burnin=burnin,
                           sim_pop=sim_pop,
                           yrs=yrs,
                           ntyr=ntyr,
                           s=s,
                           beta=exp(par[1:length(par)]),
                           sero=par.df1$sero,
                           recov=par.df1$recov,
                           age.brk=age.brk,
                           mort=par.df1$mort, #annual
                           mort_juv=par.df1$mort_juv, #annual
                           adult_fec=par.df1$adult_fec, #annual
                           wane=par.df1$wane,
                           slope.wane=slope.wane,
                           sigma=par.df1$sigma,#biwks
                           model=par.df1$model,
                           mu.sick = mu.sick,
                           mu.rec = mu.rec,
                           add.inf.mort = add.inf.mort,
                           add.rec.mort = add.rec.mort)
    
  } else if (par.df1$model =="MSIRS"){
    
    #  par = c(beta, sigma)  
    
    if ((exp(par[length(par)])>26) | (exp(par[length(par)])<1)){ #this is sigma - duration adult Abs between 1 biwk and 1 year
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
      
    }
    
    out = sim.met.MSIR.age(burnin=burnin,
                           sim_pop=sim_pop,
                           yrs=yrs,
                           ntyr=ntyr,
                           s=s,
                           beta=exp(par[1:(length(par)-1)]),
                           sero=par.df1$sero,
                           recov=par.df1$recov,
                           age.brk=age.brk,
                           mort=par.df1$mort, #annual
                           mort_juv=par.df1$mort_juv, #annual
                           adult_fec=par.df1$adult_fec, #annual
                           wane=par.df1$wane,
                           slope.wane=slope.wane,
                           sigma=exp(par[length(par)]),
                           model=par.df1$model,
                           mu.sick = mu.sick,
                           mu.rec = mu.rec,
                           add.inf.mort = add.inf.mort,
                           add.rec.mort = add.rec.mort)
    
  }else if (par.df1$model =="MSIRN"){
    
    #  par = c(beta, sigma)  
    
    if ((exp(par[length(par)])>26) | (exp(par[length(par)])<1)){ #this is sigma - duration adult Abs between 1 biwk and 1 year
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
      
    }
    
    out = sim.met.MSIRN.age(burnin=burnin,
                            sim_pop=sim_pop,
                            yrs=yrs,
                            ntyr=ntyr,
                            s=s,
                            beta=exp(par[1:(length(par)-1)]),
                            recov=par.df1$recov,
                            mort=par.df1$mort, #annual
                            mort_juv=par.df1$mort_juv, #annual
                            adult_fec=par.df1$adult_fec, #annual
                            wane=par.df1$wane,
                            age.brk=par.df1$age.brk,
                            slope.wane=slope.wane,
                            sigma=exp(par[length(par)]),
                            recru=par.df1$recru,
                            model=par.df1$model,
                            mu.sick = mu.sick,
                            mu.rec = mu.rec,
                            add.inf.mort = add.inf.mort,
                            add.rec.mort = add.rec.mort,
                            N_stat=N_stat)
  }else if (par.df1$model =="MSIRNI"){
    
    #  par = c(beta, sigma, recru)
    
    if ((exp(par[(length(par)-1)])>26) | (exp(par[(length(par)-1)])<1) | (exp(par[length(par)])>26)) {
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
      
    }
    
    out = sim.met.MSIRN.age(burnin=burnin,
                            sim_pop=sim_pop,
                            yrs=yrs,
                            ntyr=ntyr,
                            s=s,
                            beta=exp(par[1:(length(par)-2)]),
                            recov=par.df1$recov,
                            mort=par.df1$mort, #annual
                            mort_juv=par.df1$mort_juv, #annual
                            adult_fec=par.df1$adult_fec, #annual
                            wane=par.df1$wane,
                            age.brk=age.brk,
                            slope.wane=slope.wane,
                            sigma=exp(par[(length(par)-1)]),
                            model=par.df1$model,
                            recru = exp(par[length(par)]),
                            mu.sick = mu.sick,
                            mu.rec = mu.rec,
                            add.inf.mort = add.inf.mort,
                            add.rec.mort = add.rec.mort,
                            N_stat=N_stat)
    
  }else if (par.df1$model =="MSIRNS"){
    
    #  par = c(beta, sigma, recru)
    
    if ((exp(par[(length(par)-1)])>26) | (exp(par[(length(par)-1)])<1) | (exp(par[length(par)])>26)) {
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
      
    }
    
    out = sim.met.MSIRNS.age(burnin=burnin,
                             sim_pop=sim_pop,
                             yrs=yrs,
                             ntyr=ntyr,
                             s=s,
                             beta=exp(par[1:(length(par)-2)]),
                             recov=par.df1$recov,
                             mort=par.df1$mort, #annual
                             mort_juv=par.df1$mort_juv, #annual
                             adult_fec=par.df1$adult_fec, #annual
                             wane=par.df1$wane,
                             age.brk=age.brk,
                             slope.wane=slope.wane,
                             sigma=exp(par[(length(par)-1)]),
                             model=par.df1$model,
                             recru = exp(par[length(par)]),
                             mu.sick = mu.sick,
                             mu.rec = mu.rec,
                             add.inf.mort = add.inf.mort,
                             add.rec.mort = add.rec.mort,
                             N_stat=N_stat)
  }else if (par.df1$model =="MSIRNNS"){
    
    #     par = c(beta, sigma, recru, epsilon)
    
    if ((exp(par[(length(par)-2)])>26) | (exp(par[(length(par)-2)])<1) | (exp(par[(length(par)-1)])>26) | (exp(par[length(par)])>26) ) {
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
      
    }
    
    out = sim.met.MSIRNNS.age(burnin=burnin,
                              sim_pop=sim_pop,
                              yrs=yrs,
                              ntyr=ntyr,
                              s=s,
                              beta=exp(par[1:(length(par)-3)]),
                              recov=par.df1$recov,
                              mort=par.df1$mort, #annual
                              mort_juv=par.df1$mort_juv, #annual
                              adult_fec=par.df1$adult_fec, #annual
                              wane=par.df1$wane,
                              age.brk=age.brk,
                              slope.wane=slope.wane,
                              sigma=exp(par[(length(par)-2)]),
                              model=par.df1$model,
                              recru = exp(par[(length(par)-1)]),
                              epsilon = exp(par[length(par)]),
                              mu.sick = mu.sick,
                              mu.rec = mu.rec,
                              add.inf.mort = add.inf.mort,
                              add.rec.mort = add.rec.mort,
                              N_stat=N_stat)
  }else if (model =="MSIRNIRN"){
    
    #par = c(beta, sigma, recru, sigma2)
    
    if ((exp(par[(length(par)-2)])>26) | (exp(par[(length(par)-2)])<1) | (exp(par[(length(par)-1)])>26) | (exp(par[length(par)])>26) ) {
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
      
    }
    
    
    out  = sim.met.MSIRNIRN.age(burnin=burnin,
                                sim_pop=sim_pop,
                                yrs=yrs,
                                ntyr=ntyr,
                                s=s,
                                beta=exp(par[1:(length(par)-3)]),
                                age.brk = age.brk,
                                recov=par.df1$recov,
                                mort=par.df1$mort, 
                                mort_juv=par.df1$mort_juv,
                                adult_fec=par.df1$adult_fec, 
                                wane=par.df1$wane,
                                recru=exp(par[(length(par)-1)]),
                                slope.wane=slope.wane,
                                sigma=exp(par[(length(par)-2)]),
                                sigma2=exp(par[length(par)]),
                                model=par.df1$model,
                                mu.sick = par.df1$mu_sick,
                                mu.rec = par.df1$mu_rec,
                                add.inf.mort = add.inf.mort,
                                N_stat=N_stat,
                                add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRNIRN-LT"){
    
    
    
    if ((exp(par[(length(par)-3)])>26) | (exp(par[(length(par)-3)])<1) | (exp(par[(length(par)-2)])>26) | (exp(par[(length(par)-1)])>26) | (exp(par[length(par)])>26) ) {
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
      
    }
    #else if(model=="MSIRNIRN-LT"){
    # par = c(beta, sigma, recru, sigma2, epsilon)
    
    out  = sim.met.MSIRNIRN.LT.age(burnin=burnin,
                                   sim_pop=sim_pop,
                                   yrs=yrs,
                                   ntyr=ntyr,
                                   s=s,
                                   beta=exp(par[1:(length(par)-4)]),
                                   age.brk = age.brk,
                                   recov=par.df1$recov,
                                   mort=par.df1$mort, 
                                   mort_juv=par.df1$mort_juv,
                                   adult_fec=par.df1$adult_fec, 
                                   wane=par.df1$wane,
                                   recru=exp(par[(length(par)-2)]),
                                   slope.wane=slope.wane,
                                   sigma=exp(par[(length(par)-3)]),
                                   sigma2=exp(par[(length(par)-1)]),
                                   epsilon=exp(par[length(par)]),
                                   model=par.df1$model,
                                   mu.sick = par.df1$mu_sick,
                                   mu.rec = par.df1$mu_rec,
                                   add.inf.mort = add.inf.mort,
                                   N_stat=N_stat,
                                   add.rec.mort = add.rec.mort) 
    
  } else if (model =="MSIRNSIRN"){
    
    #par = c(beta, sigma, recru, sigma2, beta2)
    
    if ((exp(par[(length(par)-3)])>26) | (exp(par[(length(par)-3)])<1) | (exp(par[(length(par)-2)])>26) | (exp(par[(length(par)-1)])>26) ) {
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
      
    }
    
    out  = sim.met.MSIRNSIRN.age(burnin=burnin,
                                 sim_pop=sim_pop,
                                 yrs=yrs,
                                 ntyr=ntyr,
                                 s=s,
                                 beta=exp(par[1:(length(par)-4)]),
                                 age.brk = age.brk,
                                 recov=par.df1$recov,
                                 mort=par.df1$mort, 
                                 mort_juv=par.df1$mort_juv,
                                 adult_fec=par.df1$adult_fec, 
                                 wane=par.df1$wane,
                                 recru=exp(par[(length(par)-2)]),
                                 slope.wane=slope.wane,
                                 sigma=exp(par[(length(par)-3)]),
                                 sigma2=exp(par[(length(par)-1)]),
                                 beta2=exp(par[length(par)]),
                                 model=model,
                                 mu.sick = par.df1$mu_sick,
                                 mu.rec = par.df1$mu_rec,
                                 add.inf.mort = add.inf.mort,
                                 N_stat=N_stat,
                                 add.rec.mort = add.rec.mort) 
    
  } else if (model =="MSIRNSIRN-LT"){
    
    #par = c(beta, sigma, recru, sigma2, beta2, epsilon)
    
    if ((exp(par[(length(par)-4)])>26) | (exp(par[(length(par)-4)])<1) | (exp(par[(length(par)-3)])>26) | (exp(par[(length(par)-2)])>26)|(exp(par[(length(par))])>26)) {
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
      
    }
    
    out  = sim.met.MSIRNSIRN.LT.age(burnin=burnin,
                                    sim_pop=sim_pop,
                                    yrs=yrs,
                                    ntyr=ntyr,
                                    s=s,
                                    beta=exp(par[1:(length(par)-5)]),
                                    age.brk = age.brk,
                                    recov=par.df1$recov,
                                    mort=par.df1$mort, 
                                    mort_juv=par.df1$mort_juv,
                                    adult_fec=par.df1$adult_fec, 
                                    wane=par.df1$wane,
                                    recru=exp(par[(length(par)-3)]),
                                    slope.wane=slope.wane,
                                    sigma=exp(par[(length(par)-4)]),
                                    sigma2=exp(par[(length(par)-2)]),
                                    beta2=exp(par[(length(par)-1)]),
                                    epsilon=exp(par[length(par)]),
                                    model=model,
                                    mu.sick = par.df1$mu_sick,
                                    mu.rec = par.df1$mu_rec,
                                    add.inf.mort = add.inf.mort,
                                    N_stat=N_stat,
                                    add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRSIR"){
    
    # par = c(beta, sigma, beta2)
    
    if ((exp(par[(length(par)-1)])<1) | (exp(par[(length(par)-1)])>26)) {
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
      
    }
    
    out  = sim.met.MSIRSIR.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=exp(par[1:(length(par)-2)]),
                               age.brk = age.brk,
                               recov=par.df1$recov,
                               mort=par.df1$mort, 
                               mort_juv=par.df1$mort_juv,
                               adult_fec=par.df1$adult_fec, 
                               wane=par.df1$wane,
                               slope.wane=slope.wane,
                               sigma=exp(par[(length(par)-1)]),
                               beta2=exp(par[length(par)]),
                               model=par.df1$model,
                               mu.sick = par.df1$mu_sick,
                               mu.rec = par.df1$mu_rec,
                               add.inf.mort = add.inf.mort,
                               add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRNIR"){
    
    #par = c(beta, sigma, recru, sigma2)
    
    if ((exp(par[(length(par)-2)])>26) | (exp(par[(length(par)-2)])<1) | (exp(par[(length(par)-1)])>26) | (exp(par[length(par)])>26) ) {
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
      
    }
    
    
    out  = sim.met.MSIRNIR.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=exp(par[1:(length(par)-3)]),
                               age.brk = age.brk,
                               recov=par.df1$recov,
                               mort=par.df1$mort, 
                               mort_juv=par.df1$mort_juv,
                               adult_fec=par.df1$adult_fec, 
                               wane=par.df1$wane,
                               recru=exp(par[(length(par)-1)]),
                               slope.wane=slope.wane,
                               sigma=exp(par[(length(par)-2)]),
                               sigma2=exp(par[length(par)]),
                               model=par.df1$model,
                               mu.sick = par.df1$mu_sick,
                               mu.rec = par.df1$mu_rec,
                               add.inf.mort = add.inf.mort,
                               N_stat=N_stat,
                               add.rec.mort = add.rec.mort) 
    
  }else if (model =="MSIRNSIRNS"){
    
    #par = c(beta, sigma, recru, sigma2, beta2, recru2)
    
    if ((exp(par[(length(par)-4)])>26) | (exp(par[(length(par)-4)])<1) | (exp(par[(length(par)-3)])>26) | (exp(par[(length(par)-2)])>26) | (exp(par[(length(par))])>26) ) {
      ll=-1*10^9#write over with real number that is huge
      print(-ll)
      return(-ll)
      
    }
    
    out  = sim.met.MSIRNSIRNS.age(burnin=burnin,
                                  sim_pop=sim_pop,
                                  yrs=yrs,
                                  ntyr=ntyr,
                                  s=s,
                                  beta=exp(par[1:(length(par)-5)]),
                                  age.brk = age.brk,
                                  recov=par.df1$recov,
                                  mort=par.df1$mort, 
                                  mort_juv=par.df1$mort_juv,
                                  adult_fec=par.df1$adult_fec, 
                                  wane=par.df1$wane,
                                  recru=exp(par[(length(par)-3)]),
                                  slope.wane=slope.wane,
                                  sigma=exp(par[(length(par)-4)]),
                                  sigma2=exp(par[(length(par)-2)]),
                                  beta2=exp(par[(length(par)-1)]),
                                  recru2=exp(par[length(par)]),
                                  model=model,
                                  mu.sick = par.df1$mu_sick,
                                  mu.rec = par.df1$mu_rec,
                                  add.inf.mort = add.inf.mort,
                                  N_stat=N_stat,
                                  add.rec.mort = add.rec.mort) 
    
  } 
  #first, run your model at the desired parameters for estimation (par)
  #sometimes, they don't produce anything at all - in which case, just return a giant ll here and move on
  
  #print the guess par
  print(exp(par))
  
  if(is.na(unique(out$seropos)[1])){
    ll=-1*10^9#write over with real number that is huge
    print(-ll)
    return(-ll)
  }else{
    
    #your data should have ages by biweek
    #ggplot(out) + geom_point(aes(x=age, y=seroprevalence, color=biwk))
    # 
    # #fit a line  to these data to get the hazard by age:
    pred.fit <- with(out,smooth.spline(x=age, y=seroprevalence))
    #separate age data
    
    #data = subset(data, !is.na(age)) #already done
    #these were also done by biweek
    data = arrange(data, age) 
    
    #length(unique(sort(out$age)))#520
    #sort data by age and make a prediction for each age by biweek
    ages.new = sort(unique(data$age))
    
    out.comp <- cbind.data.frame(age=ages.new, seroprev=predict(pred.fit,x=ages.new)$y)
    #ggplot(out) + geom_point(aes(x=age, y=seroprevalence, color=biwk)) + geom_line(data=out.comp, aes(x=age, y=seroprev))
    out.comp$seroprev[out.comp$seroprev<0] =0
    out.comp$seroprev[out.comp$seroprev>1] =1
    
    
    #now, compare model and data
    
    ll=0 #sets initial log-likelihood to 0
    
    #then, compare with data by age and biweek
    for(a in 1:length(data$age)){ 
      #now run through data and compare with model
      #for-loops over every individual in the dataset, starting above the lower bound - 
      
      
      pprev <- out.comp$seroprev[out.comp$age==data$age[a]]  
      
      
      #print(pprev)
      if(cutoff=="mean"){
        ll=ll+dbinom(data$prev[a],1,pprev,log=T)    
        print(ll)
      }else if (cutoff=="lci"){
        ll=ll+dbinom(data$prev_lci[a],1,pprev,log=T)    
      }else if(cutoff=="uci"){
        ll=ll+dbinom(data$prev_uci[a],1,pprev,log=T)    
      }
      
    }
    print(-ll)
    #calculate and return the negative log-likelihood 
    #force mu.sick to be 1 or more
    
    return(-ll)
  }}


#need to have it specify where the age-breaks are and where the beta changes
#need to have it specify where the age-breaks are and where the beta changes
wrap.fit.CI <- function(model, beta, sigma, recru, epsilon, sigma2, beta2, recru2, mu.sick, mu.rec, dat1, par.df, cutoff1, species1, antigen1, s, sim_pop, yrs, ntyr,  N_stat,
                     age.brk, slope.wane, burnin, do.plot, do.save, filename){
  
  model = sapply(strsplit(par.df$model, "-mu"),'[',1)
  
  #first pull out your data...
  dat.tmp <- subset(dat1, species==species1 & antigen==antigen1)
  dat.tmp <- subset(dat.tmp, !is.na(age))
  
  #feed in the right par.dat 
  
  #dat.tmp <- subset(dat.tmp, fit==1). do you use all the data or just the longitudinal resamples?
  
  #identify your parameter(s) for fitting based on your model of choice
  if(model=="MSIR"){
    par = c(beta)  
  }else if(model=="MSIRS" | model=="MSIRN"){
    par = c(beta, sigma)  
  }else if (model=="MSIRNI" | model=="MSIRNS"){
    par = c(beta, sigma, recru)
  }else if(model == "MSIRNNS"){
    par = c(beta, sigma, recru, epsilon)
  }else if (model == "MSIRNIRN" | model == "MSIRNIR"){
    par = c(beta, sigma, recru, sigma2)
  }else if (model =="MSIRNSIRN"){
    par = c(beta, sigma, recru, sigma2, beta2)
  }else if (model == "MSIRSIR"){
    par = c(beta, sigma, beta2)
  }else if(model=="MSIRNIRN-LT"){
    par = c(beta, sigma, recru, sigma2, epsilon)
  }else if(model=="MSIRNSIRN-LT"){
    par = c(beta, sigma, recru, sigma2, beta2, epsilon)
  }else if (model =="MSIRNSIRNS"){
    par = c(beta, sigma, recru, sigma2, beta2, recru2)
  }
  
  if(mu.sick >1 & mu.rec==1){
    par <- c(par, mu.sick)
  }
  if(mu.sick ==1 & mu.rec >1){
    par <- c(par, mu.rec)
  }
  if(mu.sick >1 & mu.rec >1){
    par <- c(par, mu.sick, mu.rec)
  }
  
  print(par)
  #estimate parameters
  
  out.lik <- optim(par=log(par),fn=log.lik.age, method="Nelder-Mead", data=dat.tmp, model=model, par.df1 =par.df, burnin=burnin, 
                   cutoff=cutoff1, s=s, sim_pop=sim_pop, ntyr=ntyr, yrs=yrs, age.brk=age.brk, 
                   N_stat=N_stat, control=list(maxit=1000), hessian = TRUE)
  
  if (is.positive.definite(out.lik$hessian)){
    #save the raw CI values from the hessian (not including variance)
    hess <- solve(out.lik$hessian)
    prop_sigma <-sqrt(diag(hess))
    upper<-exp(out.lik$par+1.96*prop_sigma)
    lower<-exp(out.lik$par-1.96*prop_sigma)
    ConfidenceInterval<-data.frame(value=exp(out.lik$par), upper=upper, lower=lower)
    
    ConfidenceInterval[ConfidenceInterval<0] <- 0
    ConfidenceInterval$par <- NA
    
  }else{ 
    #make it positive definite if it is not
    hess <- solve(make.positive.definite(out.lik$hessian))
    prop_sigma <-sqrt(diag(hess))
    upper<-exp(out.lik$par+1.96*prop_sigma)
    lower<-exp(out.lik$par-1.96*prop_sigma)
    
    ConfidenceInterval<-data.frame(value=exp(out.lik$par), upper=upper, lower=lower)
    #ConfidenceInterval$par <- NA
  }
  
 
  #solve the hessian to get CIs using mvnorm
  #sigma is the covariance matrix which is the inverse of the hessian 
  #when you are minimizing the neg log liklihood in your function
  par.dist <- rmvnorm(n=500, exp(out.lik$par), sigma=hess)
  par.dist[ par.dist < 0] <- 0
  
  #now have 500 guesses for each par, all calculated in tandem.
  
  ## get predictions across this span - including variance covariance
  ## span should account for biwk of year and age across 1 year and the full age range of the model
  ## now have a matrix to hold output for all ages and biweeks
  age.span = rep(0:(s-1), each = ntyr)
  biwk.span = rep(1:ntyr, each = s)
  
  #drop the last term in each
  age.span <- age.span[1:(length(age.span)-1)]
  #biwk.span <- biwk.span[1:(length(biwk.span))]
  
  preds <- matrix(NA,nrow=length(biwk.span),ncol=length(par.dist[,1]))
  
  if(model=="MSIR"){
    
    if(mu.sick >1 & mu.rec==1){
      par.df$mu_sick = ConfidenceInterval$value[length(out.lik$par)] 
      par.df$mu_sick_lci = ConfidenceInterval$lower[length(out.lik$par)] 
      par.df$mu_sick_uci = ConfidenceInterval$upper[length(out.lik$par)] 
      
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-1)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-1)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-1)]
      
      
      out1 = sim.met.MSIR.age(burnin=burnin,
                             sim_pop=sim_pop,
                             yrs=yrs,
                             ntyr=ntyr,
                             s=s,
                             beta=par.df$beta,
                             sero=par.df$sero,
                             recov=par.df$recov,
                             age.brk=age.brk,
                             mort=par.df$mort, #annual
                             mort_juv=par.df$mort_juv, #annual
                             adult_fec=par.df$adult_fec, #annual
                             wane=par.df$wane,
                             slope.wane=slope.wane,
                             sigma=par.df$sigma,#biwks
                             model=par.df$model,
                             mu.sick = par.df$mu_sick,
                             mu.rec = par.df$mu_rec,
                             add.inf.mort = TRUE,
                             add.rec.mort = FALSE)
      

      
      for (j in 1:length(par.dist[,1])){
        
        out = sim.met.MSIR.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.dist[j,1],
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.df$sigma,#biwks
                               model=par.df$model,
                               mu.sick = par.dist[j,2],
                               mu.rec = par.df$mu_rec,
                               add.inf.mort = TRUE,
                               add.rec.mort = FALSE)
        
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
      
    }else if(mu.sick ==1 & mu.rec >1){
      
      par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
      par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
      par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
      
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-1)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-1)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-1)]
      
      out1 = sim.met.MSIR.age(burnin=burnin,
                              sim_pop=sim_pop,
                              yrs=yrs,
                              ntyr=ntyr,
                              s=s,
                              beta=par.df$beta,
                              sero=par.df$sero,
                              recov=par.df$recov,
                              age.brk=age.brk,
                              mort=par.df$mort, #annual
                              mort_juv=par.df$mort_juv, #annual
                              adult_fec=par.df$adult_fec, #annual
                              wane=par.df$wane,
                              slope.wane=slope.wane,
                              sigma=par.df$sigma,#biwks
                              model=par.df$model,
                              mu.sick = par.df$mu_sick,
                              mu.rec = par.df$mu_rec,
                              add.inf.mort = FALSE,
                              add.rec.mort = TRUE)
      
      
      for (j in 1:length(par.dist[,1])){
        
        out = sim.met.MSIR.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.dist[j,1],
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.df$sigma,#biwks
                               model=par.df$model,
                               mu.sick = par.df$mu_sick,
                               mu.rec = par.dist[j,2],
                               add.inf.mort = FALSE,
                               add.rec.mort = TRUE)
        
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
    }else if(mu.sick >1 & mu.rec >1){
      
      par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      
      par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
      par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
      par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
      
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-2)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-2)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-2)]
      
      out1 = sim.met.MSIR.age(burnin=burnin,
                              sim_pop=sim_pop,
                              yrs=yrs,
                              ntyr=ntyr,
                              s=s,
                              beta=par.df$beta,
                              sero=par.df$sero,
                              recov=par.df$recov,
                              age.brk=age.brk,
                              mort=par.df$mort, #annual
                              mort_juv=par.df$mort_juv, #annual
                              adult_fec=par.df$adult_fec, #annual
                              wane=par.df$wane,
                              slope.wane=slope.wane,
                              sigma=par.df$sigma,#biwks
                              model=par.df$model,
                              mu.sick = par.df$mu_sick,
                              mu.rec = par.df$mu_rec,
                              add.inf.mort = TRUE,
                              add.rec.mort = TRUE)
      
      
      for (j in 1:length(par.dist[,1])){
        
        out = sim.met.MSIR.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.dist[j,1],
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.df$sigma,#biwks
                               model=par.df$model,
                               mu.sick = par.dist[j,2],
                               mu.rec = par.dist[j,3],
                               add.inf.mort = TRUE,
                               add.rec.mort = TRUE)
        
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
    } else{
      
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par))] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par))]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par))]
      
      out1 = sim.met.MSIR.age(burnin=burnin,
                              sim_pop=sim_pop,
                              yrs=yrs,
                              ntyr=ntyr,
                              s=s,
                              beta=par.df$beta,
                              sero=par.df$sero,
                              recov=par.df$recov,
                              age.brk=age.brk,
                              mort=par.df$mort, #annual
                              mort_juv=par.df$mort_juv, #annual
                              adult_fec=par.df$adult_fec, #annual
                              wane=par.df$wane,
                              slope.wane=slope.wane,
                              sigma=par.df$sigma,#biwks
                              model=par.df$model,
                              mu.sick = par.df$mu_sick,
                              mu.rec = par.df$mu_rec,
                              add.inf.mort = FALSE,
                              add.rec.mort = FALSE)
      
      
      for (j in 1:length(par.dist[,1])){
        
        out = sim.met.MSIR.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.dist[j,1],
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.df$sigma,#biwks
                               model=par.df$model,
                               mu.sick = par.df$mu_sick,
                               mu.rec = par.df$mu_rec,
                               add.inf.mort = FALSE,
                               add.rec.mort = FALSE)
        
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
    }
  }else if(model == "MSIRS"){
    if(mu.sick >1 & mu.rec==1){
      par.df$mu_sick = ConfidenceInterval$value[length(out.lik$par)] 
      par.df$mu_sick_lci = ConfidenceInterval$lower[length(out.lik$par)] 
      par.df$mu_sick_uci = ConfidenceInterval$upper[length(out.lik$par)] 
      
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-2)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-2)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-2)]
      
      out1 = sim.met.MSIR.age(burnin=burnin,
                             sim_pop=sim_pop,
                             yrs=yrs,
                             ntyr=ntyr,
                             s=s,
                             beta=par.df$beta,
                             sero=par.df$sero,
                             recov=par.df$recov,
                             age.brk=age.brk,
                             mort=par.df$mort, #annual
                             mort_juv=par.df$mort_juv, #annual
                             adult_fec=par.df$adult_fec, #annual
                             wane=par.df$wane,
                             slope.wane=slope.wane,
                             sigma=par.df$sigma,
                             model=par.df$model,
                             mu.sick = par.df$mu_sick,
                             mu.rec = par.df$mu_rec,
                             add.inf.mort = TRUE,
                             add.rec.mort = FALSE)
      
    for (j in 1:length(par.dist[,1])){
      out = sim.met.MSIR.age(burnin=burnin,
                             sim_pop=sim_pop,
                             yrs=yrs,
                             ntyr=ntyr,
                             s=s,
                             beta=par.dist[j,1],
                             sero=par.df$sero,
                             recov=par.df$recov,
                             age.brk=age.brk,
                             mort=par.df$mort, #annual
                             mort_juv=par.df$mort_juv, #annual
                             adult_fec=par.df$adult_fec, #annual
                             wane=par.df$wane,
                             slope.wane=slope.wane,
                             sigma=par.dist[j,2],
                             model=par.df$model,
                             mu.sick = par.dist[j,3],
                             mu.rec = par.df$mu_rec,
                             add.inf.mort = TRUE,
                             add.rec.mort = FALSE)
      #each column is a new guess. each row is output by age and biweek
      preds[,j] <- out$seroprevalence
      
      }
      
      
    }else if(mu.sick ==1 & mu.rec >1){
      
      par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
      par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
      par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
      
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-2)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-2)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-2)]
      
      out1 = sim.met.MSIR.age(burnin=burnin,
                              sim_pop=sim_pop,
                              yrs=yrs,
                              ntyr=ntyr,
                              s=s,
                              beta=par.df$beta,
                              sero=par.df$sero,
                              recov=par.df$recov,
                              age.brk=age.brk,
                              mort=par.df$mort, #annual
                              mort_juv=par.df$mort_juv, #annual
                              adult_fec=par.df$adult_fec, #annual
                              wane=par.df$wane,
                              slope.wane=slope.wane,
                              sigma=par.df$sigma,
                              model=par.df$model,
                              mu.sick = par.df$mu_sick,
                              mu.rec = par.df$mu_rec,
                              add.inf.mort = FALSE,
                              add.rec.mort = TRUE)
      
      for (j in 1:length(par.dist[,1])){
        out = sim.met.MSIR.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.dist[j,1],
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.dist[j,2],
                               model=par.df$model,
                               mu.sick = par.df$mu_sick,
                               mu.rec = par.dist[j,3],
                               add.inf.mort = FALSE,
                               add.rec.mort = TRUE)
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
    }else if (mu.sick >1 & mu.rec >1){
      
      par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
      par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
      par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
      
      par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
      par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
      par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
      
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
      
      out1 = sim.met.MSIR.age(burnin=burnin,
                              sim_pop=sim_pop,
                              yrs=yrs,
                              ntyr=ntyr,
                              s=s,
                              beta=par.df$beta,
                              sero=par.df$sero,
                              recov=par.df$recov,
                              age.brk=age.brk,
                              mort=par.df$mort, #annual
                              mort_juv=par.df$mort_juv, #annual
                              adult_fec=par.df$adult_fec, #annual
                              wane=par.df$wane,
                              slope.wane=slope.wane,
                              sigma=par.df$sigma,
                              model=par.df$model,
                              mu.sick = par.df$mu_sick,
                              mu.rec = par.df$mu_rec,
                              add.inf.mort = TRUE,
                              add.rec.mort = TRUE)
      
      for (j in 1:length(par.dist[,1])){
        out = sim.met.MSIR.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.dist[j,1],
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.dist[j,2],
                               model=par.df$model,
                               mu.sick = par.dist[j,3],
                               mu.rec = par.dist[j,4],
                               add.inf.mort = TRUE,
                               add.rec.mort = TRUE)
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
      
    }else{
      
      
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par))] 
      par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par))] 
      par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par))] 
      
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-1)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-1)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-1)]
      
      out1 = sim.met.MSIR.age(burnin=burnin,
                              sim_pop=sim_pop,
                              yrs=yrs,
                              ntyr=ntyr,
                              s=s,
                              beta=par.df$beta,
                              sero=par.df$sero,
                              recov=par.df$recov,
                              age.brk=age.brk,
                              mort=par.df$mort, #annual
                              mort_juv=par.df$mort_juv, #annual
                              adult_fec=par.df$adult_fec, #annual
                              wane=par.df$wane,
                              slope.wane=slope.wane,
                              sigma=par.df$sigma,
                              model=par.df$model,
                              mu.sick = par.df$mu_sick,
                              mu.rec = par.df$mu_rec,
                              add.inf.mort = FALSE,
                              add.rec.mort = FALSE)
      
      for (j in 1:length(par.dist[,1])){
        out = sim.met.MSIR.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.dist[j,1],
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.dist[j,2],
                               model=par.df$model,
                               mu.sick = par.df$mu_sick,
                               mu.rec = par.df$mu_rec,
                               add.inf.mort = FALSE,
                               add.rec.mort = FALSE)
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
      
    }
  }else if(model == "MSIRN"){
    if(mu.sick >1 & mu.rec==1){
      par.df$mu_sick = ConfidenceInterval$value[length(out.lik$par)] 
      par.df$mu_sick_lci = ConfidenceInterval$lower[length(out.lik$par)] 
      par.df$mu_sick_uci = ConfidenceInterval$upper[length(out.lik$par)] 
      
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-2)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-2)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-2)]
      
      out1 = sim.met.MSIRN.age(burnin=burnin,
                              sim_pop=sim_pop,
                              yrs=yrs,
                              ntyr=ntyr,
                              s=s,
                              beta=par.df$beta,
                              sero=par.df$sero,
                              recov=par.df$recov,
                              age.brk=age.brk,
                              mort=par.df$mort, #annual
                              mort_juv=par.df$mort_juv, #annual
                              adult_fec=par.df$adult_fec, #annual
                              wane=par.df$wane,
                              slope.wane=slope.wane,
                              sigma=par.df$sigma,
                              model=par.df$model,
                              mu.sick = par.df$mu_sick,
                              mu.rec = par.df$mu_rec,
                              add.inf.mort = TRUE,
                              add.rec.mort = FALSE,
                              N_stat = "matAB")
      
      for (j in 1:length(par.dist[,1])){
        out = sim.met.MSIRN.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.dist[j,1],
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.dist[j,2],
                               model=par.df$model,
                               mu.sick = par.dist[j,3],
                               mu.rec = par.df$mu_rec,
                               add.inf.mort = TRUE,
                               add.rec.mort = FALSE,
                               N_stat = "matAB")
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
      
    }else if(mu.sick ==1 & mu.rec >1){
      
      par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
      par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
      par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
      
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-2)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-2)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-2)]
      
      out1 = sim.met.MSIRN.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.df$beta,
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.df$sigma,
                               model=par.df$model,
                               mu.sick = par.df$mu_sick,
                               mu.rec = par.df$mu_rec,
                               add.inf.mort = FALSE,
                               add.rec.mort = TRUE,
                               N_stat = "matAB")
      
      for (j in 1:length(par.dist[,1])){
        out = sim.met.MSIRN.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.dist[j,1],
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.dist[j,2],
                               model=par.df$model,
                               mu.sick = par.df$mu_sick,
                               mu.rec = par.dist[j,3],
                               add.inf.mort = FALSE,
                               add.rec.mort = TRUE,
                               N_stat = "matAB")
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
    }else if (mu.sick >1 & mu.rec >1){
      
      par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
      par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
      par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
      
      par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
      par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
      par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
      
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
      
      out1 = sim.met.MSIRN.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.df$beta,
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.df$sigma,
                               model=par.df$model,
                               mu.sick = par.df$mu_sick,
                               mu.rec = par.df$mu_rec,
                               add.inf.mort = TRUE,
                               add.rec.mort = TRUE,
                               N_stat = "matAB")
      
      for (j in 1:length(par.dist[,1])){
        out = sim.met.MSIRN.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.dist[j,1],
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.dist[j,2],
                               model=par.df$model,
                               mu.sick = par.dist[j,3],
                               mu.rec = par.dist[j,4],
                               add.inf.mort = TRUE,
                               add.rec.mort = TRUE,
                               N_stat = "matAB")
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
      
    }else{
      
      
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par))] 
      par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par))] 
      par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par))] 
      
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-1)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-1)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-1)]
      
      out1 = sim.met.MSIRN.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.df$beta,
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.df$sigma,
                               model=par.df$model,
                               mu.sick = par.df$mu_sick,
                               mu.rec = par.df$mu_rec,
                               add.inf.mort = FALSE,
                               add.rec.mort = FALSE,
                               N_stat = "matAB")
      
      for (j in 1:length(par.dist[,1])){
        out = sim.met.MSIRN.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.dist[j,1],
                               sero=par.df$sero,
                               recov=par.df$recov,
                               age.brk=age.brk,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               slope.wane=slope.wane,
                               sigma=par.dist[j,2],
                               model=par.df$model,
                               mu.sick = par.df$mu_sick,
                               mu.rec = par.df$mu_rec,
                               add.inf.mort = FALSE,
                               add.rec.mort = FALSE,
                               N_stat = "matAB")
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
      
    }
  }else if (model == "MSIRNS"){
    
    if(mu.sick >1 & mu.rec==1){
      par.df$mu_sick = ConfidenceInterval$value[length(out.lik$par)] 
      par.df$mu_sick_lci = ConfidenceInterval$lower[length(out.lik$par)] 
      par.df$mu_sick_uci = ConfidenceInterval$upper[length(out.lik$par)] 
      
      par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
      par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
      par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
      
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
      
      out1 = sim.met.MSIRNS.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.df$beta,
                               recov=par.df$recov,
                               mort=par.df$mort, #annual
                               mort_juv=par.df$mort_juv, #annual
                               adult_fec=par.df$adult_fec, #annual
                               wane=par.df$wane,
                               age.brk=age.brk,
                               slope.wane=slope.wane,
                               sigma=par.df$sigma,
                               model=par.df$model,
                               recru = par.df$recru,
                               mu.sick = par.df$mu_sick,
                               mu.rec = par.df$mu_rec,
                               add.inf.mort = TRUE,
                               add.rec.mort = FALSE,
                               N_stat="matAB")
      
      for (j in 1:length(par.dist[,1])){
        
        out = sim.met.MSIRNS.age(burnin=burnin,
                                 sim_pop=sim_pop,
                                 yrs=yrs,
                                 ntyr=ntyr,
                                 s=s,
                                 beta=par.dist[j,1],
                                 recov=par.df$recov,
                                 mort=par.df$mort, #annual
                                 mort_juv=par.df$mort_juv, #annual
                                 adult_fec=par.df$adult_fec, #annual
                                 wane=par.df$wane,
                                 age.brk=age.brk,
                                 slope.wane=slope.wane,
                                 sigma=par.dist[j,2],
                                 model=par.df$model,
                                 recru = par.dist[j,3],
                                 mu.sick = par.dist[j,4],
                                 mu.rec = par.df$mu_rec,
                                 add.inf.mort = TRUE,
                                 add.rec.mort = FALSE,
                                 N_stat="matAB")
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
    }else if(mu.sick ==1 & mu.rec >1){
      
      par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
      par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
      par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
      
      par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
      par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
      par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
      
      out1 = sim.met.MSIRNS.age(burnin=burnin,
                                sim_pop=sim_pop,
                                yrs=yrs,
                                ntyr=ntyr,
                                s=s,
                                beta=par.df$beta,
                                recov=par.df$recov,
                                mort=par.df$mort, #annual
                                mort_juv=par.df$mort_juv, #annual
                                adult_fec=par.df$adult_fec, #annual
                                wane=par.df$wane,
                                age.brk=age.brk,
                                slope.wane=slope.wane,
                                sigma=par.df$sigma,
                                model=par.df$model,
                                recru = par.df$recru,
                                mu.sick = par.df$mu_sick,
                                mu.rec = par.df$mu_rec,
                                add.inf.mort = FALSE,
                                add.rec.mort = TRUE,
                                N_stat="matAB")
      
      for (j in 1:length(par.dist[,1])){
        
        out = sim.met.MSIRNS.age(burnin=burnin,
                                 sim_pop=sim_pop,
                                 yrs=yrs,
                                 ntyr=ntyr,
                                 s=s,
                                 beta=par.dist[j,1],
                                 recov=par.df$recov,
                                 mort=par.df$mort, #annual
                                 mort_juv=par.df$mort_juv, #annual
                                 adult_fec=par.df$adult_fec, #annual
                                 wane=par.df$wane,
                                 age.brk=age.brk,
                                 slope.wane=slope.wane,
                                 sigma=par.dist[j,2],
                                 model=par.df$model,
                                 recru = par.dist[j,3],
                                 mu.sick = par.df$mu_sick,
                                 mu.rec = par.dist[j,4],
                                 add.inf.mort = FALSE,
                                 add.rec.mort = TRUE,
                                 N_stat="matAB")
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
      
    }else if (mu.sick >1 & mu.rec >1){
      par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
      par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
      par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
      
      par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      
      par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-2)] 
      par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
      par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
      
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
      par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
      par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
      
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
      
      out1 = sim.met.MSIRNS.age(burnin=burnin,
                                sim_pop=sim_pop,
                                yrs=yrs,
                                ntyr=ntyr,
                                s=s,
                                beta=par.df$beta,
                                recov=par.df$recov,
                                mort=par.df$mort, #annual
                                mort_juv=par.df$mort_juv, #annual
                                adult_fec=par.df$adult_fec, #annual
                                wane=par.df$wane,
                                age.brk=age.brk,
                                slope.wane=slope.wane,
                                sigma=par.df$sigma,
                                model=par.df$model,
                                recru = par.df$recru,
                                mu.sick = par.df$mu_sick,
                                mu.rec = par.df$mu_rec,
                                add.inf.mort = TRUE,
                                add.rec.mort = TRUE,
                                N_stat="matAB")
      
      
      for (j in 1:length(par.dist[,1])){
        
        out = sim.met.MSIRNS.age(burnin=burnin,
                                 sim_pop=sim_pop,
                                 yrs=yrs,
                                 ntyr=ntyr,
                                 s=s,
                                 beta=par.dist[j,1],
                                 recov=par.df$recov,
                                 mort=par.df$mort, #annual
                                 mort_juv=par.df$mort_juv, #annual
                                 adult_fec=par.df$adult_fec, #annual
                                 wane=par.df$wane,
                                 age.brk=age.brk,
                                 slope.wane=slope.wane,
                                 sigma=par.dist[j,2],
                                 model=par.df$model,
                                 recru = par.dist[j,3],
                                 mu.sick = par.dist[j,4],
                                 mu.rec = par.dist[j,5],
                                 add.inf.mort = TRUE,
                                 add.rec.mort = TRUE,
                                 N_stat="matAB")
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
      
    }else{
      par.df$recru = ConfidenceInterval$value[(length(out.lik$par))] 
      par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par))] 
      par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par))] 
      
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-2)] 
      par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-2)]
      par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-2)]
      
      out1 = sim.met.MSIRNS.age(burnin=burnin,
                                sim_pop=sim_pop,
                                yrs=yrs,
                                ntyr=ntyr,
                                s=s,
                                beta=par.df$beta,
                                recov=par.df$recov,
                                mort=par.df$mort, #annual
                                mort_juv=par.df$mort_juv, #annual
                                adult_fec=par.df$adult_fec, #annual
                                wane=par.df$wane,
                                age.brk=age.brk,
                                slope.wane=slope.wane,
                                sigma=par.df$sigma,
                                model=par.df$model,
                                recru = par.df$recru,
                                mu.sick = par.df$mu_sick,
                                mu.rec = par.df$mu_rec,
                                add.inf.mort = FALSE,
                                add.rec.mort = FALSE,
                                N_stat="matAB")
      
      for (j in 1:length(par.dist[,1])){
        
        out = sim.met.MSIRNS.age(burnin=burnin,
                                 sim_pop=sim_pop,
                                 yrs=yrs,
                                 ntyr=ntyr,
                                 s=s,
                                 beta=par.dist[j,1],
                                 recov=par.df$recov,
                                 mort=par.df$mort, #annual
                                 mort_juv=par.df$mort_juv, #annual
                                 adult_fec=par.df$adult_fec, #annual
                                 wane=par.df$wane,
                                 age.brk=age.brk,
                                 slope.wane=slope.wane,
                                 sigma=par.dist[j,2],
                                 model=par.df$model,
                                 recru = par.dist[j,3],
                                 mu.sick = par.df$mu_sick,
                                 mu.rec = par.df$mu_rec,
                                 add.inf.mort = FALSE,
                                 add.rec.mort = FALSE,
                                 N_stat="matAB")
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
      }
      
      
    }
    }else if (model == "MSIRNI"){
      
      if(mu.sick >1 & mu.rec==1){
        par.df$mu_sick = ConfidenceInterval$value[length(out.lik$par)] 
        par.df$mu_sick_lci = ConfidenceInterval$lower[length(out.lik$par)] 
        par.df$mu_sick_uci = ConfidenceInterval$upper[length(out.lik$par)] 
        
        par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-1)] 
        par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
        par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
        
        par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
        par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
        par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
        
        
        
        par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
        par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
        par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
        
        out1 = sim.met.MSIRN.age(burnin=burnin,
                                sim_pop=sim_pop,
                                yrs=yrs,
                                ntyr=ntyr,
                                s=s,
                                beta=par.df$beta,
                                recov=par.df$recov,
                                mort=par.df$mort, #annual
                                mort_juv=par.df$mort_juv, #annual
                                adult_fec=par.df$adult_fec, #annual
                                wane=par.df$wane,
                                age.brk=age.brk,
                                slope.wane=slope.wane,
                                sigma=par.df$sigma,
                                model=par.df$model,
                                recru = par.df$recru,
                                mu.sick = par.df$mu_sick,
                                mu.rec = par.df$mu_rec,
                                add.inf.mort = TRUE,
                                add.rec.mort = FALSE,
                                N_stat="matAB")
        
        for (j in 1:length(par.dist[,1])){
          
          out = sim.met.MSIRN.age(burnin=burnin,
                                   sim_pop=sim_pop,
                                   yrs=yrs,
                                   ntyr=ntyr,
                                   s=s,
                                   beta=par.dist[j,1],
                                   recov=par.df$recov,
                                   mort=par.df$mort, #annual
                                   mort_juv=par.df$mort_juv, #annual
                                   adult_fec=par.df$adult_fec, #annual
                                   wane=par.df$wane,
                                   age.brk=age.brk,
                                   slope.wane=slope.wane,
                                   sigma=par.dist[j,2],
                                   model=par.df$model,
                                   recru = par.dist[j,3],
                                   mu.sick = par.dist[j,4],
                                   mu.rec = par.df$mu_rec,
                                   add.inf.mort = TRUE,
                                   add.rec.mort = FALSE,
                                   N_stat="matAB")
          #each column is a new guess. each row is output by age and biweek
          preds[,j] <- out$seroprevalence
          
        }
        
      }else if(mu.sick ==1 & mu.rec >1){
        
        par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
        par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
        par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
        
        par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-1)] 
        par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
        par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
        
        par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
        par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
        par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
        
        
        par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
        par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
        par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
        
        out1 = sim.met.MSIRN.age(burnin=burnin,
                                 sim_pop=sim_pop,
                                 yrs=yrs,
                                 ntyr=ntyr,
                                 s=s,
                                 beta=par.df$beta,
                                 recov=par.df$recov,
                                 mort=par.df$mort, #annual
                                 mort_juv=par.df$mort_juv, #annual
                                 adult_fec=par.df$adult_fec, #annual
                                 wane=par.df$wane,
                                 age.brk=age.brk,
                                 slope.wane=slope.wane,
                                 sigma=par.df$sigma,
                                 model=par.df$model,
                                 recru = par.df$recru,
                                 mu.sick = par.df$mu_sick,
                                 mu.rec = par.df$mu_rec,
                                 add.inf.mort = FALSE,
                                 add.rec.mort = TRUE,
                                 N_stat="matAB")
        
        for (j in 1:length(par.dist[,1])){
          
          out = sim.met.MSIRN.age(burnin=burnin,
                                   sim_pop=sim_pop,
                                   yrs=yrs,
                                   ntyr=ntyr,
                                   s=s,
                                   beta=par.dist[j,1],
                                   recov=par.df$recov,
                                   mort=par.df$mort, #annual
                                   mort_juv=par.df$mort_juv, #annual
                                   adult_fec=par.df$adult_fec, #annual
                                   wane=par.df$wane,
                                   age.brk=age.brk,
                                   slope.wane=slope.wane,
                                   sigma=par.dist[j,2],
                                   model=par.df$model,
                                   recru = par.dist[j,3],
                                   mu.sick = par.df$mu_sick,
                                   mu.rec = par.dist[j,4],
                                   add.inf.mort = FALSE,
                                   add.rec.mort = TRUE,
                                   N_stat="matAB")
          #each column is a new guess. each row is output by age and biweek
          preds[,j] <- out$seroprevalence
          
        }
        
        
      }else if (mu.sick >1 & mu.rec >1){
        par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
        par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
        par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
        
        par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
        par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
        par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
        
        par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-2)] 
        par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
        par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
        
        par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
        par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
        par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
        
        
        
        par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
        par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
        par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
        
        out1 = sim.met.MSIRN.age(burnin=burnin,
                                 sim_pop=sim_pop,
                                 yrs=yrs,
                                 ntyr=ntyr,
                                 s=s,
                                 beta=par.df$beta,
                                 recov=par.df$recov,
                                 mort=par.df$mort, #annual
                                 mort_juv=par.df$mort_juv, #annual
                                 adult_fec=par.df$adult_fec, #annual
                                 wane=par.df$wane,
                                 age.brk=age.brk,
                                 slope.wane=slope.wane,
                                 sigma=par.df$sigma,
                                 model=par.df$model,
                                 recru = par.df$recru,
                                 mu.sick = par.df$mu_sick,
                                 mu.rec = par.df$mu_rec,
                                 add.inf.mort = TRUE,
                                 add.rec.mort = TRUE,
                                 N_stat="matAB")
        
        
        for (j in 1:length(par.dist[,1])){
          
          out = sim.met.MSIRN.age(burnin=burnin,
                                   sim_pop=sim_pop,
                                   yrs=yrs,
                                   ntyr=ntyr,
                                   s=s,
                                   beta=par.dist[j,1],
                                   recov=par.df$recov,
                                   mort=par.df$mort, #annual
                                   mort_juv=par.df$mort_juv, #annual
                                   adult_fec=par.df$adult_fec, #annual
                                   wane=par.df$wane,
                                   age.brk=age.brk,
                                   slope.wane=slope.wane,
                                   sigma=par.dist[j,2],
                                   model=par.df$model,
                                   recru = par.dist[j,3],
                                   mu.sick = par.dist[j,4],
                                   mu.rec = par.dist[j,5],
                                   add.inf.mort = TRUE,
                                   add.rec.mort = TRUE,
                                   N_stat="matAB")
          #each column is a new guess. each row is output by age and biweek
          preds[,j] <- out$seroprevalence
          
        }
        
        
      }else{
        par.df$recru = ConfidenceInterval$value[(length(out.lik$par))] 
        par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par))] 
        par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par))] 
        
        par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-1)] 
        par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
        par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
        
        
        
        par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-2)] 
        par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-2)]
        par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-2)]
        
        out1 = sim.met.MSIRN.age(burnin=burnin,
                                 sim_pop=sim_pop,
                                 yrs=yrs,
                                 ntyr=ntyr,
                                 s=s,
                                 beta=par.df$beta,
                                 recov=par.df$recov,
                                 mort=par.df$mort, #annual
                                 mort_juv=par.df$mort_juv, #annual
                                 adult_fec=par.df$adult_fec, #annual
                                 wane=par.df$wane,
                                 age.brk=age.brk,
                                 slope.wane=slope.wane,
                                 sigma=par.df$sigma,
                                 model=par.df$model,
                                 recru = par.df$recru,
                                 mu.sick = par.df$mu_sick,
                                 mu.rec = par.df$mu_rec,
                                 add.inf.mort = FALSE,
                                 add.rec.mort = FALSE,
                                 N_stat="matAB")
        
        for (j in 1:length(par.dist[,1])){
          
          out = sim.met.MSIRN.age(burnin=burnin,
                                   sim_pop=sim_pop,
                                   yrs=yrs,
                                   ntyr=ntyr,
                                   s=s,
                                   beta=par.dist[j,1],
                                   recov=par.df$recov,
                                   mort=par.df$mort, #annual
                                   mort_juv=par.df$mort_juv, #annual
                                   adult_fec=par.df$adult_fec, #annual
                                   wane=par.df$wane,
                                   age.brk=age.brk,
                                   slope.wane=slope.wane,
                                   sigma=par.dist[j,2],
                                   model=par.df$model,
                                   recru = par.dist[j,3],
                                   mu.sick = par.df$mu_sick,
                                   mu.rec = par.df$mu_rec,
                                   add.inf.mort = FALSE,
                                   add.rec.mort = FALSE,
                                   N_stat="matAB")
          #each column is a new guess. each row is output by age and biweek
          preds[,j] <- out$seroprevalence
          
        }
        
        
      }
      }else if ( model == "MSIRNNS"){
      
      if(mu.sick >1 & mu.rec==1){
        par.df$mu_sick = ConfidenceInterval$value[length(out.lik$par)] 
        par.df$mu_sick_lci = ConfidenceInterval$lower[length(out.lik$par)] 
        par.df$mu_sick_uci = ConfidenceInterval$upper[length(out.lik$par)] 
        
        par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-1)] 
        par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
        par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
        
        par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-2)] 
        par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
        par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
        
        par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
        par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
        par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
        
        
        
        par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
        par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
        par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
        
        out1 = sim.met.MSIRNNS.age(burnin=burnin,
                                  sim_pop=sim_pop,
                                  yrs=yrs,
                                  ntyr=ntyr,
                                  s=s,
                                  beta=par.df$beta,
                                  recov=par.df$recov,
                                  mort=par.df$mort, #annual
                                  mort_juv=par.df$mort_juv, #annual
                                  adult_fec=par.df$adult_fec, #annual
                                  wane=par.df$wane,
                                  age.brk=age.brk,
                                  slope.wane=slope.wane,
                                  sigma=par.df$sigma,
                                  model=par.df$model,
                                  recru = par.df$recru,
                                  epsilon = par.df$epsilon,
                                  mu.sick = par.df$mu_sick,
                                  mu.rec = par.df$mu_rec,
                                  add.inf.mort = TRUE,
                                  add.rec.mort = FALSE,
                                  N_stat="matAB")
        
        for (j in 1:length(par.dist[,1])){
        
        out = sim.met.MSIRNNS.age(burnin=burnin,
                                  sim_pop=sim_pop,
                                  yrs=yrs,
                                  ntyr=ntyr,
                                  s=s,
                                  beta=par.dist[j,1],
                                  recov=par.df$recov,
                                  mort=par.df$mort, #annual
                                  mort_juv=par.df$mort_juv, #annual
                                  adult_fec=par.df$adult_fec, #annual
                                  wane=par.df$wane,
                                  age.brk=age.brk,
                                  slope.wane=slope.wane,
                                  sigma=par.dist[j,2],
                                  model=par.df$model,
                                  recru = par.dist[j,3],
                                  epsilon = par.dist[j,4],
                                  mu.sick = par.dist[j,5],
                                  mu.rec = par.df$mu_rec,
                                  add.inf.mort = TRUE,
                                  add.rec.mort = FALSE,
                                  N_stat="matAB")
        
        #each column is a new guess. each row is output by age and biweek
        preds[,j] <- out$seroprevalence
        
        }
        
      }else if(mu.sick ==1 & mu.rec >1){
        
        par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
        par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
        par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
        
        par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-1)] 
        par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
        par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
        
        par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-2)] 
        par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
        par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
        
        par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
        par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
        par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
        
        
        
        par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
        par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
        par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
        
        out1 = sim.met.MSIRNNS.age(burnin=burnin,
                                   sim_pop=sim_pop,
                                   yrs=yrs,
                                   ntyr=ntyr,
                                   s=s,
                                   beta=par.df$beta,
                                   recov=par.df$recov,
                                   mort=par.df$mort, #annual
                                   mort_juv=par.df$mort_juv, #annual
                                   adult_fec=par.df$adult_fec, #annual
                                   wane=par.df$wane,
                                   age.brk=age.brk,
                                   slope.wane=slope.wane,
                                   sigma=par.df$sigma,
                                   model=par.df$model,
                                   recru = par.df$recru,
                                   epsilon = par.df$epsilon,
                                   mu.sick = par.df$mu_sick,
                                   mu.rec = par.df$mu_rec,
                                   add.inf.mort = FALSE,
                                   add.rec.mort = TRUE,
                                   N_stat="matAB")
        
        for (j in 1:length(par.dist[,1])){
          
          out = sim.met.MSIRNNS.age(burnin=burnin,
                                    sim_pop=sim_pop,
                                    yrs=yrs,
                                    ntyr=ntyr,
                                    s=s,
                                    beta=par.dist[j,1],
                                    recov=par.df$recov,
                                    mort=par.df$mort, #annual
                                    mort_juv=par.df$mort_juv, #annual
                                    adult_fec=par.df$adult_fec, #annual
                                    wane=par.df$wane,
                                    age.brk=age.brk,
                                    slope.wane=slope.wane,
                                    sigma=par.dist[j,2],
                                    model=par.df$model,
                                    recru = par.dist[j,3],
                                    epsilon = par.dist[j,4],
                                    mu.sick = par.df$mu_sick,
                                    mu.rec = par.dist[j,5],
                                    add.inf.mort = FALSE,
                                    add.rec.mort = TRUE,
                                    N_stat="matAB")
          
          #each column is a new guess. each row is output by age and biweek
          preds[,j] <- out$seroprevalence
          
        }
        
      }else if (mu.sick >1 & mu.rec >1){
        par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
        par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
        par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
        
        par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
        par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
        par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
        
        par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-2)] 
        par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
        par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
        
        par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-3)] 
        par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
        par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
        
        par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
        par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
        par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
        
        
        
        par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
        par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
        par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
        
        out1 = sim.met.MSIRNNS.age(burnin=burnin,
                                   sim_pop=sim_pop,
                                   yrs=yrs,
                                   ntyr=ntyr,
                                   s=s,
                                   beta=par.df$beta,
                                   recov=par.df$recov,
                                   mort=par.df$mort, #annual
                                   mort_juv=par.df$mort_juv, #annual
                                   adult_fec=par.df$adult_fec, #annual
                                   wane=par.df$wane,
                                   age.brk=age.brk,
                                   slope.wane=slope.wane,
                                   sigma=par.df$sigma,
                                   model=par.df$model,
                                   recru = par.df$recru,
                                   epsilon = par.df$epsilon,
                                   mu.sick = par.df$mu_sick,
                                   mu.rec = par.df$mu_rec,
                                   add.inf.mort = TRUE,
                                   add.rec.mort = TRUE,
                                   N_stat="matAB")
        
        
        for (j in 1:length(par.dist[,1])){
          
          out = sim.met.MSIRNNS.age(burnin=burnin,
                                    sim_pop=sim_pop,
                                    yrs=yrs,
                                    ntyr=ntyr,
                                    s=s,
                                    beta=par.dist[j,1],
                                    recov=par.df$recov,
                                    mort=par.df$mort, #annual
                                    mort_juv=par.df$mort_juv, #annual
                                    adult_fec=par.df$adult_fec, #annual
                                    wane=par.df$wane,
                                    age.brk=age.brk,
                                    slope.wane=slope.wane,
                                    sigma=par.dist[j,2],
                                    model=par.df$model,
                                    recru = par.dist[j,3],
                                    epsilon = par.dist[j,4],
                                    mu.sick = par.dist[j,5],
                                    mu.rec = par.dist[j,6],
                                    add.inf.mort = TRUE,
                                    add.rec.mort = TRUE,
                                    N_stat="matAB")
          
          #each column is a new guess. each row is output by age and biweek
          preds[,j] <- out$seroprevalence
          
        }
        
      }else{
        
        par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par))] 
        par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par))] 
        par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par))] 
        
        par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-1)] 
        par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
        par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
        
        
        
        par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
        par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
        par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
        
        
        
        par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
        par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
        par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
        
        out1 = sim.met.MSIRNNS.age(burnin=burnin,
                                   sim_pop=sim_pop,
                                   yrs=yrs,
                                   ntyr=ntyr,
                                   s=s,
                                   beta=par.df$beta,
                                   recov=par.df$recov,
                                   mort=par.df$mort, #annual
                                   mort_juv=par.df$mort_juv, #annual
                                   adult_fec=par.df$adult_fec, #annual
                                   wane=par.df$wane,
                                   age.brk=age.brk,
                                   slope.wane=slope.wane,
                                   sigma=par.df$sigma,
                                   model=par.df$model,
                                   recru = par.df$recru,
                                   epsilon = par.df$epsilon,
                                   mu.sick = par.df$mu_sick,
                                   mu.rec = par.df$mu_rec,
                                   add.inf.mort = FALSE,
                                   add.rec.mort = FALSE,
                                   N_stat="matAB")
        
        for (j in 1:length(par.dist[,1])){
          
          out = sim.met.MSIRNNS.age(burnin=burnin,
                                    sim_pop=sim_pop,
                                    yrs=yrs,
                                    ntyr=ntyr,
                                    s=s,
                                    beta=par.dist[j,1],
                                    recov=par.df$recov,
                                    mort=par.df$mort, #annual
                                    mort_juv=par.df$mort_juv, #annual
                                    adult_fec=par.df$adult_fec, #annual
                                    wane=par.df$wane,
                                    age.brk=age.brk,
                                    slope.wane=slope.wane,
                                    sigma=par.dist[j,2],
                                    model=par.df$model,
                                    recru = par.dist[j,3],
                                    epsilon = par.dist[j,4],
                                    mu.sick = par.df$mu_sick,
                                    mu.rec = par.df$mu_rec,
                                    add.inf.mort = FALSE,
                                    add.rec.mort = FALSE,
                                    N_stat="matAB")
          
          #each column is a new guess. each row is output by age and biweek
          preds[,j] <- out$seroprevalence
          
        }
        
        
      }
      }else if (model == "MSIRNIRN" ){
        
        if(mu.sick >1 & mu.rec==1){
          par.df$mu_sick = ConfidenceInterval$value[length(out.lik$par)] 
          par.df$mu_sick_lci = ConfidenceInterval$lower[length(out.lik$par)] 
          par.df$mu_sick_uci = ConfidenceInterval$upper[length(out.lik$par)] 
          
          par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
          par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
          par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
          
          par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-2)] 
          par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
          par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
          
          par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
          par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
          par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
          
          
          
          par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
          par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
          par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
          
          out1  = sim.met.MSIRNIRN.age(burnin=burnin,
                                      sim_pop=sim_pop,
                                      yrs=yrs,
                                      ntyr=ntyr,
                                      s=s,
                                      beta=par.df$beta,
                                      age.brk = age.brk,
                                      recov=par.df$recov,
                                      mort=par.df$mort, 
                                      mort_juv=par.df$mort_juv,
                                      adult_fec=par.df$adult_fec, 
                                      wane=par.df$wane,
                                      recru=par.df$recru,
                                      slope.wane=slope.wane,
                                      sigma=par.df$sigma,
                                      sigma2=par.df$sigma2,
                                      model=par.df$model,
                                      mu.sick = par.df$mu_sick,
                                      mu.rec = par.df$mu_rec,
                                      add.inf.mort = TRUE,
                                      N_stat="matAB",
                                      add.rec.mort = FALSE) 
          
          for (j in 1:length(par.dist[,1])){
            
            out  = sim.met.MSIRNIRN.age(burnin=burnin,
                                        sim_pop=sim_pop,
                                        yrs=yrs,
                                        ntyr=ntyr,
                                        s=s,
                                        beta=par.dist[j,1],
                                        age.brk = age.brk,
                                        recov=par.df$recov,
                                        mort=par.df$mort, 
                                        mort_juv=par.df$mort_juv,
                                        adult_fec=par.df$adult_fec, 
                                        wane=par.df$wane,
                                        recru=par.dist[j,3],
                                        slope.wane=slope.wane,
                                        sigma=par.dist[j,2],
                                        sigma2=par.dist[j,4],
                                        model=par.df$model,
                                        mu.sick = par.dist[j,5],
                                        mu.rec = par.df$mu_rec,
                                        add.inf.mort = TRUE,
                                        N_stat="matAB",
                                        add.rec.mort = FALSE) 
            
            #each column is a new guess. each row is output by age and biweek
            preds[,j] <- out$seroprevalence
            
          }
          
        }else if(mu.sick ==1 & mu.rec >1){
          
          par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
          par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
          par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
          
          par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
          par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
          par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
          
          par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-2)] 
          par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
          par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
          
          par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
          par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
          par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
          
          
          
          par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
          par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
          par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
          
          out1  = sim.met.MSIRNIRN.age(burnin=burnin,
                                       sim_pop=sim_pop,
                                       yrs=yrs,
                                       ntyr=ntyr,
                                       s=s,
                                       beta=par.df$beta,
                                       age.brk = age.brk,
                                       recov=par.df$recov,
                                       mort=par.df$mort, 
                                       mort_juv=par.df$mort_juv,
                                       adult_fec=par.df$adult_fec, 
                                       wane=par.df$wane,
                                       recru=par.df$recru,
                                       slope.wane=slope.wane,
                                       sigma=par.df$sigma,
                                       sigma2=par.df$sigma2,
                                       model=par.df$model,
                                       mu.sick = par.df$mu_sick,
                                       mu.rec = par.df$mu_rec,
                                       add.inf.mort = FALSE,
                                       N_stat="matAB",
                                       add.rec.mort = TRUE) 
          
          for (j in 1:length(par.dist[,1])){
            
            out  = sim.met.MSIRNIRN.age(burnin=burnin,
                                        sim_pop=sim_pop,
                                        yrs=yrs,
                                        ntyr=ntyr,
                                        s=s,
                                        beta=par.dist[j,1],
                                        age.brk = age.brk,
                                        recov=par.df$recov,
                                        mort=par.df$mort, 
                                        mort_juv=par.df$mort_juv,
                                        adult_fec=par.df$adult_fec, 
                                        wane=par.df$wane,
                                        recru=par.dist[j,3],
                                        slope.wane=slope.wane,
                                        sigma=par.dist[j,2],
                                        sigma2=par.dist[j,4],
                                        model=par.df$model,
                                        mu.sick = par.df$mu_sick,
                                        mu.rec = par.dist[j,5],
                                        add.inf.mort = FALSE,
                                        N_stat="matAB",
                                        add.rec.mort = TRUE) 
            
            #each column is a new guess. each row is output by age and biweek
            preds[,j] <- out$seroprevalence
            
          }
          
        }else if (mu.sick >1 & mu.rec >1){
          par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
          par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
          par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
          
          par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
          par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
          par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
          
          par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
          par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
          par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
          
          par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-3)] 
          par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
          par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
          
          par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
          par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
          par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
          
          
          
          par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
          par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
          par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
          
          out1  = sim.met.MSIRNIRN.age(burnin=burnin,
                                       sim_pop=sim_pop,
                                       yrs=yrs,
                                       ntyr=ntyr,
                                       s=s,
                                       beta=par.df$beta,
                                       age.brk = age.brk,
                                       recov=par.df$recov,
                                       mort=par.df$mort, 
                                       mort_juv=par.df$mort_juv,
                                       adult_fec=par.df$adult_fec, 
                                       wane=par.df$wane,
                                       recru=par.df$recru,
                                       slope.wane=slope.wane,
                                       sigma=par.df$sigma,
                                       sigma2=par.df$sigma2,
                                       model=par.df$model,
                                       mu.sick = par.df$mu_sick,
                                       mu.rec = par.df$mu_rec,
                                       add.inf.mort = TRUE,
                                       N_stat="matAB",
                                       add.rec.mort = TRUE) 
          
          for (j in 1:length(par.dist[,1])){
            
            out  = sim.met.MSIRNIRN.age(burnin=burnin,
                                        sim_pop=sim_pop,
                                        yrs=yrs,
                                        ntyr=ntyr,
                                        s=s,
                                        beta=par.dist[j,1],
                                        age.brk = age.brk,
                                        recov=par.df$recov,
                                        mort=par.df$mort, 
                                        mort_juv=par.df$mort_juv,
                                        adult_fec=par.df$adult_fec, 
                                        wane=par.df$wane,
                                        recru=par.dist[j,3],
                                        slope.wane=slope.wane,
                                        sigma=par.dist[j,2],
                                        sigma2=par.dist[j,4],
                                        model=par.df$model,
                                        mu.sick =  par.dist[j,5],
                                        mu.rec = par.dist[j,6],
                                        add.inf.mort = TRUE,
                                        N_stat="matAB",
                                        add.rec.mort = TRUE) 
            
            #each column is a new guess. each row is output by age and biweek
            preds[,j] <- out$seroprevalence
            
          }
          
        }else{
          
          par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par))] 
          par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par))] 
          par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par))] 
          
          par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-1)] 
          par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
          par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
          
          
          
          par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
          par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
          par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
          
          
          
          par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
          par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
          par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
          
          out1  = sim.met.MSIRNIRN.age(burnin=burnin,
                                       sim_pop=sim_pop,
                                       yrs=yrs,
                                       ntyr=ntyr,
                                       s=s,
                                       beta=par.df$beta,
                                       age.brk = age.brk,
                                       recov=par.df$recov,
                                       mort=par.df$mort, 
                                       mort_juv=par.df$mort_juv,
                                       adult_fec=par.df$adult_fec, 
                                       wane=par.df$wane,
                                       recru=par.df$recru,
                                       slope.wane=slope.wane,
                                       sigma=par.df$sigma,
                                       sigma2=par.df$sigma2,
                                       model=par.df$model,
                                       mu.sick = par.df$mu_sick,
                                       mu.rec = par.df$mu_rec,
                                       add.inf.mort = FALSE,
                                       N_stat="matAB",
                                       add.rec.mort = FALSE) 
          
          for (j in 1:length(par.dist[,1])){
            
            out  = sim.met.MSIRNIRN.age(burnin=burnin,
                                        sim_pop=sim_pop,
                                        yrs=yrs,
                                        ntyr=ntyr,
                                        s=s,
                                        beta=par.dist[j,1],
                                        age.brk = age.brk,
                                        recov=par.df$recov,
                                        mort=par.df$mort, 
                                        mort_juv=par.df$mort_juv,
                                        adult_fec=par.df$adult_fec, 
                                        wane=par.df$wane,
                                        recru=par.dist[j,3],
                                        slope.wane=slope.wane,
                                        sigma=par.dist[j,2],
                                        sigma2=par.dist[j,4],
                                        model=par.df$model,
                                        mu.sick =  par.df$mu_sick,
                                        mu.rec = par.df$mu_rec,
                                        add.inf.mort = FALSE,
                                        N_stat="matAB",
                                        add.rec.mort = FALSE) 
            
            #each column is a new guess. each row is output by age and biweek
            preds[,j] <- out$seroprevalence
            
          }
          
          
        }}else if (model == "MSIRNIR"){
          
          if(mu.sick >1 & mu.rec==1){
            par.df$mu_sick = ConfidenceInterval$value[length(out.lik$par)] 
            par.df$mu_sick_lci = ConfidenceInterval$lower[length(out.lik$par)] 
            par.df$mu_sick_uci = ConfidenceInterval$upper[length(out.lik$par)] 
            
            par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
            par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
            par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
            
            par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-2)] 
            par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
            par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
            
            par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
            par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
            par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
            
            
            
            par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
            par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
            par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
            
            out1  = sim.met.MSIRNIR.age(burnin=burnin,
                                       sim_pop=sim_pop,
                                       yrs=yrs,
                                       ntyr=ntyr,
                                       s=s,
                                       beta=par.df$beta,
                                       age.brk = age.brk,
                                       recov=par.df$recov,
                                       mort=par.df$mort, 
                                       mort_juv=par.df$mort_juv,
                                       adult_fec=par.df$adult_fec, 
                                       wane=par.df$wane,
                                       recru=par.df$recru,
                                       slope.wane=slope.wane,
                                       sigma=par.df$sigma,
                                       sigma2=par.df$sigma2,
                                       model=par.df$model,
                                       mu.sick = par.df$mu_sick,
                                       mu.rec = par.df$mu_rec,
                                       add.inf.mort = TRUE,
                                       N_stat="matAB",
                                       add.rec.mort = FALSE) 
            
            for (j in 1:length(par.dist[,1])){
              
              out  = sim.met.MSIRNIR.age(burnin=burnin,
                                         sim_pop=sim_pop,
                                         yrs=yrs,
                                         ntyr=ntyr,
                                         s=s,
                                         beta=par.dist[j,1],
                                         age.brk = age.brk,
                                         recov=par.df$recov,
                                         mort=par.df$mort, 
                                         mort_juv=par.df$mort_juv,
                                         adult_fec=par.df$adult_fec, 
                                         wane=par.df$wane,
                                         recru=par.dist[j,3],
                                         slope.wane=slope.wane,
                                         sigma=par.dist[j,2],
                                         sigma2=par.dist[j,4],
                                         model=par.df$model,
                                         mu.sick = par.dist[j,5],
                                         mu.rec = par.df$mu_rec,
                                         add.inf.mort = TRUE,
                                         N_stat="matAB",
                                         add.rec.mort = FALSE) 
              
              
              #each column is a new guess. each row is output by age and biweek
              preds[,j] <- out$seroprevalence
              
            }
            
          }else if(mu.sick ==1 & mu.rec >1){
            
            par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
            par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
            par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
            
            par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
            par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
            par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
            
            par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-2)] 
            par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
            par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
            
            par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
            par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
            par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
            
            
            
            par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
            par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
            par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
            
            out1  = sim.met.MSIRNIR.age(burnin=burnin,
                                        sim_pop=sim_pop,
                                        yrs=yrs,
                                        ntyr=ntyr,
                                        s=s,
                                        beta=par.df$beta,
                                        age.brk = age.brk,
                                        recov=par.df$recov,
                                        mort=par.df$mort, 
                                        mort_juv=par.df$mort_juv,
                                        adult_fec=par.df$adult_fec, 
                                        wane=par.df$wane,
                                        recru=par.df$recru,
                                        slope.wane=slope.wane,
                                        sigma=par.df$sigma,
                                        sigma2=par.df$sigma2,
                                        model=par.df$model,
                                        mu.sick = par.df$mu_sick,
                                        mu.rec = par.df$mu_rec,
                                        add.inf.mort = FALSE,
                                        N_stat="matAB",
                                        add.rec.mort = TRUE)
            
            for (j in 1:length(par.dist[,1])){
              
              out  = sim.met.MSIRNIR.age(burnin=burnin,
                                         sim_pop=sim_pop,
                                         yrs=yrs,
                                         ntyr=ntyr,
                                         s=s,
                                         beta=par.dist[j,1],
                                         age.brk = age.brk,
                                         recov=par.df$recov,
                                         mort=par.df$mort, 
                                         mort_juv=par.df$mort_juv,
                                         adult_fec=par.df$adult_fec, 
                                         wane=par.df$wane,
                                         recru=par.dist[j,3],
                                         slope.wane=slope.wane,
                                         sigma=par.dist[j,2],
                                         sigma2=par.dist[j,4],
                                         model=par.df$model,
                                         mu.sick = par.df$mu_sick,
                                         mu.rec = par.dist[j,5],
                                         add.inf.mort = FALSE,
                                         N_stat="matAB",
                                         add.rec.mort = TRUE) 
              
              
              #each column is a new guess. each row is output by age and biweek
              preds[,j] <- out$seroprevalence
              
            }
            
          }else if (mu.sick >1 & mu.rec >1){
            par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
            par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
            par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
            
            par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
            par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
            par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
            
            par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
            par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
            par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
            
            par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-3)] 
            par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
            par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
            
            par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
            par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
            par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
            
            
            
            par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
            par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
            par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
            
            out1  = sim.met.MSIRNIR.age(burnin=burnin,
                                        sim_pop=sim_pop,
                                        yrs=yrs,
                                        ntyr=ntyr,
                                        s=s,
                                        beta=par.df$beta,
                                        age.brk = age.brk,
                                        recov=par.df$recov,
                                        mort=par.df$mort, 
                                        mort_juv=par.df$mort_juv,
                                        adult_fec=par.df$adult_fec, 
                                        wane=par.df$wane,
                                        recru=par.df$recru,
                                        slope.wane=slope.wane,
                                        sigma=par.df$sigma,
                                        sigma2=par.df$sigma2,
                                        model=par.df$model,
                                        mu.sick = par.df$mu_sick,
                                        mu.rec = par.df$mu_rec,
                                        add.inf.mort = TRUE,
                                        N_stat="matAB",
                                        add.rec.mort = TRUE)
            
            for (j in 1:length(par.dist[,1])){
              
              out  = sim.met.MSIRNIR.age(burnin=burnin,
                                         sim_pop=sim_pop,
                                         yrs=yrs,
                                         ntyr=ntyr,
                                         s=s,
                                         beta=par.dist[j,1],
                                         age.brk = age.brk,
                                         recov=par.df$recov,
                                         mort=par.df$mort, 
                                         mort_juv=par.df$mort_juv,
                                         adult_fec=par.df$adult_fec, 
                                         wane=par.df$wane,
                                         recru=par.dist[j,3],
                                         slope.wane=slope.wane,
                                         sigma=par.dist[j,2],
                                         sigma2=par.dist[j,4],
                                         model=par.df$model,
                                         mu.sick = par.dist[j,5],
                                         mu.rec = par.dist[j,6],
                                         add.inf.mort = TRUE,
                                         N_stat="matAB",
                                         add.rec.mort = TRUE) 
              
              
              #each column is a new guess. each row is output by age and biweek
              preds[,j] <- out$seroprevalence
              
            }
            
          }else{
            
            par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par))] 
            par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par))] 
            par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par))] 
            
            par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-1)] 
            par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
            par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
            
            
            
            par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
            par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
            par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
            
            
            
            par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
            par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
            par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
            
            out1  = sim.met.MSIRNIR.age(burnin=burnin,
                                        sim_pop=sim_pop,
                                        yrs=yrs,
                                        ntyr=ntyr,
                                        s=s,
                                        beta=par.df$beta,
                                        age.brk = age.brk,
                                        recov=par.df$recov,
                                        mort=par.df$mort, 
                                        mort_juv=par.df$mort_juv,
                                        adult_fec=par.df$adult_fec, 
                                        wane=par.df$wane,
                                        recru=par.df$recru,
                                        slope.wane=slope.wane,
                                        sigma=par.df$sigma,
                                        sigma2=par.df$sigma2,
                                        model=par.df$model,
                                        mu.sick = par.df$mu_sick,
                                        mu.rec = par.df$mu_rec,
                                        add.inf.mort = FALSE,
                                        N_stat="matAB",
                                        add.rec.mort = FALSE)
            
            for (j in 1:length(par.dist[,1])){
              
              out  = sim.met.MSIRNIR.age(burnin=burnin,
                                         sim_pop=sim_pop,
                                         yrs=yrs,
                                         ntyr=ntyr,
                                         s=s,
                                         beta=par.dist[j,1],
                                         age.brk = age.brk,
                                         recov=par.df$recov,
                                         mort=par.df$mort, 
                                         mort_juv=par.df$mort_juv,
                                         adult_fec=par.df$adult_fec, 
                                         wane=par.df$wane,
                                         recru=par.dist[j,3],
                                         slope.wane=slope.wane,
                                         sigma=par.dist[j,2],
                                         sigma2=par.dist[j,4],
                                         model=par.df$model,
                                         mu.sick = par.df$mu_sick,
                                         mu.rec = par.df$mu_rec,
                                         add.inf.mort = FALSE,
                                         N_stat="matAB",
                                         add.rec.mort = FALSE) 
              
              
              #each column is a new guess. each row is output by age and biweek
              preds[,j] <- out$seroprevalence
              
            }
            
            
          }
        }else if (model == "MSIRNSIRN"){
          
          if(mu.sick >1 & mu.rec==1){
            par.df$mu_sick = ConfidenceInterval$value[length(out.lik$par)] 
            par.df$mu_sick_lci = ConfidenceInterval$lower[length(out.lik$par)] 
            par.df$mu_sick_uci = ConfidenceInterval$upper[length(out.lik$par)] 
            
            par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
            par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
            par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
            
            
            par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
            par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
            par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
            
            par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-3)] 
            par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
            par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
            
            par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
            par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
            par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
            
            
            
            par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
            par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
            par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
            
            out1  = sim.met.MSIRNSIRN.age(burnin=burnin,
                                         sim_pop=sim_pop,
                                         yrs=yrs,
                                         ntyr=ntyr,
                                         s=s,
                                         beta=par.df$beta,
                                         age.brk = age.brk,
                                         recov=par.df$recov,
                                         mort=par.df$mort, 
                                         mort_juv=par.df$mort_juv,
                                         adult_fec=par.df$adult_fec, 
                                         wane=par.df$wane,
                                         recru=par.df$recru,
                                         slope.wane=slope.wane,
                                         sigma=par.df$sigma,
                                         sigma2=par.df$sigma2,
                                         beta2=par.df$beta2,
                                         model=par.df$model,
                                         mu.sick = par.df$mu_sick,
                                         mu.rec = par.df$mu_rec,
                                         add.inf.mort = TRUE,
                                         N_stat="matAB",
                                         add.rec.mort = FALSE) 
            
            for (j in 1:length(par.dist[,1])){
              
              
              out  = sim.met.MSIRNSIRN.age(burnin=burnin,
                                          sim_pop=sim_pop,
                                          yrs=yrs,
                                          ntyr=ntyr,
                                          s=s,
                                          beta=par.dist[j,1],
                                          age.brk = age.brk,
                                          recov=par.df$recov,
                                          mort=par.df$mort, 
                                          mort_juv=par.df$mort_juv,
                                          adult_fec=par.df$adult_fec, 
                                          wane=par.df$wane,
                                          recru=par.dist[j,3],
                                          slope.wane=slope.wane,
                                          sigma=par.dist[j,2],
                                          sigma2=par.dist[j,4],
                                          beta2=par.dist[j,5],
                                          model=par.df$model,
                                          mu.sick = par.dist[j,6],
                                          mu.rec = par.df$mu_rec,
                                          add.inf.mort = TRUE,
                                          N_stat="matAB",
                                          add.rec.mort = FALSE) 
              
              #each column is a new guess. each row is output by age and biweek
              preds[,j] <- out$seroprevalence
              
            }
            
          }else if(mu.sick ==1 & mu.rec >1){
            
            par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
            par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
            par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
            
            par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
            par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
            par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
            
            
            par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
            par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
            par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
            
            par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-3)] 
            par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
            par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
            
            par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
            par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
            par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
            
            
            
            par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
            par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
            par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
            
            out1  = sim.met.MSIRNSIRN.age(burnin=burnin,
                                          sim_pop=sim_pop,
                                          yrs=yrs,
                                          ntyr=ntyr,
                                          s=s,
                                          beta=par.df$beta,
                                          age.brk = age.brk,
                                          recov=par.df$recov,
                                          mort=par.df$mort, 
                                          mort_juv=par.df$mort_juv,
                                          adult_fec=par.df$adult_fec, 
                                          wane=par.df$wane,
                                          recru=par.df$recru,
                                          slope.wane=slope.wane,
                                          sigma=par.df$sigma,
                                          sigma2=par.df$sigma2,
                                          beta2=par.df$beta2,
                                          model=par.df$model,
                                          mu.sick = par.df$mu_sick,
                                          mu.rec = par.df$mu_rec,
                                          add.inf.mort = FALSE,
                                          N_stat="matAB",
                                          add.rec.mort = TRUE)
            
            for (j in 1:length(par.dist[,1])){
              
              
              out  = sim.met.MSIRNSIRN.age(burnin=burnin,
                                           sim_pop=sim_pop,
                                           yrs=yrs,
                                           ntyr=ntyr,
                                           s=s,
                                           beta=par.dist[j,1],
                                           age.brk = age.brk,
                                           recov=par.df$recov,
                                           mort=par.df$mort, 
                                           mort_juv=par.df$mort_juv,
                                           adult_fec=par.df$adult_fec, 
                                           wane=par.df$wane,
                                           recru=par.dist[j,3],
                                           slope.wane=slope.wane,
                                           sigma=par.dist[j,2],
                                           sigma2=par.dist[j,4],
                                           beta2=par.dist[j,5],
                                           model=par.df$model,
                                           mu.sick = par.df$mu_sick,
                                           mu.rec = par.dist[j,6],
                                           add.inf.mort = FALSE,
                                           N_stat="matAB",
                                           add.rec.mort = TRUE) 
              
              #each column is a new guess. each row is output by age and biweek
              preds[,j] <- out$seroprevalence
              
            }
            
          }else if (mu.sick >1 & mu.rec >1){
            par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
            par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
            par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
            
            par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
            par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
            par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
            
            par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
            par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
            par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
            
            
            par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
            par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
            par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
            
            par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-4)] 
            par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
            par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
            
            par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
            par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
            par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
            
            
            
            par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
            par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
            par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
            
            out1  = sim.met.MSIRNSIRN.age(burnin=burnin,
                                          sim_pop=sim_pop,
                                          yrs=yrs,
                                          ntyr=ntyr,
                                          s=s,
                                          beta=par.df$beta,
                                          age.brk = age.brk,
                                          recov=par.df$recov,
                                          mort=par.df$mort, 
                                          mort_juv=par.df$mort_juv,
                                          adult_fec=par.df$adult_fec, 
                                          wane=par.df$wane,
                                          recru=par.df$recru,
                                          slope.wane=slope.wane,
                                          sigma=par.df$sigma,
                                          sigma2=par.df$sigma2,
                                          beta2=par.df$beta2,
                                          model=par.df$model,
                                          mu.sick = par.df$mu_sick,
                                          mu.rec = par.df$mu_rec,
                                          add.inf.mort = TRUE,
                                          N_stat="matAB",
                                          add.rec.mort = TRUE)
            
            for (j in 1:length(par.dist[,1])){
              
              
              out  = sim.met.MSIRNSIRN.age(burnin=burnin,
                                           sim_pop=sim_pop,
                                           yrs=yrs,
                                           ntyr=ntyr,
                                           s=s,
                                           beta=par.dist[j,1],
                                           age.brk = age.brk,
                                           recov=par.df$recov,
                                           mort=par.df$mort, 
                                           mort_juv=par.df$mort_juv,
                                           adult_fec=par.df$adult_fec, 
                                           wane=par.df$wane,
                                           recru=par.dist[j,3],
                                           slope.wane=slope.wane,
                                           sigma=par.dist[j,2],
                                           sigma2=par.dist[j,4],
                                           beta2=par.dist[j,5],
                                           model=par.df$model,
                                           mu.sick = par.dist[j,6],
                                           mu.rec = par.dist[j,7],
                                           add.inf.mort = TRUE,
                                           N_stat="matAB",
                                           add.rec.mort = TRUE) 
              
              #each column is a new guess. each row is output by age and biweek
              preds[,j] <- out$seroprevalence
              
            }
            
          }else{
            
            par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par))] 
            par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par))] 
            par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par))] 
            
            
            par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
            par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
            par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
            
            par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-2)] 
            par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
            par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
            
            par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
            par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
            par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
            
            
            
            par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
            par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
            par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
            
            out1  = sim.met.MSIRNSIRN.age(burnin=burnin,
                                          sim_pop=sim_pop,
                                          yrs=yrs,
                                          ntyr=ntyr,
                                          s=s,
                                          beta=par.df$beta,
                                          age.brk = age.brk,
                                          recov=par.df$recov,
                                          mort=par.df$mort, 
                                          mort_juv=par.df$mort_juv,
                                          adult_fec=par.df$adult_fec, 
                                          wane=par.df$wane,
                                          recru=par.df$recru,
                                          slope.wane=slope.wane,
                                          sigma=par.df$sigma,
                                          sigma2=par.df$sigma2,
                                          beta2=par.df$beta2,
                                          model=par.df$model,
                                          mu.sick = par.df$mu_sick,
                                          mu.rec = par.df$mu_rec,
                                          add.inf.mort = FALSE,
                                          N_stat="matAB",
                                          add.rec.mort = FALSE)
            
            for (j in 1:length(par.dist[,1])){
              
              
              out  = sim.met.MSIRNSIRN.age(burnin=burnin,
                                           sim_pop=sim_pop,
                                           yrs=yrs,
                                           ntyr=ntyr,
                                           s=s,
                                           beta=par.dist[j,1],
                                           age.brk = age.brk,
                                           recov=par.df$recov,
                                           mort=par.df$mort, 
                                           mort_juv=par.df$mort_juv,
                                           adult_fec=par.df$adult_fec, 
                                           wane=par.df$wane,
                                           recru=par.dist[j,3],
                                           slope.wane=slope.wane,
                                           sigma=par.dist[j,2],
                                           sigma2=par.dist[j,4],
                                           beta2=par.dist[j,5],
                                           model=par.df$model,
                                           mu.sick = par.df$mu_sick,
                                           mu.rec = par.df$mu_rec,
                                           add.inf.mort = FALSE,
                                           N_stat="matAB",
                                           add.rec.mort = FALSE) 
              
              #each column is a new guess. each row is output by age and biweek
              preds[,j] <- out$seroprevalence
              
            }
            
            
          }
          }else if (model == "MSIRSIR"){
            
            if(mu.sick >1 & mu.rec==1){
              par.df$mu_sick = ConfidenceInterval$value[length(out.lik$par)] 
              par.df$mu_sick_lci = ConfidenceInterval$lower[length(out.lik$par)] 
              par.df$mu_sick_uci = ConfidenceInterval$upper[length(out.lik$par)] 
              
              par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
              par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
              par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
              
              par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
              par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
              par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
              
              
              
              par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
              par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
              par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
              
              out1  = sim.met.MSIRSIR.age(burnin=burnin,
                                         sim_pop=sim_pop,
                                         yrs=yrs,
                                         ntyr=ntyr,
                                         s=s,
                                         beta=par.df$beta,
                                         age.brk = age.brk,
                                         recov=par.df$recov,
                                         mort=par.df$mort, 
                                         mort_juv=par.df$mort_juv,
                                         adult_fec=par.df$adult_fec, 
                                         wane=par.df$wane,
                                         slope.wane=slope.wane,
                                         sigma=par.df$sigma,
                                         beta2=par.df$beta2,
                                         model=par.df$model,
                                         mu.sick = par.df$mu_sick,
                                         mu.rec = par.df$mu_rec,
                                         add.inf.mort = TRUE,
                                         add.rec.mort = FALSE) 
              
              for (j in 1:length(par.dist[,1])){
                
                out  = sim.met.MSIRSIR.age(burnin=burnin,
                                           sim_pop=sim_pop,
                                           yrs=yrs,
                                           ntyr=ntyr,
                                           s=s,
                                           beta=par.dist[j,1],
                                           age.brk = age.brk,
                                           recov=par.df$recov,
                                           mort=par.df$mort, 
                                           mort_juv=par.df$mort_juv,
                                           adult_fec=par.df$adult_fec, 
                                           wane=par.df$wane,
                                           slope.wane=slope.wane,
                                           sigma=par.dist[j,2],
                                           beta2=par.dist[j,3],
                                           model=par.df$model,
                                           mu.sick = par.dist[j,4],
                                           mu.rec = par.df$mu_rec,
                                           add.inf.mort = TRUE,
                                           add.rec.mort = FALSE) 
                
                
                
                #each column is a new guess. each row is output by age and biweek
                preds[,j] <- out$seroprevalence
                
              }
              
            }else if(mu.sick ==1 & mu.rec >1){
              
              par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
              par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
              par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
              
              par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
              par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
              par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
              
              par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
              par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
              par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
              
              
              par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
              par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
              par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
              
              out1  = sim.met.MSIRSIR.age(burnin=burnin,
                                          sim_pop=sim_pop,
                                          yrs=yrs,
                                          ntyr=ntyr,
                                          s=s,
                                          beta=par.df$beta,
                                          age.brk = age.brk,
                                          recov=par.df$recov,
                                          mort=par.df$mort, 
                                          mort_juv=par.df$mort_juv,
                                          adult_fec=par.df$adult_fec, 
                                          wane=par.df$wane,
                                          slope.wane=slope.wane,
                                          sigma=par.df$sigma,
                                          beta2=par.df$beta2,
                                          model=par.df$model,
                                          mu.sick = par.df$mu_sick,
                                          mu.rec = par.df$mu_rec,
                                          add.inf.mort = FALSE,
                                          add.rec.mort = TRUE)
              
              for (j in 1:length(par.dist[,1])){
                
                out  = sim.met.MSIRSIR.age(burnin=burnin,
                                           sim_pop=sim_pop,
                                           yrs=yrs,
                                           ntyr=ntyr,
                                           s=s,
                                           beta=par.dist[j,1],
                                           age.brk = age.brk,
                                           recov=par.df$recov,
                                           mort=par.df$mort, 
                                           mort_juv=par.df$mort_juv,
                                           adult_fec=par.df$adult_fec, 
                                           wane=par.df$wane,
                                           slope.wane=slope.wane,
                                           sigma=par.dist[j,2],
                                           beta2=par.dist[j,3],
                                           model=par.df$model,
                                           mu.sick = par.df$mu_sick,
                                           mu.rec = par.dist[j,4],
                                           add.inf.mort = FALSE,
                                           add.rec.mort = TRUE) 
                
                
                
                #each column is a new guess. each row is output by age and biweek
                preds[,j] <- out$seroprevalence
                
              }
              
              
            }else if (mu.sick >1 & mu.rec >1){
              par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
              par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
              par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
              
              par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
              par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
              par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
              
              par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
              par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
              par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
              
              par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
              par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
              par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
              
              
              
              par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
              par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
              par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
              
              out1  = sim.met.MSIRSIR.age(burnin=burnin,
                                          sim_pop=sim_pop,
                                          yrs=yrs,
                                          ntyr=ntyr,
                                          s=s,
                                          beta=par.df$beta,
                                          age.brk = age.brk,
                                          recov=par.df$recov,
                                          mort=par.df$mort, 
                                          mort_juv=par.df$mort_juv,
                                          adult_fec=par.df$adult_fec, 
                                          wane=par.df$wane,
                                          slope.wane=slope.wane,
                                          sigma=par.df$sigma,
                                          beta2=par.df$beta2,
                                          model=par.df$model,
                                          mu.sick = par.df$mu_sick,
                                          mu.rec = par.df$mu_rec,
                                          add.inf.mort = TRUE,
                                          add.rec.mort = TRUE)
              
              for (j in 1:length(par.dist[,1])){
                
                out  = sim.met.MSIRSIR.age(burnin=burnin,
                                           sim_pop=sim_pop,
                                           yrs=yrs,
                                           ntyr=ntyr,
                                           s=s,
                                           beta=par.dist[j,1],
                                           age.brk = age.brk,
                                           recov=par.df$recov,
                                           mort=par.df$mort, 
                                           mort_juv=par.df$mort_juv,
                                           adult_fec=par.df$adult_fec, 
                                           wane=par.df$wane,
                                           slope.wane=slope.wane,
                                           sigma=par.dist[j,2],
                                           beta2=par.dist[j,3],
                                           model=par.df$model,
                                           mu.sick = par.dist[j,4],
                                           mu.rec = par.dist[j,5],
                                           add.inf.mort = TRUE,
                                           add.rec.mort = TRUE) 
                
                
                
                #each column is a new guess. each row is output by age and biweek
                preds[,j] <- out$seroprevalence
                
              }
              
            }else{
              par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par))] 
              par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par))] 
              par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par))] 
              
              par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-1)] 
              par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
              par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
              
              
              
              par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-2)] 
              par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-2)]
              par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-2)]
              
              out1  = sim.met.MSIRSIR.age(burnin=burnin,
                                          sim_pop=sim_pop,
                                          yrs=yrs,
                                          ntyr=ntyr,
                                          s=s,
                                          beta=par.df$beta,
                                          age.brk = age.brk,
                                          recov=par.df$recov,
                                          mort=par.df$mort, 
                                          mort_juv=par.df$mort_juv,
                                          adult_fec=par.df$adult_fec, 
                                          wane=par.df$wane,
                                          slope.wane=slope.wane,
                                          sigma=par.df$sigma,
                                          beta2=par.df$beta2,
                                          model=par.df$model,
                                          mu.sick = par.df$mu_sick,
                                          mu.rec = par.df$mu_rec,
                                          add.inf.mort = FALSE,
                                          add.rec.mort = FALSE)
              
              for (j in 1:length(par.dist[,1])){
                
                out  = sim.met.MSIRSIR.age(burnin=burnin,
                                           sim_pop=sim_pop,
                                           yrs=yrs,
                                           ntyr=ntyr,
                                           s=s,
                                           beta=par.dist[j,1],
                                           age.brk = age.brk,
                                           recov=par.df$recov,
                                           mort=par.df$mort, 
                                           mort_juv=par.df$mort_juv,
                                           adult_fec=par.df$adult_fec, 
                                           wane=par.df$wane,
                                           slope.wane=slope.wane,
                                           sigma=par.dist[j,2],
                                           beta2=par.dist[j,3],
                                           model=par.df$model,
                                           mu.sick = par.df$mu_sick,
                                           mu.rec = par.df$mu_rec,
                                           add.inf.mort = FALSE,
                                           add.rec.mort = FALSE) 
                
                
                
                #each column is a new guess. each row is output by age and biweek
                preds[,j] <- out$seroprevalence
                
              }
              
              
            }
            }else if(model=="MSIRNIRN-LT"){
              
              
              if(mu.sick >1 & mu.rec==1){
                par.df$mu_sick = ConfidenceInterval$value[length(out.lik$par)] 
                par.df$mu_sick_lci = ConfidenceInterval$lower[length(out.lik$par)] 
                par.df$mu_sick_uci = ConfidenceInterval$upper[length(out.lik$par)] 
                
                par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                
                par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                
                par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
                par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
                par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
                
                out1  = sim.met.MSIRNIRN.LT.age(burnin=burnin,
                                               sim_pop=sim_pop,
                                               yrs=yrs,
                                               ntyr=ntyr,
                                               s=s,
                                               beta=par.df$beta,
                                               age.brk = age.brk,
                                               recov=par.df$recov,
                                               mort=par.df$mort, 
                                               mort_juv=par.df$mort_juv,
                                               adult_fec=par.df$adult_fec, 
                                               wane=par.df$wane,
                                               recru=par.df$recru,
                                               slope.wane=slope.wane,
                                               sigma=par.df$sigma,
                                               sigma2=par.df$sigma2,
                                               epsilon=par.df$epsilon,
                                               model=par.df$model,
                                               mu.sick = par.df$mu_sick,
                                               mu.rec = par.df$mu_rec,
                                               add.inf.mort = TRUE,
                                               add.rec.mort = FALSE,
                                               N_stat="matAB") 
                
                for (j in 1:length(par.dist[,1])){
                  
                  out  = sim.met.MSIRNIRN.LT.age(burnin=burnin,
                                                 sim_pop=sim_pop,
                                                 yrs=yrs,
                                                 ntyr=ntyr,
                                                 s=s,
                                                 beta=par.dist[j,1],
                                                 age.brk = age.brk,
                                                 recov=par.df$recov,
                                                 mort=par.df$mort, 
                                                 mort_juv=par.df$mort_juv,
                                                 adult_fec=par.df$adult_fec, 
                                                 wane=par.df$wane,
                                                 recru=par.dist[j,3],
                                                 slope.wane=slope.wane,
                                                 sigma=par.dist[j,2],
                                                 sigma2=par.dist[j,4],
                                                 epsilon=par.dist[j,5],
                                                 model=par.df$model,
                                                 mu.sick = par.dist[j,6],
                                                 mu.rec = par.df$mu_rec,
                                                 add.inf.mort = TRUE,
                                                 add.rec.mort = FALSE,
                                                 N_stat="matAB") 
                  
                 
                  
                  #each column is a new guess. each row is output by age and biweek
                  preds[,j] <- out$seroprevalence
                  
                }
                
              }else if(mu.sick ==1 & mu.rec >1){
                
                par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
                par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
                par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
                
                par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                
                par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                
                par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
                par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
                par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
                
                out1  = sim.met.MSIRNIRN.LT.age(burnin=burnin,
                                                sim_pop=sim_pop,
                                                yrs=yrs,
                                                ntyr=ntyr,
                                                s=s,
                                                beta=par.df$beta,
                                                age.brk = age.brk,
                                                recov=par.df$recov,
                                                mort=par.df$mort, 
                                                mort_juv=par.df$mort_juv,
                                                adult_fec=par.df$adult_fec, 
                                                wane=par.df$wane,
                                                recru=par.df$recru,
                                                slope.wane=slope.wane,
                                                sigma=par.df$sigma,
                                                sigma2=par.df$sigma2,
                                                epsilon=par.df$epsilon,
                                                model=par.df$model,
                                                mu.sick = par.df$mu_sick,
                                                mu.rec = par.df$mu_rec,
                                                add.inf.mort = FALSE,
                                                add.rec.mort = TRUE,
                                                N_stat="matAB")
                
                for (j in 1:length(par.dist[,1])){
                  
                  out  = sim.met.MSIRNIRN.LT.age(burnin=burnin,
                                                 sim_pop=sim_pop,
                                                 yrs=yrs,
                                                 ntyr=ntyr,
                                                 s=s,
                                                 beta=par.dist[j,1],
                                                 age.brk = age.brk,
                                                 recov=par.df$recov,
                                                 mort=par.df$mort, 
                                                 mort_juv=par.df$mort_juv,
                                                 adult_fec=par.df$adult_fec, 
                                                 wane=par.df$wane,
                                                 recru=par.dist[j,3],
                                                 slope.wane=slope.wane,
                                                 sigma=par.dist[j,2],
                                                 sigma2=par.dist[j,4],
                                                 epsilon=par.dist[j,5],
                                                 model=par.df$model,
                                                 mu.sick = par.df$mu_sick,
                                                 mu.rec = par.dist[j,6],
                                                 add.inf.mort = FALSE,
                                                 add.rec.mort = TRUE,
                                                 N_stat="matAB") 
                  
                  
                  
                  #each column is a new guess. each row is output by age and biweek
                  preds[,j] <- out$seroprevalence
                  
                }
                
              }else if (mu.sick >1 & mu.rec >1){
                par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
                par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
                par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
                
                par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                
                par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                
                par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                
                par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                
                
                
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
                par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
                par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
                
                out1  = sim.met.MSIRNIRN.LT.age(burnin=burnin,
                                                sim_pop=sim_pop,
                                                yrs=yrs,
                                                ntyr=ntyr,
                                                s=s,
                                                beta=par.df$beta,
                                                age.brk = age.brk,
                                                recov=par.df$recov,
                                                mort=par.df$mort, 
                                                mort_juv=par.df$mort_juv,
                                                adult_fec=par.df$adult_fec, 
                                                wane=par.df$wane,
                                                recru=par.df$recru,
                                                slope.wane=slope.wane,
                                                sigma=par.df$sigma,
                                                sigma2=par.df$sigma2,
                                                epsilon=par.df$epsilon,
                                                model=par.df$model,
                                                mu.sick = par.df$mu_sick,
                                                mu.rec = par.df$mu_rec,
                                                add.inf.mort = TRUE,
                                                add.rec.mort = TRUE,
                                                N_stat="matAB")
                
                for (j in 1:length(par.dist[,1])){
                  
                  out  = sim.met.MSIRNIRN.LT.age(burnin=burnin,
                                                 sim_pop=sim_pop,
                                                 yrs=yrs,
                                                 ntyr=ntyr,
                                                 s=s,
                                                 beta=par.dist[j,1],
                                                 age.brk = age.brk,
                                                 recov=par.df$recov,
                                                 mort=par.df$mort, 
                                                 mort_juv=par.df$mort_juv,
                                                 adult_fec=par.df$adult_fec, 
                                                 wane=par.df$wane,
                                                 recru=par.dist[j,3],
                                                 slope.wane=slope.wane,
                                                 sigma=par.dist[j,2],
                                                 sigma2=par.dist[j,4],
                                                 epsilon=par.dist[j,5],
                                                 model=par.df$model,
                                                 mu.sick = par.dist[j,6],
                                                 mu.rec = par.dist[j,7],
                                                 add.inf.mort = TRUE,
                                                 add.rec.mort = TRUE,
                                                 N_stat="matAB") 
                  
                  
                  
                  #each column is a new guess. each row is output by age and biweek
                  preds[,j] <- out$seroprevalence
                  
                }
                
              }else{
                
                par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par))] 
                par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par))] 
                par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par))] 
                
                par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                
                par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                
                
                
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                
                
                
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
                par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
                par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
                
                out1  = sim.met.MSIRNIRN.LT.age(burnin=burnin,
                                                sim_pop=sim_pop,
                                                yrs=yrs,
                                                ntyr=ntyr,
                                                s=s,
                                                beta=par.df$beta,
                                                age.brk = age.brk,
                                                recov=par.df$recov,
                                                mort=par.df$mort, 
                                                mort_juv=par.df$mort_juv,
                                                adult_fec=par.df$adult_fec, 
                                                wane=par.df$wane,
                                                recru=par.df$recru,
                                                slope.wane=slope.wane,
                                                sigma=par.df$sigma,
                                                sigma2=par.df$sigma2,
                                                epsilon=par.df$epsilon,
                                                model=par.df$model,
                                                mu.sick = par.df$mu_sick,
                                                mu.rec = par.df$mu_rec,
                                                add.inf.mort = FALSE,
                                                add.rec.mort = FALSE,
                                                N_stat="matAB")
                
                for (j in 1:length(par.dist[,1])){
                  
                  out  = sim.met.MSIRNIRN.LT.age(burnin=burnin,
                                                 sim_pop=sim_pop,
                                                 yrs=yrs,
                                                 ntyr=ntyr,
                                                 s=s,
                                                 beta=par.dist[j,1],
                                                 age.brk = age.brk,
                                                 recov=par.df$recov,
                                                 mort=par.df$mort, 
                                                 mort_juv=par.df$mort_juv,
                                                 adult_fec=par.df$adult_fec, 
                                                 wane=par.df$wane,
                                                 recru=par.dist[j,3],
                                                 slope.wane=slope.wane,
                                                 sigma=par.dist[j,2],
                                                 sigma2=par.dist[j,4],
                                                 epsilon=par.dist[j,5],
                                                 model=par.df$model,
                                                 mu.sick = par.df$mu_sick,
                                                 mu.rec = par.df$mu_rec,
                                                 add.inf.mort = FALSE,
                                                 add.rec.mort = FALSE,
                                                 N_stat="matAB") 
                  
                  
                  
                  #each column is a new guess. each row is output by age and biweek
                  preds[,j] <- out$seroprevalence
                  
                }
              }
              
            }else if(model=="MSIRNSIRN-LT"){
              
              if(mu.sick >1 & mu.rec==1){
                par.df$mu_sick = ConfidenceInterval$value[length(out.lik$par)] 
                par.df$mu_sick_lci = ConfidenceInterval$lower[length(out.lik$par)] 
                par.df$mu_sick_uci = ConfidenceInterval$upper[length(out.lik$par)] 
                
                par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                
                par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                
                
                par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                
                par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                
                
                
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
                par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
                par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
                
                
                out1  = sim.met.MSIRNSIRN.LT.age(burnin=burnin,
                                                sim_pop=sim_pop,
                                                yrs=yrs,
                                                ntyr=ntyr,
                                                s=s,
                                                beta=par.df$beta,
                                                age.brk = age.brk,
                                                recov=par.df$recov,
                                                mort=par.df$mort, 
                                                mort_juv=par.df$mort_juv,
                                                adult_fec=par.df$adult_fec, 
                                                wane=par.df$wane,
                                                recru=par.df$recru,
                                                slope.wane=slope.wane,
                                                sigma=par.df$sigma,
                                                sigma2=par.df$sigma2,
                                                beta2=par.df$beta2,
                                                epsilon=par.df$epsilon,
                                                model=model,
                                                mu.sick = par.df$mu_sick,
                                                mu.rec = par.df$mu_rec,
                                                add.inf.mort = TRUE,
                                                add.rec.mort = FALSE,
                                                N_stat="matAB") 
                
                for (j in 1:length(par.dist[,1])){
                  
                  out  = sim.met.MSIRNSIRN.LT.age(burnin=burnin,
                                                  sim_pop=sim_pop,
                                                  yrs=yrs,
                                                  ntyr=ntyr,
                                                  s=s,
                                                  beta=par.dist[j,1],
                                                  age.brk = age.brk,
                                                  recov=par.df$recov,
                                                  mort=par.df$mort, 
                                                  mort_juv=par.df$mort_juv,
                                                  adult_fec=par.df$adult_fec, 
                                                  wane=par.df$wane,
                                                  recru=par.dist[j,3],
                                                  slope.wane=slope.wane,
                                                  sigma=par.dist[j,2],
                                                  sigma2=par.dist[j,4],
                                                  beta2=par.dist[j,5],
                                                  epsilon=par.dist[j,6],
                                                  model=model,
                                                  mu.sick = par.dist[j,7],
                                                  mu.rec = par.df$mu_rec,
                                                  add.inf.mort = TRUE,
                                                  add.rec.mort = FALSE,
                                                  N_stat="matAB") 
                  
                  
                  
                  #each column is a new guess. each row is output by age and biweek
                  preds[,j] <- out$seroprevalence
                  
                }
                
              }else if(mu.sick ==1 & mu.rec >1){
                
                par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
                par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
                par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
                
                par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                
                par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                
                
                par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                
                par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                
                
                
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
                par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
                par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
                
                out1  = sim.met.MSIRNSIRN.LT.age(burnin=burnin,
                                                 sim_pop=sim_pop,
                                                 yrs=yrs,
                                                 ntyr=ntyr,
                                                 s=s,
                                                 beta=par.df$beta,
                                                 age.brk = age.brk,
                                                 recov=par.df$recov,
                                                 mort=par.df$mort, 
                                                 mort_juv=par.df$mort_juv,
                                                 adult_fec=par.df$adult_fec, 
                                                 wane=par.df$wane,
                                                 recru=par.df$recru,
                                                 slope.wane=slope.wane,
                                                 sigma=par.df$sigma,
                                                 sigma2=par.df$sigma2,
                                                 beta2=par.df$beta2,
                                                 epsilon=par.df$epsilon,
                                                 model=model,
                                                 mu.sick = par.df$mu_sick,
                                                 mu.rec = par.df$mu_rec,
                                                 add.inf.mort = FALSE,
                                                 add.rec.mort = TRUE,
                                                 N_stat="matAB")
                
                for (j in 1:length(par.dist[,1])){
                  
                  out  = sim.met.MSIRNSIRN.LT.age(burnin=burnin,
                                                  sim_pop=sim_pop,
                                                  yrs=yrs,
                                                  ntyr=ntyr,
                                                  s=s,
                                                  beta=par.dist[j,1],
                                                  age.brk = age.brk,
                                                  recov=par.df$recov,
                                                  mort=par.df$mort, 
                                                  mort_juv=par.df$mort_juv,
                                                  adult_fec=par.df$adult_fec, 
                                                  wane=par.df$wane,
                                                  recru=par.dist[j,3],
                                                  slope.wane=slope.wane,
                                                  sigma=par.dist[j,2],
                                                  sigma2=par.dist[j,4],
                                                  beta2=par.dist[j,5],
                                                  epsilon=par.dist[j,6],
                                                  model=model,
                                                  mu.sick = par.df$mu_sick,
                                                  mu.rec = par.dist[j,7],
                                                  add.inf.mort = FALSE,
                                                  add.rec.mort = TRUE,
                                                  N_stat="matAB") 
                  
                  
                  
                  #each column is a new guess. each row is output by age and biweek
                  preds[,j] <- out$seroprevalence
                  
                }
                
              }else if (mu.sick >1 & mu.rec >1){
                par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
                par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
                par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
                
                par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                
                par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)]
                
                par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                
                
                par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                
                par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-6)] 
                par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-6)] 
                par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-6)] 
                
                
                
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-7)] 
                par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-7)]
                par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-7)]
                
                out1  = sim.met.MSIRNSIRN.LT.age(burnin=burnin,
                                                 sim_pop=sim_pop,
                                                 yrs=yrs,
                                                 ntyr=ntyr,
                                                 s=s,
                                                 beta=par.df$beta,
                                                 age.brk = age.brk,
                                                 recov=par.df$recov,
                                                 mort=par.df$mort, 
                                                 mort_juv=par.df$mort_juv,
                                                 adult_fec=par.df$adult_fec, 
                                                 wane=par.df$wane,
                                                 recru=par.df$recru,
                                                 slope.wane=slope.wane,
                                                 sigma=par.df$sigma,
                                                 sigma2=par.df$sigma2,
                                                 beta2=par.df$beta2,
                                                 epsilon=par.df$epsilon,
                                                 model=model,
                                                 mu.sick = par.df$mu_sick,
                                                 mu.rec = par.df$mu_rec,
                                                 add.inf.mort = TRUE,
                                                 add.rec.mort = TRUE,
                                                 N_stat="matAB")
                
                for (j in 1:length(par.dist[,1])){
                  
                  out  = sim.met.MSIRNSIRN.LT.age(burnin=burnin,
                                                  sim_pop=sim_pop,
                                                  yrs=yrs,
                                                  ntyr=ntyr,
                                                  s=s,
                                                  beta=par.dist[j,1],
                                                  age.brk = age.brk,
                                                  recov=par.df$recov,
                                                  mort=par.df$mort, 
                                                  mort_juv=par.df$mort_juv,
                                                  adult_fec=par.df$adult_fec, 
                                                  wane=par.df$wane,
                                                  recru=par.dist[j,3],
                                                  slope.wane=slope.wane,
                                                  sigma=par.dist[j,2],
                                                  sigma2=par.dist[j,4],
                                                  beta2=par.dist[j,5],
                                                  epsilon=par.dist[j,6],
                                                  model=model,
                                                  mu.sick = par.dist[j,7],
                                                  mu.rec = par.dist[j,8],
                                                  add.inf.mort = TRUE,
                                                  add.rec.mort = TRUE,
                                                  N_stat="matAB") 
                  
                  
                  
                  #each column is a new guess. each row is output by age and biweek
                  preds[,j] <- out$seroprevalence
                  
                }
                
              }else{
                
                par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par))] 
                par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par))] 
                par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par))]
                
                
                par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                
                
                par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                
                par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                
                
                
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
                par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
                par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
                
                out1  = sim.met.MSIRNSIRN.LT.age(burnin=burnin,
                                                 sim_pop=sim_pop,
                                                 yrs=yrs,
                                                 ntyr=ntyr,
                                                 s=s,
                                                 beta=par.df$beta,
                                                 age.brk = age.brk,
                                                 recov=par.df$recov,
                                                 mort=par.df$mort, 
                                                 mort_juv=par.df$mort_juv,
                                                 adult_fec=par.df$adult_fec, 
                                                 wane=par.df$wane,
                                                 recru=par.df$recru,
                                                 slope.wane=slope.wane,
                                                 sigma=par.df$sigma,
                                                 sigma2=par.df$sigma2,
                                                 beta2=par.df$beta2,
                                                 epsilon=par.df$epsilon,
                                                 model=model,
                                                 mu.sick = par.df$mu_sick,
                                                 mu.rec = par.df$mu_rec,
                                                 add.inf.mort = FALSE,
                                                 add.rec.mort = FALSE,
                                                 N_stat="matAB")
                
                for (j in 1:length(par.dist[,1])){
                  
                  out  = sim.met.MSIRNSIRN.LT.age(burnin=burnin,
                                                  sim_pop=sim_pop,
                                                  yrs=yrs,
                                                  ntyr=ntyr,
                                                  s=s,
                                                  beta=par.dist[j,1],
                                                  age.brk = age.brk,
                                                  recov=par.df$recov,
                                                  mort=par.df$mort, 
                                                  mort_juv=par.df$mort_juv,
                                                  adult_fec=par.df$adult_fec, 
                                                  wane=par.df$wane,
                                                  recru=par.dist[j,3],
                                                  slope.wane=slope.wane,
                                                  sigma=par.dist[j,2],
                                                  sigma2=par.dist[j,4],
                                                  beta2=par.dist[j,5],
                                                  epsilon=par.dist[j,6],
                                                  model=model,
                                                  mu.sick = par.df$mu_sick,
                                                  mu.rec = par.df$mu_rec,
                                                  add.inf.mort = FALSE,
                                                  add.rec.mort = FALSE,
                                                  N_stat="matAB") 
                  
                  
                  
                  #each column is a new guess. each row is output by age and biweek
                  preds[,j] <- out$seroprevalence
                  
                }
                
                
              }
            }else if (model =="MSIRNSIRNS"){
              if(mu.sick >1 & mu.rec==1){
                par.df$mu_sick = ConfidenceInterval$value[length(out.lik$par)] 
                par.df$mu_sick_lci = ConfidenceInterval$lower[length(out.lik$par)] 
                par.df$mu_sick_uci = ConfidenceInterval$upper[length(out.lik$par)] 
                
                par.df$recru2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                par.df$recru2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                par.df$recru2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                
                par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                
                
                par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                
                par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                
                
                
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
                par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
                par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
                
                
                out1  = sim.met.MSIRNSIRNS.age(burnin=burnin,
                                              sim_pop=sim_pop,
                                              yrs=yrs,
                                              ntyr=ntyr,
                                              s=s,
                                              beta=par.df$beta,
                                              age.brk = age.brk,
                                              recov=par.df$recov,
                                              mort=par.df$mort, 
                                              mort_juv=par.df$mort_juv,
                                              adult_fec=par.df$adult_fec, 
                                              wane=par.df$wane,
                                              recru=par.df$recru,
                                              slope.wane=slope.wane,
                                              sigma=par.df$sigma,
                                              sigma2=par.df$sigma2,
                                              beta2=par.df$beta2,
                                              recru2=par.df$recru2,
                                              model=model,
                                              mu.sick = par.df$mu_sick,
                                              mu.rec = par.df$mu_rec,
                                              add.inf.mort = TRUE,
                                              add.rec.mort = FALSE,
                                              N_stat="matAB") 
                
                for (j in 1:length(par.dist[,1])){
                  
                  out  = sim.met.MSIRNSIRNS.age(burnin=burnin,
                                                sim_pop=sim_pop,
                                                yrs=yrs,
                                                ntyr=ntyr,
                                                s=s,
                                                beta=par.dist[j,1],
                                                age.brk = age.brk,
                                                recov=par.df$recov,
                                                mort=par.df$mort, 
                                                mort_juv=par.df$mort_juv,
                                                adult_fec=par.df$adult_fec, 
                                                wane=par.df$wane,
                                                recru=par.dist[j,3],
                                                slope.wane=slope.wane,
                                                sigma=par.dist[j,2],
                                                sigma2=par.dist[j,4],
                                                beta2=par.dist[j,5],
                                                recru2=par.dist[j,6],
                                                model=model,
                                                mu.sick = par.dist[j,7],
                                                mu.rec = par.df$mu_rec,
                                                add.inf.mort = TRUE,
                                                add.rec.mort = FALSE,
                                                N_stat="matAB") 
                  
                  
                  #each column is a new guess. each row is output by age and biweek
                  preds[,j] <- out$seroprevalence
                  
                }
                
              }else if(mu.sick ==1 & mu.rec >1){
                
                par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
                par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
                par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
                
                par.df$recru2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                par.df$recru2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                par.df$recru2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                
                par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                
                
                par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                
                par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                
                
                
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
                par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
                par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
                
                out1  = sim.met.MSIRNSIRNS.age(burnin=burnin,
                                               sim_pop=sim_pop,
                                               yrs=yrs,
                                               ntyr=ntyr,
                                               s=s,
                                               beta=par.df$beta,
                                               age.brk = age.brk,
                                               recov=par.df$recov,
                                               mort=par.df$mort, 
                                               mort_juv=par.df$mort_juv,
                                               adult_fec=par.df$adult_fec, 
                                               wane=par.df$wane,
                                               recru=par.df$recru,
                                               slope.wane=slope.wane,
                                               sigma=par.df$sigma,
                                               sigma2=par.df$sigma2,
                                               beta2=par.df$beta2,
                                               recru2=par.df$recru2,
                                               model=model,
                                               mu.sick = par.df$mu_sick,
                                               mu.rec = par.df$mu_rec,
                                               add.inf.mort = FALSE,
                                               add.rec.mort = TRUE,
                                               N_stat="matAB") 
                
                for (j in 1:length(par.dist[,1])){
                  
                  out  = sim.met.MSIRNSIRNS.age(burnin=burnin,
                                                sim_pop=sim_pop,
                                                yrs=yrs,
                                                ntyr=ntyr,
                                                s=s,
                                                beta=par.dist[j,1],
                                                age.brk = age.brk,
                                                recov=par.df$recov,
                                                mort=par.df$mort, 
                                                mort_juv=par.df$mort_juv,
                                                adult_fec=par.df$adult_fec, 
                                                wane=par.df$wane,
                                                recru=par.dist[j,3],
                                                slope.wane=slope.wane,
                                                sigma=par.dist[j,2],
                                                sigma2=par.dist[j,4],
                                                beta2=par.dist[j,5],
                                                recru2=par.dist[j,6],
                                                model=model,
                                                mu.sick = par.df$mu_sick,
                                                mu.rec = par.dist[j,7],
                                                add.inf.mort = FALSE,
                                                add.rec.mort = TRUE,
                                                N_stat="matAB") 
                  
                  
                  #each column is a new guess. each row is output by age and biweek
                  preds[,j] <- out$seroprevalence
                  
                }
                
              }else if (mu.sick >1 & mu.rec >1){
                par.df$mu_rec = ConfidenceInterval$value[length(out.lik$par)] 
                par.df$mu_rec_lci = ConfidenceInterval$lower[length(out.lik$par)] 
                par.df$mu_rec_uci = ConfidenceInterval$upper[length(out.lik$par)] 
                
                par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                
                par.df$recru2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                par.df$recru2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                par.df$recru2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)]
                
                par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                
                
                par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                
                par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-6)] 
                par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-6)] 
                par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-6)] 
                
                
                
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-7)] 
                par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-7)]
                par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-7)]
                
                out1  = sim.met.MSIRNSIRNS.age(burnin=burnin,
                                               sim_pop=sim_pop,
                                               yrs=yrs,
                                               ntyr=ntyr,
                                               s=s,
                                               beta=par.df$beta,
                                               age.brk = age.brk,
                                               recov=par.df$recov,
                                               mort=par.df$mort, 
                                               mort_juv=par.df$mort_juv,
                                               adult_fec=par.df$adult_fec, 
                                               wane=par.df$wane,
                                               recru=par.df$recru,
                                               slope.wane=slope.wane,
                                               sigma=par.df$sigma,
                                               sigma2=par.df$sigma2,
                                               beta2=par.df$beta2,
                                               recru2=par.df$recru2,
                                               model=model,
                                               mu.sick = par.df$mu_sick,
                                               mu.rec = par.df$mu_rec,
                                               add.inf.mort = TRUE,
                                               add.rec.mort = TRUE,
                                               N_stat="matAB") 
                
                for (j in 1:length(par.dist[,1])){
                  
                  out  = sim.met.MSIRNSIRNS.age(burnin=burnin,
                                                sim_pop=sim_pop,
                                                yrs=yrs,
                                                ntyr=ntyr,
                                                s=s,
                                                beta=par.dist[j,1],
                                                age.brk = age.brk,
                                                recov=par.df$recov,
                                                mort=par.df$mort, 
                                                mort_juv=par.df$mort_juv,
                                                adult_fec=par.df$adult_fec, 
                                                wane=par.df$wane,
                                                recru=par.dist[j,3],
                                                slope.wane=slope.wane,
                                                sigma=par.dist[j,2],
                                                sigma2=par.dist[j,4],
                                                beta2=par.dist[j,5],
                                                recru2=par.dist[j,6],
                                                model=model,
                                                mu.sick = par.dist[j,7],
                                                mu.rec = par.dist[j,8],
                                                add.inf.mort = TRUE,
                                                add.rec.mort = TRUE,
                                                N_stat="matAB") 
                  
                  
                  #each column is a new guess. each row is output by age and biweek
                  preds[,j] <- out$seroprevalence
                  
                }
                
              }else{
                
                par.df$recru2 = ConfidenceInterval$value[(length(out.lik$par))] 
                par.df$recru2_lci = ConfidenceInterval$lower[(length(out.lik$par))] 
                par.df$recru2_uci = ConfidenceInterval$upper[(length(out.lik$par))]
                
                
                par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                
                
                par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                
                par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                
                
                
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
                par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
                par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
                
                out1  = sim.met.MSIRNSIRNS.age(burnin=burnin,
                                               sim_pop=sim_pop,
                                               yrs=yrs,
                                               ntyr=ntyr,
                                               s=s,
                                               beta=par.df$beta,
                                               age.brk = age.brk,
                                               recov=par.df$recov,
                                               mort=par.df$mort, 
                                               mort_juv=par.df$mort_juv,
                                               adult_fec=par.df$adult_fec, 
                                               wane=par.df$wane,
                                               recru=par.df$recru,
                                               slope.wane=slope.wane,
                                               sigma=par.df$sigma,
                                               sigma2=par.df$sigma2,
                                               beta2=par.df$beta2,
                                               recru2=par.df$recru2,
                                               model=model,
                                               mu.sick = par.df$mu_sick,
                                               mu.rec = par.df$mu_rec,
                                               add.inf.mort = FALSE,
                                               add.rec.mort = FALSE,
                                               N_stat="matAB") 
                
                for (j in 1:length(par.dist[,1])){
                  
                  out  = sim.met.MSIRNSIRNS.age(burnin=burnin,
                                                sim_pop=sim_pop,
                                                yrs=yrs,
                                                ntyr=ntyr,
                                                s=s,
                                                beta=par.dist[j,1],
                                                age.brk = age.brk,
                                                recov=par.df$recov,
                                                mort=par.df$mort, 
                                                mort_juv=par.df$mort_juv,
                                                adult_fec=par.df$adult_fec, 
                                                wane=par.df$wane,
                                                recru=par.dist[j,3],
                                                slope.wane=slope.wane,
                                                sigma=par.dist[j,2],
                                                sigma2=par.dist[j,4],
                                                beta2=par.dist[j,5],
                                                recru2=par.dist[j,6],
                                                model=model,
                                                mu.sick = par.df$mu_sick,
                                                mu.rec = par.df$mu_rec,
                                                add.inf.mort = FALSE,
                                                add.rec.mort = FALSE,
                                                N_stat="matAB") 
                  
                  
                  #each column is a new guess. each row is output by age and biweek
                  preds[,j] <- out$seroprevalence
                  
                }
                
                
              }
            }
  
  
  
  preds[preds<0] <- 0
  
  pred.qs  <- apply(X=preds,MARGIN=1,FUN=quantile,probs=c(0.025,0.975), na.rm=TRUE) 
  pred.qs.iqr  <- apply(X=preds,MARGIN=1,FUN=quantile,probs=c(0.25,0.75), na.rm=TRUE) 
  
  #and combine to produce the time series:
  out1$seroprev_lci = pred.qs[1,]
  out1$seroprev_uci = pred.qs[2,]
  
  out1$seroprev_lci_25 = pred.qs.iqr[1,]
  out1$seroprev_uci_75 = pred.qs.iqr[2,]
  
  out1$model <- par.df$model
  out1$cutoff <-  cutoff1
  out1$antigen <- antigen1
  out1$species <- species1
  
  return(out1)  
  
  
  
}


slurm_array_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

par.fit <- read.csv(file =  paste0("fit-pars-", slurm_array_id, ".csv"), header = T, stringsAsFactors = F)




out.fit <- wrap.fit.CI(model=par.fit$model,
                              dat1=dat,
                              species1 = "Eidolon dupreanum",
                              antigen1=par.fit$antigen,
                              cutoff1 = par.fit$cutoff,
                              burnin=20,
                              sim_pop=10000,
                              yrs=50,
                              ntyr=26,
                              s=20,
                              beta = par.fit$beta,
                              par.df = par.fit,
                              age.brk = 20,
                              slope.wane=1,
                              sigma=par.fit$sigma, #annual rate of loss of adult Abs
                              mu.sick = par.fit$mu_sick, #multiplier x probability annual mortality in infected individuals
                              mu.rec = par.fit$mu_rec, #multiplier x probability annual mortality in recovered individuals
                              recru = par.fit$recru, #annual rate of return to infectiousness
                              epsilon = par.fit$epsilon, #annual rate of return to longterm AB- immunity
                              sigma2 = par.fit$sigma2,
                              beta2 = par.fit$beta2,
                              recru2 = par.fit$recru2,
                              N_stat="matAB",
                              do.plot=FALSE,
                              do.save=FALSE,
                              filename = NA)

filename = paste0(unique(par.fit$model), "-", unique(par.fit$antigen), "-", unique(par.fit$cutoff), slurm_array_id, ".Rdata")
save(out.fit, file=filename)
