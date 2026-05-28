rm(list=ls())
#set wd
#.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4")

library(plyr)
library(dplyr)   
library(epitools)
library(deSolve)
library(lubridate)
library(matrixcalc)
library(Matrix)
library(mvtnorm)
library(tidyr)
library(scales)
library(corpcor)

#load data
#homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/model_fitting/GhV-fits/run-CIs-boot"
#setwd(homewd)
#dat <- read.csv(file =  paste0(homewd,"/data/age-seroprevalence-data.csv"), header = T, stringsAsFactors = F)
dat <- read.csv(file =  "age-seroprevalence-data.csv", header = T, stringsAsFactors = F)
sero.dat <- read.csv(file =  "biweekly-serology-data.csv", header = T, stringsAsFactors = F)
shed.dat <- read.csv(file =  "biweekly-shedding-data.csv", header = T, stringsAsFactors = F)
head(dat)
head(sero.dat)
head(shed.dat)

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
  #age.mat.dat <- c( c(1:length(age.mat)), unlist(age.mat))
  #names(age.mat.dat) <- c("age", "pop")
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
make.prop <- function(df){
  tot = df[,ncol(df)]
  df[,2:(ncol(df)-1)] <- df[,2:(ncol(df)-1)]/tot
  return(df)
}

#build functions
buildTMat_all_age <- function(c, model, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	beta, beta2, recru, recov, sigma, sigma2, loss, loss2, wane,wane.fix,sigma.fix,epsilon, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
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
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  if (length(wane.fix)==1) wane.fix <- rep(wane.fix,nage)
  if (length(sigma2)==1) sigma2 <- rep(sigma2,nage)
  if (length(loss2)==1) loss2 <- rep(loss2,nage)
  if (length(loss)==1) loss <- rep(loss,nage)
  if (length(epsilon)==1) epsilon <- rep(epsilon,nage)
  
  
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
  
  
  
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  
  #I1 is in the 5th place on all
  #I2 is in 10 or 11 depending on the model
  
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,5] #always. for all model forms, I1 comes #5
  if(model=="MSIRSIR"){
    I_m2 = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,10] #I2 at #10
  }else if (model == "MSIRNSIRN" | model == "MSIRNSIRN-LT" | model == "MSIRNSIRNS"){
    I_m2 = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,11] #I2 at #11
  }else{
    I_m2 = 0
  }
  
  foi = 1-exp(-beta*((sum(I_m) + sum(I_m2))/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  foi2 = 1-exp(-beta2*((sum(I_m) + sum(I_m2))/sum(Npop_epi))) #uses different beta because susceptibles in the second class might be less suscetible
  
  mat1 <- matrix(0,c,c) # 
  
  Tmat <- matrix(0,c*nage,c*nage) #MSIRN for every age cohort
  
  for (j in 1:nage) { 
    
  if(model =="MSIRS"){
      
      #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
      mat1[] <- 0
      mat1[1,1] <- 1- wane.fix[j] 
      mat1[2,1] <- wane.fix[j] #M1 to M2
      mat1[2,2] <- 1- wane.fix[j]  
      mat1[3,2] <- wane.fix[j] #M2 to M3
      mat1[3,3] <- 1 - waning.maternal[j]
      mat1[4,3] <- waning.maternal[j] #M3 to S
      mat1[4,4] <- 1-foi[j]
      mat1[5,4] <- foi[j] #S to I
      mat1[5,5] <- 1-recov
      mat1[6,5] <- recov #I to R1
      mat1[6,6] <- 1-sigma.fix   
      mat1[7,6] <- sigma.fix #R1 to R2     
      mat1[7,7] <- 1-sigma.fix   
      mat1[8,7] <- sigma.fix #R2 to R3
      mat1[8,8] <- 1-sigma[j]
      mat1[4,8] <- sigma[j] #R3 to S   
      
      
    }else if(model =="MSIRN" | model=="MSIRNI"){
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- wane.fix[j] 
    mat1[2,1] <- wane.fix[j] #M1 to M2
    mat1[2,2] <- 1- wane.fix[j]  
    mat1[3,2] <- wane.fix[j] #M2 to M3
    mat1[3,3] <- 1 - waning.maternal[j]
    mat1[4,3] <- waning.maternal[j] #M3 to S
    mat1[4,4] <- 1-foi[j]
    mat1[5,4] <- foi[j] #S to I
    mat1[5,5] <- 1-recov
    mat1[6,5] <- recov #I to R1
    mat1[6,6] <- 1-sigma.fix   
    mat1[7,6] <- sigma.fix #R1 to R2     
    mat1[7,7] <- 1-sigma.fix   
    mat1[8,7] <- sigma.fix #R2 to R3
    mat1[8,8] <- 1-sigma[j]
    mat1[9,8] <- sigma[j] #R3 to N   
    mat1[9,9] <- 1- recru[j] #once in N, stay there, unless you recrudesce   
    mat1[5,9] <- recru[j] #recrudescence. moving from N back to I
      
    } else if (model=="MSIRNS"){
    
      mat1[] <- 0
      mat1[1,1] <- 1- wane.fix[j] 
      mat1[2,1] <- wane.fix[j] #M1 to M2
      mat1[2,2] <- 1- wane.fix[j]  
      mat1[3,2] <- wane.fix[j] #M2 to M3
      mat1[3,3] <- 1 - waning.maternal[j]
      mat1[4,3] <- waning.maternal[j] #M3 to S
      mat1[4,4] <- 1-foi[j]
      mat1[5,4] <- foi[j] #S to I
      mat1[5,5] <- 1-recov
      mat1[6,5] <- recov #I to R1
      mat1[6,6] <- 1-sigma.fix   
      mat1[7,6] <- sigma.fix #R1 to R2     
      mat1[7,7] <- 1-sigma.fix   
      mat1[8,7] <- sigma.fix #R2 to R3
      mat1[8,8] <- 1-sigma[j]
      mat1[9,8] <- sigma[j] #R3 to N   
      mat1[9,9] <- 1- loss[j] #once in N, stay there, unless you lose immunity
      mat1[4,9] <- loss[j] #loss of immunity. moving from N back to S
      
    } else if (model=="MSIRNNI"){
      
      mat1[] <- 0
      mat1[1,1] <- 1- wane.fix[j] 
      mat1[2,1] <- wane.fix[j] #M1 to M2
      mat1[2,2] <- 1- wane.fix[j]  
      mat1[3,2] <- wane.fix[j] #M2 to M3
      mat1[3,3] <- 1 - waning.maternal[j]
      mat1[4,3] <- waning.maternal[j] #M3 to S
      mat1[4,4] <- 1-foi[j]
      mat1[5,4] <- foi[j] #S to I
      mat1[5,5] <- 1-recov
      mat1[6,5] <- recov #I to R1
      mat1[6,6] <- 1-sigma.fix   
      mat1[7,6] <- sigma.fix #R1 to R2     
      mat1[7,7] <- 1-sigma.fix   
      mat1[8,7] <- sigma.fix #R2 to R3
      mat1[8,8] <- 1-sigma[j]
      mat1[9,8] <- sigma[j] #R3 to N   
      mat1[9,9] <- 1- recru[j] - epsilon[j] #once in N, stay there, unless you recrudesce to I or advance to N2  
      mat1[10,9] <- epsilon[j] #N1 to N2
      mat1[10,10] <- 1 #once in N2 you stay there
      mat1[5,9] <- recru[j] #recrudescence. moving from N1 back to I
      
      
    } else if(model == "MSIRNNS"){
    
      mat1[] <- 0
      mat1[1,1] <- 1- wane.fix[j] 
      mat1[2,1] <- wane.fix[j] #M1 to M2
      mat1[2,2] <- 1- wane.fix[j]  
      mat1[3,2] <- wane.fix[j] #M2 to M3
      mat1[3,3] <- 1 - waning.maternal[j]
      mat1[4,3] <- waning.maternal[j] #M3 to S
      mat1[4,4] <- 1-foi[j]
      mat1[5,4] <- foi[j] #S to I
      mat1[5,5] <- 1-recov
      mat1[6,5] <- recov #I to R1
      mat1[6,6] <- 1-sigma.fix   
      mat1[7,6] <- sigma.fix #R1 to R2     
      mat1[7,7] <- 1-sigma.fix   
      mat1[8,7] <- sigma.fix #R2 to R3
      mat1[8,8] <- 1-sigma[j]
      mat1[9,8] <- sigma[j] #R3 to N   
      mat1[9,9] <- 1- loss[j] - epsilon[j] #once in N, stay there, unless you lose immunity or advance to N2
      mat1[10,9] <- epsilon[j] #N1 to N2
      mat1[10,10] <- 1 #once in N2 you stay
      mat1[4,9] <- loss[j] #loss of immunity. moving from N1 back to S
      
    } else if(model =="MSIRSIR"){ #11 classes 
    
      mat1[] <- 0
      mat1[1,1] <- 1- wane.fix[j] 
      mat1[2,1] <- wane.fix[j] #M1 to M2
      mat1[2,2] <- 1- wane.fix[j]  
      mat1[3,2] <- wane.fix[j] #M2 to M3
      mat1[3,3] <- 1 - waning.maternal[j]
      mat1[4,3] <- waning.maternal[j] #M3 to S1
      mat1[4,4] <- 1-foi[j]
      mat1[5,4] <- foi[j] #S1 to I
      mat1[5,5] <- 1-recov
      mat1[6,5] <- recov #I to R1.1
      mat1[6,6] <- 1-sigma.fix   
      mat1[7,6] <- sigma.fix #R1.1 to R1.2     
      mat1[7,7] <- 1-sigma.fix   
      mat1[8,7] <- sigma.fix #R1.2 to R1.3
      mat1[8,8] <- 1-sigma[j]
      mat1[9,8] <- sigma[j] #R1.3 to S2   
      mat1[9,9] <- 1- foi2[j] #S2
      mat1[10,9] <- foi2[j] #S2 to I2
      mat1[10,10] <- 1-recov #I2
      mat1[11,10] <- recov #I2 to R2.1
      mat1[11,11] <- 1 #R2.1 stay here, don't lose titers
      
      
    }else if(model=="MSIRNIR"){ #13 classes
    
      mat1[] <- 0
      mat1[1,1] <- 1- wane.fix[j] 
      mat1[2,1] <- wane.fix[j] #M1 to M2
      mat1[2,2] <- 1- wane.fix[j]  
      mat1[3,2] <- wane.fix[j] #M2 to M3
      mat1[3,3] <- 1 - waning.maternal[j]
      mat1[4,3] <- waning.maternal[j] #M3 to S
      mat1[4,4] <- 1-foi[j]
      mat1[5,4] <- foi[j] #S to I
      mat1[5,5] <- 1-recov
      mat1[6,5] <- recov #I to R1.1
      mat1[6,6] <- 1-sigma.fix   
      mat1[7,6] <- sigma.fix #R1.1 to R1.2     
      mat1[7,7] <- 1-sigma.fix   
      mat1[8,7] <- sigma.fix #R1.2 to R1.3
      mat1[8,8] <- 1-sigma[j]
      mat1[9,8] <- sigma[j] #R1.3 to N   
      mat1[9,9] <- 1- recru[j] #once in N, stay there, unless you recrudesce to I2  
      mat1[10,9] <- recru[j] #recrudescence. moving from N to I2
      mat1[10,10] <- 1-recov #stay in I2
      mat1[11,10] <- recov# I2 to R2.1
      mat1[11,11] <- 1-sigma.fix # stay in R2.1
      mat1[12,11] <- sigma.fix # R2.1 to R2.2
      mat1[12,12] <- 1-sigma.fix # stay in R2.2
      mat1[13,12] <- sigma.fix # R2.2 to R2.3
      mat1[13,13] <- 1-sigma2[j] #stay in R2.3
      mat1[9,13] <- sigma2[j]# R2.3 to N
      
      
    } else if (model=="MSIRNIRN" | model=="MSIRNIRN-LT"){ 
    
      mat1[] <- 0
      mat1[1,1] <- 1- wane.fix[j] 
      mat1[2,1] <- wane.fix[j] #M1 to M2
      mat1[2,2] <- 1- wane.fix[j]  
      mat1[3,2] <- wane.fix[j] #M2 to M3
      mat1[3,3] <- 1 - waning.maternal[j]
      mat1[4,3] <- waning.maternal[j] #M3 to S
      mat1[4,4] <- 1-foi[j]
      mat1[5,4] <- foi[j] #S to I
      mat1[5,5] <- 1-recov
      mat1[6,5] <- recov #I to R1.1
      mat1[6,6] <- 1-sigma.fix   
      mat1[7,6] <- sigma.fix #R1.1 to R1.2     
      mat1[7,7] <- 1-sigma.fix   
      mat1[8,7] <- sigma.fix #R1.2 to R1.3
      mat1[8,8] <- 1-sigma[j]
      mat1[9,8] <- sigma[j] #R1.3 to N1   
      mat1[9,9] <- 1- recru[j] - epsilon[j] #once in N1, stay there, unless you recrudesce to I2 or advance to N2
      mat1[10,9] <- recru[j] #recrudescence. moving from N1 to I2
      mat1[10,10] <- 1-recov #stay in I2
      mat1[11,10] <- recov# I2 to R2.1
      mat1[11,11] <- 1-sigma.fix # stay in R2.1
      mat1[12,11] <- sigma.fix # R2.1 to R2.2
      mat1[12,12] <- 1-sigma.fix # stay in R2.2
      mat1[13,12] <- sigma.fix # R2.2 to R2.3
      mat1[13,13] <- 1-sigma2[j] #stay in R2.3
      mat1[14,13] <- sigma2[j]# R2.3 to N2 where you stay
      mat1[14,14] <- 1# N2 where you stay
      mat1[14,9] <- epsilon[j]# moving from N1 direct to N2
      
    } else if (model=="MSIRNSIRN" | model=="MSIRNSIRN-LT" ){ 
    
    mat1[] <- 0
    mat1[1,1] <- 1- wane.fix[j] 
    mat1[2,1] <- wane.fix[j] #M1 to M2
    mat1[2,2] <- 1- wane.fix[j]  
    mat1[3,2] <- wane.fix[j] #M2 to M3
    mat1[3,3] <- 1 - waning.maternal[j]
    mat1[4,3] <- waning.maternal[j] #M3 to S
    mat1[4,4] <- 1-foi[j]
    mat1[5,4] <- foi[j] #S to I
    mat1[5,5] <- 1-recov
    mat1[6,5] <- recov #I to R1.1
    mat1[6,6] <- 1-sigma.fix   
    mat1[7,6] <- sigma.fix #R1.1 to R1.2     
    mat1[7,7] <- 1-sigma.fix   
    mat1[8,7] <- sigma.fix #R1.2 to R1.3
    mat1[8,8] <- 1-sigma[j]
    mat1[9,8] <- sigma[j] #R1.3 to N1   
    mat1[9,9] <- 1- loss[j] - epsilon[j] #once in N1, stay there, unless you lose immunity to S2 or advance to N2
    mat1[10,9] <- loss[j] #loss of immunity. moving from N1 to S2
    mat1[10,10] <- 1- foi2[j] #staying in S2
    mat1[11,10] <- foi2[j] #from S2 to I2
    mat1[11,11] <- 1-recov #stay in I2
    mat1[12,11] <- recov# I2 to R2.1
    mat1[12,12] <- 1-sigma.fix # stay in R2.1
    mat1[13,12] <- sigma.fix # R2.1 to R2.2
    mat1[13,13] <- 1-sigma.fix # stay in R2.2
    mat1[14,13] <- sigma.fix # R2.2 to R2.3
    mat1[14,14] <- 1-sigma2[j] #stay in R2.3
    mat1[15,14] <- sigma2[j]# R2.3 to N2 where you stay
    mat1[15,15] <- 1# N2 where you stay
    mat1[15,9] <- epsilon[j]# moving from N1 direct to N2 (LT)
    
    }else if (model=="MSIRNSIRNS" ){
      
      mat1[] <- 0
      mat1[1,1] <- 1- wane.fix[j] 
      mat1[2,1] <- wane.fix[j] #M1 to M2
      mat1[2,2] <- 1- wane.fix[j]  
      mat1[3,2] <- wane.fix[j] #M2 to M3
      mat1[3,3] <- 1 - waning.maternal[j]
      mat1[4,3] <- waning.maternal[j] #M3 to S
      mat1[4,4] <- 1-foi[j]
      mat1[5,4] <- foi[j] #S to I
      mat1[5,5] <- 1-recov
      mat1[6,5] <- recov #I to R1.1
      mat1[6,6] <- 1-sigma.fix   
      mat1[7,6] <- sigma.fix #R1.1 to R1.2     
      mat1[7,7] <- 1-sigma.fix   
      mat1[8,7] <- sigma.fix #R1.2 to R1.3
      mat1[8,8] <- 1-sigma[j]
      mat1[9,8] <- sigma[j] #R1.3 to N1   
      mat1[9,9] <- 1- loss[j]  #once in N1, stay there, unless you lose immunity to S2 
      mat1[10,9] <- loss[j] #loss of immunity. moving from N1 to S2
      mat1[10,10] <- 1- foi2[j] #staying in S2
      mat1[11,10] <- foi2[j] #from S2 to I2
      mat1[11,11] <- 1-recov #stay in I2
      mat1[12,11] <- recov# I2 to R2.1
      mat1[12,12] <- 1-sigma.fix # stay in R2.1
      mat1[13,12] <- sigma.fix # R2.1 to R2.2
      mat1[13,13] <- 1-sigma.fix # stay in R2.2
      mat1[14,13] <- sigma.fix # R2.2 to R2.3
      mat1[14,14] <- 1-sigma2[j] #stay in R2.3
      mat1[15,14] <- sigma2[j]# R2.3 to N2 where you stay
      mat1[15,15] <- 1-loss2[j]# N2 where you stay unless you lose immunity back to S2
      mat1[10,15] <- loss2[j] # loss of immunity. moving from from N2 to S2
      
      
    }
    
    
    #put in surv. we first say it is equal across all epidemiological classes
    surv <- rep(1-mort.vect[j],c);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    
    if (add.inf.mort==TRUE){ #we'll let infecteds die more too
      
      #5th place in all models
      surv[5] <- 1-mort.vect.sick[j]
      
      #but the second infected class will vary
      
      if(model=="MSIRSIR"){
        surv[10] <- 1-mort.vect.sick[j] #I2 at #10
        
      }else if (model == "MSIRNSIRN" | model == "MSIRNSIRN-LT" | model == "MSIRNSIRNS"){
        surv[11] <- 1-mort.vect.sick[j] #I2 at #11
      }
      
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too - location will vary for all models
      
      #all models will have R1.1/R1.2/R1.3 in places 6,7,8
      surv[6] <- 1-mort.vect.rec[j]
      surv[7] <- 1-mort.vect.rec[j]
      surv[8] <- 1-mort.vect.rec[j]
      
      if ((model == "MSIRNIR") | (model=="MSIRSIR" ) | (model=="MSIRNIRN") | (model=="MSIRNIRN-LT") ){
        surv[11] <-  1-mort.vect.rec[j]
        surv[12] <-  1-mort.vect.rec[j]
        surv[13] <-  1-mort.vect.rec[j]
      } else if ((model=="MSIRNSIRNS") | (model=="MSIRNSIRN") | (model=="MSIRNSIRN-LT")){
        surv[12] <-  1-mort.vect.rec[j]
        surv[13] <-  1-mort.vect.rec[j]
        surv[14] <-  1-mort.vect.rec[j]
      }
      
   } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*c+1):(j*c+c),((j-1)*c+1):(j*c)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*c+1):(j*c),((j-1)*c+1):(j*c)] <- mat1*surv*(1-age.rate[j])
      
      
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
      Tmat[((j-1)*c+1):(j*c),((j-1)*c+1):(j*c)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildFMatrix_all <- function(c, model, age.classes, adult_fec, surv.biwk, biwk, N_stat){ 	#one for every age class
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
  Fmat <- matrix(0,c*nage,c*nage)
  
  for (j in 1:nage) { #no fertility in first age class (mat immune)
    if(model=="MSIRS"){
      Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],4))  #the mat immune offpring (so mom = I and R)
      Fmat[2,((j-1)*c+1):(j*c)]<- c(rep(fert.biwk[j],4),rep(0,4)) #the susceptible offspring (so mom= M and S)
      
    }else if(model=="MSIRN"|model=="MSIRNI"| model=="MSIRNS"){
    if(N_stat=="matAB"){
      #try it assuming N moms produce maternally immune pups
      Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],(c-4))) #the mat immune offspring (so mom = I, R, and N)
      Fmat[2,((j-1)*c+1):(j*c)]<-  c(rep(fert.biwk[j],4),rep(0,c-4)) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S)
      
    }else if (N_stat=="matSus"){
      Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],(c-5)), 0)  #the mat immune offpring (so mom = I and R but not N)
      Fmat[2,((j-1)*c+1):(j*c)]<- c(rep(fert.biwk[j],4),rep(0,(c-5)), fert.biwk[j]) #the susceptible offspring (so mom= M, S, and N)
    }
    }else if(model=="MSIRNNI"|model=="MSIRNNS"){
      if(N_stat=="matAB"){
        #try it assuming N moms produce maternally immune pups
        Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],(c-4))) #the mat immune offspring (so mom = I, R, and N)
        Fmat[2,((j-1)*c+1):(j*c)]<-  c(rep(fert.biwk[j],4),rep(0,c-4)) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S)
        
      }else if (N_stat=="matSus"){
        Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],(c-6)), rep(0,2))  #the mat immune offpring (so mom = I and R but not N)
        Fmat[2,((j-1)*c+1):(j*c)]<- c(rep(fert.biwk[j],4),rep(0,(c-6)), rep(fert.biwk[j],2)) #the susceptible offspring (so mom= M, S, and N)
      }
    } else  if(model=="MSIRSIR"){
      if(N_stat=="matAB"){
        #try it assuming N moms produce maternally immune pups
        Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],(c-4))) #the mat immune offspring (so mom = I, R, and N and S2)
        Fmat[2,((j-1)*c+1):(j*c)]<-  c(rep(fert.biwk[j],4),rep(0,c-4)) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S)
        
      }else if (N_stat=="matSus"){
        Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],(c-7)), rep(0,3))  #the mat immune offpring (so mom = I and R but not N or S2)
        Fmat[2,((j-1)*c+1):(j*c)]<- c(rep(fert.biwk[j],4),rep(0,(c-7)), rep(fert.biwk[j],3)) #the susceptible offspring (so mom= M, S, and N)
      }
    }else if (model=="MSIRNIR"){
     
      if(N_stat=="matAB"){
        #try it assuming N moms produce maternally immune pups
        Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],(c-4))) #the mat immune offspring (so mom = I, R, and N and S2)
        Fmat[2,((j-1)*c+1):(j*c)]<-  c(rep(fert.biwk[j],4),rep(0,c-4)) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S)
        
      }else if (N_stat=="matSus"){
        Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],(4)), 0,0, rep(fert.biwk[j],3))  #the mat immune offpring (so mom = I and R but not N or S2)
        Fmat[2,((j-1)*c+1):(j*c)]<- c(rep(fert.biwk[j],4),rep(0,(4)), rep(fert.biwk[j],5)) #the susceptible offspring (so mom= M, S, and N)
      } 
    }else if (model=="MSIRNIRN" | model=="MSIRNIRN-LT"){
      if(N_stat=="matAB"){
        #try it assuming N moms produce maternally immune pups
        Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],(c-4))) #the mat immune offspring (so mom = I, R, and N and S2)
        Fmat[2,((j-1)*c+1):(j*c)]<-  c(rep(fert.biwk[j],4),rep(0,c-4)) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S)
        
      }else if (N_stat=="matSus"){
        Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],(4)), 0,0, rep(fert.biwk[j],3),0)  #the mat immune offpring (so mom = I and R but not N or S2)
        Fmat[2,((j-1)*c+1):(j*c)]<- c(rep(fert.biwk[j],4),rep(0,(4)), rep(fert.biwk[j],2), rep(0,3), fert.biwk[j]) #the susceptible offspring (so mom= M, S, and N)
      } 
    }else if (model=="MSIRNSIRN" | model=="MSIRNSIRN-LT"|model=="MSIRNSIRNS"){
      if(N_stat=="matAB"){
        #try it assuming N moms produce maternally immune pups
        Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],(c-4))) #the mat immune offspring (so mom = I, R, and N and S2)
        Fmat[2,((j-1)*c+1):(j*c)]<-  c(rep(fert.biwk[j],4),rep(0,c-4)) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S)
        
      }else if (N_stat=="matSus"){
        Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],(4)), 0,0, rep(fert.biwk[j],4),0)  #the mat immune offpring (so mom = I and R but not N or S2)
        Fmat[2,((j-1)*c+1):(j*c)]<- c(rep(fert.biwk[j],4),rep(0,(4)), rep(fert.biwk[j],2), rep(0,4), fert.biwk[j]) #the susceptible offspring (so mom= M, S, and N)
      } 
        
    }
    }
    
    
  
  
  return(Fmat)
}


#and sim models
sim.met.age.all <- function(c,  sim_pop, yrs, ntyr, age.brk, s, beta,beta2, recov, mort, mort_juv, adult_fec, wane,wane.fix,loss, loss2, epsilon, sigma, sigma2, recru, slope.wane,mu.sick, mu.rec, sigma.fix, add.inf.mort,add.rec.mort, model, N_stat){
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  
  #first pull out your data...
  
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
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  
  
  M1_init = rep(0, s)
  M2_init = rep(0, s)
  M3_init = rep(0, s)
  
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  R1_init = rep(0, s)
  R2_init = rep(0, s)
  R3_init = rep(0, s)
  
  S_init = bat.mat - I_init 
  
  if(model =="MSIRS"){
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init)
    
  }else if(model =="MSIRN" | model =="MSIRNI" | model =="MSIRNS"){
    N_init = rep(0, s)
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init, N_init)
    
  }else if(model =="MSIRNNI"| model == "MSIRNNS"){
    N1_init = rep(0, s)
    N2_init = rep(0, s)
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init, N1_init, N2_init)
    
  }else if (model == "MSIRSIR"){
    S2_init = rep(0, s)
    I2_init = rep(0, s)
    R2_1_init = rep(0, s) #just stays in here
    
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init, S2_init, I2_init,  R2_1_init)
    
  }else if (model =="MSIRNIR"){
    N_init = rep(0, s)
    I2_init = rep(0, s)
    R2_1_init = rep(0, s)
    R2_2_init = rep(0, s)
    R2_3_init = rep(0, s)
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init, N_init, I2_init,  R2_1_init,  R2_2_init,  R2_3_init)
    
  }else if (model == "MSIRNIRN" | model == "MSIRNIRN-LT"){
    
    N1_init = rep(0, s)
    N2_init = rep(0, s)
    I2_init = rep(0, s)
    R2_1_init = rep(0, s)
    R2_2_init = rep(0, s)
    R2_3_init = rep(0, s)
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init, N1_init, I2_init,  R2_1_init,  R2_2_init,  R2_3_init, N2_init)
    
  }else if (model == "MSIRNSIRN" | model == "MSIRNSIRN-LT"| model == "MSIRNSIRNS"){
    
    N1_init = rep(0, s)
    N2_init = rep(0, s)
    I2_init = rep(0, s)
    S2_init = rep(0, s)
    R2_1_init = rep(0, s)
    R2_2_init = rep(0, s)
    R2_3_init = rep(0, s)
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init, N1_init, S2_init, I2_init,  R2_1_init,  R2_2_init,  R2_3_init, N2_init)
    
  }
  
  
  
  
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
    Tmat <- buildTMat_all_age(c=c, model = model, Npop= N_pop_ts[,i], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort)^(1/ntyr), 
                              surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr),
                              surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr),
                              surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, beta2=beta2, sigma=(1-exp(-(sigma*(1/ntyr)))),
                              sigma2=(1-exp(-(sigma2*(1/ntyr)))), loss=(1-exp(-(loss*(1/ntyr)))), loss2=(1-exp(-(loss2*(1/ntyr)))), 
                              recov=recov, wane.fix=wane.fix, sigma.fix=sigma.fix, epsilon = (1-exp(-(epsilon*(1/ntyr)))),  recru=(1-exp(-(recru*(1/ntyr)))),
                              wane=(1-exp(-(wane*(1/ntyr)))), add.rec.mort=add.rec.mort,  add.inf.mort=add.inf.mort, 
                              age.rate=(1/ntyr), slope.wane=slope.wane)
    
    
    #print(i)
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #print(biwk1)
    
    #feed into fertility matrix:
    Fmat <- buildFMatrix_all(c=c, model = model, age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
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
  
  #this will depend on the model
  if(model == "MSIRS"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1", "R1.1", "R1.2", "R1.3")  
  }else if(model =="MSIRN" | model =="MSIRNI" | model =="MSIRNS"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1",  "R1.1", "R1.2", "R1.3", "N1")  
  }else if(model =="MSIRNNI"| model == "MSIRNNS"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1",  "R1.1", "R1.2", "R1.3", "N1", "N2")  
  }else if(model == "MSIRSIR"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1",  "R1.1", "R1.2", "R1.3", "S2", "I2", "R2.1")  
  }else if(model =="MSIRNIR"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1",  "R1.1", "R1.2", "R1.3", "N1", "I2", "R2.1", "R2.2", "R2.3")  
  }else if (model == "MSIRNIRN" | model == "MSIRNIRN-LT"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1",  "R1.1", "R1.2", "R1.3", "N1", "I2", "R2.1", "R2.2", "R2.3", "N2")  
  }else if (model == "MSIRNSIRN" | model == "MSIRNSIRN-LT"| model == "MSIRNSIRNS"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1",  "R1.1", "R1.2", "R1.3", "N1", "S2", "I2", "R2.1", "R2.2", "R2.3", "N2")  
  }
  
  
  
  dat.tot$N = rowSums((dat.tot[,2:ncol(dat.tot)]))
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = make.prop(dat.tot)
  prop.tot <- dplyr::select(prop.tot, -(N))
  dat.tot <- dplyr::select(dat.tot, -(N))
  
  dat.tot <- dat.tot %>%
    pivot_longer(
      cols = -time,          # all columns except time
      names_to = "class",    # column for state names
      values_to = "count"    # column for counts
    )
  
  prop.tot <- prop.tot %>%
    pivot_longer(
      cols = -time,          # all columns except time
      names_to = "class",    # column for state names
      values_to = "proportion"    # column for counts
    )
  
  #
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) 
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) 
  
  dat.tot$proportion = dat.tot$count/sim_pop
  
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  
  age.dat.tot$class = rep(unique(dat.tot$class), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  
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
  #ggplot(eq.dat.sero) + geom_line(aes(x=times, y=count, color = class)) + facet_wrap(~class, scales="free")
  
  eq.dat.sero$sero[eq.dat.sero$class=="M1" | eq.dat.sero$class=="M2"| eq.dat.sero$class=="M3" | eq.dat.sero$class=="R1.1" |eq.dat.sero$class=="R1.2" | eq.dat.sero$class=="R1.3"| eq.dat.sero$class=="R2.1" |eq.dat.sero$class=="R2.2" | eq.dat.sero$class=="R2.3" ] <- eq.dat.sero$count[eq.dat.sero$class=="M1" | eq.dat.sero$class=="M2"| eq.dat.sero$class=="M3" | eq.dat.sero$class=="R1.1" |eq.dat.sero$class=="R1.2" | eq.dat.sero$class=="R1.3"| eq.dat.sero$class=="R2.1" |eq.dat.sero$class=="R2.2" | eq.dat.sero$class=="R2.3" ]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  #eq.sero.all <- ddply(eq.dat.sero, .(times, class, age), summarize, N= sum(count))
  
  #eq.sero.all  %>% mutate_if(is.factor, as.character) ->  eq.sero.all
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  #eq.sero.all$age = as.numeric(eq.sero.all$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  
  #eq.sum <- ddply(eq.sero.all,.(times, age), summarise, Ntot=sum(N))
  #eq.sero.all <- merge(eq.sero.all, eq.sum, by = c("times", "age"))
  #eq.sero.all$proportion = as.numeric(eq.sero.all$N)/as.numeric(eq.sero.all$Ntot)
  #out.prop.state <- eq.sero.all
  
  
  
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
  
  # out.prop.state$biwk=NA
  # for (i in 1:length(biwk.seq)){
  #   out.prop.state$biwk[out.prop.state$times==uni.times[i]] <- biwk.seq[i]  
  # }
  # out.prop.state$biwk <- as.factor(out.prop.state$biwk)
  
  prop.eq = subset(dat.tot, time >(yrs-1))
  prop.eq = dplyr::select(prop.eq,-(proportion))
  
  prop.eq$sero <- 0
  prop.eq$sero[prop.eq$class=="M1" | prop.eq$class=="M2"| prop.eq$class=="M3" | prop.eq$class=="R1.1" |prop.eq$class=="R1.2" | prop.eq$class=="R1.3"| prop.eq$class=="R2.1" |prop.eq$class=="R2.2" | prop.eq$class=="R2.3" ] <- prop.eq$count[prop.eq$class=="M1" | prop.eq$class=="M2"| prop.eq$class=="M3" | prop.eq$class=="R1.1" |prop.eq$class=="R1.2" | prop.eq$class=="R1.3"| prop.eq$class=="R2.1" |prop.eq$class=="R2.2" | prop.eq$class=="R2.3" ] 
  prop.eq$shed <- 0
  prop.eq$shed[prop.eq$class=="I1" | prop.eq$class=="I2"] <- prop.eq$count[prop.eq$class=="I1" | prop.eq$class=="I2"] 
  prop.eq$time <- round(prop.eq$time, 6)
  prop.sum <- ddply(prop.eq,.(time), summarize, seropos = sum(sero), prevpos = sum(shed), N=sum(count))
  prop.sum$seroprevalence <- prop.sum$seropos/prop.sum$N
  prop.sum$prevalence <- prop.sum$prevpos/prop.sum$N
  #head(prop.sum)
  #ggplot(data = prop.sum) + geom_line(aes(x=time, y=seroprevalence))
  #ggplot(data = prop.sum) + geom_line(aes(x=time, y=prevalence))
  prop.sum$biwk=NA
  
  for (i in 1:length(biwk.seq)){
    prop.sum$biwk[prop.sum$time==uni.times[i]] <- biwk.seq[i]  
  }
  prop.sum$biwk <- as.factor(prop.sum$biwk)
  
  #ggplot(data = prop.eq) + geom_line(aes(x=time, y=proportion, color=as.factor(biwk))) + facet_wrap(~class, scales="free")
  #and also get 
  
  # if (mort==.207){
  #   month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
  #   month.list <- rep(month.list, each=2)
  #   month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  # }else if (mort==.456){
  #   month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
  #   month.list <- rep(month.list, each=2)
  #   month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  # }
  # #need to make month.list same length as biwk 
  # uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  # out.mat$month <- NA
  # for (i in 1:length(uni.biwk)){
  #   out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  # }
  # 
  # if (mort==.207){
  #   out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  # }else if (mort==.456){
  #   out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  # }
  #ggplot(data=out.mat) + geom_line(aes(x=age, y=seroprevalence))
  #ggplot(data=out.mat) + geom_line(aes(x=age, y=seroprevalence)) + facet_wrap(~biwk, scales="free")
  
  return(list(out.mat, prop.sum))#, seas.prev))
}
sim.met.age.simple <- function(c,  sim_pop, yrs, ntyr, age.brk, s, par.df1, N_stat){
  
  
  model = par.df1$model
  beta=par.df1$beta
  loss=par.df1$loss
  loss2 = par.df1$loss2
  epsilon=par.df1$epsilon
  beta2=par.df1$beta2 
  recru=par.df1$recru
  sigma=par.df1$sigma
  sigma2=par.df1$sigma2
  recov=par.df1$recov
  mort=par.df1$mort
  mort_juv=par.df1$mort_juv
  adult_fec=par.df1$adult_fec
  wane=par.df1$wane
  wane.fix = par.df1$wane_fix
  slope.wane=1
  mu.sick=par.df1$mu_sick
  mu.rec=par.df1$mu_rec
  sigma.fix = par.df1$sigma_fix
  
  if(mu.sick>1 & mu.rec==1){
    add.inf.mort=TRUE
    add.rec.mort = FALSE
  }else if(mu.sick==1 & mu.rec>1){
    add.inf.mort=FALSE
    add.rec.mort = TRUE
  }else if(mu.sick==1 & mu.rec==1){
    add.inf.mort=FALSE
    add.rec.mort = FALSE
  }else if(mu.sick>1 & mu.rec>1){
    add.inf.mort=TRUE
    add.rec.mort = TRUE
  }
  
  mu.sick.juv = mu.sick
  mu.rec.juv = mu.rec
  if(mu.sick.juv>2){
    mu.sick.juv <- 2
  }
  if(mu.rec.juv>2){
    mu.rec.juv <- 2
  }
  
  #first pull out your data...
  
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
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  
  
  M1_init = rep(0, s)
  M2_init = rep(0, s)
  M3_init = rep(0, s)
  
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  R1_init = rep(0, s)
  R2_init = rep(0, s)
  R3_init = rep(0, s)
  
  S_init = bat.mat - I_init 
  
  if(model =="MSIRS"){
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init)
    
  }else if(model =="MSIRN" | model =="MSIRNI" | model =="MSIRNS"){
    N_init = rep(0, s)
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init, N_init)
    
  }else if(model =="MSIRNNI"| model == "MSIRNNS"){
    N1_init = rep(0, s)
    N2_init = rep(0, s)
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init, N1_init, N2_init)
    
  }else if (model == "MSIRSIR"){
    S2_init = rep(0, s)
    I2_init = rep(0, s)
    R2_1_init = rep(0, s) #just stays in here
    
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init, S2_init, I2_init,  R2_1_init)
    
  }else if (model =="MSIRNIR"){
    N_init = rep(0, s)
    I2_init = rep(0, s)
    R2_1_init = rep(0, s)
    R2_2_init = rep(0, s)
    R2_3_init = rep(0, s)
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init, N_init, I2_init,  R2_1_init,  R2_2_init,  R2_3_init)
    
  }else if (model == "MSIRNIRN" | model == "MSIRNIRN-LT"){
    
    N1_init = rep(0, s)
    N2_init = rep(0, s)
    I2_init = rep(0, s)
    R2_1_init = rep(0, s)
    R2_2_init = rep(0, s)
    R2_3_init = rep(0, s)
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init, N1_init, I2_init,  R2_1_init,  R2_2_init,  R2_3_init, N2_init)
    
  }else if (model == "MSIRNSIRN" | model == "MSIRNSIRN-LT"| model == "MSIRNSIRNS"){
    
    N1_init = rep(0, s)
    N2_init = rep(0, s)
    I2_init = rep(0, s)
    S2_init = rep(0, s)
    R2_1_init = rep(0, s)
    R2_2_init = rep(0, s)
    R2_3_init = rep(0, s)
    
    N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init, N1_init, S2_init, I2_init,  R2_1_init,  R2_2_init,  R2_3_init, N2_init)
    
  }
  
  
  
  
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
    Tmat <- buildTMat_all_age(c=c, model = model, Npop= N_pop_ts[,i], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort)^(1/ntyr), 
                              surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr),
                              surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr),
                              surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta, beta2=beta2, sigma=(1-exp(-(sigma*(1/ntyr)))),
                              sigma2=(1-exp(-(sigma2*(1/ntyr)))), loss=(1-exp(-(loss*(1/ntyr)))), loss2=(1-exp(-(loss2*(1/ntyr)))), 
                              recov=recov, wane.fix=wane.fix, sigma.fix=sigma.fix, epsilon = (1-exp(-(epsilon*(1/ntyr)))),  recru=(1-exp(-(recru*(1/ntyr)))),
                              wane=(1-exp(-(wane*(1/ntyr)))), add.rec.mort=add.rec.mort,  add.inf.mort=add.inf.mort, 
                              age.rate=(1/ntyr), slope.wane=slope.wane)
    
    
    #print(i)
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #print(biwk1)
    
    #feed into fertility matrix:
    Fmat <- buildFMatrix_all(c=c, model = model, age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
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
  
  #this will depend on the model
  if(model == "MSIRS"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1", "R1.1", "R1.2", "R1.3")  
  }else if(model =="MSIRN" | model =="MSIRNI" | model =="MSIRNS"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1",  "R1.1", "R1.2", "R1.3", "N1")  
  }else if(model =="MSIRNNI"| model == "MSIRNNS"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1",  "R1.1", "R1.2", "R1.3", "N1", "N2")  
  }else if(model == "MSIRSIR"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1",  "R1.1", "R1.2", "R1.3", "S2", "I2", "R2.1")  
  }else if(model =="MSIRNIR"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1",  "R1.1", "R1.2", "R1.3", "N1", "I2", "R2.1", "R2.2", "R2.3")  
  }else if (model == "MSIRNIRN" | model == "MSIRNIRN-LT"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1",  "R1.1", "R1.2", "R1.3", "N1", "I2", "R2.1", "R2.2", "R2.3", "N2")  
  }else if (model == "MSIRNSIRN" | model == "MSIRNSIRN-LT"| model == "MSIRNSIRNS"){
    names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1",  "R1.1", "R1.2", "R1.3", "N1", "S2", "I2", "R2.1", "R2.2", "R2.3", "N2")  
  }
  
  
  
  dat.tot$N = rowSums((dat.tot[,2:ncol(dat.tot)]))
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = make.prop(dat.tot)
  prop.tot <- dplyr::select(prop.tot, -(N))
  dat.tot <- dplyr::select(dat.tot, -(N))
  
  dat.tot <- dat.tot %>%
    pivot_longer(
      cols = -time,          # all columns except time
      names_to = "class",    # column for state names
      values_to = "count"    # column for counts
    )
  
  prop.tot <- prop.tot %>%
    pivot_longer(
      cols = -time,          # all columns except time
      names_to = "class",    # column for state names
      values_to = "proportion"    # column for counts
    )
  
  #
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) 
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) 
  
  dat.tot$proportion = dat.tot$count/sim_pop
  
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  
  age.dat.tot$class = rep(unique(dat.tot$class), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  
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
  #ggplot(eq.dat.sero) + geom_line(aes(x=times, y=count, color = class)) + facet_wrap(~class, scales="free")
  
  eq.dat.sero$sero[eq.dat.sero$class=="M1" | eq.dat.sero$class=="M2"| eq.dat.sero$class=="M3" | eq.dat.sero$class=="R1.1" |eq.dat.sero$class=="R1.2" | eq.dat.sero$class=="R1.3"| eq.dat.sero$class=="R2.1" |eq.dat.sero$class=="R2.2" | eq.dat.sero$class=="R2.3" ] <- eq.dat.sero$count[eq.dat.sero$class=="M1" | eq.dat.sero$class=="M2"| eq.dat.sero$class=="M3" | eq.dat.sero$class=="R1.1" |eq.dat.sero$class=="R1.2" | eq.dat.sero$class=="R1.3"| eq.dat.sero$class=="R2.1" |eq.dat.sero$class=="R2.2" | eq.dat.sero$class=="R2.3" ]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  #eq.sero.all <- ddply(eq.dat.sero, .(times, class, age), summarize, N= sum(count))
  
  #eq.sero.all  %>% mutate_if(is.factor, as.character) ->  eq.sero.all
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  #eq.sero.all$age = as.numeric(eq.sero.all$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  
  #eq.sum <- ddply(eq.sero.all,.(times, age), summarise, Ntot=sum(N))
  #eq.sero.all <- merge(eq.sero.all, eq.sum, by = c("times", "age"))
  #eq.sero.all$proportion = as.numeric(eq.sero.all$N)/as.numeric(eq.sero.all$Ntot)
  #out.prop.state <- eq.sero.all
  
  
  
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
  
  # out.prop.state$biwk=NA
  # for (i in 1:length(biwk.seq)){
  #   out.prop.state$biwk[out.prop.state$times==uni.times[i]] <- biwk.seq[i]  
  # }
  # out.prop.state$biwk <- as.factor(out.prop.state$biwk)
  
  prop.eq = subset(dat.tot, time >(yrs-1))
  prop.eq = dplyr::select(prop.eq,-(proportion))
  
  prop.eq$sero <- 0
  prop.eq$sero[prop.eq$class=="M1" | prop.eq$class=="M2"| prop.eq$class=="M3" | prop.eq$class=="R1.1" |prop.eq$class=="R1.2" | prop.eq$class=="R1.3"| prop.eq$class=="R2.1" |prop.eq$class=="R2.2" | prop.eq$class=="R2.3" ] <- prop.eq$count[prop.eq$class=="M1" | prop.eq$class=="M2"| prop.eq$class=="M3" | prop.eq$class=="R1.1" |prop.eq$class=="R1.2" | prop.eq$class=="R1.3"| prop.eq$class=="R2.1" |prop.eq$class=="R2.2" | prop.eq$class=="R2.3" ] 
  prop.eq$shed <- 0
  prop.eq$shed[prop.eq$class=="I1" | prop.eq$class=="I2"] <- prop.eq$count[prop.eq$class=="I1" | prop.eq$class=="I2"] 
  prop.eq$time <- round(prop.eq$time, 6)
  prop.sum <- ddply(prop.eq,.(time), summarize, seropos = sum(sero), prevpos = sum(shed), N=sum(count))
  prop.sum$seroprevalence <- prop.sum$seropos/prop.sum$N
  prop.sum$prevalence <- prop.sum$prevpos/prop.sum$N
  #head(prop.sum)
  #ggplot(data = prop.sum) + geom_line(aes(x=time, y=seroprevalence))
  #ggplot(data = prop.sum) + geom_line(aes(x=time, y=prevalence))
  prop.sum$biwk=NA
  
  for (i in 1:length(biwk.seq)){
    prop.sum$biwk[prop.sum$time==uni.times[i]] <- biwk.seq[i]  
  }
  prop.sum$biwk <- as.factor(prop.sum$biwk)
  
  #adjust for PCR weight
  prop.sum$prevalence =  prop.sum$prevalence*par.df1$weightPCR
  
  #ggplot(data = prop.eq) + geom_line(aes(x=time, y=proportion, color=as.factor(biwk))) + facet_wrap(~class, scales="free")
  #and also get 
  
  # if (mort==.207){
  #   month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
  #   month.list <- rep(month.list, each=2)
  #   month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  # }else if (mort==.456){
  #   month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
  #   month.list <- rep(month.list, each=2)
  #   month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  # }
  # #need to make month.list same length as biwk 
  # uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  # out.mat$month <- NA
  # for (i in 1:length(uni.biwk)){
  #   out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  # }
  # 
  # if (mort==.207){
  #   out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  # }else if (mort==.456){
  #   out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  # }
  #ggplot(data=out.mat) + geom_line(aes(x=age, y=seroprevalence))
  #ggplot(data=out.mat) + geom_line(aes(x=age, y=seroprevalence)) + facet_wrap(~biwk, scales="free")
  
  return(list(out.mat, prop.sum))#, seas.prev))
}



expand_serodata <- function(df){
  
  line_list <- do.call(
    rbind,
    lapply(1:nrow(df), function(i){
      
      n_pos <- df$Npos[i]
      n_neg <- df$Ntot[i] - n_pos
      
      data.frame(
        age = rep(df$age[i], df$Ntot[i]),
        serostatus = c(
          rep(1, n_pos),
          rep(0, n_neg)
        )
      )
    })
  )
  
  rownames(line_list) <- NULL
  
  return(line_list)
}
expand_serodata_biwk <- function(df){
  
  line_list <- do.call(
    rbind,
    lapply(1:nrow(df), function(i){
      
      n_pos <- df$Npos[i]
      n_neg <- df$Ntot[i] - n_pos
      
      data.frame(
        biwk = rep(df$biwk[i], df$Ntot[i]),
        status = c(
          rep(1, n_pos),
          rep(0, n_neg)
        )
      )
    })
  )
  
  rownames(line_list) <- NULL
  
  return(line_list)
}


#likelihood function
#this for a biweekly observation process
#this for a biweekly observation process
log.lik.sim <- function(par, c, data, sero.dat, shed.dat, cutoff, sim_pop, yrs, ntyr, age.brk, s, par.df1, model, N_stat){
  
  wPCR = exp(par[length(par)])
  par = par[1:(length(par)-1)]
  
  if (par.df1$mu_sick>1 & par.df1$mu_rec==1){
    add.inf.mort <- TRUE
    add.rec.mort <- FALSE
    mu.sick = exp(par[length(par)])
    par <- par[1:(length(par)-1)]
    mu.rec = 1
    
    #mu sick can't be <1 or much >4 because you'd have 0 survivorship
    if(mu.sick<1 | mu.sick > 4){
      ll=-1*10^9#write over with real number that is huge
      #print(-ll)
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
      #print(-ll)
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
      #print(-ll)
      return(-ll)
    }
    
  }else{
    
    add.inf.mort <- FALSE
    add.rec.mort <- FALSE
    mu.rec = 1
    mu.sick = 1
  }
  
  if (wPCR>1){
    ll=-1*10^9#write over with real number that is huge
    #print(-ll)
    return(-ll)
  }
  
  
  
  if(model=="MSIRS" | model=="MSIRN"){
    
    beta = exp(par[1])
    sigma = exp(par[2])
    
    recru = par.df1$recru
    sigma2 = par.df1$sigma2
    loss2 = par.df1$loss2
    loss = par.df1$loss
    epsilon = par.df1$epsilon
    beta2 = par.df1$beta2
     
  }else if (model=="MSIRNI"){
    beta = exp(par[1])
    sigma = exp(par[2])
    recru = exp(par[3])
    
    sigma2 = par.df1$sigma2
    loss2 = par.df1$loss2
    loss = par.df1$loss
    epsilon = par.df1$epsilon
    beta2 = par.df1$beta2
    
  }else if (model=="MSIRNS"){
    
    beta = exp(par[1])
    sigma = exp(par[2])
    loss = exp(par[3])
    
    sigma2 = par.df1$sigma2
    loss2 = par.df1$loss2
    recru = par.df1$recru
    epsilon = par.df1$epsilon
    beta2 = par.df1$beta2
    
  }else if(model == "MSIRNNS"){
    beta = exp(par[1])
    sigma = exp(par[2])
    loss = exp(par[3])
    epsilon = exp(par[4])
    
    recru = par.df1$recru
    sigma2 = par.df1$sigma2
    loss2 = par.df1$loss2
    beta2 = par.df1$beta2
    
    
  }else if(model == "MSIRNNI"){
    beta = exp(par[1])
    sigma = exp(par[2])
    recru = exp(par[3])
    epsilon = exp(par[4])
    
    sigma2 = par.df1$sigma2
    loss2 = par.df1$loss2
    loss = par.df1$loss
    beta2 = par.df1$beta2
    
  }else if (model == "MSIRNIRN" | model == "MSIRNIR"){
    
    beta = exp(par[1])
    sigma = exp(par[2])
    recru = exp(par[3])
    sigma2 = exp(par[4])
    
    
    loss2 = par.df1$loss2
    loss = par.df1$loss
    epsilon = par.df1$epsilon
    beta2 = par.df1$beta2
    
  }else if (model =="MSIRNSIRN"){
    beta = exp(par[1])
    sigma = exp(par[2])
    loss = exp(par[3])
    sigma2 = exp(par[4])
    beta2 = exp(par[5])
    
    
    loss2 = par.df1$loss2
    recru = par.df1$recru
    epsilon = par.df1$epsilon
    
    
  }else if (model == "MSIRSIR"){
    beta = exp(par[1])
    sigma = exp(par[2])
    beta2 = exp(par[3])
    
    sigma2 = par.df1$sigma2
    loss2 = par.df1$loss2
    loss = par.df1$loss
    epsilon = par.df1$epsilon
    recru = par.df1$recru
    
  }else if(model=="MSIRNIRN-LT"){
    beta = exp(par[1])
    sigma = exp(par[2])
    recru = exp(par[3])
    sigma2 = exp(par[4])
    epsilon = exp(par[5])
    
    
    loss2 = par.df1$loss2
    loss = par.df1$loss
    beta2 = par.df1$beta2
    
    
  }else if(model=="MSIRNSIRN-LT"){
    
    beta = exp(par[1])
    sigma = exp(par[2])
    loss = exp(par[3])
    sigma2 = exp(par[4])
    beta2 = exp(par[5])
    epsilon = exp(par[6])
    
    
    loss2 = par.df1$loss2
    recru = par.df1$recru
    
  }else if (model =="MSIRNSIRNS"){
    beta = exp(par[1])
    sigma = exp(par[2])
    loss = exp(par[3])
    sigma2 = exp(par[4])
    beta2 = exp(par[5])
    loss2 = exp(par[6])
    
    
    recru = par.df1$recru
    epsilon = par.df1$epsilon
    
  } 
  
  if(sum((1 - exp(-recru/26)), (1 - exp(-epsilon/26)))>1){
    
    ll=-1*10^9#write over with real number that is huge
    #print(-ll)
    return(-ll)
    
  }else if(sum((1 - exp(-loss/26)), (1 - exp(-epsilon/26)))>1){
    ll=-1*10^9#write over with real number that is huge
    #print(-ll)
    return(-ll)
    
  }else{
  
  
  out = sim.met.age.all(c=c,  sim_pop=sim_pop, yrs=yrs, ntyr=ntyr, age.brk=age.brk, s=s, 
                        beta=beta, loss=loss, loss2 = loss2,epsilon=epsilon,
                        beta2=beta2, recru=recru,sigma=sigma, sigma2=sigma2, 
                        recov=par.df1$recov, mort=par.df1$mort, mort_juv=par.df1$mort_juv, 
                        adult_fec=par.df1$adult_fec, wane=par.df1$wane,
                        wane.fix = par.df1$wane_fix, slope.wane=1, 
                        mu.sick=mu.sick, mu.rec=mu.rec,sigma.fix = par.df1$sigma_fix,
                        add.inf.mort=add.inf.mort,add.rec.mort=add.rec.mort, model=par.df1$model, N_stat = N_stat)
  age.seroprev <- out[[1]]
  prev.ts <- out[[2]]
  
  #ggplot(prev.ts) + geom_line(aes(x=as.numeric(biwk), y=prevalence))
  #ggplot(prev.ts) + geom_line(aes(x=as.numeric(biwk), y=seroprevalence))
  prev.ts$prevalence <- prev.ts$prevalence*wPCR
  
  
  
 
  #first, run your model at the desired parameters for estimation (par)
  #sometimes, they don't produce anything at all - in which case, just return a giant ll here and move on
  
  #print the guess par
  #print(exp(par))
  
  if(is.na(unique(age.seroprev$seropos)[1])){
    ll=-1*10^9#write over with real number that is huge
    #print(-ll)
    return(-ll)
  }else{
    
    #your data should have ages by biweek
    #ggplot(age.seroprev) + geom_point(aes(x=age, y=seroprevalence, color=biwk))
    # 
    # #fit a line  to these data to get the hazard by age:
    pred.fit <- with(age.seroprev,smooth.spline(x=age, y=seroprevalence))
    #separate age data
    
    #data = subset(data, !is.na(age)) #already done
    #these were also done by biweek
    data = arrange(data, age) #one data source gives you age seroprevalence
    sero.dat = arrange(sero.dat, biwk)
    shed.dat = arrange(shed.dat, biwk)
    
    #length(unique(sort(out$age)))#520
    #sort data by age and make a prediction for each age by biweek
    ages.new = sort(unique(data$age))
    
    out.comp <- cbind.data.frame(age=ages.new, seroprev=predict(pred.fit,x=ages.new)$y)
    #ggplot(age.seroprev) + geom_point(aes(x=age, y=seroprevalence, color=biwk)) + geom_line(data=out.comp, aes(x=age, y=seroprev))
    out.comp$seroprev[out.comp$seroprev<0] =0
    out.comp$seroprev[out.comp$seroprev>1] =1
    
    
    #now, compare model and data
    
    ll=0 #sets initial log-likelihood to 0
    ll_age = 0
    ll_sero = 0
    ll_shed = 0
    #then, compare with data by age and biweek
    for(a in 1:length(data$age)){ 
      #now run through data and compare with model
      #for-loops over every individual in the dataset, starting above the lower bound - 
      
      
      pprev <- out.comp$seroprev[out.comp$age==data$age[a]]  
      
      
      #evaluate age-seroprevalence 
      if(cutoff=="mean"){
        ll_age=ll_age+dbinom(data$serostatus[a],1,pprev,log=T)    
       # print(ll_age)
      }else if (cutoff=="lci"){
        ll_age=ll_age+dbinom(data$serostatus_lci[a],1,pprev,log=T)    
      }else if(cutoff=="uci"){
        ll_age=ll_age+dbinom(data$serostatus_uci[a],1,pprev,log=T)    
      }
      
    }
    
    #evaluate seasonal seroprevalence within a year
    for(a in 1:nrow(sero.dat)){ 
      #ggplot(prev.ts) + geom_line(aes(x=as.numeric(biwk), y=seroprevalence))
      seroprev <-prev.ts$seroprevalence[prev.ts$biwk==sero.dat$biwk[a]]  
      
      if(cutoff=="mean"){
        ll_sero=ll_sero+dbinom(sero.dat$status[a],1,seroprev,log=T)    
        #print(ll_sero)
      }else if (cutoff=="lci"){
        ll_sero=ll_sero+dbinom(sero.dat$status_lci[a],1,seroprev,log=T)    
      }else if(cutoff=="uci"){
        ll_sero=ll_sero+dbinom(sero.dat$status_uci[a],1,seroprev,log=T)    
      }
    }
    
    #evaluate shedding prevalence within a year
    for(a in 1:nrow(shed.dat)){ 
      #ggplot(prev.ts) + geom_line(aes(x=as.numeric(biwk), y=prevalence))
      #ggplot(shed.dat) + geom_jitter(aes(x=as.numeric(biwk), y=shed_pos), height=0.1, width =0.1)
      prevseas <-prev.ts$prevalence[prev.ts$biwk==shed.dat$biwk[a]]  
      
      
      #now run through data and compare with model
      #for-loops over every individual in the dataset, starting above the lower bound - 
      ll_shed=ll_shed+dbinom(shed.dat$status[a],1,prevseas,log=T)    
      #print(ll_shed)
     
    }
   
    if(-ll_age==Inf){
      ll_age=-1*10^9#write over with real number that is huge
      
    }
    if (-ll_sero == Inf){
      ll_sero=-1*10^9#write over with real number that is huge
    }
    if(-ll_shed==Inf){
      ll_shed=-1*10^9#write over with real number that is huge
    }
    
    
    ll = ll_age + ll_sero + ll_shed
    
    #print(-ll)
    #calculate and return the negative log-likelihood 
    #force mu.sick to be 1 or more
    
    return(-ll)
  }}}
sim.observations <- function(par, c, data, sero.dat, shed.dat, cutoff, sim_pop, yrs, ntyr, age.brk, s, par.df1, model, N_stat){
  
  wPCR = exp(par[length(par)])
  par = par[1:(length(par)-1)]
  
  if (par.df1$mu_sick>1 & par.df1$mu_rec==1){
    add.inf.mort <- TRUE
    add.rec.mort <- FALSE
    mu.sick = exp(par[length(par)])
    par <- par[1:(length(par)-1)]
    mu.rec = 1
    
    
  }else if (par.df1$mu_sick ==1 & par.df1$mu_rec>1){
    add.rec.mort <- TRUE
    add.inf.mort <- FALSE
    mu.rec = exp(par[length(par)])
    par <- par[1:(length(par)-1)]
    mu.sick = 1
    
    
    
  }else if (par.df1$mu_sick > 1 & par.df1$mu_rec>1){
    add.inf.mort <- TRUE
    add.rec.mort <- TRUE
    
    mu.rec = exp(par[length(par)])# last
    par <- par[1:(length(par)-1)]
    mu.sick = exp(par[length(par)]) #second to last
    par <- par[1:(length(par)-1)]
    
    
  }else{
    
    add.inf.mort <- FALSE
    add.rec.mort <- FALSE
    mu.rec = 1
    mu.sick = 1
  }
  
  
  
  if(model=="MSIRS" | model=="MSIRN"){
    
    beta = exp(par[1])
    sigma = exp(par[2])
    
    recru = par.df1$recru
    sigma2 = par.df1$sigma2
    loss2 = par.df1$loss2
    loss = par.df1$loss
    epsilon = par.df1$epsilon
    beta2 = par.df1$beta2
    
  }else if (model=="MSIRNI"){
    beta = exp(par[1])
    sigma = exp(par[2])
    recru = exp(par[3])
    
    sigma2 = par.df1$sigma2
    loss2 = par.df1$loss2
    loss = par.df1$loss
    epsilon = par.df1$epsilon
    beta2 = par.df1$beta2
    
  }else if (model=="MSIRNS"){
    
    beta = exp(par[1])
    sigma = exp(par[2])
    loss = exp(par[3])
    
    sigma2 = par.df1$sigma2
    loss2 = par.df1$loss2
    recru = par.df1$recru
    epsilon = par.df1$epsilon
    beta2 = par.df1$beta2
    
  }else if(model == "MSIRNNS"){
    beta = exp(par[1])
    sigma = exp(par[2])
    loss = exp(par[3])
    epsilon = exp(par[4])
    
    recru = par.df1$recru
    sigma2 = par.df1$sigma2
    loss2 = par.df1$loss2
    beta2 = par.df1$beta2
    
    
  }else if(model == "MSIRNNI"){
    beta = exp(par[1])
    sigma = exp(par[2])
    recru = exp(par[3])
    epsilon = exp(par[4])
    
    sigma2 = par.df1$sigma2
    loss2 = par.df1$loss2
    loss = par.df1$loss
    beta2 = par.df1$beta2
    
  }else if (model == "MSIRNIRN" | model == "MSIRNIR"){
    
    beta = exp(par[1])
    sigma = exp(par[2])
    recru = exp(par[3])
    sigma2 = exp(par[4])
    
    
    loss2 = par.df1$loss2
    loss = par.df1$loss
    epsilon = par.df1$epsilon
    beta2 = par.df1$beta2
    
  }else if (model =="MSIRNSIRN"){
    beta = exp(par[1])
    sigma = exp(par[2])
    loss = exp(par[3])
    sigma2 = exp(par[4])
    beta2 = exp(par[5])
    
    
    loss2 = par.df1$loss2
    recru = par.df1$recru
    epsilon = par.df1$epsilon
    
    
  }else if (model == "MSIRSIR"){
    beta = exp(par[1])
    sigma = exp(par[2])
    beta2 = exp(par[3])
    
    sigma2 = par.df1$sigma2
    loss2 = par.df1$loss2
    loss = par.df1$loss
    epsilon = par.df1$epsilon
    recru = par.df1$recru
    
  }else if(model=="MSIRNIRN-LT"){
    beta = exp(par[1])
    sigma = exp(par[2])
    recru = exp(par[3])
    sigma2 = exp(par[4])
    epsilon = exp(par[5])
    
    
    loss2 = par.df1$loss2
    loss = par.df1$loss
    beta2 = par.df1$beta2
    
    
  }else if(model=="MSIRNSIRN-LT"){
    
    beta = exp(par[1])
    sigma = exp(par[2])
    loss = exp(par[3])
    sigma2 = exp(par[4])
    beta2 = exp(par[5])
    epsilon = exp(par[6])
    
    
    loss2 = par.df1$loss2
    recru = par.df1$recru
    
  }else if (model =="MSIRNSIRNS"){
    beta = exp(par[1])
    sigma = exp(par[2])
    loss = exp(par[3])
    sigma2 = exp(par[4])
    beta2 = exp(par[5])
    loss2 = exp(par[6])
    
    
    recru = par.df1$recru
    epsilon = par.df1$epsilon
    
  } 
  
  
  #first, run your model at the desired parameters for estimation (par)
  #sometimes, they don't produce anything at all - in which case, just return a giant ll here and move on
  out = sim.met.age.all(c=c,  sim_pop=sim_pop, yrs=yrs, ntyr=ntyr, age.brk=age.brk, s=s, 
                          beta=beta, loss=loss, loss2 = loss2,epsilon=epsilon,
                          beta2=beta2, recru=recru,sigma=sigma, sigma2=sigma2, 
                          recov=par.df1$recov, mort=par.df1$mort, mort_juv=par.df1$mort_juv, 
                          adult_fec=par.df1$adult_fec, wane=par.df1$wane,
                          wane.fix = par.df1$wane_fix, slope.wane=1, 
                          mu.sick=mu.sick, mu.rec=mu.rec,sigma.fix = par.df1$sigma_fix,
                          add.inf.mort=add.inf.mort,add.rec.mort=add.rec.mort, model=par.df1$model, N_stat = N_stat)
    age.seroprev <- out[[1]]
    prev.ts <- out[[2]]
    
    #ggplot(prev.ts) + geom_line(aes(x=as.numeric(biwk), y=prevalence))
    #ggplot(prev.ts) + geom_line(aes(x=as.numeric(biwk), y=seroprevalence))
    prev.ts$prevalence <- prev.ts$prevalence*wPCR
    
    #now add the binomial noise for each estimate based on the counts from the data
    
      
      #your data should have ages by biweek
      #ggplot(age.seroprev) + geom_point(aes(x=age, y=seroprevalence, color=biwk))
      # 
      # #fit a line  to these data to get the hazard by age:
      pred.fit <- with(age.seroprev,smooth.spline(x=age, y=seroprevalence))
      #separate age data
      
      #data = subset(data, !is.na(age)) #already done
      #these were also done by biweek
      data = arrange(data, age) #one data source gives you age seroprevalence
      sero.dat = arrange(sero.dat, biwk)
      shed.dat = arrange(shed.dat, biwk)
      
      #length(unique(sort(out$age)))#520
      #sort data by age and make a prediction for each age by biweek
      ages.new = sort(unique(data$age))
      
      data.sum <- ddply(data,.(age), summarize, Npos=sum(prev), Ntot=length(prev))
      
      out.comp <- cbind.data.frame(age=ages.new, seroprev=predict(pred.fit,x=ages.new)$y)
      #ggplot(age.seroprev) + geom_point(aes(x=age, y=seroprevalence, color=biwk)) + geom_line(data=out.comp, aes(x=age, y=seroprev))
      out.comp$seroprev[out.comp$seroprev<0] =0
      out.comp$seroprev[out.comp$seroprev>1] =1
      
      
      #now, compare adjust model output based on sample size of the data
      y_age_sim <- rbinom(n = length(out.comp$seroprev), size = data.sum$Ntot, prob = out.comp$seroprev)
      #this gives a simulated trajectory of the number of seropositives per age class, using
      #the original sample size to estimate uncertainty
      
      #new simulated dataset 
      sim.age.sero <- data.sum
      sim.age.sero$Npos <- y_age_sim
      
      #and expand into linelist
      sim_agesero_linelist <- expand_serodata(sim.age.sero)
      
      #and do the same for the biweekly prevalence and seroprevalence data
      #we happen to know that the sero.dat data is missing some biweeks of sampling.
      #need to fill these in with zeros
      
      biwk.unique <- sort(unique(sero.dat$biwk))
      if(length(biwk.unique)<26){
        add.biwk <- setdiff(1:26, biwk.unique)
        add.tmp <- sero.dat[1:length(add.biwk),]
        add.tmp$seropos <- 0
        add.tmp$biwk <- add.biwk
        sero.dat <- rbind(sero.dat, add.tmp)
        sero.dat <- arrange(sero.dat, biwk)
      }
      
      sero.biwk.sum <- ddply(sero.dat,.(biwk), summarize, Npos=sum(seropos), Ntot=length(seropos))
      
      #generate new simulated prediction and write over
      y_sero_sum <- rbinom(n = length(prev.ts$seroprevalence), size = sero.biwk.sum$Ntot, prob = prev.ts$seroprevalence)
      sero.biwk.sum$Npos <- y_sero_sum
      
      #now expand into linelist
      sim_biwk_sero_linelist <- expand_serodata_biwk(sero.biwk.sum)
      
      
      #and lastly the prevalence data
      
      biwk.unique.prev <- sort(unique(shed.dat$biwk))
      if(length(biwk.unique.prev)<26){
        add.biwk.prev <- setdiff(1:26, biwk.unique.prev)
        add.tmp.prev <- shed.dat[1:length(add.biwk.prev),]
        add.tmp.prev$shed_pos <- 0
        add.tmp.prev$biwk <- add.biwk.prev
        shed.dat <- rbind(shed.dat, add.tmp.prev)
        shed.dat <- arrange(shed.dat, biwk)
      }
      
      shed.biwk.sum <- ddply(shed.dat,.(biwk), summarize, Npos=sum(shed_pos), Ntot=length(shed_pos))
      
      #generate new simulated prediction and write over
      y_shed_sum <- rbinom(n = length(prev.ts$prevalence), size = shed.biwk.sum$Ntot, prob = prev.ts$prevalence)
      shed.biwk.sum$Npos <- y_shed_sum
      
      sim_biwk_prev_linelist <- expand_serodata_biwk(shed.biwk.sum)
      
      
      #now return these three simulated datasets
      return(list(sim_agesero_linelist, sim_biwk_sero_linelist, sim_biwk_prev_linelist))
      
    }

fit.boot.data <- function(par, num_states, dat.tmp, model, par.df, sero.dat2, shed.dat2, cutoff1, s, sim_pop, ntyr, yrs, age.brk, N_stat){
  
  out.lik <- optim(par=log(par),fn=log.lik.sim, method="Nelder-Mead", c = num_states, data=as.data.frame(dat.tmp),
                   model=model, par.df1 =par.df, sero.dat = as.data.frame(sero.dat2), shed.dat=as.data.frame(shed.dat2),
                   cutoff=cutoff1, s=s, sim_pop=sim_pop, ntyr=ntyr, yrs=yrs, age.brk=age.brk, 
                   N_stat=N_stat, control=list(maxit=1000), hessian = FALSE)
  
  #then collect fitted parameters and return
  
  ConfidenceInterval <- data.frame(
    value = exp(out.lik$par)
  )
  
  # Optional: truncate negatives (though exp() should already prevent this)
  ConfidenceInterval[ConfidenceInterval < 0] <- 0
  
  
  
  if (par.df$mu_sick>1 & par.df$mu_rec==1){
    add.inf.mort <- TRUE
    add.rec.mort <- FALSE
    
  }else if (par.df$mu_sick ==1 & par.df$mu_rec>1){
    
    add.inf.mort <- FALSE
    add.rec.mort <- TRUE
    
  }else if (par.df$mu_sick > 1 & par.df$mu_rec>1){
    add.inf.mort <- TRUE
    add.rec.mort <- TRUE
    
  }else{
    
    add.inf.mort <- FALSE
    add.rec.mort <- FALSE
  }
  
  mu.sick = par.df$mu_sick
  mu.rec = par.df$mu_rec
  
  
  #collect the fitted values into the par.df dataframe
  
  if(model == "MSIRS" | model == "MSIRN"){
    if(mu.sick >1 & mu.rec==1){
      
      
      par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
      # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
      # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
      # 
      
      par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      # 
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
      # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
      # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
      # 
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
      # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
      # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
      # 
      
    }else if(mu.sick ==1 & mu.rec >1){
      
      par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
      # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
      # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
      # 
      par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      # 
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
      # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
      # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
      # 
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
      # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
      # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
      # 
    }else if (mu.sick >1 & mu.rec >1){
      
      par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
      # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
      # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
      # 
      
      par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      # 
      par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-2)] 
      # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
      # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
      # 
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
      # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
      # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
      # 
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
      # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
      # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
      # 
      
    }else{
      
      par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
      # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
      # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
      # 
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      # 
      
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-2)] 
      # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-2)]
      # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-2)]
      # 
      
    }
  }else if (model == "MSIRNI" ){
    
    if(mu.sick >1 & mu.rec==1){
      
      par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
      # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
      # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
      # 
      
      par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      # 
      par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-2)] 
      # par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
      # par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
      # 
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
      # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
      # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
      # 
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
      # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
      # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
      # 
    }else if(mu.sick ==1 & mu.rec >1){
      
      par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
      # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
      # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
      # 
      par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      # 
      par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-2)] 
      # par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
      # par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
      # 
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
      # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
      # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
      # 
      
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
      # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
      # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
      # 
      
    }else if (mu.sick >1 & mu.rec >1){
      
      par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
      # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
      # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
      # 
      par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
      # 
      par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-2)] 
      # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
      # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
      # 
      par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-3)] 
      # par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
      # par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
      # 
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
      # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
      # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
      # 
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
      # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
      # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
      # 
    }else{
      
      par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
      # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
      # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
      # 
      par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-1)] 
      # par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
      # par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
      # 
      par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
      # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
      # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
      # 
      par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
      # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
      # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
      # 
      
    }}else if (model == "MSIRNS"){
      
      if(mu.sick >1 & mu.rec==1){
        
        par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
        # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
        # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
        # 
        par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
        # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
        # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
        # 
        par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-2)] 
        # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
        # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
        # 
        par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
        # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
        # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
        # 
        par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
        # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
        # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
        # 
      }else if(mu.sick ==1 & mu.rec >1){
        
        par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
        # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
        # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
        # 
        par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
        # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
        # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
        # 
        par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-2)] 
        # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
        # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
        # 
        par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
        # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
        # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
        # 
        par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
        # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
        # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
        # 
      }else if (mu.sick >1 & mu.rec >1){
        
        par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
        # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
        # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
        # 
        par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
        # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
        # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
        # 
        par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-2)] 
        # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
        # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
        # 
        par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-3)] 
        # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
        # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
        # 
        par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
        # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
        # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
        # 
        par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
        # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
        # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
        # 
      }else{
        
        par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
        # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
        # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
        # 
        par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-1)] 
        # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
        # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
        # 
        par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
        # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
        # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
        # 
        par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
        # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
        # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
        # 
        # 
      }}else if ( model == "MSIRNNS"){
        
        if(mu.sick >1 & mu.rec==1){
          
          par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
          # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
          # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
          # 
          par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
          # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
          # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
          # 
          par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-2)] 
          # par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
          # par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
          # 
          par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-3)] 
          # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
          # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
          # 
          par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
          # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
          # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
          # 
          par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
          # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
          # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
          # 
        }else if(mu.sick ==1 & mu.rec >1){
          
          par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
          # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
          # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
          # 
          par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
          # par.df$mu_rec_lci =  ConfidenceInterval$lower[(length(out.lik$par)-1)] 
          # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
          # 
          par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-2)] 
          # par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
          # par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
          # 
          par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-3)] 
          # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
          # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
          # 
          par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
          # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
          # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
          # 
          par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
          # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
          # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
          # 
        }else if (mu.sick >1 & mu.rec >1){
          
          par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
          # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
          # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)] 
          # 
          par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
          # par.df$mu_rec_lci =  ConfidenceInterval$lower[(length(out.lik$par)-1)] 
          # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
          # 
          par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-2)] 
          # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
          # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
          # 
          par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-3)] 
          # par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
          # par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
          # 
          par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-4)] 
          # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
          # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
          # 
          par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
          # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
          # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
          # 
          par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
          # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
          # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
          # 
        }else{
          
          par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
          # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
          # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
          # 
          par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-1)] 
          # par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
          # par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
          # 
          par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-2)] 
          # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
          # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
          # 
          par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
          # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
          # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
          # 
          par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
          # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
          # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
          # 
          
        }}else if (model == "MSIRNIRN" | model == "MSIRNIR"){
          
          if(mu.sick >1 & mu.rec==1){
            
            par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
            # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
            # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
            # 
            par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
            # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
            # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
            # 
            par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
            # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
            # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
            # 
            par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-3)] 
            # par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
            # par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
            # 
            par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
            # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
            # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
            # 
            par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
            # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
            # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
            # 
          }else if(mu.sick ==1 & mu.rec >1){
            
            par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
            # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
            # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
            # 
            par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
            # par.df$mu_rec_lci =  ConfidenceInterval$lower[(length(out.lik$par)-1)] 
            # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
            # 
            par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
            # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
            # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)]
            # 
            par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-3)] 
            # par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
            # par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
            # 
            par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
            # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
            # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
            # 
            par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
            # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
            # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
            # 
          }else if (mu.sick >1 & mu.rec >1){
            
            par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
            # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
            # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
            # 
            par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
            # par.df$mu_rec_lci =  ConfidenceInterval$lower[(length(out.lik$par)-1)] 
            # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
            # 
            par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-2)] 
            # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
            # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
            # 
            par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
            # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
            # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
            # 
            par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-4)] 
            # par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
            # par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
            # 
            par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
            # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
            # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
            # 
            par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
            # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
            # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
            # 
          }else{
            
            par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
            # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
            # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
            # 
            par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
            # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
            # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
            # 
            par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-2)] 
            # par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
            # par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
            # 
            par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
            # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
            # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
            # 
            par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
            # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
            # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
            # 
            
          }}else if (model == "MSIRNSIRN"){
            
            if(mu.sick >1 & mu.rec==1){
              
              par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
              # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
              # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
              # 
              par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
              # par.df$mu_sick_lci =  ConfidenceInterval$lower[(length(out.lik$par)-1)] 
              # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
              # 
              par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
              # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
              # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
              # 
              par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
              # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
              # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
              # 
              par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-4)] 
              # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
              # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
              # 
              par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
              # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
              # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
              # 
              par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
              # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
              # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
              # 
            }else if(mu.sick ==1 & mu.rec >1){
              
              par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
              # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
              # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
              # 
              par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
              # par.df$mu_rec_lci =  ConfidenceInterval$lower[(length(out.lik$par)-1)] 
              # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
              # 
              par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
              # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
              # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
              # 
              par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
              # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
              # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
              # 
              par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-4)] 
              # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
              # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
              # 
              par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
              # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
              # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
              # 
              par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
              # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
              # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
              # 
            }else if (mu.sick >1 & mu.rec >1){
              
              par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
              # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
              # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
              # 
              par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
              # par.df$mu_rec_lci =  ConfidenceInterval$lower[(length(out.lik$par)-1)] 
              # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
              # 
              par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-2)] 
              # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
              # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
              # 
              par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
              # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
              # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
              # 
              
              par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-4)] 
              # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
              # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
              # 
              par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-5)] 
              # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
              # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
              # 
              par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-6)] 
              # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-6)] 
              # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-6)] 
              # 
              par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-7)] 
              # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-7)]
              # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-7)]
              # 
            }else{
              
              par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
              # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
              # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
              # 
              par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
              # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
              # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
              # 
              # 
              par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
              # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
              # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
              # 
              par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-3)] 
              # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
              # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
              # 
              par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
              # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
              # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
              # 
              par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
              # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
              # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
              # 
              
            }}else if (model == "MSIRSIR"){
              
              if(mu.sick >1 & mu.rec==1){
                
                par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                # 
                par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                # par.df$mu_sick_lci =  ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
                # 
                par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                # 
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                # 
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
                # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
                # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
                # 
              }else if(mu.sick ==1 & mu.rec >1){
                
                par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                # 
                par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
                # 
                par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                # 
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                # 
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-4)] 
                # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-4)]
                # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-4)]
                # 
              }else if (mu.sick >1 & mu.rec >1){
                
                par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                # 
                par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
                # 
                par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                # 
                par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                # 
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                # 
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
                # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
                # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
                # 
              }else{
                
                par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                # 
                
                par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                # par.df$beta2_lci =  ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
                # 
                par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                # 
                
                
                par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-3)] 
                # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-3)]
                # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-3)]
                # 
                
              }}else if(model=="MSIRNIRN-LT"){
                
                
                
                if(mu.sick >1 & mu.rec==1){
                  
                  par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                  # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                  # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                  # 
                  par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                  # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                  # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                  # 
                  par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                  # par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                  # par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                  # 
                  par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                  # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                  # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                  # 
                  par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                  # par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                  # par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                  # 
                  par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                  # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                  # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                  # 
                  par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
                  # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
                  # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
                  # 
                }else if(mu.sick ==1 & mu.rec >1){
                  
                  par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                  # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                  # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                  # 
                  par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                  # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                  # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                  # 
                  par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                  # par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                  # par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                  # 
                  par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                  # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                  # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                  # 
                  par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                  # par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                  # par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                  # 
                  par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                  # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                  # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                  # 
                  par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
                  # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
                  # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
                  # 
                }else if (mu.sick >1 & mu.rec >1){
                  
                  par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                  # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                  # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                  # 
                  par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                  # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                  # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                  # 
                  par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                  # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                  # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                  # 
                  par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                  # par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                  # par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                  # 
                  par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                  # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                  # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                  # 
                  par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                  # par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                  # par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                  # 
                  par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-6)] 
                  # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-6)] 
                  # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-6)] 
                  # 
                  par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-7)] 
                  # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-7)]
                  # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-7)]
                  # 
                }else{
                  
                  par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                  # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                  # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                  # 
                  par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                  # par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                  # par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                  # 
                  par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                  # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                  # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                  # 
                  par.df$recru = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                  # par.df$recru_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                  # par.df$recru_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                  # 
                  par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                  # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                  # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                  # 
                  par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-5)] 
                  par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-5)]
                  par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-5)]
                }
                
              }else if(model=="MSIRNSIRN-LT"){
                
                if(mu.sick >1 & mu.rec==1){
                  
                  par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                  # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                  # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                  # 
                  par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                  # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                  # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                  # 
                  par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                  # par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                  # par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                  # 
                  par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                  # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                  # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                  # 
                  par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                  # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                  # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                  # 
                  par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                  # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                  # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                  # 
                  par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-6)] 
                  # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-6)] 
                  # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-6)] 
                  # 
                  par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-7)] 
                  # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-7)]
                  # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-7)]
                  # 
                }else if(mu.sick ==1 & mu.rec >1){
                  
                  par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                  # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                  # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                  # 
                  par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                  # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                  # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                  # 
                  par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                  # par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                  # par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                  # 
                  par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                  # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                  # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                  # 
                  par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                  # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                  # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                  # 
                  par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                  # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                  # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                  # 
                  par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-6)] 
                  # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-6)] 
                  # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-6)] 
                  # 
                  par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-7)] 
                  # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-7)]
                  # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-7)]
                  # 
                }else if (mu.sick >1 & mu.rec >1){
                  
                  par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                  # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                  # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                  # 
                  par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                  # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                  # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                  # 
                  par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                  # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                  # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                  # 
                  par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                  # par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                  # par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)]
                  # 
                  par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                  # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                  # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                  # 
                  par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                  # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                  # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                  # 
                  par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-6)] 
                  # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-6)] 
                  # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-6)] 
                  # 
                  par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-7)] 
                  # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-7)] 
                  # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-7)] 
                  # 
                  par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-8)] 
                  # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-8)]
                  # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-8)]
                  # 
                }else{
                  
                  par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                  # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                  # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                  # 
                  par.df$epsilon = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                  # par.df$epsilon_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                  # par.df$epsilon_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                  # 
                  par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                  # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                  # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                  # 
                  par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                  # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                  # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                  # 
                  par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                  # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                  # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                  # 
                  par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                  # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                  # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                  # 
                  par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
                  # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
                  # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
                  # 
                  
                }
              }else if (model =="MSIRNSIRNS"){
                if(mu.sick >1 & mu.rec==1){
                  
                  par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                  # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                  # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                  # 
                  par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                  # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                  # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                  # 
                  par.df$loss2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                  # par.df$loss2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                  # par.df$loss2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                  # 
                  par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                  # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                  # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                  # 
                  par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                  # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                  # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                  # 
                  par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                  # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                  # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                  # 
                  par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-6)] 
                  # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-6)] 
                  # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-6)] 
                  # 
                  par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-7)] 
                  # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-7)]
                  # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-7)]
                  # 
                }else if(mu.sick ==1 & mu.rec >1){
                  
                  par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                  # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                  # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                  # 
                  par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                  # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                  # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                  # 
                  par.df$loss2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                  # par.df$loss2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                  # par.df$loss2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                  # 
                  par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                  # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                  # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                  # 
                  par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                  # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                  # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                  # 
                  par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                  # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                  # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                  # 
                  par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-6)] 
                  # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-6)] 
                  # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-6)] 
                  # 
                  par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-7)] 
                  # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-7)]
                  # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-7)]
                  # 
                }else if (mu.sick >1 & mu.rec >1){
                  
                  par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                  # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                  # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                  # 
                  par.df$mu_rec = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                  # par.df$mu_rec_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                  # par.df$mu_rec_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)] 
                  # 
                  par.df$mu_sick = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                  # par.df$mu_sick_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                  # par.df$mu_sick_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                  # 
                  par.df$loss2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                  # par.df$loss2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                  # par.df$loss2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)]
                  # 
                  par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                  # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                  # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                  # 
                  par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                  # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                  # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                  # 
                  par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-6)] 
                  # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-6)] 
                  # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-6)] 
                  # 
                  par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-7)] 
                  # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-7)] 
                  # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-7)] 
                  # 
                  par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-8)] 
                  # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-8)]
                  # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-8)]
                  # 
                }else{
                  
                  par.df$weightPCR <- ConfidenceInterval$value[length(out.lik$par)] 
                  # par.df$weightPCR_lci <- ConfidenceInterval$lower[length(out.lik$par)] 
                  # par.df$weightPCR_uci <- ConfidenceInterval$upper[length(out.lik$par)]
                  # 
                  par.df$loss2 = ConfidenceInterval$value[(length(out.lik$par)-1)] 
                  # par.df$loss2_lci = ConfidenceInterval$lower[(length(out.lik$par)-1)] 
                  # par.df$loss2_uci = ConfidenceInterval$upper[(length(out.lik$par)-1)]
                  # 
                  par.df$beta2 = ConfidenceInterval$value[(length(out.lik$par)-2)] 
                  # par.df$beta2_lci = ConfidenceInterval$lower[(length(out.lik$par)-2)] 
                  # par.df$beta2_uci = ConfidenceInterval$upper[(length(out.lik$par)-2)] 
                  # 
                  par.df$sigma2 = ConfidenceInterval$value[(length(out.lik$par)-3)] 
                  # par.df$sigma2_lci = ConfidenceInterval$lower[(length(out.lik$par)-3)] 
                  # par.df$sigma2_uci = ConfidenceInterval$upper[(length(out.lik$par)-3)] 
                  # 
                  par.df$loss = ConfidenceInterval$value[(length(out.lik$par)-4)] 
                  # par.df$loss_lci = ConfidenceInterval$lower[(length(out.lik$par)-4)] 
                  # par.df$loss_uci = ConfidenceInterval$upper[(length(out.lik$par)-4)] 
                  # 
                  par.df$sigma = ConfidenceInterval$value[(length(out.lik$par)-5)] 
                  # par.df$sigma_lci = ConfidenceInterval$lower[(length(out.lik$par)-5)] 
                  # par.df$sigma_uci = ConfidenceInterval$upper[(length(out.lik$par)-5)] 
                  # 
                  par.df$beta = ConfidenceInterval$value[1:(length(out.lik$par)-6)] 
                  # par.df$beta_lci = ConfidenceInterval$lower[1:(length(out.lik$par)-6)]
                  # par.df$beta_uci = ConfidenceInterval$upper[1:(length(out.lik$par)-6)]
                  # 
                  # 
                }
              }
  
  
  
  
  #save you likelihood and convergence history
  par.df$llik=out.lik$value
  par.df$convergence = out.lik$convergence
  return(par.df)
}

#need to have it specify where the age-breaks are and where the beta changes
wrap.fit.CI.boot <- function(N_boot, model, num_states, beta, sigma, recru, epsilon, loss, loss2, sigma2, beta2, mu.sick, mu.rec, dat1, sero.dat1, shed.dat1, par.df, cutoff1, species1, antigen1, s, sim_pop, yrs, ntyr,  N_stat, weight.PCR,
                     age.brk, slope.wane, do.plot, do.save, filename1, filename2, filename3){
  
  #first pull out your data...
  dat.tmp <- subset(dat1, species==species1 & antigen==antigen1)
  sero.dat2 = subset(sero.dat1, bat_species==species1 & antigen==antigen1)
  shed.dat2 = subset(shed.dat1, bat_species==species1 & antigen==antigen1)
  dat.tmp <- subset(dat.tmp, !is.na(age))
  
  #feed in the right par.dat 
  
  #edit the model name
  model <- sapply(strsplit(model, "-mu"),'[',1)
  par.df$model <- sapply(strsplit(par.df$model, "-mu"),'[',1)
  
  #dat.tmp <- subset(dat.tmp, fit==1). do you use all the data or just the longitudinal resamples?
  
  #identify your parameter(s) for fitting based on your model of choice
  if(model=="MSIRN" | model== "MSIRS"){
    par = c(beta, sigma)  
  }else if (model=="MSIRNI"){
    par = c(beta, sigma, recru)
  }else if (model=="MSIRNS"){
    par = c(beta, sigma, loss)
  }else if(model == "MSIRNNS"){
    par = c(beta, sigma, loss, epsilon)
  }else if(model == "MSIRNNI"){
    par = c(beta, sigma, recru, epsilon)
  }else if (model == "MSIRNIRN" | model == "MSIRNIR"){
    par = c(beta, sigma, recru, sigma2)
  }else if (model =="MSIRNSIRN"){
    par = c(beta, sigma, loss, sigma2, beta2)
  }else if (model == "MSIRSIR"){
    par = c(beta, sigma, beta2)
  }else if(model=="MSIRNIRN-LT"){
    par = c(beta, sigma, recru, sigma2, epsilon)
  }else if(model=="MSIRNSIRN-LT"){
    par = c(beta, sigma, loss, sigma2, beta2, epsilon)
  }else if (model =="MSIRNSIRNS"){
    par = c(beta, sigma, loss, sigma2, beta2, loss2)
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
  
  par <- c(par, weight.PCR)
  #print(par)
  #estimate parameters
  
  #no need to optimize parameters as these are already optimized
  # out.lik <- optim(par=log(par),fn=log.lik.age, method="Nelder-Mead", c = num_states, data=dat.tmp,
  #                  model=model, par.df1 =par.df, sero.dat = sero.dat2, shed.dat=shed.dat2,
  #                  cutoff=cutoff1, s=s, sim_pop=sim_pop, ntyr=ntyr, yrs=yrs, age.brk=age.brk, 
  #                  N_stat=N_stat, control=list(maxit=1000), hessian = FALSE)
  # 
  
  #and save these optimized parameters
  out.par <- data.frame(
    value = par
  )
  
  
  # now simulate the model using the best fit parameters as input
  # for as many iterations as you would like to bootstrap, 
  # eventually, after producing simulated data, you will refit to each simulation
  
  #here are all the simulated datasets
  out.sim <- replicate(n=N_boot, sim.observations(par=log(out.par$value), c = num_states, data=dat.tmp,
                              model=model, par.df1 =par.df, sero.dat = sero.dat2, shed.dat=shed.dat2,
                              cutoff=cutoff1, s=s, sim_pop=sim_pop, ntyr=ntyr, yrs=yrs, age.brk=age.brk, 
                              N_stat=N_stat), simplify=FALSE)
  
  save(out.sim, file = filename1)
  
  #now separate these into three lists of the three types of data
  all.age.sero <- lapply(out.sim, '[', 1)
  all.sero.biwk <- lapply(out.sim, '[', 2)
  all.shed.biwk <- lapply(out.sim, '[', 3)
  
  fit.sims <- mapply(FUN=fit.boot.data, dat.tmp=all.age.sero, sero.dat2=all.sero.biwk, shed.dat2=all.shed.biwk, 
                     MoreArgs=list(par=out.par$value, num_states=num_states, model=model,
                                   par.df=par.df, cutoff1=cutoff1, s=s, sim_pop=sim_pop,
                                   ntyr=ntyr, yrs=yrs, age.brk=age.brk, N_stat=N_stat), SIMPLIFY = FALSE)

  save(fit.sims, file = filename2)
  #then run the models again using the fitted parameters and save the output seroprevalence
  #to estimate the CIs on the seroprevalence
  
  
  run.fitted.sims = mapply(FUN=sim.met.age.simple, par.df1=fit.sims, 
                           MoreArgs=list(c=num_states,  sim_pop=sim_pop, yrs=yrs, ntyr=ntyr, 
                                         age.brk=age.brk, s=s,N_stat = N_stat), 
                           SIMPLIFY = FALSE)
   
  save(run.fitted.sims, file = filename3)
  #collect the age-seroprev data and the biweekly-seroprev/prev data
  all.age.sero.fitted <- lapply(run.fitted.sims, '[', 1)
  all.biwk.sero.fitted <- lapply(run.fitted.sims, '[', 2)
  

  #then estimate CIs on seroprevalence and parameters
  
  #seroprevalence by age
  seroprev_list <- do.call(rbind,lapply(all.age.sero.fitted, extract.seroprevalence))
  
  #and get the CIs
  seroprev_lci <- apply(seroprev_list, 2, quantile, 0.025)
  seroprev_uci <- apply(seroprev_list, 2, quantile, 0.975)
  seroprev_median <- apply(seroprev_list, 2, median)
  seroprev_lci_25 <- apply(seroprev_list, 2, quantile, 0.25)
  seroprev_uci_75 <- apply(seroprev_list, 2, quantile, 0.75)
  
  
  #seroprevalence by biweek
  seroprev_list_biwk <- do.call(rbind,lapply(all.biwk.sero.fitted, extract.seroprevalence))
  
  #and get the CIs
  seroprev_lci_biwk <- apply(seroprev_list_biwk, 2, quantile, 0.025)
  seroprev_uci_biwk <- apply(seroprev_list_biwk, 2, quantile, 0.975)
  seroprev_median_biwk <- apply(seroprev_list_biwk, 2, median)
  seroprev_lci_25_biwk <- apply(seroprev_list_biwk, 2, quantile, 0.25)
  seroprev_uci_75_biwk <- apply(seroprev_list_biwk, 2, quantile, 0.75)
  
  #prevalence by biweek
  prev_list_biwk <- do.call(rbind,lapply(all.biwk.sero.fitted, extract.prevalence))
  
  #and get the CIs
  prev_lci_biwk <- apply(prev_list_biwk, 2, quantile, 0.025)
  prev_uci_biwk <- apply(prev_list_biwk, 2, quantile, 0.975)
  prev_median_biwk <- apply(prev_list_biwk, 2, median)
  prev_lci_25_biwk <- apply(prev_list_biwk, 2, quantile, 0.25)
  prev_uci_75_biwk <- apply(prev_list_biwk, 2, quantile, 0.75)
  
  add.inf.mort = FALSE
  add.rec.mort = FALSE
  if(par.df$mu_sick>1){
    add.inf.mort = TRUE
  }
  if(par.df$mu_rec>1){
    add.rec.mort = TRUE
  }
  #run the deterministic model with the best fit parameters
  #first, the point estimate sim, using fitted parameters
  out1 = sim.met.age.all(c=par.df$num_states,  sim_pop=sim_pop, yrs=yrs, ntyr=ntyr, age.brk=age.brk, s=s, 
                         beta=par.df$beta, loss=par.df$loss, loss2 = par.df$loss2,
                         epsilon=par.df$epsilon, beta2=par.df$beta2, recru=par.df$recru,
                         sigma=par.df$sigma, sigma2=par.df$sigma2, 
                         recov=par.df$recov, mort=par.df$mort, mort_juv=par.df$mort_juv, 
                         adult_fec=par.df$adult_fec, wane=par.df$wane,
                         wane.fix = par.df$wane_fix, slope.wane=1, 
                         mu.sick=par.df$mu_sick, mu.rec=par.df$mu_rec,sigma.fix = par.df$sigma_fix,
                         add.inf.mort=add.inf.mort,add.rec.mort=add.rec.mort, model=par.df$model, N_stat = N_stat)
  
  #attach the prevalence and seroprevalence to the output
  age.sero.out <- out1[[1]]
  
  age.sero.out$seroprev_lci = seroprev_lci
  age.sero.out$seroprev_uci = seroprev_uci
  
  age.sero.out$seroprev_lci_25 = seroprev_lci_25
  age.sero.out$seroprev_uci_75 = seroprev_uci_75
  
  age.sero.out$seroprev_median = seroprev_median
  
  age.sero.out$cutoff <-  cutoff1
  age.sero.out$antigen <- antigen1
  age.sero.out$species <- species1
  
  biwk.ts.out <- out1[[2]]
  biwk.ts.out$prevalence <- biwk.ts.out$prevalence*par.df$weightPCR
  
  biwk.ts.out$seroprev_lci = seroprev_lci_biwk
  biwk.ts.out$seroprev_uci = seroprev_uci_biwk
  
  biwk.ts.out$seroprev_lci_25 = seroprev_lci_25_biwk
  biwk.ts.out$seroprev_uci_75 = seroprev_uci_75_biwk
  
  biwk.ts.out$seroprev_median = seroprev_median_biwk
  
  biwk.ts.out$prev_lci = prev_lci_biwk
  biwk.ts.out$prev_uci = prev_uci_biwk
  
  biwk.ts.out$prev_lci_25 = prev_lci_25_biwk
  biwk.ts.out$prev_uci_75 = prev_uci_75_biwk
  
  biwk.ts.out$prev_median = prev_median_biwk
  
  
  biwk.ts.out$cutoff <-  cutoff1
  biwk.ts.out$antigen <- antigen1
  biwk.ts.out$species <- species1
  
  
  #now collect the CIs on all the parameters --
  #these will come from collecting the quantiles on the fit.sims list
  boot_par <- do.call(rbind, lapply(fit.sims, get.parameter.matrix, model=model))
  
  par_ci <- apply(
    boot_par,
    2,
    quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  par_ci <- as.data.frame(par_ci)
  #now, attach these back to the fitted par.df
  #depending on the model in question
  
  
  #now, extract the fitted parameter matrix based on the model type
  if(model=="MSIRN" | model== "MSIRS"){
    par.df$beta_lci <- par_ci$beta[1]
    par.df$beta_uci <- par_ci$beta[2]
    par.df$sigma_lci <- par_ci$sigma[1]
    par.df$sigma_uci <- par_ci$sigma[2]
    
  }else if (model=="MSIRNI"){
    par.df$beta_lci <- par_ci$beta[1]
    par.df$beta_uci <- par_ci$beta[2]
    par.df$sigma_lci <- par_ci$sigma[1]
    par.df$sigma_uci <- par_ci$sigma[2]
    par.df$recru_lci <- par_ci$recru[1]
    par.df$recru_uci <- par_ci$recru[2]
    
  }else if (model=="MSIRNS"){
    par.df$beta_lci <- par_ci$beta[1]
    par.df$beta_uci <- par_ci$beta[2]
    par.df$sigma_lci <- par_ci$sigma[1]
    par.df$sigma_uci <- par_ci$sigma[2]
    par.df$loss_lci <- par_ci$loss[1]
    par.df$loss_uci <- par_ci$loss[2]
  }else if(model == "MSIRNNS"){
    par.df$beta_lci <- par_ci$beta[1]
    par.df$beta_uci <- par_ci$beta[2]
    par.df$sigma_lci <- par_ci$sigma[1]
    par.df$sigma_uci <- par_ci$sigma[2]
    par.df$loss_lci <- par_ci$loss[1]
    par.df$loss_uci <- par_ci$loss[2]
    par.df$epsilon_lci <- par_ci$epsilon[1]
    par.df$epsilon_uci <- par_ci$epsilon[2]
  }else if(model == "MSIRNNI"){
    
    par.df$beta_lci <- par_ci$beta[1]
    par.df$beta_uci <- par_ci$beta[2]
    par.df$sigma_lci <- par_ci$sigma[1]
    par.df$sigma_uci <- par_ci$sigma[2]
    par.df$recru_lci <- par_ci$recru[1]
    par.df$recru_uci <- par_ci$recru[2]
    par.df$epsilon_lci <- par_ci$epsilon[1]
    par.df$epsilon_uci <- par_ci$epsilon[2]
    
  }else if (model == "MSIRNIRN" | model == "MSIRNIR"){
    par.df$beta_lci <- par_ci$beta[1]
    par.df$beta_uci <- par_ci$beta[2]
    par.df$sigma_lci <- par_ci$sigma[1]
    par.df$sigma_uci <- par_ci$sigma[2]
    par.df$recru_lci <- par_ci$recru[1]
    par.df$recru_uci <- par_ci$recru[2]
    par.df$sigma2_lci <- par_ci$sigma2[1]
    par.df$sigma2_uci <- par_ci$sigma2[2]
    
  }else if (model =="MSIRNSIRN"){
    par.df$beta_lci <- par_ci$beta[1]
    par.df$beta_uci <- par_ci$beta[2]
    par.df$sigma_lci <- par_ci$sigma[1]
    par.df$sigma_uci <- par_ci$sigma[2]
    par.df$loss_lci <- par_ci$loss[1]
    par.df$loss_uci <- par_ci$loss[2]
    par.df$sigma2_lci <- par_ci$sigma2[1]
    par.df$sigma2_uci <- par_ci$sigma2[2]
    par.df$beta2_lci <- par_ci$beta2[1]
    par.df$beta2_uci <- par_ci$beta2[2]
  }else if (model == "MSIRSIR"){
    par.df$beta_lci <- par_ci$beta[1]
    par.df$beta_uci <- par_ci$beta[2]
    par.df$sigma_lci <- par_ci$sigma[1]
    par.df$sigma_uci <- par_ci$sigma[2]
    par.df$beta2_lci <- par_ci$beta2[1]
    par.df$beta2_uci <- par_ci$beta2[2]
    
  }else if(model=="MSIRNIRN-LT"){
    par.df$beta_lci <- par_ci$beta[1]
    par.df$beta_uci <- par_ci$beta[2]
    par.df$sigma_lci <- par_ci$sigma[1]
    par.df$sigma_uci <- par_ci$sigma[2]
    par.df$recru_lci <- par_ci$recru[1]
    par.df$recru_uci <- par_ci$recru[2]
    par.df$sigma2_lci <- par_ci$sigma2[1]
    par.df$sigma2_uci <- par_ci$sigma2[2]
    par.df$epsilon_lci <- par_ci$epsilon[1]
    par.df$epsilon_uci <- par_ci$epsilon[2]
    
  }else if(model=="MSIRNSIRN-LT"){
    par.df$beta_lci <- par_ci$beta[1]
    par.df$beta_uci <- par_ci$beta[2]
    par.df$sigma_lci <- par_ci$sigma[1]
    par.df$sigma_uci <- par_ci$sigma[2]
    par.df$loss_lci <- par_ci$loss[1]
    par.df$loss_uci <- par_ci$loss[2]
    par.df$sigma2_lci <- par_ci$sigma2[1]
    par.df$sigma2_uci <- par_ci$sigma2[2]
    par.df$beta2_lci <- par_ci$beta2[1]
    par.df$beta2_uci <- par_ci$beta2[2]
    par.df$epsilon_lci <- par_ci$epsilon[1]
    par.df$epsilon_uci <- par_ci$epsilon[2]
    
  }else if (model =="MSIRNSIRNS"){
    
    par.df$beta_lci <- par_ci$beta[1]
    par.df$beta_uci <- par_ci$beta[2]
    par.df$sigma_lci <- par_ci$sigma[1]
    par.df$sigma_uci <- par_ci$sigma[2]
    par.df$loss_lci <- par_ci$loss[1]
    par.df$loss_uci <- par_ci$loss[2]
    par.df$sigma2_lci <- par_ci$sigma2[1]
    par.df$sigma2_uci <- par_ci$sigma2[2]
    par.df$beta2_lci <- par_ci$beta2[1]
    par.df$beta2_uci <- par_ci$beta2[2]
    par.df$loss2_lci <- par_ci$loss2[1]
    par.df$loss2_uci <- par_ci$loss2[2]
    
  }
  
  #and the mu_sick, mu-rec
  
  if(par.df$mu_sick>1){
    par.df$mu_sick_lci <- par_ci$mu_sick[1]
    par.df$mu_sick_uci <- par_ci$mu_sick[2]
  }
  
  if(par.df$mu_rec>1){
    par.df$mu_rec_lci <- par_ci$mu_rec[1]
    par.df$mu_rec_uci <- par_ci$mu_rec[2]
  }
  
  par.df$weightPCR_lci <- par_ci$weightPCR[1]
  par.df$weightPCR_uci <- par_ci$weightPCR[2]
  
  
  #we penalize only for the estimated parameters in AIC
  if(par.df$mu_sick >1 & par.df$mu_rec==1){
    model <- paste0(model, "-mu-sick")
  }
  if(par.df$mu_sick ==1 & par.df$mu_rec >1){
    model <- paste0(model, "-mu-rec")
  }
  if(par.df$mu_sick >1 & par.df$mu_rec >1){
    model <- paste0(model, "-mu-sick-mu-rec")
  }
  
  biwk.ts.out$model <-  model
  age.sero.out$model <- model
  
  par.df$model <- model
  par.df$cutoff <-  cutoff1
  par.df$antigen <- antigen1
  par.df$species <- species1
  
  
  return(list(par.df, age.sero.out, biwk.ts.out)) 
  
  
  
  
}

extract.seroprevalence <- function(df){
  df <- as.data.frame(df)
  seroprev <- df$seroprevalence
  return(seroprev)
}

extract.prevalence <- function(df){
  df <- as.data.frame(df)
  prev <- df$prevalence
  return(prev)
}
get.parameter.matrix <- function(df, model){
  df = as.data.frame(df)
  
  #now, extract the fitted parameter matrix based on the model type
  if(model=="MSIRN" | model== "MSIRS"){
    par = c(df$beta, df$sigma)  
  }else if (model=="MSIRNI"){
    par = c(df$beta, df$sigma, df$recru)
  }else if (model=="MSIRNS"){
    par = c(df$beta, df$sigma, df$loss)
    names(par) <- c("beta", "sigma", "loss")
  }else if(model == "MSIRNNS"){
    par = c(df$beta, df$sigma, df$loss, df$epsilon)
    names(par) <- c("beta", "sigma", "loss", "epsilon")
  }else if(model == "MSIRNNI"){
    par = c(df$beta, df$sigma, df$recru, df$epsilon)
  }else if (model == "MSIRNIRN" | model == "MSIRNIR"){
    par = c(df$beta, df$sigma, df$recru, df$sigma2)
  }else if (model =="MSIRNSIRN"){
    par = c(df$beta, df$sigma, df$loss, df$sigma2, df$beta2)
  }else if (model == "MSIRSIR"){
    par = c(df$beta, df$sigma, df$beta2)
  }else if(model=="MSIRNIRN-LT"){
    par = c(df$beta, df$sigma, df$recru, df$sigma2, df$epsilon)
  }else if(model=="MSIRNSIRN-LT"){
    par = c(df$beta, df$sigma, df$loss, df$sigma2, df$beta2, df$epsilon)
  }else if (model =="MSIRNSIRNS"){
    par = c(df$beta, df$sigma, df$loss, df$sigma2, df$beta2, df$loss2)
  }
  
  if(df$mu_sick>1){
    par1 = c(par, df$mu_sick)
    names(par1) <- c(names(par), "mu_sick")
    par = par1
  }
  if(df$mu_rec>1){
    par1 = c(par, df$mu_rec)
    names(par1) <- c(names(par), "mu_rec")
    par = par1
  }
  
  par1 <- c(par, df$weightPCR)
  names(par1) <- c(names(par), "weightPCR")
  par = par1

  
  return(par)
  
}

slurm_array_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

par.fit <- read.csv(file =  paste0("fit-pars-", slurm_array_id, ".csv"), header = T, stringsAsFactors = F)


out.fit <- wrap.fit.CI.boot(N_boot=5,
                              model=par.fit$model,
                              dat1=dat,
                              sero.dat1 = sero.dat,
                              shed.dat1 = shed.dat,
                              species1 = "Eidolon dupreanum",
                              antigen1=par.fit$antigen,
                              cutoff1 = par.fit$cutoff,
                              num_states = par.fit$num_states,
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
                              loss = par.fit$loss,
                              loss2 = par.fit$loss2,
                              weight.PCR = par.fit$weightPCR,
                              N_stat="matAB",
                              do.plot=FALSE,
                              do.save=FALSE,
                              filename1 = paste0(unique(par.fit$model), "-", unique(par.fit$antigen), "-", unique(par.fit$cutoff), "-sim-data.Rdata"),
                              filename2 = paste0(unique(par.fit$model), "-", unique(par.fit$antigen), "-", unique(par.fit$cutoff), "-fitted-sims.Rdata"),
                              filename3 = paste0(unique(par.fit$model), "-", unique(par.fit$antigen), "-", unique(par.fit$cutoff), "-runs-w-fits.Rdata"))

filename = paste0(unique(out.fit$model), "-", unique(par.fit$antigen), "-", unique(par.fit$cutoff), slurm_array_id, ".Rdata")
save(out.fit, file=filename)
