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
library(numDeriv)

#this is the version that builds a comparable matrix to fit the top five model structures
#load data
#homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/model_fitting/GhV-fits/run-CIs-matrix"
#setwd(homewd)
#dat <- read.csv(file =  paste0(homewd,"/data/age-seroprevalence-data.csv"), header = T, stringsAsFactors = F)
dat <- read.csv(file =  "age-seroprevalence-data.csv", header = T, stringsAsFactors = F)
sero.dat <- read.csv(file =  "biweekly-serology-data.csv", header = T, stringsAsFactors = F)
shed.dat <- read.csv(file =  "biweekly-shedding-data.csv", header = T, stringsAsFactors = F)
head(dat)
head(sero.dat)
head(shed.dat)


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

#build functions - here this is for one matrix structure only and applicable only to MSIRS/MSIRNS/MSIRNNS/MSIRNI/MSIRNNI models
buildTMat_all_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, surv.biwk.sick, surv.juv.biwk.sick,surv.biwk.rec, surv.juv.biwk.rec,	beta, recru, recru2, recov, sigma1, sigma2, loss, wane,wane.fix,sigma.fix,epsilon, add.inf.mort, add.rec.mort, age.rate, slope.wane){
  
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
  
  # if (length(beta2)==1) beta2 <- rep(beta2,nage) # 
  # if (length(beta2)==s) beta2 <- beta2
  # if (length(beta2) > 1 & length(beta2)<s){
  #   beta2.list <- list()
  #   for (i in 1:length(beta2)){
  #     beta2.list[[i]] = rep(beta2[i], age.brk[i])
  #   }
  #   beta2 = c(unlist(beta2.list))
  # } 
  # 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma1)==1) sigma1 <- rep(sigma1,nage)
  if (length(sigma2)==1) sigma2 <- rep(sigma2,nage)
  if (length(recru)==1) recru <- rep(recru,nage)
  if (length(recru2)==1) recru2 <- rep(recru2,nage)
  #if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  #if (length(mu.rec)==1) mu.rec <- rep(mu.rec, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  if (length(wane.fix)==1) wane.fix <- rep(wane.fix,nage)
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
  
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,5] #always. for all model forms, I1 comes #5
  foi = 1-exp(-beta*((sum(I_m))/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  
  mat1 <- matrix(0,c,c) # 
  
  Tmat <- matrix(0,c*nage,c*nage) #
  
  for (j in 1:nage) { 
    
    
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
    mat1[8,8] <- 1-sigma1[j] -sigma2[j] - recru2[j] #R3 stays R3
    mat1[9,8] <- sigma2[j] #R3 to N1
    mat1[4,8] <-  sigma1[j]#R3 to S in MSIRS
    mat1[5,8] <-  recru2[j]#R3 to I in MSIRI
    mat1[9,9] <- 1- recru[j] - loss[j] - epsilon[j] #once in N1, stay there, unless you recrudesce to I or lose immunity (loss) to S or advance (epsilon) to N2  
    mat1[4,9] <- loss[j] #N1 to S
    mat1[10,9] <- epsilon[j] #N1 to N2
    mat1[10,10] <- 1 #once in N2 you stay there
    mat1[5,9] <- recru[j] #recrudescence. moving from N1 back to I
    
    
    #put in surv. we first say it is equal across all epidemiological classes
    surv <- rep(1-mort.vect[j],c);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    
    if (add.inf.mort==TRUE){ #we'll let infecteds die more too
      
      #5th place in all models
      surv[5] <- 1-mort.vect.sick[j]
      
      #but the second infected class will vary
      # 
      # if(model=="MSIRSIR"){
      #   surv[10] <- 1-mort.vect.sick[j] #I2 at #10
      #   
      # }else if (model == "MSIRNSIRN" | model == "MSIRNSIRN-LT" | model == "MSIRNSIRNS"){
      #   surv[11] <- 1-mort.vect.sick[j] #I2 at #11
      # }
      
    }
    if (add.rec.mort==TRUE){ #we'll let seropositives die more too - location will vary for all models
      
      #all models will have R1.1/R1.2/R1.3 in places 6,7,8
      surv[6] <- 1-mort.vect.rec[j]
      surv[7] <- 1-mort.vect.rec[j]
      surv[8] <- 1-mort.vect.rec[j]
      
      # if ((model == "MSIRNIR") | (model=="MSIRSIR" ) | (model=="MSIRNIRN") | (model=="MSIRNIRN-LT") ){
      #   surv[11] <-  1-mort.vect.rec[j]
      #   surv[12] <-  1-mort.vect.rec[j]
      #   surv[13] <-  1-mort.vect.rec[j]
      # } else if ((model=="MSIRNSIRNS") | (model=="MSIRNSIRN") | (model=="MSIRNSIRN-LT")){
      #   surv[12] <-  1-mort.vect.rec[j]
      #   surv[13] <-  1-mort.vect.rec[j]
      #   surv[14] <-  1-mort.vect.rec[j]
      # }
      
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
buildFMatrix_all <- function(c,  age.classes, adult_fec, surv.biwk, biwk, N_stat){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3 #using adult survival rate here assumes a post-breeding census: mom must survive to the end of the previous timestep, then give birth
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
    
    if(N_stat=="matAB"){
      #try it assuming N moms produce maternally immune pups
      Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],(c-4))) #the mat immune offspring (so mom = I, R, and N)
      Fmat[4,((j-1)*c+1):(j*c)]<-  c(rep(fert.biwk[j],4),rep(0,c-4)) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S)
      
    }else if (N_stat=="matSus"){
      Fmat[1,((j-1)*c+1):(j*c)] <- c(rep(0, 4),rep(fert.biwk[j],4), rep(0,2))  #the mat immune offpring (so mom = I and R but not N)
      Fmat[4,((j-1)*c+1):(j*c)]<- c(rep(fert.biwk[j],4),rep(0,4), rep(fert.biwk[j],2)) #the susceptible offspring (so mom= M, S, and N)
    }
  } 
  
  
  
  return(Fmat)
}


#and sim models
sim.met.age.all <- function(c,  sim_pop, yrs, ntyr, age.brk, s, beta,recov, mort, mort_juv, adult_fec, wane,wane.fix,loss,  epsilon, sigma1, sigma2, recru, recru2, slope.wane,mu.sick, mu.rec, sigma.fix, add.inf.mort,add.rec.mort, model, N_stat){
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
  
  N1_init = rep(0, s)
  N2_init = rep(0, s)
  
  N_tot = cbind(M1_init, M2_init, M3_init, S_init, I_init, R1_init, R2_init, R3_init, N1_init, N2_init)
  
  
  
  
  
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
    Tmat <- buildTMat_all_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort)^(1/ntyr), 
                              surv.juv.biwk = (1-mort_juv)^(1/ntyr), surv.juv.biwk.sick=(1-(mort_juv*mu.sick.juv))^(1/ntyr),
                              surv.biwk.sick=(1-(mort*mu.sick))^(1/ntyr), surv.juv.biwk.rec=(1-(mort_juv*mu.rec.juv))^(1/ntyr),
                              surv.biwk.rec=(1-(mort*mu.rec))^(1/ntyr),	beta=beta,  sigma1=(1-exp(-(sigma1*(1/ntyr)))),
                              sigma2=(1-exp(-(sigma2*(1/ntyr)))), loss=(1-exp(-(loss*(1/ntyr)))), recru2=(1-exp(-(recru2*(1/ntyr)))),
                              recov=recov, wane.fix=wane.fix, sigma.fix=sigma.fix, epsilon = (1-exp(-(epsilon*(1/ntyr)))),  recru=(1-exp(-(recru*(1/ntyr)))),
                              wane=(1-exp(-(wane*(1/ntyr)))), add.rec.mort=add.rec.mort,  add.inf.mort=add.inf.mort, 
                              age.rate=(1/ntyr), slope.wane=slope.wane)
    
    
    #print(i)
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #print(biwk1)
    
    #feed into fertility matrix:
    Fmat <- buildFMatrix_all(c=c, age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
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
  
  names(dat.tot) = c("time", "M1", "M2", "M3", "S1", "I1",  "R1.1", "R1.2", "R1.3", "N1", "N2")  
  
  
  
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



#likelihood function
#this for a biweekly observation process
log.lik.age <- function(par, name.par, c, data, sero.dat, shed.dat, cutoff, sim_pop, yrs, ntyr, age.brk, s, par.df1, model, N_stat){
  
  names(par) <- name.par
  wPCR = exp(par["wPCR"])
  
  
  if (par.df1$mu_sick>1 & par.df1$mu_rec==1){
    add.inf.mort <- TRUE
    add.rec.mort <- FALSE
    mu.sick = exp(par["mu_sick"])
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
    mu.rec = exp(par["mu_rec"])
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
    
    mu.rec = exp(par["mu_rec"])
    mu.sick = exp(par["mu_sick"])
    
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
  
  
  
  if(model=="MSIRS" ){
    
    beta = exp(par["beta"])
    sigma1 = exp(par["sigma1"])
    
    recru = par.df1$recru
    sigma2 = par.df1$sigma2
    loss = par.df1$loss
    epsilon = par.df1$epsilon
    recru2 = par.df1$recru2
    
  }else if (model=="MSIRI"){
    beta = exp(par["beta"])
    recru2 = exp(par["recru2"])
    
    recru = par.df1$recru
    sigma1 = par.df1$sigma1
    loss = par.df1$loss
    epsilon = par.df1$epsilon
    sigma2 = par.df1$sigma2
    
  }else if (model=="MSIRSI"){
    beta = exp(par["beta"])
    recru2 = exp(par["recru2"])
    sigma1 = exp(par["sigma1"])
    
    recru = par.df1$recru
    loss = par.df1$loss
    epsilon = par.df1$epsilon
    sigma2 = par.df1$sigma2
    
  }else if (model=="MSIRN"){
    beta = exp(par["beta"])
    sigma2 = exp(par["sigma2"])
    
    recru = par.df1$recru
    sigma1 = par.df1$sigma1
    loss = par.df1$loss
    epsilon = par.df1$epsilon
    recru2 = par.df1$recru2
    
  }else if (model=="MSIRNI"){
    beta = exp(par["beta"])
    sigma2 = exp(par["sigma2"])
    recru = exp(par["recru"])
    
    sigma1 = par.df1$sigma1
    loss = par.df1$loss
    epsilon = par.df1$epsilon
    recru2 = par.df1$recru2
    
  }else if (model=="MSIRNS"){
    
    beta = exp(par["beta"])
    sigma2 = exp(par["sigma2"])
    loss = exp(par["loss"])
    
    sigma1 = par.df1$sigma1
    epsilon = par.df1$epsilon
    recru = par.df1$recru
    recru2 = par.df1$recru2
    
  }else if (model=="MSIRNSI"){
    
    beta = exp(par["beta"])
    sigma2 = exp(par["sigma2"])
    loss = exp(par["loss"])
    recru = exp(par["recru"])
    
    sigma1 = par.df1$sigma1
    epsilon = par.df1$epsilon
    recru2 = par.df1$recru2
    
  }else if(model == "MSIRNNS"){
    beta = exp(par["beta"])
    sigma2 = exp(par["sigma2"])
    loss = exp(par["loss"])
    epsilon = exp(par["epsilon"])
    
    sigma1 = par.df1$sigma1
    recru = par.df1$recru
    recru2 = par.df1$recru2
    
    
  }else if(model == "MSIRNNI"){
    beta = exp(par["beta"])
    sigma2 = exp(par["sigma2"])
    recru = exp(par["recru"])
    epsilon = exp(par["epsilon"])
    
    
    sigma1 = par.df1$sigma1
    loss = par.df1$loss
    recru2 = par.df1$recru2
    
  }else if(model == "MSIRNNSI"){
    beta = exp(par["beta"])
    sigma2 = exp(par["sigma2"])
    recru = exp(par["recru"])
    epsilon = exp(par["epsilon"])
    loss = exp(par["loss"])
    
    sigma1 = par.df1$sigma1
    recru2 = par.df1$recru2
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
                          beta=beta, loss=loss, epsilon=epsilon, recru2=recru2,
                          recru=recru,sigma1=sigma1, sigma2=sigma2, 
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
    
    if(sum(is.finite(unique(age.seroprev$seroprevalence)))>0 | sum(is.finite(unique(age.seroprev$age)))>0){
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
          ll_age=ll_age+dbinom(data$prev[a],1,pprev,log=T)    
          # print(ll_age)
        }else if (cutoff=="lci"){
          ll_age=ll_age+dbinom(data$prev_lci[a],1,pprev,log=T)    
        }else if(cutoff=="uci"){
          ll_age=ll_age+dbinom(data$prev_uci[a],1,pprev,log=T)    
        }
        
      }
      
      #evaluate seasonal seroprevalence within a year
      for(a in 1:nrow(sero.dat)){ 
        #ggplot(prev.ts) + geom_line(aes(x=as.numeric(biwk), y=seroprevalence))
        seroprev <-prev.ts$seroprevalence[prev.ts$biwk==sero.dat$biwk[a]]  
        
        if(cutoff=="mean"){
          ll_sero=ll_sero+dbinom(sero.dat$seropos[a],1,seroprev,log=T)    
          #print(ll_sero)
        }else if (cutoff=="lci"){
          ll_sero=ll_sero+dbinom(sero.dat$seropos_lci[a],1,seroprev,log=T)    
        }else if(cutoff=="uci"){
          ll_sero=ll_sero+dbinom(sero.dat$seropos_uci[a],1,seroprev,log=T)    
        }
      }
      
      #evaluate shedding prevalence within a year
      for(a in 1:nrow(shed.dat)){ 
        #ggplot(prev.ts) + geom_line(aes(x=as.numeric(biwk), y=prevalence))
        #ggplot(shed.dat) + geom_jitter(aes(x=as.numeric(biwk), y=shed_pos), height=0.1, width =0.1)
        prevseas <-prev.ts$prevalence[prev.ts$biwk==shed.dat$biwk[a]]  
        
        
        #now run through data and compare with model
        #for-loops over every individual in the dataset, starting above the lower bound - 
        ll_shed=ll_shed+dbinom(shed.dat$shed_pos[a],1,prevseas,log=T)    
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

#will need to restructure below
wrap.fit.cov <- function(model, num_states, beta, sigma1, recru, recru2, epsilon, loss,  sigma2,  mu.sick, mu.rec, dat1, sero.dat1, shed.dat1, par.df, cutoff1, species1, antigen1, s, sim_pop, yrs, ntyr,  N_stat, weight.PCR,
                     age.brk, slope.wane, do.plot, do.save, filename){
  
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
  #identify your parameter(s) for fitting based on your model of choice
  if(model== "MSIRN"){
    par = c(beta, sigma2)
    name.par <- c("beta", "sigma2")
  }else if(model== "MSIRS"){
    par = c(beta, sigma1)  
    name.par <- c("beta", "sigma1")
  }else if(model== "MSIRI"){
    par = c(beta, recru2)  
    name.par <- c("beta", "recru2")
  }else if(model== "MSIRSI"){
    par = c(beta, sigma1, recru2)  
    name.par <- c("beta", "sigma1", "recru2")
  }else if (model=="MSIRNI"){
    par = c(beta, sigma2, recru)
    name.par <- c("beta", "sigma2", "recru")
  }else if (model=="MSIRNS"){
    par = c(beta, sigma2, loss)
    name.par <- c("beta", "sigma2", "loss")
  }else if (model=="MSIRNSI"){
    par = c(beta, sigma2, loss, recru)
    name.par <- c("beta", "sigma2", "loss", "recru")
  }else if(model == "MSIRNNS"){
    par = c(beta, sigma2, loss, epsilon)
    name.par <- c("beta", "sigma2", "loss", "epsilon")
  }else if(model == "MSIRNNI"){
    par = c(beta, sigma2, recru, epsilon)
    name.par <- c("beta", "sigma2", "recru", "epsilon")
  }else if(model == "MSIRNNSI"){
    par = c(beta, sigma2, recru, loss, epsilon)
    name.par <- c("beta", "sigma2", "recru", "loss", "epsilon")
  }
  
  if(mu.sick >1 & mu.rec==1){
    par <- c(par, mu.sick)
    name.par <- c(name.par, "mu_sick")
  }
  if(mu.sick ==1 & mu.rec >1){
    par <- c(par, mu.rec)
    name.par <- c(name.par, "mu_rec")
  }
  if(mu.sick >1 & mu.rec >1){
    par <- c(par, mu.sick, mu.rec)
    name.par <- c(name.par, "mu_sick", "mu_rec")
  }
  
  par <- c(par, weight.PCR)
  name.par <- c(name.par, "wPCR")
  #print(par)
  
  #no need to estimate parameters since code brings in parameters that
  #are already best fit and then functions to estimate the CIs on
  #seroprevalence and prevalence
  
  # Vector 'par' contains your saved optimized point estimates on the natural scale
  log_par_hat <- log(par)
  
  # Calculate high-precision Hessian at your saved point estimates
  # using hessian function from numDeriv package
  
  H_accurate <- hessian(
    func = log.lik.age, 
    x = log_par_hat, 
    # Pass all extra arguments required by log.lik.age:
    name.par = name.par, c = num_states, data = dat.tmp,
    model = model, par.df1 = par.df, sero.dat = sero.dat2, shed.dat = shed.dat2,
    cutoff = cutoff1, s = s, sim_pop = sim_pop, ntyr = ntyr, yrs = yrs, age.brk = age.brk, 
    N_stat = N_stat
  )
  
  # Symmetrize
  H_sym <- (H_accurate + t(H_accurate)) / 2
  
  # Invert to get the covariance matrix
  # If solve() throws an error due to singularity, fall back to Matrix::nearPD
  cov_mat <- tryCatch({
    solve(H_sym)
  }, error = function(e) {
    message("Hessian still singular/non-PD. Applying nearPD fallback...")
    as.matrix(solve(Matrix::nearPD(H_sym)$mat))
  })
  
  # Sample parameters and proceed with prevalence/seroprevalence predictions
  # sigma is the covariance matrix which is the inverse of the hessian 
  # when you are minimizing the neg log liklihood in your function
  par.dist.log <- rmvnorm(n = 500, mean = log_par_hat, sigma = cov_mat)
  par.dist <- exp(par.dist.log)
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
  
  #predictions for seroprevalence here
  preds <- matrix(NA,nrow=length(biwk.span),ncol=length(par.dist[,1]))
  
  
  #predictions for biweekly prevalence and seroprevalence will be shorter
  biwk_preds_prev <- matrix(NA,nrow=26,ncol=length(par.dist[,1]))
  biwk_preds_seroprev <- biwk_preds_prev 
  
  
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
  
  
  
  # now run one sim with the point parameters and a series of sims across 
  # the covariance matrix to produce the confidence intervals on the seroprevalence
  # and prevalence outputs
  
  #first, the point estimate sim, using fitted parameters
  out1 = sim.met.age.all(c=par.df$num_states,  sim_pop=sim_pop, yrs=yrs, ntyr=ntyr, age.brk=age.brk, s=s, 
                         beta=par.df$beta, loss=par.df$loss, recru2=par.df$recru2,
                         epsilon=par.df$epsilon,  recru=par.df$recru,
                         sigma1=par.df$sigma1, sigma2=par.df$sigma2, 
                         recov=par.df$recov, mort=par.df$mort, mort_juv=par.df$mort_juv, 
                         adult_fec=par.df$adult_fec, wane=par.df$wane,
                         wane.fix = par.df$wane_fix, slope.wane=1, 
                         mu.sick=par.df$mu_sick, mu.rec=par.df$mu_rec,sigma.fix = par.df$sigma_fix,
                         add.inf.mort=add.inf.mort,add.rec.mort=add.rec.mort, model=par.df$model, N_stat = N_stat)
  
  # name the parameter guesses with the corresponding names from 
  names(par.dist)  <- name.par
  #and the series across the covariance matrix
  for (j in 1:length(par.dist[,1])){
    
    if(model=="MSIRS" ){
      
      beta=par.dist[j,"beta"]
      sigma1=par.dist[j,"sigma1"]
      
      recru = par.df$recru
      sigma2 = par.df$sigma2
      loss = par.df$loss
      epsilon = par.df$epsilon
      recru2 = par.df$recru2
      
      
    }else if (model=="MSIRI"){
      beta=par.dist[j,"beta"]
      recru2 = par.dist[j,"recru2"]
      
      recru = par.df$recru
      sigma2 = par.df$sigma2
      sigma1 = par.df$sigma1
      loss = par.df$loss
      epsilon = par.df$epsilon
      
    }else if (model=="MSIRSI"){
      beta=par.dist[j,"beta"]
      recru2 = par.dist[j,"recru2"]
      sigma1=par.dist[j,"sigma1"]
      
      recru = par.df$recru
      sigma2 = par.df$sigma2
      loss = par.df$loss
      epsilon = par.df$epsilon
      
    }else if (model=="MSIRNI"){
      beta=par.dist[j,"beta"]
      sigma2=par.dist[j,"sigma2"]
      recru = par.dist[j,"recru"]
      
      sigma1 = par.df$sigma1
      recru2 = par.df$recru2
      loss = par.df$loss
      epsilon = par.df$epsilon
      
    }else if (model=="MSIRNS"){
      
      beta=par.dist[j,"beta"]
      sigma2=par.dist[j,"sigma2"]
      loss = par.dist[j,"loss"]
      
      recru2 = par.df$recru2
      sigma1 = par.df$sigma1
      recru = par.df$recru
      epsilon = par.df$epsilon
      
    }else if (model=="MSIRNS"){
      
      beta=par.dist[j,"beta"]
      sigma2=par.dist[j,"sigma2"]
      loss = par.dist[j,"loss"]
      recru = par.dist[j,"recru"]
      
      recru2 = par.df$recru2
      sigma1 = par.df$sigma1
      epsilon = par.df$epsilon
      
    }else if(model == "MSIRNNS"){
      
      beta=par.dist[j,"beta"]
      sigma2=par.dist[j,"sigma2"]
      loss = par.dist[j,"loss"]
      epsilon = par.dist[j,"epsilon"]
      
      recru2 = par.df$recru2
      recru = par.df$recru
      sigma1 = par.df$sigma1
      

    }else if(model == "MSIRNNI"){
      beta=par.dist[j,"beta"]
      sigma2=par.dist[j,"sigma2"]
      recru = par.dist[j,"recru"]
      epsilon = par.dist[j,"epsilon"]
      
      recru2 = par.df$recru2
      sigma1 = par.df$sigma1
      loss = par.df$loss
      
    }else if(model == "MSIRNNSI"){
      beta=par.dist[j,"beta"]
      sigma2=par.dist[j,"sigma2"]
      recru = par.dist[j,"recru"]
      loss = par.dist[j,"loss"]
      epsilon = par.dist[j,"epsilon"]
      
      recru2 = par.df$recru2
      sigma1 = par.df$sigma1
      
    }else if(model == "MSIRN"){
      beta=par.dist[j,"beta"]
      sigma2=par.dist[j,"sigma2"]
      
      
      recru2 = par.df$recru2
      recru = par.df$recru
      sigma1 = par.df$sigma1
      loss = par.df$loss
      epsilon = par.df$epsilon
      
    }
    
    weight.PCR = par.dist[j, "wPCR"]
    
    if(mu.sick >1 & mu.rec==1){
      
      mu.sick = par.dist[j, "mu_sick"]
      mu.rec = par.df$mu_rec
      
    }else if (mu.sick ==1 & mu.rec >1){
      
      mu.sick = par.df$mu_sick
      mu.rec = par.dist[j, "mu_rec"]
      
    }else if(mu.sick >1 & mu.rec >1){
      
      mu.sick = par.dist[j, "mu_sick"]
      mu.rec = par.dist[j, "mu_rec"]
      
    }else{
      mu.sick = par.df$mu_sick
      mu.rec = par.df$mu_rec
      
    }
    
    
    
    #the parameters here do not pull from the par.df which are the best fits 
    #but instead from the selections specified in the if/else statements above 
    #that sample across par.dist, the covariance matrix.
    
    out = sim.met.age.all(c=par.df$num_states,  sim_pop=sim_pop, 
                          yrs=yrs, ntyr=ntyr, age.brk=age.brk, s=s, 
                          beta=beta, loss=loss, 
                          epsilon=epsilon, recru=recru, recru2=recru2,
                          sigma1=sigma1, sigma2=sigma2, 
                          recov=par.df$recov, mort=par.df$mort, mort_juv=par.df$mort_juv, 
                          adult_fec=par.df$adult_fec, wane=par.df$wane,
                          wane.fix = par.df$wane_fix, slope.wane=1, 
                          mu.sick=mu.sick, mu.rec=mu.rec,sigma.fix = par.df$sigma_fix,
                          add.inf.mort=add.inf.mort,add.rec.mort=add.rec.mort, 
                          model=par.df$model, N_stat = N_stat)
    
    #each column is a new guess. each row is output by age and biweek
    seroprev.age.df <- out[[1]]
    seroprev.prev.ts <- out[[2]]
    
    preds[,j] <- seroprev.age.df$seroprevalence
    biwk_preds_seroprev[,j] <- seroprev.prev.ts$seroprevalence
    biwk_preds_prev[,j] <- seroprev.prev.ts$prevalence*weight.PCR
    
    #write some script to evaluate fit of model against data at these parameters
    
    #and evaluate likelihood against that already estimate
    prev.tmp = seroprev.prev.ts$prevalence*weight.PCR
    seroprev.tmp = seroprev.prev.ts$seroprevalence
    age.sero.tmp = seroprev.age.df$seroprevalence
    
    #your data should have ages by biweek
    #ggplot(age.seroprev) + geom_point(aes(x=age, y=seroprevalence, color=biwk))
    # 
    # #fit a line  to these data to get the hazard by age:
    pred.fit <- with(seroprev.age.df,smooth.spline(x=age, y=seroprevalence))
    
    #arreanf
    data = arrange(dat.tmp, age) #one data source gives you age seroprevalence
    sero.dat = arrange(sero.dat1, biwk)
    shed.dat = arrange(shed.dat1, biwk)
    
    #length(unique(sort(out$age)))#520
    #sort data by age and make a prediction for each age by biweek
    ages.new = sort(unique(data$age))
    
    out.comp <- cbind.data.frame(age=ages.new, seroprev=predict(pred.fit,x=ages.new)$y)
    #ggplot(seroprev.age.df) + geom_point(aes(x=age, y=seroprevalence, color=biwk)) + geom_line(data=out.comp, aes(x=age, y=seroprev))
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
      if(cutoff1=="mean"){
        ll_age=ll_age+dbinom(data$prev[a],1,pprev,log=T)    
        # print(ll_age)
      }else if (cutoff1=="lci"){
        ll_age=ll_age+dbinom(data$prev_lci[a],1,pprev,log=T)    
      }else if(cutoff1=="uci"){
        ll_age=ll_age+dbinom(data$prev_uci[a],1,pprev,log=T)    
      }
      
    }
    
    #evaluate seasonal seroprevalence within a year
    for(a in 1:nrow(sero.dat)){ 
      #ggplot(prev.ts) + geom_line(aes(x=as.numeric(biwk), y=seroprevalence))
      seroprev <-seroprev.prev.ts$seroprevalence[seroprev.prev.ts$biwk==sero.dat$biwk[a]]  
      
      if(cutoff1=="mean"){
        ll_sero=ll_sero+dbinom(sero.dat$seropos[a],1,seroprev,log=T)    
        #print(ll_sero)
      }else if (cutoff1=="lci"){
        ll_sero=ll_sero+dbinom(sero.dat$seropos_lci[a],1,seroprev,log=T)    
      }else if(cutoff1=="uci"){
        ll_sero=ll_sero+dbinom(sero.dat$seropos_uci[a],1,seroprev,log=T)    
      }
    }
    
    #evaluate shedding prevalence within a year
    for(a in 1:nrow(shed.dat)){ 
      #ggplot(prev.ts) + geom_line(aes(x=as.numeric(biwk), y=prevalence))
      #ggplot(shed.dat) + geom_jitter(aes(x=as.numeric(biwk), y=shed_pos), height=0.1, width =0.1)
      prevseas <-seroprev.prev.ts$prevalence[seroprev.prev.ts$biwk==shed.dat$biwk[a]]  
      
      
      #now run through data and compare with model
      #for-loops over every individual in the dataset, starting above the lower bound - 
      ll_shed=ll_shed+dbinom(shed.dat$shed_pos[a],1,prevseas,log=T)    
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
    
    #and if this ll is better than the previous fit, replace it here
    if(par.df$ll > -ll){
      
      #replace with all of the new values
      par.df$weightPCR = weight.PCR
      
      par.df$beta = beta
      par.df$loss=loss
      par.df$epsilon=epsilon
      par.df$recru=recru
      par.df$recru2=recru2
      par.df$sigma2=sigma2
      par.df$sigma1=sigma1
      par.df$llik = -ll
      
      #replace the point estimate run
      out1 <- out
      
    }
    
    
    
  }
  
  
  
  #return your model stats, including estimated pars, the lci/uci and the output seroprev with CIs
  
  
  #we penalize only for the estimated parameters in AIC
  par.df$k = as.numeric(length(par))
  par.df$AIC = (2*par.df$k + 2*(as.numeric(par.df$llik)))
  
  #check AICc for small sample size -- this should be the combination of the three different dataset sizes
  par.df$n = length(dat.tmp$age) + nrow(sero.dat2) + nrow(shed.dat2)
  par.df$AICc = par.df$AIC + ((2*(par.df$k)*(par.df$k+1))/(par.df$n-par.df$k-1)) 
  
  
  
  if(mu.sick >1 & mu.rec==1){
    model <- paste0(model, "-mu-sick")
  }
  if(mu.sick ==1 & mu.rec >1){
    model <- paste0(model, "-mu-rec")
  }
  if(mu.sick >1 & mu.rec >1){
    model <- paste0(model, "-mu-sick-mu-rec")
  }
  
  par.df$model <- model
  par.df$cutoff <-  cutoff1
  par.df$antigen <- antigen1
  par.df$species <- species1
  
  
  
  #age -seroprev
  
  preds[preds<0] <- 0
  pred.qs  <- apply(X=preds,MARGIN=1,FUN=quantile,probs=c(0.025,0.975), na.rm=TRUE) 
  pred.qs.iqr  <- apply(X=preds,MARGIN=1,FUN=quantile,probs=c(0.25,0.75), na.rm=TRUE) 
  
  #biweekly prev/seroprev
  biwk_preds_seroprev[biwk_preds_seroprev<0] <- 0
  biwk_preds_seroprev.qs  <- apply(X=biwk_preds_seroprev,MARGIN=1,FUN=quantile,probs=c(0.025,0.975), na.rm=TRUE) 
  biwk_preds_seroprev.qs.iqr  <- apply(X=biwk_preds_seroprev,MARGIN=1,FUN=quantile,probs=c(0.25,0.75), na.rm=TRUE) 
  
  biwk_preds_prev[biwk_preds_prev<0] <- 0
  biwk_preds_prev.qs  <- apply(X=biwk_preds_prev,MARGIN=1,FUN=quantile,probs=c(0.025,0.975), na.rm=TRUE) 
  biwk_preds_prev.qs.iqr  <- apply(X=biwk_preds_prev,MARGIN=1,FUN=quantile,probs=c(0.25,0.75), na.rm=TRUE) 
  
  
  #and combine to produce the time series:
  age.sero.out <- out1[[1]]
  
  age.sero.out$seroprev_lci = pred.qs[1,]
  age.sero.out$seroprev_uci = pred.qs[2,]
  
  age.sero.out$seroprev_lci_25 = pred.qs.iqr[1,]
  age.sero.out$seroprev_uci_75 = pred.qs.iqr[2,]
  
  
  age.sero.out$cutoff <-  cutoff1
  age.sero.out$antigen <- antigen1
  age.sero.out$species <- species1
  
  biwk.ts.out <- out1[[2]]
  biwk.ts.out$prevalence <- biwk.ts.out$prevalence*par.df$weightPCR
  
  biwk.ts.out$seroprev_lci = biwk_preds_seroprev.qs[1,]
  biwk.ts.out$seroprev_uci = biwk_preds_seroprev.qs[2,]
  
  biwk.ts.out$seroprev_lci_25 = biwk_preds_seroprev.qs.iqr[1,]
  biwk.ts.out$seroprev_uci_75 = biwk_preds_seroprev.qs.iqr[2,]
  
  biwk.ts.out$prev_lci = biwk_preds_prev.qs[1,]
  biwk.ts.out$prev_uci = biwk_preds_prev.qs[2,]
  
  biwk.ts.out$prev_lci_25 = biwk_preds_prev.qs.iqr[1,]
  biwk.ts.out$prev_uci_75 = biwk_preds_prev.qs.iqr[2,]
  
  
  biwk.ts.out$cutoff <-  cutoff1
  biwk.ts.out$antigen <- antigen1
  biwk.ts.out$species <- species1
  
  
  biwk.ts.out$model <-  model
  age.sero.out$model <- model
  
  
  
  return(list(par.df, age.sero.out, biwk.ts.out)) 
  
  
  
}


slurm_array_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

par.fit <- read.csv(file =  paste0("fit-pars-", slurm_array_id, ".csv"), header = T, stringsAsFactors = F)


out.fit <- wrap.fit.cov(model=par.fit$model,
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
                        sigma1=par.fit$sigma1, #annual rate of loss of adult Abs in MSIRS context (direct from R to S)
                        mu.sick = par.fit$mu_sick, #multiplier x probability annual mortality in infected individuals
                        mu.rec = par.fit$mu_rec, #multiplier x probability annual mortality in recovered individuals
                        recru = par.fit$recru, #annual rate of return to infectiousness
                        recru2 = par.fit$recru2, #annual rate of return to infectiousness from R class with no N class (MSIRI context)
                        epsilon = par.fit$epsilon, #annual rate of return to longterm AB- immunity
                        sigma2 = par.fit$sigma2, #annual rate of loss of adult Abs in all other contexts (from N1 back to S)
                        loss = par.fit$loss,
                        weight.PCR = par.fit$weightPCR,
                        N_stat="matAB",
                        do.plot=FALSE,
                        do.save=FALSE,
                        filename = NA)

filename = paste0(unique(out.fit$model), "-", unique(par.fit$antigen), "-", unique(par.fit$cutoff), slurm_array_id, ".Rdata")
save(out.fit, file=filename)


