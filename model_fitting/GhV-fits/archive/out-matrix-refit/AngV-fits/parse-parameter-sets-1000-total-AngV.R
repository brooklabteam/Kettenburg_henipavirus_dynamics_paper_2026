rm(list=ls())

library(plyr)
library(dplyr)

#load data
homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
setwd(homewd)

#load all the parameter inits

par.dat <- read.csv(file= paste0(homewd, "/model_fitting/AngV-fits/fit-pars-init-AngV.csv"), header = T, stringsAsFactors = F)

par.dat <- arrange(par.dat, antigen, cutoff, model, mu_sick, mu_rec)

#and split by row, as each of these is a model
par.dat.list <- split(par.dat, as.numeric(rownames(par.dat)))


save.par.list <- function(df1, rownum, folder_prefix){
  write.csv(df1, file = paste0(folder_prefix,"/fit-pars-",rownum, ".csv"), row.names = FALSE)
}
rep.starts.AngV <- function(df, n.times){
  
  #get your mu.sick and mu.rec status
  if(df$mu_sick >1){
    add.inf.mort= TRUE
  }else{
    add.inf.mort= FALSE
  }
  
  if(df$mu_rec >1){
    add.rec.mort= TRUE
  }else{
    add.rec.mort= FALSE
  }
  #replicate the parameters 1000 times
  big.df <- data.table::rbindlist(replicate(n=n.times, expr = df, simplify = FALSE))
  
  # iterate the starting parameters relevant to this model (some models have more parameters that are estimated)
  # there are a few rules in our parameters
  # sigma cannot be > 26 or <1
  # recru, recru2, epsilon, sigma2 cannot be > 26
  # wane is not useful greater than 450. but can constrain at 200 for starting values
  
  big.df$wane <- sample(x=seq(.01,200, length.out=100000), size=n.times, replace = TRUE)
  
  if (df$model=="MSIR"){
    big.df$beta <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    
  }else if(df$model=="MSIRS" | df$model=="MSIRN"){
    big.df$beta <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma <- sample(x=seq(1,26, length.out=100000), size=n.times, replace = TRUE)
    
  }else if (df$model=="MSIRNI"){
    big.df$beta <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma <- sample(x=seq(1,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$recru <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
  }else if (df$model=="MSIRNS"){
    big.df$beta <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma <- sample(x=seq(1,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$loss <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
  }else if(df$model == "MSIRNNS"){
    big.df$beta <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma <- sample(x=seq(1,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$loss <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$epsilon <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
  }else if(df$model == "MSIRNNI"){
    big.df$beta <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma <- sample(x=seq(1,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$recru <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$epsilon <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    
  }else if (df$model == "MSIRNIRN" | df$model == "MSIRNIR"){
    big.df$beta <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma <- sample(x=seq(1,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$recru <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma2 <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    
  }else if (df$model =="MSIRNSIRN"){
    big.df$beta <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma <- sample(x=seq(1,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$loss <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma2 <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$beta2 <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    
    
  }else if (df$model == "MSIRSIR"){
    big.df$beta <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma <- sample(x=seq(1,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$beta2 <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    
    
  }else if(df$model=="MSIRNIRN-LT"){
    big.df$beta <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma <- sample(x=seq(1,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$recru <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma2 <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$epsilon <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    
  }else if(df$model=="MSIRNSIRN-LT"){
    big.df$beta <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma <- sample(x=seq(1,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$loss <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma2 <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$beta2 <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$epsilon <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    
    
  }else if (df$model =="MSIRNSIRNS"){
    big.df$beta <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma <- sample(x=seq(1,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$loss <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$sigma2 <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$beta2 <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
    big.df$loss2 <- sample(x=seq(.01,26, length.out=100000), size=n.times, replace = TRUE)
  }
  
  big.df$weightPCR <- sample(x=seq(0,1, length.out=100000), size=n.times, replace = TRUE)
  
  
  #mu rec and mu sick can't be <1 and can't be >4 because you'd have 0 survivorship
  if(add.inf.mort ==TRUE){
    big.df$mu_sick <- sample(x=seq(1,4, length.out=1000), size=n.times, replace = TRUE)
  }
  if(add.rec.mort==TRUE){
    big.df$mu_rec <- sample(x=seq(1,4, length.out=1000), size=n.times, replace = TRUE)
  }
  
  #and return big.df
  #compile into one really big df and save to the same folder
  
  return(big.df)
  
}

rep.par.df <- data.table::rbindlist(lapply(par.dat.list, rep.starts.AngV, n.times=19)) #52 lines, each repeated 19 times
rep.list <- split(rep.par.df, as.numeric(rownames(rep.par.df)))

#now save each line as its own csv for fitting
#for 1
mapply(df1=rep.list,rownum=as.list(1:length(rep.list)), FUN=save.par.list, MoreArgs = list(folder_prefix= paste0(homewd,"/model_fitting/AngV-fits/fit-1-AngV")))

#run these, then run more iterations of the best fit models