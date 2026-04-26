rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

#set wd
setwd(paste0(homewd,"/model_fitting/GhV-fits/out-CIs/"))
#load all the data and add to a list


function.load <- function(prefix, list, i){
  load(paste0(prefix,i,".Rdata")) 
  list[[i]] <- out.fit
  print(list)
  
  return(list)
}

dat.fit <- list()

for (i in 1:20){
  skip_to_next <- FALSE
  tryCatch({dat.fit <- function.load(prefix = "mean", list=dat.fit, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}
dat.fit.df <- data.table::rbindlist(dat.fit, fill = TRUE)



write.csv(dat.fit.df, file=paste0(homewd,"/model_fitting/GhV-fits/GhV_fits_per_modeltype_CIs.csv"), row.names = F)
