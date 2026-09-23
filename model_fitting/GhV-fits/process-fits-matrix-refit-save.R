rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(data.table)

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

#set wd
setwd(paste0(homewd,"/model_fitting/GhV-fits/out-matrix-refit/"))
#load all the data and add to a list


function.load <- function(prefix, list, i){
  load(paste0(prefix,i,".Rdata")) 
  list[[i]] <- out.fit
  print(list)
  
  return(list)
}

dat.fit <- list()

for (i in 1:1000){
  skip_to_next <- FALSE
  tryCatch({dat.fit <- function.load(prefix = "mean", list=dat.fit, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}
dat.fit.df <- data.table:: rbindlist(map(dat.fit, 1))
age.sero.df <- data.table:: rbindlist(map(dat.fit, 2))
biwk.df <- data.table:: rbindlist(map(dat.fit, 3))


write.csv(dat.fit.df, file=paste0(homewd,"/model_fitting/GhV-fits/all_fits_matrix_refit.csv"), row.names = F)
write.csv(age.sero.df, file=paste0(homewd,"/model_fitting/GhV-fits/all_fits_matrix_age_sero_df_refit.csv"), row.names = F)
write.csv(biwk.df, file=paste0(homewd,"/model_fitting/GhV-fits/all_fits_matrix_biwk_df_refit.csv"), row.names = F)

