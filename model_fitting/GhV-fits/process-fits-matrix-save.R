rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(data.table)

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

#set wd
setwd(paste0(homewd,"/model_fitting/GhV-fits/out-matrix-fit-1a/"))
#load all the data and add to a list


function.load <- function(prefix, list, i){
  load(paste0(prefix,i,".Rdata")) 
  list[[i]] <- out.fit
  print(list)
  
  return(list)
}

dat.fit.1a <- list()

for (i in 1:1000){
  skip_to_next <- FALSE
  tryCatch({dat.fit.1a <- function.load(prefix = "mean", list=dat.fit.1a, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}
dat.fit.df.1a <- data.table:: rbindlist(map(dat.fit.1a, 1))
age.sero.df.1a <- data.table:: rbindlist(map(dat.fit.1a, 2))
biwk.df.1a <- data.table:: rbindlist(map(dat.fit.1a, 3))

#for now, skip the data that came with (elements 2 and 3 of each big list)

write.csv(dat.fit.df.1a, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1a.csv"), row.names = F)
write.csv(age.sero.df.1a, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_1a.csv"), row.names = F)
write.csv(biwk.df.1a, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_1a.csv"), row.names = F)



setwd(paste0(homewd,"/model_fitting/GhV-fits/out-matrix-fit-1b/"))

dat.fit.1b <- list()

for (i in 201:1000){
  skip_to_next <- FALSE
  tryCatch({dat.fit.1b <- function.load(prefix = "mean", list=dat.fit.1b, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}
dat.fit.df.1b <- data.table:: rbindlist(map(dat.fit.1b, 1))
age.sero.df.1b <- data.table:: rbindlist(map(dat.fit.1b, 2))
biwk.df.1b <- data.table:: rbindlist(map(dat.fit.1b, 3))

#for now, skip the data that came with (elements 2 and 3 of each big list)

write.csv(dat.fit.df.1b, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1b.csv"), row.names = F)
write.csv(age.sero.df.1b, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_1b.csv"), row.names = F)
write.csv(biwk.df.1b, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_1b.csv"), row.names = F)



setwd(paste0(homewd,"/model_fitting/GhV-fits/out-matrix-fit-2a/"))

dat.fit.2a <- list()

for (i in 1:1000){
  skip_to_next <- FALSE
  tryCatch({dat.fit.2a <- function.load(prefix = "mean", list=dat.fit.2a, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}
dat.fit.df.2a <- data.table:: rbindlist(map(dat.fit.2a, 1))
age.sero.df.2a <- data.table:: rbindlist(map(dat.fit.2a, 2))
biwk.df.2a <- data.table:: rbindlist(map(dat.fit.2a, 3))

#for now, skip the data that came with (elements 2 and 3 of each big list)

write.csv(dat.fit.df.2a, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2a.csv"), row.names = F)
write.csv(age.sero.df.2a, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_2a.csv"), row.names = F)
write.csv(biwk.df.2a, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_2a.csv"), row.names = F)



setwd(paste0(homewd,"/model_fitting/GhV-fits/out-matrix-fit-2b/"))

dat.fit.2b <- list()

for (i in 201:1000){
  skip_to_next <- FALSE
  tryCatch({dat.fit.2b <- function.load(prefix = "mean", list=dat.fit.2b, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}
dat.fit.df.2b <- data.table:: rbindlist(map(dat.fit.2b, 1))
age.sero.df.2b <- data.table:: rbindlist(map(dat.fit.2b, 2))
biwk.df.2b <- data.table:: rbindlist(map(dat.fit.2b, 3))

#for now, skip the data that came with (elements 2 and 3 of each big list)

write.csv(dat.fit.df.2b, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2b.csv"), row.names = F)
write.csv(age.sero.df.2b, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_2b.csv"), row.names = F)
write.csv(biwk.df.2b, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_2b.csv"), row.names = F)


 