rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

#set wd
setwd(paste0(homewd,"/model_fitting/GhV-fits/out-CIs-matrix/"))
#load all the data and add to a list


function.load <- function(prefix, list, i){
  load(paste0(prefix,i,".Rdata")) 
  list[[i]] <- out.fit
  print(list)
  
  return(list)
}

dat.fit <- list()

for (i in 1:40){
  skip_to_next <- FALSE
  tryCatch({dat.fit <- function.load(prefix = "mean", list=dat.fit, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}

# #now remove any that are blank
# 
# for(i in length(dat.fit):1){
#   if(is.null(dat.fit[[i]])){
#     dat.fit[[i]] <- c()
#   }
# }




#extract first, second, third, item in list
dat.par <- data.table::rbindlist(sapply(dat.fit, '[', 1))
dat.age.sero <- data.table::rbindlist(sapply(dat.fit, '[', 2))
dat.biweekly.ts <- data.table::rbindlist(sapply(dat.fit, '[', 3))

head(dat.par)
head(dat.age.sero)
head(dat.biweekly.ts)

#and save these 
write.csv(dat.par, file=paste0(homewd,"/model_fitting/GhV-fits/par_fits_GhV_CI_final.csv"), row.names = F)
write.csv(dat.age.sero, file=paste0(homewd,"/model_fitting/GhV-fits/age_sero_df_CI_final.csv"), row.names = F)
write.csv(dat.biweekly.ts, file=paste0(homewd,"/model_fitting/GhV-fits/biweekly_ts_GhV_CI_final.csv"), row.names = F)

