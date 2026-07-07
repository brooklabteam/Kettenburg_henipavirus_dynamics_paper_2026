rm(list=ls())

library(plyr)
library(dplyr)

#load data
homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
setwd(homewd)

#load all the parameter inits -- these are the top hits from the first round of matrix fitting

par.dat <- read.csv(file= paste0(homewd,"/model_fitting/GhV-fits/top-fits-matrix-refit.csv"), header = T, stringsAsFactors = F)

par.dat <- arrange(par.dat, antigen, cutoff, model, mu_sick, mu_rec)

#and split by row, as each of these is a model
par.dat.list <- split(par.dat, as.numeric(rownames(par.dat)))


save.par.list <- function(df1, rownum, folder_prefix){
  write.csv(df1, file = paste0(folder_prefix,"/fit-pars-",rownum, ".csv"), row.names = FALSE)
}

#now save each line as its own csv for fitting
#for 1
mapply(df1=par.dat.list,rownum=as.list(1:length(par.dat.list)), FUN=save.par.list, MoreArgs = list(folder_prefix= paste0(homewd,"/model_fitting/GhV-fits/run-CIs-matrix")))

