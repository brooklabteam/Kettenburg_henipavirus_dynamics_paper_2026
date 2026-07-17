rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

#set wd
setwd(paste0(homewd,"/model_fitting/GhV-fits/out-CIs/"))


par.dat <- read.csv(file= paste0(homewd,"/model_fitting/GhV-fits/GhV_fits_per_modeltype_CIs.csv"), header = T, stringsAsFactors = F)
par.dat$convergence # all converged

#add in those missing
par.dat2 <- read.csv(file= paste0(homewd,"/model_fitting/GhV-fits/top_refits_per_modeltype.csv"), header = T, stringsAsFactors = F)
tmp <- setdiff(par.dat2$model, par.dat$model)

sub <- list()
for(i in 1:length(tmp)){
  sub[[i]] <- par.dat2[par.dat2$model==tmp[i],]
}
sub <- data.table::rbindlist(sub)

par.dat <- rbind(par.dat, sub)

par.dat <- arrange(par.dat, antigen, cutoff, model, mu_sick, mu_rec)

#and split by row, as each of these is a model
par.dat.list <- split(par.dat, as.numeric(rownames(par.dat)))

# we don't need to shuffle these any more, but we want to split them into 20 separate 
# csvs for fitting the confidence intervals as well

save.par.list <- function(df1, rownum, folder_prefix){
  write.csv(df1, file = paste0(folder_prefix,"/fit-pars-",rownum, ".csv"), row.names = FALSE)
}

#now save each line as its own csv for fitting of CIs
mapply(df1=par.dat.list,rownum=as.list(1:length(par.dat.list)), FUN=save.par.list, MoreArgs = list(folder_prefix= paste0(homewd,"/model_fitting/GhV-fits/run-CIs")))
