
rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(data.table)

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

setwd(homewd)

dat.fit.df.1a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1a.csv"), stringsAsFactors = F, header = T)
dat.fit.df.1b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1b.csv"), stringsAsFactors = F, header = T)
dat.fit.df.2a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2a.csv"), stringsAsFactors = F, header = T)
dat.fit.df.2b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2b.csv"), stringsAsFactors = F, header = T)

age.sero.df.1a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_1a.csv"), stringsAsFactors = F, header = T)
age.sero.df.1b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_1b.csv"), stringsAsFactors = F, header = T)
age.sero.df.2a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_2a.csv"), stringsAsFactors = F, header = T)
age.sero.df.2b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_2b.csv"), stringsAsFactors = F, header = T)


biwk.df.1a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_1a.csv"), stringsAsFactors = F, header = T)
biwk.df.1b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_1b.csv"), stringsAsFactors = F, header = T)
biwk.df.2a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_2a.csv"), stringsAsFactors = F, header = T)
biwk.df.2b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_2b.csv"), stringsAsFactors = F, header = T)



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

for (i in 601:800){
  skip_to_next <- FALSE
  tryCatch({dat.fit.1a <- function.load(prefix = "mean", list=dat.fit.1a, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}
dat.fit.df.1a.add <- data.table:: rbindlist(map(dat.fit.1a, 1))
age.sero.df.1a.add <- data.table:: rbindlist(map(dat.fit.1a, 2))
biwk.df.1a.add <- data.table:: rbindlist(map(dat.fit.1a, 3))

#dat.fit.df.1a <- rbind(dat.fit.1a, dat.fit.df.1a.add)
#age.sero.df.1a <- rbind(age.sero.df.1a, age.sero.df.1a.add)
#biwk.df.1a <- rbind(biwk.df.1a, biwk.df.1a.add)



write.csv(dat.fit.df.1a.add, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1a_add.csv"), row.names = F)
write.csv(age.sero.df.1a.add, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_1a_add.csv"), row.names = F)
write.csv(biwk.df.1a.add, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_1a_add.csv"), row.names = F)


setwd(paste0(homewd,"/model_fitting/GhV-fits/out-matrix-fit-1b/"))

dat.fit.1b <- list()

for (i in 1:200){
  skip_to_next <- FALSE
  tryCatch({dat.fit.1b <- function.load(prefix = "mean", list=dat.fit.1b, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}
dat.fit.df.1b.add <- data.table:: rbindlist(map(dat.fit.1b, 1))
age.sero.df.1b.add <- data.table:: rbindlist(map(dat.fit.1b, 2))
biwk.df.1b.add <- data.table:: rbindlist(map(dat.fit.1b, 3))

# dat.fit.df.1b <- rbind(dat.fit.1b, dat.fit.df.1b.add)
# age.sero.df.1b <- rbind(age.sero.df.1b, age.sero.df.1b.add)
# biwk.df.1b <- rbind(biwk.df.1b, biwk.df.1b.add)
# 

write.csv(dat.fit.df.1b.add, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1b_add.csv"), row.names = F)
write.csv(age.sero.df.1b.add, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_1b_add.csv"), row.names = F)
write.csv(biwk.df.1b.add, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_1b_add.csv"), row.names = F)


setwd(paste0(homewd,"/model_fitting/GhV-fits/out-matrix-fit-2a/"))

dat.fit.2a <- list()

for (i in 601:800){
  skip_to_next <- FALSE
  tryCatch({dat.fit.2a <- function.load(prefix = "mean", list=dat.fit.2a, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}
dat.fit.df.2a.add <- data.table:: rbindlist(map(dat.fit.2a, 1))
age.sero.df.2a.add <- data.table:: rbindlist(map(dat.fit.2a, 2))
biwk.df.2a.add <- data.table:: rbindlist(map(dat.fit.2a, 3))

# dat.fit.df.2a <- rbind(dat.fit.2a, dat.fit.df.2a.add)
# age.sero.df.2a <- rbind(age.sero.df.2a, age.sero.df.2a.add)
# biwk.df.2a <- rbind(biwk.df.2a, biwk.df.2a.add)

write.csv(dat.fit.df.2a.add, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2a_add.csv"), row.names = F)
write.csv(age.sero.df.2a.add, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_2a_add.csv"), row.names = F)
write.csv(biwk.df.2a.add, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_2a_add.csv"), row.names = F)


setwd(paste0(homewd,"/model_fitting/GhV-fits/out-matrix-fit-2b/"))

dat.fit.2b <- list()

for (i in 1:200){
  skip_to_next <- FALSE
  tryCatch({dat.fit.2b <- function.load(prefix = "mean", list=dat.fit.2b, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}
dat.fit.df.2b.add <- data.table:: rbindlist(map(dat.fit.2b, 1))
age.sero.df.2b.add <- data.table:: rbindlist(map(dat.fit.2b, 2))
biwk.df.2b.add <- data.table:: rbindlist(map(dat.fit.2b, 3))

# dat.fit.df.2b <- rbind(dat.fit.2b, dat.fit.df.2b.add)
# age.sero.df.2b <- rbind(age.sero.df.2b, age.sero.df.2b.add)
# biwk.df.2b <- rbind(biwk.df.2b, biwk.df.2b.add)



write.csv(dat.fit.df.2b.add, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2b_add.csv"), row.names = F)
write.csv(age.sero.df.2b.add, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_2b_add.csv"), row.names = F)
write.csv(biwk.df.2b.add, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_2b_add.csv"), row.names = F)

rm(list=ls())


homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

setwd(homewd)

dat.fit.df.1a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1a.csv"), stringsAsFactors = F, header = T)
dat.fit.df.1b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1b.csv"), stringsAsFactors = F, header = T)
dat.fit.df.2a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2a.csv"), stringsAsFactors = F, header = T)
dat.fit.df.2b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2b.csv"), stringsAsFactors = F, header = T)

age.sero.df.1a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_1a.csv"), stringsAsFactors = F, header = T)
age.sero.df.1b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_1b.csv"), stringsAsFactors = F, header = T)
age.sero.df.2a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_2a.csv"), stringsAsFactors = F, header = T)
age.sero.df.2b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_2b.csv"), stringsAsFactors = F, header = T)


biwk.df.1a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_1a.csv"), stringsAsFactors = F, header = T)
biwk.df.1b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_1b.csv"), stringsAsFactors = F, header = T)
biwk.df.2a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_2a.csv"), stringsAsFactors = F, header = T)
biwk.df.2b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_2b.csv"), stringsAsFactors = F, header = T)




dat.fit.df.1a.add <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1a_add.csv"), stringsAsFactors = F, header = T)
age.sero.df.1a.add <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_1a_add.csv"), stringsAsFactors = F, header = T)
biwk.df.1a.add <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_1a_add.csv"), stringsAsFactors = F, header = T)


dat.fit.df.1b.add <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1b_add.csv"), stringsAsFactors = F, header = T)
age.sero.df.1b.add <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_1b_add.csv"), stringsAsFactors = F, header = T)
biwk.df.1b.add <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_1b_add.csv"), stringsAsFactors = F, header = T)


dat.fit.df.2a.add <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2a_add.csv"), stringsAsFactors = F, header = T)
age.sero.df.2a.add <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_2a_add.csv"), stringsAsFactors = F, header = T)
biwk.df.2a.add <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_2a_add.csv"), stringsAsFactors = F, header = T)

dat.fit.df.2b.add <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2b_add.csv"), stringsAsFactors = F, header = T)
age.sero.df.2b.add <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_2b_add.csv"), stringsAsFactors = F, header = T)
biwk.df.2b.add <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_2b_add.csv"), stringsAsFactors = F, header = T)


dat.fit.df.1a <- rbind(dat.fit.df.1a, dat.fit.df.1a.add)
age.sero.df.1a <- rbind(age.sero.df.1a, age.sero.df.1a.add)
biwk.df.1a <- rbind(biwk.df.1a, biwk.df.1a.add)



dat.fit.df.1b <- rbind(dat.fit.df.1b, dat.fit.df.1b.add)
age.sero.df.1b <- rbind(age.sero.df.1b, age.sero.df.1b.add)
biwk.df.1b <- rbind(biwk.df.1b, biwk.df.1b.add)



dat.fit.df.2a <- rbind(dat.fit.df.2a, dat.fit.df.2a.add)
age.sero.df.2a <- rbind(age.sero.df.2a, age.sero.df.2a.add)
biwk.df.2a <- rbind(biwk.df.2a, biwk.df.2a.add)



dat.fit.df.2b <- rbind(dat.fit.df.2b, dat.fit.df.2b.add)
age.sero.df.2b <- rbind(age.sero.df.2b, age.sero.df.2b.add)
biwk.df.2b <- rbind(biwk.df.2b, biwk.df.2b.add)



#save

write.csv(dat.fit.df.1a, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1a.csv"), row.names = F)
write.csv(age.sero.df.1a, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_1a.csv"), row.names = F)
write.csv(biwk.df.1a, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_1a.csv"), row.names = F)



write.csv(dat.fit.df.1b, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1b.csv"), row.names = F)
write.csv(age.sero.df.1b, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_1b.csv"), row.names = F)
write.csv(biwk.df.1b, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_1b.csv"), row.names = F)



write.csv(dat.fit.df.2a, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2a.csv"), row.names = F)
write.csv(age.sero.df.2a, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_2a.csv"), row.names = F)
write.csv(biwk.df.2a, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_2a.csv"), row.names = F)



write.csv(dat.fit.df.2b, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2b.csv"), row.names = F)
write.csv(age.sero.df.2b, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_age_sero_df_2b.csv"), row.names = F)
write.csv(biwk.df.2b, file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_biwk_df_2b.csv"), row.names = F)


