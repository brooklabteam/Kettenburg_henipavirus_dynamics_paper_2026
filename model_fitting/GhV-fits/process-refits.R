rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

#set wd
setwd(paste0(homewd,"/model_fitting/GhV-fits/refit-out/"))
#load all the data and add to a list


function.load <- function(prefix, list, i){
  load(paste0(prefix,i,".Rdata")) 
  list[[i]] <- out.fit
  print(list)
  
  return(list)
}

dat.fit <- list()

for (i in 1:1240){
  skip_to_next <- FALSE
  tryCatch({dat.fit <- function.load(prefix = "mean", list=dat.fit, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}
dat.fit.df <- data.table::rbindlist(dat.fit, fill = TRUE)



write.csv(dat.fit.df, file=paste0(homewd,"/model_fitting/GhV-fits/all_fits_refit1.csv"), row.names = F)

rm(list=ls())


homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

setwd(homewd)

dat <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/all_fits_refit1.csv"), stringsAsFactors = F, header = T)

dat = subset(dat, convergence == 0 & !is.na(model))
head(dat)

unique(dat$model)
length(unique(dat$model)) #20 - that is every model represented (5 base models x 4 iterations of each)

#plot the number that completed per model - 19 possible
dat.sum <- ddply(dat,.(model), summarise, N=length(model))


ggplot(dat.sum) + geom_point(aes(model, y=N))

arrange(dat.sum, N) 
# there are 2 models that had < 20 completions but still had 14 and 17 each, 
# so these are thoroughly vetted. 7 models had >50 iterations each


#now bank the top 50 model fits overall by AIC

get.deltaAIC <- function(df){
  minAIC=min(df$AIC)
  df$delta_AIC = df$AIC - minAIC
  return(df)
}
get.deltallik <- function(df){
  minllik=min(df$llik)
  df$delta_llik= df$llik - minllik
  return(df)
}

dat <- get.deltallik(dat)
dat <- get.deltalAIC(dat)

dat <- arrange(dat, delta_AIC)
# and save:
top.fits = dat[1:50,]
write.csv(top.fits, file =paste0(homewd,"/model_fitting/GhV-fits/top-refits-AIC.csv"), row.names = F)

dat <- arrange(dat, delta_llik)
top.fits = dat[1:50,]
write.csv(top.fits, file =paste0(homewd,"/model_fitting/GhV-fits/top-refits-llik.csv"), row.names = F)

dat <- dplyr::select(dat, -(delta_AIC), -(delta_llik))
#rank by delta AIC within a species/antigen/model combination
dat.split <- dlply(dat,.(species, antigen, model))

#select top 5 of fits from each model
get.top5<- function(df){
  df <- arrange(df, llik) #arrange by llik, take the top 10% of rows
  df.out = df[1:5,]
  return(df.out)
}

dat5 <- data.table::rbindlist(lapply(dat.split,get.top5))

#select top from each model

#rank by delta AIC within a species/antigen/model combination
dat.split <- dlply(dat,.(species, antigen, model))

#select top 5 of fits from each model
get.top.mod<- function(df){
  df <- arrange(df, llik) #arrange by llik, take the top model
  df.out = df[1,]
  return(df.out)
}


# #then, compare them
# 
# #rank by delta AIC within a species/antigen combination
dat.split <- dlply(dat5,.(species, antigen))


dat <- data.table::rbindlist(lapply(dat.split, get.deltaAIC))
head(dat)

dat.split <- dlply(dat,.(species, antigen))
dat <- data.table::rbindlist(lapply(dat.split, get.deltallik))
head(dat)

dat.split.mod <- dlply(dat,.(species, antigen, model))
dat.top <- data.table::rbindlist(lapply(dat.split.mod,get.top.mod))
dat.top <- subset(dat.top,!is.na(model))
write.csv(dat.top, paste0(homewd,"/model_fitting/GhV-fits/top_refits_per_modeltype.csv"), row.names = F)

#delta llik:
#first, no mu
pA <- ggplot(data=subset(dat, mu_rec==1 & mu_sick==1 & antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip() +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik),width=.09, size=1, alpha=.4, color="gray35");pA
#best fit:  MSIRNS

#mu-rec
pB<- ggplot(data=subset(dat, mu_rec>1 & mu_sick==1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip() +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik, color=mu_rec),width=.09, size=1, alpha=.4);pB
 #MSIRNS

#mu-sick
pC<- ggplot(data=subset(dat, mu_rec==1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip() +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik, color=mu_sick),width=.09, size=1, alpha=.4);pC
#MSIRS, MSIRNS, MSIRNI

#both mu
pD<- ggplot(data=subset(dat, mu_rec>1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip() +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  geom_jitter(aes(x=model, y=delta_llik, size=mu_sick, color=mu_rec),width=.09, alpha=.4) +
  scale_size_continuous(range = c(.01,1));pD
#MSIRS,  MSIRNI


pGhV <- cowplot::plot_grid(pA,pB,pC,pD, nrow = 2,ncol = 2, align = "hv")



ggsave(file = paste0(homewd,"/model_fitting/GhV-fits/refit-out-plots/GhV-fits-llik10.png"),
       plot = pGhV,
       units="mm",  
       width=100, 
       height=70, 
       scale=3, 
       dpi=300)


#delta AIC
#no mu
pA <- ggplot(data=subset(dat, mu_rec==1 & mu_sick==1 & antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip() +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_AIC, group=model),  
   #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_AIC),width=.09, size=1, alpha=.4, color="gray35");pA
# MSIRNS

pB <- ggplot(data=subset(dat, mu_rec>1 & mu_sick==1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip() +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_AIC, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_AIC, color=mu_rec),width=.09, size=1, alpha=.4);pB
# MSIRNS

pC<- ggplot(data=subset(dat, mu_rec==1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip() +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_AIC, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_AIC, color=mu_sick),width=.09, size=1, alpha=.4);pC
#MSIRS

pD<- ggplot(data=subset(dat, mu_rec>1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip() +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  geom_jitter(aes(x=model, y=delta_AIC, size=mu_sick, color=mu_rec),width=.09, alpha=.4) +
  scale_size_continuous(range = c(.01,1));pD
#MSIRS, MSIRNI

pGhV <- cowplot::plot_grid(pA,pB,pC,pD, nrow = 2,ncol = 2, align = "hv")



ggsave(file = paste0(homewd,"/model_fitting/GhV-fits/refit-out-plots/GhV-fits-AIC10.pdf"),
       plot = pGhV,
       units="mm",  
       width=100, 
       height=70, 
       scale=3, 
       dpi=300)


# 
# #sort and save by species antigen
# dat <- arrange(dat, species, antigen, delta_AIC)
# 
# 
# #now take the best fit pars for each model
# 
# dat.split.mod <- dlply(dat, .(species, antigen, model))
# 
# get.top.mod<- function(df){
#   df <- arrange(df, llik) #arrange by llik, take the top 10% of rows
#   df.out = df[1,]
#   return(df.out)
# }
# 
# top.fits <- data.table::rbindlist(lapply(dat.split.mod, get.top.mod))
