rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

#set wd
setwd(paste0(homewd,"/model_fitting/GhV-fits/out-1/"))
#load all the data and add to a list


function.load <- function(prefix, list, i){
  load(paste0(prefix,i,".Rdata")) 
  list[[i]] <- out.fit
  print(list)
  
  return(list)
}

dat.fit <- list()

for (i in 1:988){
  skip_to_next <- FALSE
  tryCatch({dat.fit <- function.load(prefix = "mean", list=dat.fit, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}
dat.fit.df <- data.table::rbindlist(dat.fit, fill = TRUE)



write.csv(dat.fit.df, file=paste0(homewd,"/model_fitting/GhV-fits/all_fits_round1.csv"), row.names = F)

rm(list=ls())


homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

setwd(homewd)

dat <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/all_fits_round1.csv"), stringsAsFactors = F, header = T)

dat = subset(dat, convergence == 0 & !is.na(model))
head(dat)

unique(dat$model)
length(unique(dat$model)) #52 - that is every model represented

#plot the number that completed per model - 19 possible
dat.sum <- ddply(dat,.(model), summarise, N=length(model))


ggplot(dat.sum) + geom_point(aes(model, y=N))

arrange(dat.sum, N) 
# there are 4 models that had < 5 completions (2 only had 1). 
# all 4 were the most complex model "MSIRNSIRNS".
# we tried submitting some re-runs on this model with a longer run time on savio, 
# but it did not improve much, and the model was poor on the whole (on average,
# the worst llik and AIC scores overall)


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
  df <- arrange(df, llik) #arrange by llik, take the top 10% of rows
  df.out = df[1,]
  return(df.out)
}


# #then, compare them
# 
# #rank by delta AIC within a species/antigen combination
 dat.split <- dlply(dat5,.(species, antigen))

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

dat <- data.table::rbindlist(lapply(dat.split, get.deltaAIC))
head(dat)

dat.split <- dlply(dat,.(species, antigen))
dat <- data.table::rbindlist(lapply(dat.split, get.deltallik))
head(dat)

dat.split.mod <- dlply(dat,.(species, antigen, model))
dat.top <- data.table::rbindlist(lapply(dat.split.mod,get.top.mod))
dat.top <- subset(dat.top,!is.na(model))
write.csv(dat.top, paste0(homewd,"/model_fitting/GhV-fits/top_fits_round1.csv"), row.names = F)

#delta llik:
#first, no mu
pA <- ggplot(data=subset(dat, mu_rec==1 & mu_sick==1 & antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,100)) +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik),width=.09, size=1, alpha=.4, color="gray35");pA
#best fit: MSIRS, MSIRNS, MSIRNNS, MSIRNNI,MSIRNI

#mu-rec
pB<- ggplot(data=subset(dat, mu_rec>1 & mu_sick==1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,100)) +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik, color=mu_rec),width=.09, size=1, alpha=.4);pB
#MSIRS, MSIRNS,  MSIRNNS, MSIRNNI,  MSIRNI

#mu-sick
pC<- ggplot(data=subset(dat, mu_rec==1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,100)) +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik, color=mu_sick),width=.09, size=1, alpha=.4);pC
#MSIRS, MSIRNS,  MSIRNNS, MSIRNNI,  MSIRNI

#both mu
pD<- ggplot(data=subset(dat, mu_rec>1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,10)) +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  geom_jitter(aes(x=model, y=delta_llik, size=mu_sick, color=mu_rec),width=.09, alpha=.4) +
  scale_size_continuous(range = c(.01,1));pD
#MSIRS, MSIRNS,   MSIRNNI,  MSIRNI


pGhV <- cowplot::plot_grid(pA,pB,pC,pD, nrow = 2,ncol = 2, align = "hv")



ggsave(file = paste0(homewd,"/model_fitting/GhV-fits/out-1-plots/GhV-fits-llik10.png"),
       plot = pGhV,
       units="mm",  
       width=100, 
       height=70, 
       scale=3, 
       dpi=300)


#delta AIC
#no mu
pA <- ggplot(data=subset(dat, mu_rec==1 & mu_sick==1 & antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,300)) +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_AIC, group=model),  
   #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_AIC),width=.09, size=1, alpha=.4, color="gray35");pA
#MSIRS, MSIRNS,  MSIRNNS, MSIRNNI,  MSIRNI

pB <- ggplot(data=subset(dat, mu_rec>1 & mu_sick==1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,300)) +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_AIC, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_AIC, color=mu_rec),width=.09, size=1, alpha=.4);pB
#MSIRS, MSIRNS,  MSIRNNS, MSIRNNI,  MSIRNI

pC<- ggplot(data=subset(dat, mu_rec==1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,300)) +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_AIC, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_AIC, color=mu_sick),width=.09, size=1, alpha=.4);pC
#MSIRS, MSIRNS,  MSIRNNS, MSIRNNI,  MSIRNI

pD<- ggplot(data=subset(dat, mu_rec>1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,300)) +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  geom_jitter(aes(x=model, y=delta_AIC, size=mu_sick, color=mu_rec),width=.09, alpha=.4) +
  scale_size_continuous(range = c(.01,1));pD
#MSIRS, MSIRNS,  MSIRNNS, MSIRNNI,  MSIRNI

pGhV <- cowplot::plot_grid(pA,pB,pC,pD, nrow = 2,ncol = 2, align = "hv")



ggsave(file = paste0(homewd,"/model_fitting/GhV-fits/out-1-plots/GhV-fits-AIC10.png"),
       plot = pGhV,
       units="mm",  
       width=100, 
       height=70, 
       scale=3, 
       dpi=300)



#sort and save by species antigen
dat <- arrange(dat, species, antigen, delta_AIC)

#take the top performing models:MSIRS, MSIRNNI, MSIRNI, MSIRNNS, MSIRNS
dat.top.mods = subset(dat, model == "MSIRS" |model == "MSIRS-mu-sick" |model == "MSIRS-mu-rec" |model == "MSIRS-mu-sick-mu-rec" |
                           model == "MSIRNNI" |model == "MSIRNNI-mu-sick" |model == "MSIRNNI-mu-rec" |model == "MSIRNNI-mu-sick-mu-rec" |
                          model == "MSIRNI" |model == "MSIRNI-mu-sick" |model == "MSIRNI-mu-rec" |model == "MSIRNI-mu-sick-mu-rec" |
                        model == "MSIRNNS" |model == "MSIRNNS-mu-sick" |model == "MSIRNNS-mu-rec" |model == "MSIRNNS-mu-sick-mu-rec" |
                        model == "MSIRNS" |model == "MSIRNS-mu-sick" |model == "MSIRNS-mu-rec" |model == "MSIRNS-mu-sick-mu-rec" )

#now take the best fit pars for each model

dat.split.mod <- dlply(dat.top.mods, .(species, antigen, model))

get.top.mod<- function(df){
  df <- arrange(df, llik) #arrange by llik, take the top 10% of rows
  df.out = df[1,]
  return(df.out)
}

top.fits <- data.table::rbindlist(lapply(dat.split.mod, get.top.mod))

write.csv(top.fits, file =paste0(homewd,"/model_fitting/GhV-fits/top-fits-round1.csv"), row.names = F)
