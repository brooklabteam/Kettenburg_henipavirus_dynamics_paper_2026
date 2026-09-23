rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(data.table)

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

setwd(homewd)

#972
dat <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/all_fits_matrix_refit.csv"), stringsAsFactors = F, header = T)
head(dat)


#769
dat = subset(dat, convergence == 0 & !is.na(dat$model))

unique(dat$model)
length(unique(dat$model)) #20 - that is every model represented in the refit subset

#plot the number that completed per model - 50 possible
dat.sum <- ddply(dat,.(model), summarise, N=length(model)) #18 - 50 for each


ggplot(dat.sum) + geom_point(aes(model, y=N))

arrange(dat.sum, N) #most >30, one only 14, one only 23



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
write.csv(dat.top, paste0(homewd,"/model_fitting/GhV-fits/top_fits_per_model_matrix_refit.csv"), row.names = F)

arrange(dat.top, llik)
#MSIRSI-mu-rec, MSIRSI-mu-sick-mu-rec, MSIRSI-mu-sick

#delta llik:
#first, no mu
pA <- ggplot(data=subset(dat, mu_rec==1 & mu_sick==1 & antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip() +#ylim=c(0,25)
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik),width=.09, size=1, alpha=.4, color="gray35");pA
#best fit: MSIRSI, MSIRNNSI, MSIRNSI, MSIRS, MSIRNS


#mu-rec
pB<- ggplot(data=subset(dat, mu_rec>1 & mu_sick==1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip() +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik, color=mu_rec),width=.09, size=1, alpha=.4);pB
#MSIRSI, MSIRNSI, MSIRNNSI,  MSIRS, MSIRNS

#mu-sick
pC<- ggplot(data=subset(dat, mu_rec==1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip() +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik, color=mu_sick),width=.09, size=1, alpha=.4);pC
#MSIRSI, MSIRNNSI,  MSIRNSI,  MSIRS, MSIRNS

#both mu
pD<- ggplot(data=subset(dat, mu_rec>1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  geom_jitter(aes(x=model, y=delta_llik, size=mu_sick, color=mu_rec),width=.09, alpha=.4) +
  scale_size_continuous(range = c(.01,1));pD
#MSIRSI, MSIRNNSI,  MSIRNSI,  MSIRS, MSIRNS


pGhV <- cowplot::plot_grid(pA,pB,pC,pD, nrow = 2,ncol = 2, align = "hv")



# ggsave(file = paste0(homewd,"/model_fitting/GhV-fits/out-1-plots/GhV-fits-llik10.png"),
#        plot = pGhV,
#        units="mm",  
#        width=100, 
#        height=70, 
#        scale=3, 
#        dpi=300)


#delta AIC
#no mu
pA <- ggplot(data=subset(dat, mu_rec==1 & mu_sick==1 & antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_AIC, group=model),  
   #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_AIC),width=.09, size=1, alpha=.4, color="gray35");pA
#MSIRSI,  MSIRS,  MSIRNSI,  MSIRNS, MSIRNNSI, 

pB <- ggplot(data=subset(dat, mu_rec>1 & mu_sick==1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_AIC, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_AIC, color=mu_rec),width=.09, size=1, alpha=.4);pB
#MSIRSI,  MSIRS,  MSIRNSI,  MSIRNS, MSIRNNSI, 

pC<- ggplot(data=subset(dat, mu_rec==1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_AIC, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_AIC, color=mu_sick),width=.09, size=1, alpha=.4);pC
#MSIRSI,  MSIRS,  MSIRNS, MSIRNSI,   MSIRNNSI, 

pD<- ggplot(data=subset(dat, mu_rec>1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  geom_jitter(aes(x=model, y=delta_AIC, size=mu_sick, color=mu_rec),width=.09, alpha=.4) +
  scale_size_continuous(range = c(.01,1));pD
#MSIRS, MSIRNS,  MSIRNNS, MSIRNNI,  MSIRNI

pGhV <- cowplot::plot_grid(pA,pB,pC,pD, nrow = 2,ncol = 2, align = "hv")



# ggsave(file = paste0(homewd,"/model_fitting/GhV-fits/out-1-plots/GhV-fits-AIC10.pdf"),
#        plot = pGhV,
#        units="mm",  
#        width=100, 
#        height=70, 
#        scale=3, 
#        dpi=300)


#and compare with round 1 - replace if necessary
round1.dat <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/top_fits_per_model_matrix_round1.csv"), header = T, stringsAsFactors = F)

top.mods1 = subset(round1.dat, model == "MSIRS" |model == "MSIRS-mu-sick" |model == "MSIRS-mu-rec" |model == "MSIRS-mu-sick-mu-rec" |
                                       model == "MSIRSI" |model == "MSIRSI-mu-sick" |model == "MSIRSI-mu-rec" |model == "MSIRSI-mu-sick-mu-rec" |
                                       model == "MSIRNSI" |model == "MSIRNSI-mu-sick" |model == "MSIRNSI-mu-rec" |model == "MSIRNSI-mu-sick-mu-rec" |
                                       model == "MSIRNNSI" |model == "MSIRNNSI-mu-sick" |model == "MSIRNNSI-mu-rec" |model == "MSIRNNSI-mu-sick-mu-rec" |
                                       model == "MSIRNS" |model == "MSIRNS-mu-sick" |model == "MSIRNS-mu-rec" |model == "MSIRNS-mu-sick-mu-rec" )

other.mods1 = subset(round1.dat, model != "MSIRS" &model != "MSIRS-mu-sick" &model != "MSIRS-mu-rec" &model != "MSIRS-mu-sick-mu-rec" &
                     model != "MSIRSI" &model != "MSIRSI-mu-sick" &model != "MSIRSI-mu-rec" &model != "MSIRSI-mu-sick-mu-rec" &
                     model != "MSIRNSI" &model != "MSIRNSI-mu-sick" &model != "MSIRNSI-mu-rec" &model != "MSIRNSI-mu-sick-mu-rec" &
                     model != "MSIRNNSI" &model != "MSIRNNSI-mu-sick" &model != "MSIRNNSI-mu-rec" &model != "MSIRNNSI-mu-sick-mu-rec" &
                     model != "MSIRNS" &model != "MSIRNS-mu-sick" &model != "MSIRNS-mu-rec" &model != "MSIRNS-mu-sick-mu-rec" )



top.mods1 <- arrange(top.mods1, model)
top.mods1$run <- "round1"
dat.top <- arrange(dat.top, model)
dat.top$run <- "refit"

#now compare and take the best one of the top fits
top.list <- dlply(dat.top, .(model))
mod1.list <- dlply(top.mods1, .(model))
 
comp.df <- function(df1, df2){
  df.comp <- rbind(df1, df2)
  df.comp <- arrange(df.comp, llik)
  df.out <- df.comp[1,]
  return(df.out)
}

out.list <- mapply(comp.df, top.list, mod1.list, SIMPLIFY = FALSE)

out.df <- data.table::rbindlist(out.list) 
subset(out.df, run=="refit") #most were best in round 1 except MSIRS, MSIRSI, MSIRSI-mu-rec
#these three need CIs to be rerun

out.df <- dplyr::select(out.df, -(run))
#and combine with the original dataset
out.df.all <- rbind(out.df, other.mods1)

#and save
write.csv(out.df.all, file = paste0(homewd,"/model_fitting/GhV-fits/final_best_fit_pars_per_model_matrix.csv"), row.names = F)

