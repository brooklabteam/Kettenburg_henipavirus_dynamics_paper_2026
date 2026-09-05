rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(data.table)

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

#set wd
setwd(paste0(homewd,"/model_fitting/GhV-fits/"))
#load all the data and add to a list

#load the partial datasets and curate

dat.fit.df.1a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1a.csv"), stringsAsFactors = F, header = T)
dat.fit.df.1b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_1b.csv"), stringsAsFactors = F, header = T)
dat.fit.df.2a <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2a.csv"), stringsAsFactors = F, header = T)
dat.fit.df.2b <- read.csv(file=paste0(homewd,"/model_fitting/GhV-fits/fits_matrix_2b.csv"), stringsAsFactors = F, header = T)

#combine
dat <- rbind(dat.fit.df.1a, dat.fit.df.1b, dat.fit.df.2a, dat.fit.df.2b)
head(dat)
unique(dat$llik)
nrow(dat) #3929

dat = subset(dat, convergence == 0 & !is.na(model))

#save all fits
write.csv(dat, file = paste0(homewd,"/model_fitting/GhV-fits/all_fits_matrix.csv"), row.names = F)
# 3348

unique(dat$model)
length(unique(dat$model)) #40 = every model represented

#plot the number that completed per model - 62 possible
dat.sum <- ddply(dat,.(model), summarise, N=length(model)) 
arrange(dat.sum, N) #>=60 for 37 models; >43 for all

ggplot(dat.sum) + geom_point(aes(model, y=N))


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
write.csv(dat.top, paste0(homewd,"/model_fitting/GhV-fits/top_fits_per_model_matrix_round1.csv"), row.names = F)


#delta llik:
#first, no mu
pA <- ggplot(data=subset(dat, mu_rec==1 & mu_sick==1 & antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik),width=.09, size=1, alpha=.4, color="gray35");pA
#best fit: MSIRSI, MSIRNNSI, MSIRNSI, MSIRS, MSIRNS, MSIRNNS. MSIRN performs worst!

#mu-rec
pB<- ggplot(data=subset(dat, mu_rec>1 & mu_sick==1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik, color=mu_rec),width=.09, size=1, alpha=.4);pB
#MSIRSI, MSIRNSI, MSIRNNSI, MSIRS, MSIRNS, MSIRNNS. MSIRN performs worst

#mu-sick
pC<- ggplot(data=subset(dat, mu_rec==1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik, color=mu_sick),width=.09, size=1, alpha=.4);pC
#MSIRSI, MSIRNNSI, MSIRNSI, MSIRS, MSIRNS, MSIRNNS. MSIRN performs worst

#both mu
pD<- ggplot(data=subset(dat, mu_rec>1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  geom_jitter(aes(x=model, y=delta_llik, size=mu_sick, color=mu_rec),width=.09, alpha=.4) +
  scale_size_continuous(range = c(.01,1));pD
#MSIRSI, MSIRNSI, MSIRNNSI, MSIRS, MSIRNS, MSIRNI. MSIRN performs worst


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
#best fit: MSIRSI, MSIRS,  MSIRNSI, MSIRNS, MSIRNNSI, MSIRNNS. MSIRN performs worst!

pB <- ggplot(data=subset(dat, mu_rec>1 & mu_sick==1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_AIC, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_AIC, color=mu_rec),width=.09, size=1, alpha=.4);pB
#MSIRSI, MSIRS, MSIRNS, MSIRNSI, MSIRNNSI, MSIRNNS. MSIRN performs worst

pC<- ggplot(data=subset(dat, mu_rec==1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_AIC, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_AIC, color=mu_sick),width=.09, size=1, alpha=.4);pC
#MSIRSI, MSIRS, MSIRNS,MSIRNSI,MSIRNNSI, MSIRNNS. MSIRN performs worst

pD<- ggplot(data=subset(dat, mu_rec>1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  geom_jitter(aes(x=model, y=delta_AIC, size=mu_sick, color=mu_rec),width=.09, alpha=.4) +
  scale_size_continuous(range = c(.01,1));pD
#MSIRSI, MSIRS, MSIRNSI, MSIRNS, MSIRNNSI, MSIRNI. MSIRN performs worst

pGhV <- cowplot::plot_grid(pA,pB,pC,pD, nrow = 2,ncol = 2, align = "hv")


# ggsave(file = paste0(homewd,"/model_fitting/GhV-fits/out-1-plots/GhV-fits-AIC10.pdf"),
#        plot = pGhV,
#        units="mm",  
#        width=100, 
#        height=70, 
#        scale=3, 
#        dpi=300)



#sort and save by species antigen
dat <- arrange(dat, species, antigen, delta_llik)


#model ranking by llik is:
#MSIRSI-mu-sick-mu-rec, MSIRSI-mu-rec, MSIRSI-mu-sick, MSIRSI,
#next comes MSIRNSI (mu-rec, then mu-sick-mu-rec); MSIRNNSI (mu-rec, then mu-sick-mu-rec, then plain);
#then MSIRS-mu-rec, MSIRS-mu-sick-mu-rec,
#then MSIRNSI, MSIRNS-mu-rec, MSIRNNSI-mu-sick....


#take the top 5 performing models and all constructions of each: 
#MSIRSI, MSIRNSI, MSIRNNSI, MSIRS, MSIRNS

dat.top.mods = subset(dat, model == "MSIRS" |model == "MSIRS-mu-sick" |model == "MSIRS-mu-rec" |model == "MSIRS-mu-sick-mu-rec" |
                            model == "MSIRSI" |model == "MSIRSI-mu-sick" |model == "MSIRSI-mu-rec" |model == "MSIRSI-mu-sick-mu-rec" |
                           model == "MSIRNSI" |model == "MSIRNSI-mu-sick" |model == "MSIRNSI-mu-rec" |model == "MSIRNSI-mu-sick-mu-rec" |
                         model == "MSIRNNSI" |model == "MSIRNNSI-mu-sick" |model == "MSIRNNSI-mu-rec" |model == "MSIRNNSI-mu-sick-mu-rec" |
                         model == "MSIRNS" |model == "MSIRNS-mu-sick" |model == "MSIRNS-mu-rec" |model == "MSIRNS-mu-sick-mu-rec" )

#now take the best fit pars for each model

dat.split.mod <- dlply(dat.top.mods, .(species, antigen, model))

get.top.mod<- function(df){
  df <- arrange(df, llik) #arrange by llik, take the top 10% of rows
  df.out = df[1,]
  return(df.out)
}

top.fits <- data.table::rbindlist(lapply(dat.split.mod, get.top.mod))

write.csv(top.fits, file =paste0(homewd,"/model_fitting/GhV-fits/top-fits-matrix-narrow-round1.csv"), row.names = F)


#then, also plot the top 20 parameter sets for each of the top 5 models and compare
rm(list=ls())


homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

setwd(homewd)

#3348
dat <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/all_fits_matrix.csv"), stringsAsFactors = F, header = T)

#3348
dat = subset(dat, convergence == 0 & !is.na(model))

unique(dat$model)
length(unique(dat$model)) #40 - that is every model represented

#rank by llik within a species/antigen/model combination
dat.split <- dlply(dat,.(species, antigen, model))

#select top 20 of fits from each model
get.top10<- function(df){
  df <- arrange(df, llik) #arrange by llik, take the top 10% of rows
  df.out = df[1:10,]
  return(df.out)
}


dat10 <- data.table::rbindlist(lapply(dat.split,get.top10))

#now rank by model and llik and save
dat10 <- arrange(dat10, model, llik)


write.csv(dat10, file =paste0(homewd,"/model_fitting/GhV-fits/top-10-matrix-round1.csv"), row.names = F)

#and take these and plot the range of parameter values for the parameters that are fitted for each

plot.one <- function(df){
  if(unique(df$model)== "MSIRS"){
    plot.df <- dplyr::select(df, model, beta, sigma1)
  }else if(unique(df$model)== "MSIRS-mu-sick"){
    plot.df <- dplyr::select(df, model, beta, sigma1, mu_sick)
  }else if(unique(df$model)== "MSIRS-mu-rec"){
    plot.df <- dplyr::select(df, model, beta, sigma1, mu_rec)
  }else if(unique(df$model)== "MSIRS-mu-sick-mu-rec"){
    plot.df <- dplyr::select(df, model, beta, sigma1, mu_sick, mu_rec)
  }else if (unique(df$model)=="MSIRNI"){
    plot.df <- dplyr::select(df, model, beta, sigma2, recru)
  }else if (unique(df$model)=="MSIRNI-mu-sick"){
    plot.df <- dplyr::select(df, model, beta, sigma2, recru, mu_sick)
  }else if (unique(df$model)=="MSIRNI-mu-rec"){
    plot.df <- dplyr::select(df, model, beta, sigma2, recru, mu_rec)
  }else if (unique(df$model)=="MSIRNI-mu-sick-mu-rec"){
    plot.df <- dplyr::select(df, model, beta, sigma2, recru, mu_sick, mu_rec)
  }else if (unique(df$model)=="MSIRNS"){
    plot.df <- dplyr::select(df, model, beta, sigma2, loss)
  }else if (unique(df$model)=="MSIRNS-mu-sick"){
    plot.df <- dplyr::select(df, model, beta, sigma2, loss, mu_sick)
  }else if (unique(df$model)=="MSIRNS-mu-rec"){
    plot.df <- dplyr::select(df, model, beta, sigma2, loss, mu_rec)
  }else if (unique(df$model)=="MSIRNS-mu-sick-mu-rec"){
    plot.df <- dplyr::select(df, model, beta, sigma2, loss, mu_sick, mu_rec)
  }else if(unique(df$model) == "MSIRNNS"){
    plot.df <- dplyr::select(df, model, beta, sigma2, loss, epsilon)
  }else if(unique(df$model) == "MSIRNNS-mu-sick"){
    plot.df <- dplyr::select(df, model, beta, sigma2, loss, epsilon, mu_sick)
  }else if(unique(df$model) == "MSIRNNS-mu-rec"){
    plot.df <- dplyr::select(df, model, beta, sigma2, loss, epsilon, mu_rec)
  }else if(unique(df$model) == "MSIRNNS-mu-sick-mu-rec"){
    plot.df <- dplyr::select(df, model, beta, sigma2, loss, epsilon, mu_sick, mu_rec)
  }else if(unique(df$model) == "MSIRNNI"){
    plot.df <- dplyr::select(df, model,beta, sigma2, recru, epsilon)
  }else if(unique(df$model) == "MSIRNNI-mu-sick"){
    plot.df <- dplyr::select(df, model,beta, sigma2, recru, epsilon, mu_sick)
  }else if(unique(df$model) == "MSIRNNI-mu-rec"){
    plot.df <- dplyr::select(df, model,beta, sigma2, recru, epsilon, mu_rec)
  }else if(unique(df$model) == "MSIRNNI-mu-sick-mu-rec"){
    plot.df <- dplyr::select(df, model,beta, sigma2, recru, epsilon, mu_sick, mu_rec)
  }else if(unique(df$model) == "MSIRNSI"){
    plot.df <- dplyr::select(df, model,beta, sigma2, loss, recru)
  }else if(unique(df$model) == "MSIRNSI-mu-sick"){
    plot.df <- dplyr::select(df, model,beta, sigma2, loss, recru, mu_sick)
  }else if(unique(df$model) == "MSIRNSI-mu-rec"){
    plot.df <- dplyr::select(df, model,beta, sigma2, loss, recru, mu_rec)
  }else if(unique(df$model) == "MSIRNSI-mu-sick-mu-rec"){
    plot.df <- dplyr::select(df, model,beta, sigma2, loss, recru, mu_sick, mu_rec)
  }else if(unique(df$model) == "MSIRSI"){
    plot.df <- dplyr::select(df, model,beta, sigma1, recru2)
  }else if(unique(df$model) == "MSIRSI-mu-sick"){
    plot.df <- dplyr::select(df, model,beta, sigma1, recru2, mu_sick)
  }else if(unique(df$model) == "MSIRSI-mu-rec"){
    plot.df <- dplyr::select(df, model,beta, sigma1, recru2, mu_rec)
  }else if(unique(df$model) == "MSIRSI-mu-sick-mu-rec"){
    plot.df <- dplyr::select(df, model,beta, sigma1, recru2, mu_sick, mu_rec)
  }else if(unique(df$model) == "MSIRNNSI"){
    plot.df <- dplyr::select(df, model,beta, sigma2, recru, loss, epsilon)
  }else if(unique(df$model) == "MSIRNNSI-mu-sick"){
    plot.df <- dplyr::select(df, model,beta, sigma2, recru, loss, epsilon, mu_sick)
  }else if(unique(df$model) == "MSIRNNSI-mu-rec"){
    plot.df <- dplyr::select(df, model,beta, sigma2, recru, loss, epsilon, mu_rec)
  }else if(unique(df$model) == "MSIRNNSI-mu-sick-mu-rec"){
    plot.df <- dplyr::select(df, model,beta, sigma2, recru, loss, epsilon, mu_sick, mu_rec)
  }else if (unique(df$model)=="MSIRI"){
    plot.df <- dplyr::select(df, model, beta, recru2)
  }else if (unique(df$model)=="MSIRI-mu-sick"){
    plot.df <- dplyr::select(df, model, beta, recru2, mu_sick)
  }else if (unique(df$model)=="MSIRI-mu-rec"){
    plot.df <- dplyr::select(df, model, beta, recru2, mu_rec)
  }else if (unique(df$model)=="MSIRI-mu-sick-mu-rec"){
    plot.df <- dplyr::select(df, model, beta, recru2, mu_sick, mu_rec)
  }else if (unique(df$model)=="MSIRN"){
    plot.df <- dplyr::select(df, model, beta, sigma2)
  }else if (unique(df$model)=="MSIRN-mu-sick"){
    plot.df <- dplyr::select(df, model, beta, sigma2, mu_sick)
  }else if (unique(df$model)=="MSIRN-mu-rec"){
    plot.df <- dplyr::select(df, model, beta, sigma2, mu_rec)
  }else if (unique(df$model)=="MSIRN-mu-sick-mu-rec"){
    plot.df <- dplyr::select(df, model, beta, sigma2, mu_sick, mu_rec)
  }
  
  
    
  
  
  plot.long <-  reshape(plot.df, varying = c(2:ncol(plot.df)), v.names = "value", timevar = "parameter", times= c(names(plot.df)[2:ncol(plot.df)]), direction = "long")
  
  p1 <- ggplot(plot.long) +
    geom_violin(aes(x=parameter, y=value), draw_quantiles =c(.25,.5,.75)) + 
    facet_wrap(parameter~., scales="free") + geom_jitter(aes(x=parameter, y=value), width = .1) +
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
                       axis.title = element_blank(), axis.text = element_text(size=14),
                       strip.text = element_text(size=14)) + ggtitle(unique(plot.df$model))
  return(p1)
}

plot.top.10.fits <- function(df, filename1, filename2, filename3, filename4){
  df.list <- dlply(df, .(model))
  
  plot.list  <- lapply(df.list, plot.one)

  all.par1 <- cowplot::plot_grid(plotlist = plot.list[1:12], ncol = 4, nrow =3)
  all.par2 <- cowplot::plot_grid(plotlist = plot.list[13:24], ncol = 4, nrow =3)
  all.par3 <- cowplot::plot_grid(plotlist = plot.list[25:36], ncol = 4, nrow =3)
  all.par4 <- cowplot::plot_grid(plotlist = plot.list[37:40], ncol = 4, nrow =2)
  
  #and save
  ggsave(file = filename1,
              plot = all.par1,
              units="mm",  
              width=150, 
              height=150, 
              scale=3, 
              dpi=300)
  
  ggsave(file = filename2,
         plot = all.par2,
         units="mm",  
         width=150, 
         height=150, 
         scale=3, 
         dpi=300)
  
  
  ggsave(file = filename3,
         plot = all.par3,
         units="mm",  
         width=150, 
         height=150, 
         scale=3, 
         dpi=300)
  
  
  ggsave(file = filename4,
         plot = all.par4,
         units="mm",  
         width=150, 
         height=150, 
         scale=3, 
         dpi=300)
         
}

plot.top.10.fits(df = dat10, 
                 filename1 = paste0(homewd, "/model_fitting/GhV-fits/figures/best-fit-pars-1.png"), 
                 filename2 = paste0(homewd, "/model_fitting/GhV-fits/figures/best-fit-pars-2.png"),
                 filename3 = paste0(homewd, "/model_fitting/GhV-fits/figures/best-fit-pars-3.png"),
                 filename4 = paste0(homewd, "/model_fitting/GhV-fits/figures/best-fit-pars-4.png") )
