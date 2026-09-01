rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(data.table)

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

#set wd
setwd(paste0(homewd,"/model_fitting/GhV-fits/out-matrix-fit-all/"))
#load all the data and add to a list


function.load <- function(prefix, list, i){
  load(paste0(prefix,i,".Rdata")) 
  list[[i]] <- out.fit
  print(list)
  
  return(list)
}

dat.fit <- list()

for (i in 1:4000){
  skip_to_next <- FALSE
  tryCatch({dat.fit <- function.load(prefix = "mean", list=dat.fit, i=i)}, error = function(i) {skip_to_next <<- TRUE})
  if(skip_to_next) { next }    
}
dat.fit.df <- data.table:: rbindlist(map(dat.fit, 1))
age.sero.df <- data.table:: rbindlist(map(dat.fit, 2))
biwk.df <- data.table:: rbindlist(map(dat.fit, 3))

#for now, skip the data that came with (elements 2 and 3 of each big list)

write.csv(dat.fit.df, file=paste0(homewd,"/model_fitting/GhV-fits/all_fits_matrix.csv"), row.names = F)
write.csv(age.sero.df, file=paste0(homewd,"/model_fitting/GhV-fits/all_fits_matrix_age_sero_df.csv"), row.names = F)
write.csv(biwk.df, file=paste0(homewd,"/model_fitting/GhV-fits/all_fits_matrix_biwk_df.csv"), row.names = F)

rm(list=ls())


homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

setwd(homewd)

#1898
dat <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/all_fits_matrix.csv"), stringsAsFactors = F, header = T)
head(dat)
unique(dat$llik)

dat = subset(dat, convergence == 0 & !is.na(model))
# 1424 
# write.csv(dat, paste0(homewd,"/model_fitting/GhV-fits/all_converged_fits_matrix_round1.csv"), row.names = F)
# head(dat)

unique(dat$model)
length(unique(dat$model)) #20 - that is every model represented

#plot the number that completed per model - 62 possible
dat.sum <- ddply(dat,.(model), summarise, N=length(model)) #98-100 for each


ggplot(dat.sum) + geom_point(aes(model, y=N))

arrange(dat.sum, N) 
#MSIRNNS-mu-sick-mu-rec had only 7 completion, others had ~20-40 but most had >50
# all models had > 98 completions


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
#best fit: MSIRS, MSIRNS, MSIRNNS, MSIRNNI,MSIRNI

#mu-rec
pB<- ggplot(data=subset(dat, mu_rec>1 & mu_sick==1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik, color=mu_rec),width=.09, size=1, alpha=.4);pB
#MSIRS, MSIRNS,  MSIRNNS, MSIRNNI,  MSIRNI

#mu-sick
pC<- ggplot(data=subset(dat, mu_rec==1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_llik, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_llik, color=mu_sick),width=.09, size=1, alpha=.4);pC
#MSIRS,  MSIRNNI,  MSIRNI

#both mu
pD<- ggplot(data=subset(dat, mu_rec>1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_llik, group=model),  
               show.legend=F, fill=NA) +
  geom_jitter(aes(x=model, y=delta_llik, size=mu_sick, color=mu_rec),width=.09, alpha=.4) +
  scale_size_continuous(range = c(.01,1));pD
#MSIRS, MSIRNS,   MSIRNNI,  MSIRNI


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
#MSIRS, MSIRNS,  MSIRNNS, MSIRNNI,  MSIRNI

pB <- ggplot(data=subset(dat, mu_rec>1 & mu_sick==1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_AIC, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_AIC, color=mu_rec),width=.09, size=1, alpha=.4);pB
#MSIRS, MSIRNS,  MSIRNNS, MSIRNNI,  MSIRNI

pC<- ggplot(data=subset(dat, mu_rec==1 & mu_sick>1& antigen=="GhV") ) + theme_bw() +
  theme(axis.text.x = element_text(angle=90))+ facet_grid(~antigen) + coord_flip(ylim=c(0,25)) +
  geom_boxplot(aes(x=model, y=delta_AIC, group=model),  
               show.legend=F, fill=NA) +
  #geom_violin(aes(x=model, y=delta_AIC, group=model),  
  #           draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+
  geom_jitter(aes(x=model, y=delta_AIC, color=mu_sick),width=.09, size=1, alpha=.4);pC
#MSIRS, MSIRNS,  MSIRNNS, MSIRNNI,  MSIRNI

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

write.csv(top.fits, file =paste0(homewd,"/model_fitting/GhV-fits/top-fits-matrix-round1.csv"), row.names = F)


#then, also plot the top 20 parameter sets for each of the models and compare
rm(list=ls())


homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

setwd(homewd)

#1898
dat <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/all_fits_matrix.csv"), stringsAsFactors = F, header = T)

#1424
dat = subset(dat, convergence == 0 & !is.na(model))

unique(dat$model)
length(unique(dat$model)) #20 - that is every model represented

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


write.csv(dat20, file =paste0(homewd,"/model_fitting/GhV-fits/top-20-matrix-round1.csv"), row.names = F)

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

plot.top.10.fits <- function(df, filename1, filename2){
  df.list <- dlply(df, .(model))
  
  plot.list  <- lapply(df.list, plot.one)

  all.par1 <- cowplot::plot_grid(plotlist = plot.list[1:12], ncol = 4, nrow =3)
  all.par2 <- cowplot::plot_grid(plotlist = plot.list[13:20], ncol = 4, nrow =2)
  
  #and save
  ggsave(file = filename1,
              plot = all.par1,
              units="mm",  
              width=150, 
              height=150, 
              scale=3, 
              dpi=300)
  
  ggsave(file = filename2,
         plot = all.par1,
         units="mm",  
         width=150, 
         height=150, 
         scale=3, 
         dpi=300)
         
}

plot.top.10.fits(df = dat10, filename1 = paste0(homewd, "/model_fitting/GhV-fits/figures/best-fit-pars-1.png"), filename2 = paste0(homewd, "/model_fitting/GhV-fits/figures/best-fit-pars-2.png") )
