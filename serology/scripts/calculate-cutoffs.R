rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(mclust)
library(ggh4x)
library(ggcorrplot)
library(factoextra)
library(paran)
library(hopkins)
library(cluster)
library(corrr) 
library(FactoMineR)
library(ggpmisc)
library(reshape2)

#use this script to calculate cutoffs for the luminex


#first record your wd for this repo (home) 
#cara adds hers here too - just comment out whichever you are not using
#homewd = "/Users/carabrook/Developer/luminex-serology-madabat"
homewd = "/Users/gwenddolenkettenburg/Desktop/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

#now load in the merged luminex data from the "data" folder

#make sure this is latest catching file
dat <- read.csv(paste0(homewd, "/serology/data/lum_catch_merge.csv"), header = T, stringsAsFactors = F)
head(dat)


length(unique(dat$sampleid))
dat$sampleid[duplicated(dat$sampleid)]

test <- ddply(dat, .(bat_species, bat_sex), summarise, Ntot=length(unique(sampleid)))
test

#now for cutoffs, just slim to antigen data only and convert to long format
lumdat <- dplyr::select(dat, sampleid, bat_species, names(dat)[2:7])
head(lumdat)
lumdat.long <- melt(lumdat)
head(lumdat.long)
names(lumdat.long)[names(lumdat.long)=="variable"] <- "antigen"
names(lumdat.long)[names(lumdat.long)=="value"] <- "MFI"
unique(lumdat.long$antigen)
sort(unique(lumdat.long$MFI))
class(lumdat.long$MFI)
lumdat.long <- subset(lumdat.long, !is.na(MFI))


# # look at these really high mock binding individuals
# ggplot(subset(lumdat.long, MFI<0)) + 
#   geom_histogram(aes(MFI)) + facet_wrap(bat_species~antigen, scales = "free")
# 
# ddply(lumdat.long,.(bat_species, antigen), summarise, minMFI=min(MFI))
# unique(lumdat.long$sampleid[lumdat.long$MFI< -1000]) # there are just a few (11) Rousettus with really wild mock values
# arrange(subset(lumdat.long, MFI< -1000), sampleid)

# convert any negative values to 1
lumdat.long$MFI[lumdat.long$MFI<1] <- 1
lumdat.long$log10MFI <- log10(lumdat.long$MFI)

#add family
lumdat.long$virus_family[lumdat.long$antigen== "NiV"| lumdat.long$antigen=="HeV"| lumdat.long$antigen=="CedV"| lumdat.long$antigen=="GhV"| lumdat.long$antigen=="AngV"| lumdat.long$antigen=="MojV"] <- "Paramyxovirus"
unique(lumdat.long$virus_family)

#and genus information
lumdat.long$virus_genus <- "Henipavirus"
unique(lumdat.long$virus_genus)

#factor antigens by relatedness
lumdat.long$antigen <- factor(lumdat.long$antigen, levels=c( "NiV", "HeV","CedV", "GhV", "AngV", "MojV"))

lumdat.long$bat_species <- factor(lumdat.long$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))


# Now apply mixture models on these subsets of the data as selected in the "choose-antigens.R" script
bat.slim = subset(lumdat.long, bat_species=="Eidolon dupreanum" & antigen =="HeV" |
                               bat_species=="Eidolon dupreanum" & antigen =="NiV" |
                               bat_species=="Eidolon dupreanum" & antigen =="GhV" |
                               bat_species=="Eidolon dupreanum" & antigen =="AngV" |
                               bat_species=="Rousettus madagascariensis" & antigen =="GhV" )

# plot the distribution of MFI by each species - antigen pair
ggplot(data=bat.slim, aes(MFI)) + geom_histogram() + facet_wrap(bat_species~antigen) + scale_x_log10()
 
#now, fit the mixture models to the data and calculate cutoffs 
fit.mixtures <- function(df, return.dat,  mclust_thresh, mclust_thresh_lci, mclust_thresh_uci){
  name_bat= unique(df$bat_species)
  name_antigen =unique(df$antigen)
  
  print(name_bat)
  print(name_antigen)
  
  # in the case of RavV, we are finding that it selects a 
  # threshold too low (<500) with "E" assumption, so we force it to the
  # same number of clusters optimized for BdbV, which is 8
  if(name_bat=="Rousettus madagascariensis" & name_antigen=="RavV"){
    # it chooses 4 clusters with 1 empty cluster when allowed all
    # possibilities and 5 clusters with 1 empty cluster when allowed that too
    # so we will force it to 6 or more where it does not make any empty clusters
    out_clust <-  Mclust(df$log10MFI, 6:9, c("E")) 
    #when given this range, it chooses 8 which is the optimal cluster number for BdBV
 # }else if(name_bat=="Eidolon dupreanum" & name_antigen=="AngV"){
  #  out_clust <-  Mclust(df$log10MFI, 1:5, c("E"))
  
  }else{
    #otherwise, we allow the data to decide and it always chooses many
    #apply the cluster analysis, querying a possibility of  total clusters
    out_clust <-  Mclust(df$log10MFI, 1:9, c("E")) #changed this to 9  
  
  }
  
  #9 possible clusters, V = uneven variance, E = even variance
  # originally tested both distribution types but best support in all cases
  # and most conservative estimates for "E" so we will go with that.
  
  #can also plot this
  #plot(out_clust, what = "BIC", xlab = "log10_MFI") 
  #and can plot what distribution each data point is classified based on above
  #bottom row gives all data
  #higher rows give what distribution each data point is classed in
  #plot(out_clust, what = "classification", xlab = "log10_MFI") 
  #and plot the degree of uncertainty for which each data value is classified to a 
  #certain classification
  #plot(out_clust, what = "uncertainty", xlab = "log10_MFI")
  #colored by classification with height corresponding to uncertainty
  
  #select the number of distributions that give the best BIC
  n_dist = length(unique(out_clust$classification))
  
  # sometimes, mclust erroneously plots an empty cluster with the 
  # same mean as one with data. This is an example of overfitting. 
  # In these cases,make sure you have the right number of parameters
  if (n_dist==as.numeric(names(which(out_clust$BIC[,out_clust$modelName] == max(out_clust$BIC, na.rm = T))))){
    #now pull the parameters from the best fit
    # prop is the mixing proportion (lambda from Brook et al 2019)
    # mean is mu 
    # variance is sigma
    param.df <- cbind.data.frame(bat_species= unique(df$bat_species),
                                 antigen=unique(df$antigen),
                                 n_clust = n_dist,
                                 id_clust = seq(1, n_dist, 1),
                                 prop = out_clust$parameters$pro,
                                 mean = out_clust$parameters$mean,
                                 sd = sqrt(out_clust$parameters$variance$sigmasq))
    
  }else{ #here's where it overfit an empty cluster
    if(length(out_clust$parameters$variance$sigmasq)<1){
      param.df <- cbind.data.frame(bat_species= unique(df$bat_species),
                                   antigen=unique(df$antigen),
                                   n_clust = n_dist,
                                   id_clust = seq(1, n_dist, 1),
                                   prop = out_clust$parameters$pro[sort(unique(out_clust$classification))],
                                   mean = out_clust$parameters$mean[sort(unique(out_clust$classification))],
                                   sd = sqrt(out_clust$parameters$variance$sigmasq[sort(unique(out_clust$classification))&!is.na(out_clust$parameters$variance$sigmasq[sort(unique(out_clust$classification))])]))#this should be a single value, except in the "V" case above, and that does not overfit the distributions
      
    }else{
      param.df <- cbind.data.frame(bat_species= unique(df$bat_species),
                                   antigen=unique(df$antigen),
                                   n_clust = n_dist,
                                   id_clust = seq(1, n_dist, 1),
                                   prop = out_clust$parameters$pro[sort(unique(out_clust$classification))],
                                   mean = out_clust$parameters$mean[sort(unique(out_clust$classification))],
                                   sd = sqrt(out_clust$parameters$variance$sigmasq))#this should be a single value, except in the "V" case above, and that does not overfit the distributions
      
    }
    
    
  }
  
  
  param.df$optimal_mod = out_clust$modelName
  
  
  
  # For here, we select the cutoff certainty - we'll make this an 
  # input into the function so we can change it at various points 
  #log10MFI_threshold_low <- max(out_clust$data[out_clust$z[, 1] > mclust_thresh])
  # this gives the maximum value in the first distribution (negative distribution) 
  # that is negative with the certainty we specified. if we had more than 2 distributions, 
  # it currently is fixed at still giving the value for the first cluster
  
  
  if(ncol(out_clust$z)>1){
    
    # here we select the maximum value of the distribution that is 1 below the highest distribution (so if there are 3 clusters, this
    # is the max of cluster 2, or if there are 2 clusters, this is the max of cluster 1)
    # occassionally, the lower cluster gets split on either side of the upper, so we need to correct for that and make
    # sure the high values get assigned to the highest cluster 
    log10MFI_max_second_to_highest_dist <- max(out_clust$data[out_clust$classification == max(sort(unique(out_clust$classification))[-n_dist]) & out_clust$data < max(out_clust$data[out_clust$classification == max(sort(unique(out_clust$classification)))])])  
    
    
    
  }else{
    log10MFI_max_second_to_highest_dist = 0
    
  }
  
  log10MFI_threshold_hi <- min(out_clust$data[out_clust$z[, max(unique(out_clust$classification))] > mclust_thresh & out_clust$data > log10MFI_max_second_to_highest_dist]) 
  # if (log10MFI_threshold_hi==Inf){
  #   
  #   print(paste0("using alt mean threshold for ", name_antigen, " in ", name_bat))
  #   
  #   log10MFI_threshold_hi <- min(out_clust$data[out_clust$z[, max(unique(out_clust$classification))] > mclust_thresh_alt & out_clust$data > log10MFI_max_second_to_highest_dist]) 
  # }
  log10MFI_threshold_hi_lci <- min(out_clust$data[out_clust$z[, max(unique(out_clust$classification))] > mclust_thresh_lci & out_clust$data > log10MFI_max_second_to_highest_dist]) 
  # if (log10MFI_threshold_hi_lci==Inf){
  #   
  #   print(paste0("using alt lci threshold for ", name_antigen, " in ", name_bat))
  #   
  #   log10MFI_threshold_hi_lci <- min(out_clust$data[out_clust$z[, max(unique(out_clust$classification))] > mclust_thresh_lci_alt & out_clust$data > log10MFI_max_second_to_highest_dist]) 
  # }
  log10MFI_threshold_hi_uci <- min(out_clust$data[out_clust$z[, max(unique(out_clust$classification))] > mclust_thresh_uci & out_clust$data > log10MFI_max_second_to_highest_dist]) 
  # if (log10MFI_threshold_hi_uci==Inf){
  #   
  #   print(paste0("using alt uci for ", name_antigen, " in ", name_bat))
  # 
  #   log10MFI_threshold_hi_uci <- min(out_clust$data[out_clust$z[, max(unique(out_clust$classification))] > mclust_thresh_uci_alt & out_clust$data > log10MFI_max_second_to_highest_dist]) 
  # }
  # this gives the minimum value in the upper-most distribution (positive distribution) 
  # that is positive with the certainty we specified. if we estimate two distributions, this
  # is the second distribution only. if we estimate 3, this is the third
  # it also corrects to make sure it is ABOVE the max value for the distribution below it.
  
  

  #' We can plot this 
  plot(out_clust, what = "density", xlab = paste0("log10 MFI ", name_antigen, " for ", name_bat)) # Base plot in black
  lines(
    seq(1, 6, by = 0.01),               # over a sequence
    dnorm(seq(1, 6, by = 0.01),         # generate a normal distribution
          out_clust$parameters$mean[1], # with mean of the first cluster
          sqrt(out_clust$parameters$variance$sigmasq[1]))* # and SD of the first cluster
      out_clust$parameters$pro[1],      # Times the proportion of values in the cluster
    col = "blue", lty = 2)              # colored blue and dotted
  # if (name_bat=="Rousettus madagascariensis" & name_antigen=="RavV"){
  #   lines(
  #     seq(1, 6, by = 0.01), 
  #     dnorm(seq(1, 6, by = 0.01), 
  #           out_clust$parameters$mean[ncol(out_clust$z)], 
  #           sqrt(out_clust$parameters$variance$sigmasq[ncol(out_clust$z)]))*
  #       out_clust$parameters$pro[ncol(out_clust$z)],
  #     col = "red", lty = 2)
  # }else{
    lines(
    seq(1, 6, by = 0.01), 
    dnorm(seq(1, 6, by = 0.01), 
          out_clust$parameters$mean[max(unique(out_clust$classification))], 
          sqrt(out_clust$parameters$variance$sigmasq))*
      out_clust$parameters$pro[max(unique(out_clust$classification))],
    col = "red", lty = 2)
   #}
  lines(density(df$log10MFI), col = "purple") 
  # Add vertical dotted lines with threshold values.
  # These are different thresholds for what we consider the 
  # threshold value that defines a seropositive or seronegative case
  # abline(v = log10MFI_threshold_low , lty = 3, col="blue") # max value to be negative
  abline(v = log10MFI_threshold_hi, lty = 3, col="red") #min value to be positive. this will move right as your certainty value increases
  abline(v = log10MFI_threshold_hi_lci, lty = 3, col="cornflowerblue") 
  abline(v = log10MFI_threshold_hi_uci, lty = 3, col="seagreen") 
  
  
  #where does that put the cutoff?
  #store it in the parameter dataset
  param.df$cutoff_lMFI_clust <-   log10MFI_threshold_hi
  param.df$cutoff_MFI_clust <- 10^(log10MFI_threshold_hi)
  param.df$cutoff_lMFI_clust_lci <-   log10MFI_threshold_hi_lci
  param.df$cutoff_MFI_clust_lci <- 10^(log10MFI_threshold_hi_lci)
  param.df$cutoff_lMFI_clust_uci <-   log10MFI_threshold_hi_uci
  param.df$cutoff_MFI_clust_uci <- 10^(log10MFI_threshold_hi_uci)
  
  
  if(return.dat==TRUE){
    
    
    df$seropos <- 0
    df$seropos[df$log10MFI>=unique(param.df$cutoff_lMFI_clust)] <- 1
    df$seropos_lci <- 0
    df$seropos_lci[df$log10MFI>=unique(param.df$cutoff_lMFI_clust_lci)] <- 1
    df$seropos_uci <- 0
    df$seropos_uci[df$log10MFI>=unique(param.df$cutoff_lMFI_clust_uci)] <- 1
    
    return(list(param.df, df))
    
  }else{
    
    #return just the pars
    return(param.df)
  }
  
} 

#slim to pathogens of interest
#lumdat.long <- subset(lumdat.long, bat_species=="Eidolon dupreanum" & antigen =="HeV" | bat_species=="Eidolon dupreanum" & antigen =="GhV" | bat_species=="Eidolon dupreanum" & antigen =="MadV" |  bat_species=="Eidolon dupreanum" & antigen =="BDBV" | bat_species=="Rousettus madagascariensis" & antigen =="GhV" | bat_species=="Rousettus madagascariensis" & antigen =="NiV" | bat_species=="Rousettus madagascariensis" & antigen =="RAVV" |bat_species=="Rousettus madagascariensis" & antigen =="BDBV")

lum.list <- dlply(bat.slim, .(bat_species, antigen))

out <- lapply(lum.list, fit.mixtures, return.dat=T,  mclust_thresh=.90, mclust_thresh_lci = .85, mclust_thresh_uci=.95) #list of lists
mixture.par <- data.table::rbindlist(lapply(out, function(x) x[[1]])) #this is the parameter estimates from the distributions.
lumdat.clust <- data.table::rbindlist(lapply(out, function(x) x[[2]])) #this is the raw data classed as seropositive or seronegative based on the chosen threshold

write.csv(mixture.par, paste0(homewd, "/serology/data/mix_sub_parameters_75_90_96.csv"), row.names = F)
write.csv(lumdat.clust, paste0(homewd, "/serology/data/mix_sub_data_75_90_96.csv"), row.names = F)

#make a data frame in wide format too to merge with the viral shedding data
df_wide <- dcast(
  lumdat.clust,
  sampleid ~ antigen,
  value.var = c("MFI", "log10MFI","seropos","seropos_lci","seropos_uci","bat_species","virus_family","virus_genus")
)
write.csv(df_wide, paste0(homewd, "/serology/data/mix_sub_data_75_90_96_wide.csv"), row.names = F)


#also run for 

#now plot the modeled best fit distributions and the cutoff for each
make.sim.dat <- function(df1){
  mod.dat <- dnorm(x = seq(-1,5,  by= 0.01), mean = df1$mean, sd=df1$sd)*df1$prop
  return(mod.dat)
  
}
normal.sim <- function(df){
  row.list <- df %>% group_split(rownames(df))
  
  mod.out <- lapply(row.list, make.sim.dat)
  
    dat <- cbind.data.frame(log10MFI = rep(seq(-1,5,  by= 0.01), nrow(df)),
                            density =  c(unlist(mod.out)),
                            id_clust = c(rep(1:nrow(df), each=length(seq(-1,5,  by= 0.01)))))
    
 
  #and some details of the model 
  dat$bat_species = unique(df$bat_species)
  dat$antigen = unique(df$antigen)
  return(dat)
}
density.dat <- function(df){
  dens.df  <- cbind.data.frame(log10MFI=density(df$log10MFI)$x,
                               density = density(df$log10MFI)$y)
  
  #plot(density(df$log10MFI))
  #with(dens.df, plot(log10MFI, density))
  
  dens.df$bat_species = unique(df$bat_species)
  dens.df$antigen = unique(df$antigen)
  return(dens.df)
}
plot.dat.model.clust.lci.uci <- function(dat, mod.df, pEidfile, pPterfile, pRoufile, plot_save){
  #first, generate simulated data for each distribution using the 
  #estimated parameters
  
  #split the parameter estimates by species and antigen
  par.split <- dlply(mod.df, .(bat_species, antigen))
  
  model.dat <- data.table::rbindlist(lapply(par.split, normal.sim))
  #head(model.dat)
  model.dat$id_clust <- as.factor(model.dat$id_clust)
  
  #and get density within each histogram
  dat.split <- dlply(dat, .(bat_species, antigen))
  dens.dat <- data.table::rbindlist(lapply(dat.split, density.dat))
  
  #and plot all together across the three species
  
  pEid <- ggplot(data= subset(dens.dat, bat_species=="Eidolon dupreanum"))+
    geom_line(aes(x= log10MFI, y=density)) + theme_bw()+
    facet_wrap(antigen~., scales = "free")    +       # Breaking up plots by species
    theme(panel.grid = element_blank(),
          plot.title = element_text(face="italic")) +
    geom_line(data= subset(model.dat, bat_species=="Eidolon dupreanum"),
              aes(x = log10MFI,
                  y=density, color=id_clust)) +  # Using MFI on the X axis...
    #geom_vline(data=subset(mod.df, bat_species=="Eidolon dupreanum"),
    #          aes(xintercept = cutoff_lMFI), color= "red", linetype=2, size=.7) +
    geom_vline(data=subset(mod.df, bat_species=="Eidolon dupreanum"),
               aes(xintercept = cutoff_lMFI_clust), color= "red", linetype=2, size=.7) +
    geom_vline(data=subset(mod.df, bat_species=="Eidolon dupreanum"),
               aes(xintercept = cutoff_lMFI_clust_lci), color= "orange", linetype=2, size=.7) +
    geom_vline(data=subset(mod.df, bat_species=="Eidolon dupreanum"),
               aes(xintercept = cutoff_lMFI_clust_uci), color= "darkred", linetype=2, size=.7) +
    
    labs(title = "E. dupreanum") +
    xlab("log10 MFI")
  
  #we can be very confident that values above the dashed line are seropositive... 
  #but I think we can probably be more lenient as well...
  
  pPter <- ggplot(data= subset(dens.dat, bat_species=="Pteropus rufus"))+
    geom_line(aes(x= log10MFI, y=density)) + theme_bw()+
    theme(panel.grid = element_blank(),
          plot.title = element_text(face="italic")) +
    geom_line(data= subset(model.dat, bat_species=="Pteropus rufus"),
              aes(x = log10MFI,
                  y=density, color=id_clust)) +  # Using MFI on the X axis...
    geom_vline(data=subset(mod.df, bat_species=="Pteropus rufus"),
               aes(xintercept = cutoff_lMFI_clust), color= "red", linetype=2, size=.7) +
    geom_vline(data=subset(mod.df, bat_species=="Pteropus rufus"),
               aes(xintercept = cutoff_lMFI_clust_lci), color= "orange", linetype=2, size=.7) +
    geom_vline(data=subset(mod.df, bat_species=="Pteropus rufus"),
               aes(xintercept = cutoff_lMFI_clust_uci), color= "darkred", linetype=2, size=.7) +
    # geom_vline(data=subset(mod.df, bat_species=="Pteropus rufus"),
    #           aes(xintercept = cutoff_lMFI), color= "red", linetype=2, size=.7) +
    facet_wrap(antigen~., scales = "free_y")    +       # Breaking up plots by species
    labs(title = "P. rufus") +
    xlab("log10 MFI")
  
  
  pRou<- ggplot(data= subset(dens.dat, bat_species=="Rousettus madagascariensis"))+
    geom_line(aes(x= log10MFI, y=density)) + theme_bw()+
    theme(panel.grid = element_blank(),
          plot.title = element_text(face="italic")) +
    geom_line(data= subset(model.dat, bat_species=="Rousettus madagascariensis" & !is.na(density)),
              aes(x = log10MFI,
                  y=density, color=id_clust)) +  # Using MFI on the X axis...
    # geom_vline(data=subset(mod.df, bat_species=="Rousettus madagascariensis"),
    #           aes(xintercept = cutoff_lMFI), color= "red", linetype=2, size=.7) +
    geom_vline(data=subset(mod.df, bat_species=="Rousettus madagascariensis"),
               aes(xintercept = cutoff_lMFI_clust), color= "red", linetype=2, size=.7) +
    geom_vline(data=subset(mod.df, bat_species=="Rousettus madagascariensis"),
               aes(xintercept = cutoff_lMFI_clust_lci), color= "orange", linetype=2, size=.7) +
    geom_vline(data=subset(mod.df, bat_species=="Rousettus madagascariensis"),
               aes(xintercept = cutoff_lMFI_clust_uci), color= "darkred", linetype=2, size=.7) +
    facet_wrap(antigen~., scales = "free_y")    +       # Breaking up plots by species
    labs(title = "R. madagascariensis") +
    xlab("log10 MFI")
  
  #and save these in turn
  
  if(plot_save == TRUE){
    
    if(!is.na(pEidfile)){
      ggsave(file = pEidfile,
             plot=pEid,
             units="mm",  
             width=90, 
             height=40, 
             scale=3)#, 
      
    }
  
    if(!is.na(pPterfile) ){
    ggsave(file = pPterfile,
           plot=pPter,
           units="mm",  
           width=90, 
           height=40, 
           scale=3)#, 
    }
    
    if(!is.na(pRoufile) ){
    ggsave(file = pRoufile,
           plot=pRou,
           units="mm",  
           width=90, 
           height=40, 
           scale=3)#, 
    }
  }
}


#now plot the data and the fitted mixture models, including the 3 different cutoffs
plot.dat.model.clust.lci.uci(dat =lumdat.clust, mod.df=mixture.par, plot_save = T,
                         pEidfile= paste0(homewd, "/serology/figures/Eidolon_mixtures_85_90_95.png"),
                         pPterfile= NA,
                         pRoufile= paste0(homewd, "/serology/figures/Rousettus_mixtures_85_90_95.png"))

#and what is the prevalence at the 90 cutoff?
lumdat.sum <- ddply(lumdat.clust, .(bat_species, antigen), summarise, Ntot=length(seropos), Npos=sum(seropos), Npos_lci=sum(seropos_lci), Npos_uci=sum(seropos_uci))
lumdat.sum$seroprevalence <- lumdat.sum$Npos/lumdat.sum$Ntot
lumdat.sum$seroprevalence_lci <- lumdat.sum$Npos_lci/lumdat.sum$Ntot
lumdat.sum$seroprevalence_uci <- lumdat.sum$Npos_uci/lumdat.sum$Ntot


print(lumdat.sum)
write.csv(lumdat.sum, paste0(homewd, "/serology/data/lumdat.sum.sub.csv"), row.names = F)

### what proportion of bats are seropositivefor multiple related antigens?
head(lumdat.clust)
lumdat.clust$keep <- 1
lumdat.clust$keep[lumdat.clust$bat_species=="Rousettus madagascariensis" & lumdat.clust$antigen=="GhV"] <- 0
lum.run = subset(lumdat.clust, keep==1)
lum.split <- dlply(lum.run,.(sampleid))

sum.sero <- function(df, cutoff){
  if(cutoff=="mean"){
 if(sum(df$seropos)==1) {
    df.new <- cbind.data.frame(stat="single-seropos", antigen = unique(df$antigen[df$seropos==1]))
  } else if (sum(df$seropos)==0) {
    df.new <- cbind.data.frame(stat="seroneg", antigen = NA)
  } else if (sum(df$seropos)>1){
    df.new <- cbind.data.frame(stat="multi-seropos", antigen = do.call(paste, c(as.list(as.character(df$antigen[df$seropos==1])), sep = "-")))
  }
  df.new$sampleid = unique(df$sampleid)
  df.new$bat_species = unique(df$bat_species)
  df.new$cutoff = cutoff
  }else if(cutoff=="lci"){
    if(sum(df$seropos_lci==1)) {
      df.new <- cbind.data.frame(stat="single-seropos", antigen = unique(df$antigen[df$seropos_lci==1]))
    } else if (sum(df$seropos_lci)==0) {
      df.new <- cbind.data.frame(stat="seroneg", antigen = NA)
    } else if (sum(df$seropos_lci)>1){
      df.new <- cbind.data.frame(stat="multi-seropos", antigen = do.call(paste, c(as.list(as.character(df$antigen[df$seropos_lci==1])), sep = "-")))
    }
    df.new$sampleid = unique(df$sampleid)
    df.new$bat_species = unique(df$bat_species)
    df.new$cutoff = cutoff
  }else if(cutoff=="uci"){
    if(sum(df$seropos_uci)==1) {
      df.new <- cbind.data.frame(stat="single-seropos", antigen = unique(df$antigen[df$seropos_uci==1]))
    } else if (sum(df$seropos_uci)==0) {
      df.new <- cbind.data.frame(stat="seroneg", antigen = NA)
    } else if (sum(df$seropos_uci)>1){
      df.new <- cbind.data.frame(stat="multi-seropos", antigen = do.call(paste, c(as.list(as.character(df$antigen[df$seropos_uci==1])), sep = "-")))
    }
    df.new$sampleid = unique(df$sampleid)
    df.new$bat_species = unique(df$bat_species)
    df.new$cutoff = cutoff
  }
   return(df.new) 
}

lum.sum <- data.table::rbindlist(lapply(lum.split, sum.sero, cutoff="mean"))
head(lum.sum)

lum.sum2 <- data.table::rbindlist(lapply(lum.split, sum.sero, cutoff="lci"))
head(lum.sum2)

lum.sum3 <- data.table::rbindlist(lapply(lum.split, sum.sero, cutoff="uci"))
head(lum.sum3)

lum.sum <- rbind(lum.sum, lum.sum2, lum.sum3)

tot.bats <- ddply(lum.sum,.(bat_species, cutoff), summarise, N=length(sampleid))
out.all <- ddply(lum.sum,.(bat_species, stat, antigen, cutoff), summarise, Npos=length(sampleid))
out.all <- merge(out.all, tot.bats, by =c("bat_species", "cutoff"))

out.all$proportion <- out.all$Npos/out.all$N*100

#and save this as a potential table
write.csv(out.all, file = paste0(homewd,"/serology/data/prop-double-seropos.csv"), row.names = F)
# 
# #what about difference in seroprevalence by site?
# dat.site <- dplyr::select(dat, sampleid, roost_site)
# dat.site$roost_site <- sapply(strsplit(dat.site$roost_site, "_"), '[', 1)
# lumdat.clust <- merge(lumdat.clust, dat.site, by="sampleid", all.x = T)
# 
# lumdat.sum.site <- ddply(lumdat.clust, .(bat_species,  roost_site, antigen), summarise, Ntot=length(seropos), Npos=sum(seropos), Npos_lci=sum(seropos_lci), Npos_uci=sum(seropos_uci))
# lumdat.sum.site$seroprevalence <- lumdat.sum.site$Npos/lumdat.sum.site$Ntot
# lumdat.sum.site$seroprevalence_lci <- lumdat.sum.site$Npos_lci/lumdat.sum.site$Ntot
# lumdat.sum.site$seroprevalence_uci <- lumdat.sum.site$Npos_uci/lumdat.sum.site$Ntot
# 

