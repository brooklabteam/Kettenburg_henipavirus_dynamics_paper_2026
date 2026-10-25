rm(list=ls())

#this script uses recapture data and aged data to put some restraints on
#the duration of both adult and maternal Ab-mediated immunity.
#this precedes the fitting of our mechanistic models

library(plyr)
library(dplyr)   
library(ggplot2)
library(stats)
library(mgcv)
library(reshape2) 
library(ggh4x)
library(minpack.lm)
library(devtools)
library(rstan)
library(ggplot2)

#fit declining maternal immunity function to data


################################################################################
################################################################################
homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"

################################################################################
################################################################################


dat <- read.csv(paste0(homewd,"/model_fitting/GhV-fits/matrix-fit-1a/age-seroprevalence-data.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)[names(dat)=="species"] <- "bat_species"
dat$log10MFI <- log10(dat$MFI)

#focus just on GhV for now
dat.mat = subset(dat, antigen=="GhV" & age <=1 & bat_species=="Eidolon dupreanum")
sort(dat.mat$log10MFI)

#first, let's see what is susceptible:
dat.mat$suscept <- dat.mat$prev-1
dat.mat$suscept[dat.mat$suscept==-1] <- 1

#in cases where things where not suscept by relaxed lci but were suscept by strict uci
# only 3 - make 50/50
dat.mat$suscept[dat.mat$prev_uci==0 & dat.mat$prev_lci==1] <- sample(c(0,1), 3, replace = T)


fit <- glm(suscept ~ log10MFI, family = binomial, data = dat.mat)
summary(fit)

alpha <- coef(fit)[1]
beta  <- coef(fit)[2]

become_susc <- function(y, alpha, beta) {
  1 / (1 + exp(-(alpha + beta * y)))
}

y= log10(seq(1:28000))


plot(dat.mat$log10MFI, dat.mat$suscept, pch=16, col=rgb(0,0,0,0.3),
     xlab="log10MFI",
     ylab="Probability susceptible",
     main="Fitted susceptibility–MFI relationship\nfor waning maternal immunity")
lines(y, become_susc(y=y, alpha = alpha, beta = beta), col="blue", lwd=2)

#new function 
become_susc <- function(y, alpha=35.50088, beta=-11.57193) {
  1 / (1 + exp(-(alpha + beta * y)))
}


#now get your titer trajectory for maternal immunity
#want to remove those that might have been born seronegative -
#remove any with age <1 week that are seroneg
dat.mat[dat.mat$age<= (1/52) & dat.mat$prev==0,] # just 1
dat.mat = dat.mat[dat.mat$sampleid!="KEL402",]

#some cutoffs were established by mixture modeling
cutoff.df <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 3)), antigen=rep(c( "GhV"),3))
cutoff.df$MFI = c(1141,1763, 1619.25)#, 2037.75, 905.25, 1667.25)
cutoff.df$type = rep(c("mean", "uci", "lci"),each=1)
cutoff.df$log10MFI<- log10(cutoff.df$MFI)


fit <- nls(log10MFI ~ alpha * exp(-lambda * age) + eta,
           data = dat.mat,
           start = list(alpha = max(dat.mat$log10MFI, na.rm=TRUE),
                        lambda = 0.5,
                        eta = min(dat.mat$log10MFI, na.rm=TRUE)))


summary(fit)
sqrt(mean(resid(fit)^2)) #standarddeviation


#exp(-coef(fit)["lambda"] * 1) 


## plot time course maternal immunity

head(dat.mat)
pred.dat <- cbind.data.frame(age = seq(0,1,.01), pred = predict(fit, newdata = data.frame(age = seq(0,1,.01))))
ggplot(dat.mat) +geom_point(aes(x=age,y=log10MFI)) + 
                 geom_hline(aes(yintercept =cutoff.df$log10MFI[cutoff.df$type=="mean"]), linetype=2) +
                geom_line(data=pred.dat, aes(x=age, y=pred), color="blue") + theme_bw()

# blue line here is the nls fit

#now make the trajectory
coef(fit)
# Parameters for maternal waning
lambda <- 3.188032    # rate of decay
eta    <- 1.188038   # asymptote (residual baseline)
sd_val <- 0.2   # variability around mean trajectory - this is less than that estimated from the data

#this function is the maternal antibody titer decay
traj.mean.mat.exp <- function(x, y, lambda, eta, sd = 0.15, dt = 1/365) {
  
  
  # Compute per-time-step decay multiplier
  r <- exp(-lambda * (dt))
  
  # Expected mean at next time step for each x
  mu <- eta + r * (x - eta)
  
  
  
  # Build the transition matrix
  u <- outer(y, mu, function(y, mu) dnorm(y, mean = mu, sd = sd))
  
  # Normalize each column to integrate to 1
  dx <- diff(c(x[1] - (x[2]-x[1])/2,
               (x[-1] + x[-length(x)])/2,
               x[length(x)] + (x[length(x)]-x[length(x)-1])/2))
  u <- sweep(u, 2, colSums(u * dx), FUN = "/")
  
  return(u)
}


set.seed(123)

# Define titer grid
x <- seq(0, max(dat.mat$log10MFI), length.out = 100)
dx <- diff(x)[1]


# Build kernel matrix
K <- traj.mean.mat.exp(x, y=x, lambda = lambda, eta=eta, sd = .1, dt=1/365)

# Function to sample next titer from K
sample_next_titer <- function(current_idx, K) {
  probs <- K[, current_idx]
  probs <- probs / sum(probs)
  sample(seq_along(probs), 1, prob = probs)
}

# Simulate N individuals over T time steps
N <- 100   # number of simulated individuals
T <- 365   # number of time steps (days)
trajectories <- matrix(NA, nrow = T, ncol = N)

# Start all individuals near the max titer
start_idx <- which.min(abs(x - max(x)))
trajectories[1, ] <- x[start_idx]

for (t in 2:T) {
  for (i in 1:N) {
    current_idx <- which.min(abs(x - trajectories[t-1, i]))
    next_idx <- sample_next_titer(current_idx, K)
    trajectories[t, i] <- x[next_idx]
  }
}

# Compute deterministic mean trajectory (black line)
mu <- eta + exp(-lambda * ((0:(T-1))/365)) * (max(x) - eta)

#convert stochastic runs to dataset
rownames(trajectories) <- 1:nrow(trajectories)
colnames(trajectories) <- paste0("sim-", 1:ncol(trajectories))

# Convert matrix to data frame
df.stoch <- melt(trajectories)   # 
colnames(df.stoch) <- c("age", "bat", "log10MFI")
df.stoch$age <- df.stoch$age/365
df.mean <- cbind.data.frame(age=(0:(T-1))/365, log10MFI=mu)

# Plot with ggplot
p1 <- ggplot(df.stoch) +
  geom_line(aes(x = age, y = log10MFI, group=bat), color="grey") +
  geom_point(data = dat.mat, aes(x=age,y=log10MFI)) + 
  geom_hline(aes(yintercept =cutoff.df$log10MFI[cutoff.df$type=="mean"]), linetype=2) +
  geom_line(data=pred.dat, aes(x=age, y=pred), color="blue", size=1) + 
  theme_bw() + theme(panel.grid = element_blank()) +
  geom_line(data=df.mean, aes(x=age, y=log10MFI), color="red", size=1); p1
#blue line is the statistical model fit from nls
#red line is the deterministic mathematic model fit of the titer trajectory

#replot with just the deterministic model
p2 <- ggplot(df.stoch) + ggtitle("Simulated maternal antibody trajectories for 100 bats") +
  geom_line(aes(x = age, y = log10MFI, group=bat, alpha=bat), color="grey", show.legend = F) +
  geom_point(data = dat.mat, aes(x=age,y=log10MFI)) + 
  geom_hline(aes(yintercept =cutoff.df$log10MFI[cutoff.df$type=="mean"]), linetype=2) +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.title = element_text(size=18),
                     plot.title = element_text(size=18),
                     axis.text = element_text(size=14)) +
  ylab("log10 MFI") +
  geom_line(data=df.mean, aes(x=age, y=log10MFI), color="red", size=1); p2

# save this -- should probably be a component of at least a supplemental
# figure at some points
ggsave(file = paste0(homewd, "/final_figures_tables/tmp/within-host-matAb-fit-GhV.png"),
       plot = p2,
       units="mm",  
       width=65, 
       height=55, 
       scale=2.5, 
       dpi=300)
  
min(df.mean$age[df.mean$log10MFI<cutoff.df$log10MFI[cutoff.df$type=="mean"]]) 
#0.1753425 years - intersection of red and dashed - duration maternal immunity

min(df.mean$age[df.mean$log10MFI<cutoff.df$log10MFI[cutoff.df$type=="mean"]])*365 
#64 days - just over 2 months (by red)
min(df.mean$age[df.mean$log10MFI<cutoff.df$log10MFI[cutoff.df$type=="mean"]])*365/14
#4.57 biweeks
min(pred.dat$age[pred.dat$pred<cutoff.df$log10MFI[cutoff.df$type=="mean"]])*365/14 
#54.75 days - just under 2 months (by blue)
#3.9 biweeks

# blue line is nls fit to the data
# red line is the deterministic modeled 
# outcome from the fitted viral
# kinetics model - it shows that mean maternal 
# immunity wanes by 64 days via red, but 55 days via blue
# this is 2x 14 day timesteps as deterministic compartments
# each with exit probability 1, and then the 3rd compartment with 
# a probablistic exit
# 55 days - 28 days (M1 + M2) = 27 days (duration M3)
# 1/p3 = 27 days; p3 = 0.0373 days-1
# expressed in years for input into the model:
# 13.51852
# because timesteps are in biweeks (with 26 biweeks in a year), the probability of exit in this timestep is
# (1-exp(-(13.51852*(1/26)))) = 0.4054456
# for reporting in the paper, we fix the duration of maternal immunity at 28 days (two biweeks)
# then estimate the remaining time in the third maternally immune class.
# In total, this results in a total mean duration of maternal immunity of 
# roughly 55 days, or roughly 3.9 biweeks


################################################################################
################################################################################

#now  fit adult Ab trajectory to data
#this functions does it all and saves the data for plotting
run.adultAB.deterministic <- function(catch.dat, lum.dat, homewd, 
                        filename_fitted_bats,  filename_deterministic_df){
  
  #load the recaptures
  recap = subset(catch.dat, recap_new=="yes")
  # head(recap)
  # names(recap)
  #slim to columns of interest
  recap <- dplyr::select(recap, barcode_id, sampleid, roost_site, processing_date, bat_species, bat_age_class, bat_sex, mass_forearm_residual)
  # head(recap)
  
  
  #just take  eidolon
  recap=subset(recap, bat_species=="Eidolon dupreanum")
  recap$processing_date <- as.Date(recap$processing_date)
  # min(recap$processing_date)#"2013-11-24"
  # max(recap$processing_date)#"2026-01-24"
  
  
  
  #and now load the new luminex data 
  
  #slim
  lum.dat <- dplyr::select(lum.dat, sampleid,bat_species, antigen, MFI, seropos, seropos_lci, seropos_uci)
  
  #head(lum.dat)
  #length(unique(lum.dat$sampleid)) #1055
  #head(recap)
  
  lum.dat = subset(lum.dat, bat_species=="Eidolon dupreanum" & antigen=="GhV")
  
  
  recap <- merge(recap, lum.dat, by=c("bat_species", "sampleid"), all.x = T)
  #head(recap)
  # min(recap$processing_date)#"2013-11-24"
  # max(recap$processing_date)#"2026-01-22"
  # 
  
  recap = subset(recap, antigen=="GhV")
  
  # split by antigen and barcode and make sure there are multiple timepoints. 
  # color as bats that stay neg, seroconvert, serorevert or do both
  id.recap <- function(df){
    if(nrow(df)>1){
      df <- arrange(df, processing_date)
      #df.new <- dplyr::select(df, -(keep))
      df.new <- df
      df.new$class <- NA
      if(sum(as.numeric(as.character(df$seropos)))==0){
        df.new$class <- "always seronegative" # all seroneg
      }else if(sum(as.numeric(as.character(df$seropos)))==length(df$seropos)){
        df.new$class <- "always seropositive" #all seropos
      }else{
        vec.sero <- unique(diff(as.numeric(as.character(df$seropos))))
        if(sum(unique(vec.sero))==0){
          df.new$class <- "seroconversion-seroreversion" #both
        }else if((sum(unique(match(vec.sero,-1)),na.rm=T)==0) & (sum(unique(match(vec.sero,1)),na.rm=T)==1)){
          df.new$class <- "seroconversion" #both
        }else if((sum(unique(match(vec.sero,1)),na.rm=T)==0) & (sum(unique(match(vec.sero,-1)),na.rm=T)==1)){
          df.new$class <- "seroreversion" #both
        }
      }
      df.new$N_cap<- nrow(df.new)
      return(df.new)
    }
    
  }
  recap.list <- dlply(recap, .(antigen, barcode_id))
  
  recap.sum <- data.table::rbindlist(lapply(recap.list, id.recap))
  
  recap.sum$seropos <- as.factor(recap.sum$seropos)
  
  #and plot
  recap.sum$is_seroneg <- 1
  recap.sum$is_seroneg[recap.sum$class!="always seronegative"] <- 0
  recap.sum$is_seroneg <- as.factor(recap.sum$is_seroneg)
  
  #unique(recap.sum$class)
  recap.sum$seropos <- as.character(recap.sum$seropos)
  recap.sum$seropos[recap.sum$seropos==0] <- "seronegative"
  recap.sum$seropos[recap.sum$seropos==1] <- "seropositive"
  recap.sum$seropos <- factor(recap.sum$seropos, levels=c("seronegative", "seropositive"))
  
  
  #get those with < 1 year bewteen sampling dates
  recap.list <- dlply(recap.sum, .(barcode_id))
  
  count_spacing <- function(df, time_days){
    df <- arrange(df, processing_date)
    diff= c(diff(df$processing_date))
    
    df$days_diff <- c(NA, diff)
    
    df$time_series <- cumsum(c(1,diff))
    if(sum(df$days_diff[!is.na(df$days_diff)]<time_days)){
      return(df)
    }
    
  }
  
  recap.short <-  data.table::rbindlist(lapply(recap.list, count_spacing, time_days=6*31))
  
  #recap.short <- subset(recap.short, N_cap>2)
  
  
  recap.short$simple_class <- recap.short$class
  recap.short$simple_class <- as.character(recap.short$simple_class)
  recap.short$simple_class[recap.short$simple_class!="always seronegative"] <- "dynamic"
  recap.short = subset(recap.short, simple_class=="dynamic")
  #alphaz2 = c("always seronegative" = .3, "dynamic" = 1 )
  
  recap.short$era <- "2014-2016"
  recap.short$era[recap.short$processing_date>as.Date("2018-01-01")] <- "2018-2020"
  
  
  
  #now try restricting to only triptics that go "low-high-low". center on the peak
  #also need the first and second point to be within 6 months of each other
  recap.list <- dlply(recap.short, .(barcode_id))
  find.tryptic <- function(df){
    if(nrow(df)>=3){
      
      if(df$MFI[1]<df$MFI[2]){
        df.out = df[1:3,]
        
        
        df.out$new_ts <- NA
        df.out$new_ts[df.out$MFI==max(df.out$MFI)] <- 0
        
        df.out$new_ts[is.na(df.out$new_ts)] <- df.out$time_series[is.na(df.out$new_ts)] + -df.out$time_series[df.out$MFI==max(df.out$MFI)]
        
        return(df.out)
      }else if(df$MFI[2]<df$MFI[3] & nrow(df)>3){
        df.out = df[2:nrow(df),]
        
        #remove anything separated by >6 months and recompute diff
        df.out =subset(df.out, days_diff<(6*31))
        df.out$days_diff <- c(NA, diff(df.out$processing_date))
        
        if(sum(df.out$days_diff[!is.na(df.out$days_diff)]<(6*31))>1){ #otherwise you lose it
          
          
          
          df.out$new_ts <- NA
          df.out$new_ts[df.out$MFI==max(df.out$MFI)] <- 0
          
          df.out$new_ts[is.na(df.out$new_ts)] <- df.out$time_series[is.na(df.out$new_ts)] + -df.out$time_series[df.out$MFI==max(df.out$MFI)]
          
          return(df.out)
        }     }
      
      
    }
    
    
  }
  
  
  triptic_df <- data.table::rbindlist(lapply(recap.list, find.tryptic))
  
  triptic_df$seropos <- as.character(triptic_df$seropos)
  
  
  pre.MFI <- lm(MFI~new_ts, data = subset(triptic_df, barcode_id=="956000006105427" & new_ts<1))
  post.MFI <- lm(MFI~new_ts, data = subset(triptic_df, barcode_id=="956000006105427" & new_ts>-1))
  # 
  #MFI = 1141 is the cutoff for seropositivity
  add.pseudo <- triptic_df[triptic_df$barcode_id=="956000006105427",]
  add.pseudo <- add.pseudo[1:2,]
  add.pseudo$MFI <- c(predict(pre.MFI, newdata = data.frame(new_ts = -90)), predict(post.MFI, newdata = data.frame(new_ts = 100)))
  add.pseudo$new_ts <- c(-90, 100)
  add.pseudo$seropos <- "seronegative"
  
  triptic_df <- rbind(triptic_df, add.pseudo)
  
  
  #names(triptic_df)
  triptic.list <- dlply(triptic_df,.(barcode_id))
  
  prep.data.stan <- function(df){
    df$t_obs <- df$new_ts + abs(min(df$new_ts))+1
    df$last_neg <- max(df$t_obs[df$new_ts<0 & df$seropos=="seronegative"])
    df$first_pos <- min(df$t_obs[df$seropos=="seropositive"])
    
    return(df)
  }
  
  triptic_df <- data.table::rbindlist(lapply(triptic.list, prep.data.stan))
  
  #now run stan
  # Map barcode_id to consecutive integers for Stan
  bat_levels <- unique(triptic_df$barcode_id)
  triptic_df <- triptic_df %>%
    mutate(bat_idx = as.numeric(factor(barcode_id, levels = bat_levels)))
  
  # Prepare Stan data
  stan_data <- list(
    N = nrow(triptic_df),
    N_bats = length(bat_levels),
    bat_id = triptic_df$bat_idx,
    t_obs = triptic_df$t_obs,
    MFI = triptic_df$MFI,
    t_exp_min = triptic_df$last_neg,
    t_exp_max = triptic_df$first_pos
  )
  
  # Fit Stan model to tryptic data get trajectory of antibody dynamics
  # in exposed adults
  fit <- stan(file = paste0(homewd,"/model_fitting/GhV-fits/latent_exposure_model.stan"),
              data = stan_data,
              iter = 10000, chains = 4, cores = 4, warmup = 2000,
              control = list(adapt_delta = 0.99))
  #print(fit, pars = c("A", "B", "sigma", "t_exp"), probs = c(0.05, 0.5, 0.95))
  
  # here is the model that is fit: y(τ)=A⋅τ⋅exp(−B⋅τ)
  #A = amplitude factor (A=PeakValue⋅B⋅e)
  #B = rate modifier (time_to_peak = 1/B) but also plays a role in the decay on the other end...
  
  posterior <- as.data.frame(fit)
  #write.csv(posterior, file = paste0(homewd, "/model_fitting/GhV-fits/adult-Ab-stan-fit-posterior.csv"), row.names = F)
  #now get the deterministic trajectory
  
  #and the deterministic
  simulate_new_bat_deterministic_df <- function(Adist, Bdist, t_seq, t_exp = 0) {
    
    A <- mean(Adist)
    B <- mean(Bdist)
    t_exp <- 0   # or mean(posterior$t_exp[bat_id]) if you want average exposure
    
    
    tau <- pmax(t_seq - t_exp, 0)
    
    # simulate raw MFI
    MFI <- A * tau * exp(-B * tau) 
    
    #and the lci and uci
    ci_A <- quantile(Adist, probs = c(0.025, 0.975))
    ci_B <- quantile(Bdist, probs = c(0.025, 0.975))
    
    MFI_uci <- ci_A['97.5%'] * tau * exp(-ci_B['2.5%'] * tau)
    MFI_lci <- ci_A['2.5%'] * tau * exp(-ci_B['97.5%'] * tau)
    
    df.MFI <- cbind.data.frame(t=t_seq,
                               MFI = MFI,
                               MFI_lci = MFI_lci,
                               MFI_uci = MFI_uci)
    
    ggplot(df.MFI) + geom_line(aes(x=t_seq, y=MFI)) + geom_ribbon(aes(x=t_seq, ymin=MFI_lci, ymax=MFI_uci), alpha=.3) + theme_bw() +
    
    return(df.MFI)
  }
  
  det.MFI <- simulate_new_bat_deterministic_df(Adist=posterior$mu_A, Bdist=posterior$mu_B, t_seq=0:150)
  
  peak.det <- det.MFI$t[det.MFI$MFI==max(det.MFI$MFI)]
  
  #and the individual bats
  stan_df <- as.data.frame(stan_data)
  
  #align stan_df bats around the peak
  stan_df$t_exp_max[stan_df$bat_id==2] <- 91
  stan_df$new_ts <- stan_df$t_obs- stan_df$t_exp_max + peak.det
  
  #and save 
  write.csv(stan_df, file = filename_fitted_bats, row.names = F)
  write.csv(det.MFI, file = filename_deterministic_df, row.names = F)
}

run.adultAB.stochastic <- function(catch.dat, lum.dat, homewd, 
                        filename_fitted_bats, filename_stochastic_sims, filename_deterministic_sim){
  
  #load the recaptures
  recap = subset(catch.dat, recap_new=="yes")
  # head(recap)
  # names(recap)
  #slim to columns of interest
  recap <- dplyr::select(recap, barcode_id, sampleid, roost_site, processing_date, bat_species, bat_age_class, bat_sex, mass_forearm_residual)
  # head(recap)
  
  
  #just take  eidolon
  recap=subset(recap, bat_species=="Eidolon dupreanum")
  recap$processing_date <- as.Date(recap$processing_date)
  # min(recap$processing_date)#"2013-11-24"
  # max(recap$processing_date)#"2026-01-24"
  
  
  
  #and now load the new luminex data 
  
  #slim
  lum.dat <- dplyr::select(lum.dat, sampleid,bat_species, antigen, MFI, seropos, seropos_lci, seropos_uci)
  
  #head(lum.dat)
  #length(unique(lum.dat$sampleid)) #1055
  #head(recap)
  
  lum.dat = subset(lum.dat, bat_species=="Eidolon dupreanum" & antigen=="GhV")
  
  
  recap <- merge(recap, lum.dat, by=c("bat_species", "sampleid"), all.x = T)
  #head(recap)
  # min(recap$processing_date)#"2013-11-24"
  # max(recap$processing_date)#"2026-01-22"
  # 
  
  recap = subset(recap, antigen=="GhV")
  
  # split by antigen and barcode and make sure there are multiple timepoints. 
  # color as bats that stay neg, seroconvert, serorevert or do both
  id.recap <- function(df){
    if(nrow(df)>1){
      df <- arrange(df, processing_date)
      #df.new <- dplyr::select(df, -(keep))
      df.new <- df
      df.new$class <- NA
      if(sum(as.numeric(as.character(df$seropos)))==0){
        df.new$class <- "always seronegative" # all seroneg
      }else if(sum(as.numeric(as.character(df$seropos)))==length(df$seropos)){
        df.new$class <- "always seropositive" #all seropos
      }else{
        vec.sero <- unique(diff(as.numeric(as.character(df$seropos))))
        if(sum(unique(vec.sero))==0){
          df.new$class <- "seroconversion-seroreversion" #both
        }else if((sum(unique(match(vec.sero,-1)),na.rm=T)==0) & (sum(unique(match(vec.sero,1)),na.rm=T)==1)){
          df.new$class <- "seroconversion" #both
        }else if((sum(unique(match(vec.sero,1)),na.rm=T)==0) & (sum(unique(match(vec.sero,-1)),na.rm=T)==1)){
          df.new$class <- "seroreversion" #both
        }
      }
      df.new$N_cap<- nrow(df.new)
      return(df.new)
    }
    
  }
  recap.list <- dlply(recap, .(antigen, barcode_id))
  
  recap.sum <- data.table::rbindlist(lapply(recap.list, id.recap))
  
  recap.sum$seropos <- as.factor(recap.sum$seropos)
  
  #and plot
  recap.sum$is_seroneg <- 1
  recap.sum$is_seroneg[recap.sum$class!="always seronegative"] <- 0
  recap.sum$is_seroneg <- as.factor(recap.sum$is_seroneg)
  
  #unique(recap.sum$class)
  recap.sum$seropos <- as.character(recap.sum$seropos)
  recap.sum$seropos[recap.sum$seropos==0] <- "seronegative"
  recap.sum$seropos[recap.sum$seropos==1] <- "seropositive"
  recap.sum$seropos <- factor(recap.sum$seropos, levels=c("seronegative", "seropositive"))
  
  
  #get those with < 1 year bewteen sampling dates
  recap.list <- dlply(recap.sum, .(barcode_id))
  
  count_spacing <- function(df, time_days){
    df <- arrange(df, processing_date)
    diff= c(diff(df$processing_date))
    
    df$days_diff <- c(NA, diff)
    
    df$time_series <- cumsum(c(1,diff))
    if(sum(df$days_diff[!is.na(df$days_diff)]<time_days)){
      return(df)
    }
    
  }
  
  recap.short <-  data.table::rbindlist(lapply(recap.list, count_spacing, time_days=6*31))
  
  #recap.short <- subset(recap.short, N_cap>2)
  
  
  recap.short$simple_class <- recap.short$class
  recap.short$simple_class <- as.character(recap.short$simple_class)
  recap.short$simple_class[recap.short$simple_class!="always seronegative"] <- "dynamic"
  recap.short = subset(recap.short, simple_class=="dynamic")
  #alphaz2 = c("always seronegative" = .3, "dynamic" = 1 )
  
  recap.short$era <- "2014-2016"
  recap.short$era[recap.short$processing_date>as.Date("2018-01-01")] <- "2018-2020"
  
  
  
  #now try restricting to only triptics that go "low-high-low". center on the peak
  #also need the first and second point to be within 6 months of each other
  recap.list <- dlply(recap.short, .(barcode_id))
  find.tryptic <- function(df){
    if(nrow(df)>=3){
      
      if(df$MFI[1]<df$MFI[2]){
        df.out = df[1:3,]
        
        
        df.out$new_ts <- NA
        df.out$new_ts[df.out$MFI==max(df.out$MFI)] <- 0
        
        df.out$new_ts[is.na(df.out$new_ts)] <- df.out$time_series[is.na(df.out$new_ts)] + -df.out$time_series[df.out$MFI==max(df.out$MFI)]
        
        return(df.out)
      }else if(df$MFI[2]<df$MFI[3] & nrow(df)>3){
        df.out = df[2:nrow(df),]
        
        #remove anything separated by >6 months and recompute diff
        df.out =subset(df.out, days_diff<(6*31))
        df.out$days_diff <- c(NA, diff(df.out$processing_date))
        
        if(sum(df.out$days_diff[!is.na(df.out$days_diff)]<(6*31))>1){ #otherwise you lose it
          
          
          
          df.out$new_ts <- NA
          df.out$new_ts[df.out$MFI==max(df.out$MFI)] <- 0
          
          df.out$new_ts[is.na(df.out$new_ts)] <- df.out$time_series[is.na(df.out$new_ts)] + -df.out$time_series[df.out$MFI==max(df.out$MFI)]
          
          return(df.out)
        }     }
      
      
    }
    
    
  }
  
  
  triptic_df <- data.table::rbindlist(lapply(recap.list, find.tryptic))
  
  triptic_df$seropos <- as.character(triptic_df$seropos)
  
  
  pre.MFI <- lm(MFI~new_ts, data = subset(triptic_df, barcode_id=="956000006105427" & new_ts<1))
  post.MFI <- lm(MFI~new_ts, data = subset(triptic_df, barcode_id=="956000006105427" & new_ts>-1))
  # 
  #MFI = 1141 is the cutoff for seropositivity
  add.pseudo <- triptic_df[triptic_df$barcode_id=="956000006105427",]
  add.pseudo <- add.pseudo[1:2,]
  add.pseudo$MFI <- c(predict(pre.MFI, newdata = data.frame(new_ts = -90)), predict(post.MFI, newdata = data.frame(new_ts = 100)))
  add.pseudo$new_ts <- c(-90, 100)
  add.pseudo$seropos <- "seronegative"
  
  triptic_df <- rbind(triptic_df, add.pseudo)
  
  
  #names(triptic_df)
  triptic.list <- dlply(triptic_df,.(barcode_id))
  
  prep.data.stan <- function(df){
    df$t_obs <- df$new_ts + abs(min(df$new_ts))+1
    df$last_neg <- max(df$t_obs[df$new_ts<0 & df$seropos=="seronegative"])
    df$first_pos <- min(df$t_obs[df$seropos=="seropositive"])
    
    return(df)
  }
  
  triptic_df <- data.table::rbindlist(lapply(triptic.list, prep.data.stan))
  
  #now run stan
  # Map barcode_id to consecutive integers for Stan
  bat_levels <- unique(triptic_df$barcode_id)
  triptic_df <- triptic_df %>%
    mutate(bat_idx = as.numeric(factor(barcode_id, levels = bat_levels)))
  
  # Prepare Stan data
  stan_data <- list(
    N = nrow(triptic_df),
    N_bats = length(bat_levels),
    bat_id = triptic_df$bat_idx,
    t_obs = triptic_df$t_obs,
    MFI = triptic_df$MFI,
    t_exp_min = triptic_df$last_neg,
    t_exp_max = triptic_df$first_pos
  )
  
  # Fit Stan model to tryptic data get trajectory of antibody dynamics
  # in exposed adults
  fit <- stan(file = paste0(homewd,"/model_fitting/GhV-fits/latent_exposure_model.stan"),
              data = stan_data,
              iter = 10000, chains = 4, cores = 4, warmup = 2000,
              control = list(adapt_delta = 0.99))
  #print(fit, pars = c("A", "B", "sigma", "t_exp"), probs = c(0.05, 0.5, 0.95))
  
  posterior <- as.data.frame(fit)
  
  #get the hyperparameters (population-level)
  hyperpars <- posterior %>%
    select(mu_A, sigma_A, mu_B, sigma_B, sigma)
  
  # quick summary
  #summary(hyperpars)
  
  #get the bat-specific parameters 
  #(A and B for each of the 4 bats included in the fitting)
  bat_params <- posterior %>%
    select(matches("^A\\[|^B\\["))  # matches A[1], A[2], ... B[1], B[2], ...
  
  # latent exposure times
  t_exp_post <- posterior %>%
    select(matches("^t_exp\\["))     # t_exp[1], t_exp[2], ...
  
  bat_summary <- bat_params %>%
    summarise(across(everything(), list(
      mean = ~mean(.),
      median = ~median(.),
      q05 = ~quantile(., 0.05),
      q95 = ~quantile(., 0.95)
    )))
  
  #head(bat_summary)
  
  # to plot the individual bat fits:
  
  sim_fun <- function(A, B, t_exp_val, sigma, t_seq) {
    tau <- pmax(t_seq - t_exp_val, 0)
    A * tau * exp(-B * tau) + rnorm(length(tau), 0, sigma)
  }
  
  # simulate trajectory for bat 1
  sim_traj <- sim_fun(A=mean(posterior$`A[1]`), B=mean(posterior$`B[1]`), t_exp_val=mean(posterior$`t_exp[1]`), sigma=mean(posterior$sigma), t_seq=0:150)
  
  #plot(0:150, sim_traj) #this is the Ab trajectory for bat 1
  t_seq = 0:150
  
  sim_list <- lapply(1:4, function(i) {
    
    A_hat <- mean(posterior[[paste0("A[", i, "]")]])
    B_hat <- mean(posterior[[paste0("B[", i, "]")]])
    sigma_hat <- mean(posterior[["sigma"]])
    t_exp_val <- mean(posterior[[paste0("t_exp[", i, "]")]])
    data.frame(
      t = t_seq,
      MFI_sim = sim_fun(A_hat, B_hat, t_exp_val, sigma_hat, t_seq),
      bat = i
    )
  })
  
  sim_df <- do.call(rbind, sim_list)
  
  stan_df <- as.data.frame(stan_data)
  
  
  #now simulate many new bats from the population-level hyperparameters
  hyperpars <- posterior %>%
    select(mu_A, sigma_A, mu_B, sigma_B, sigma)
  
  
  
  simulate_new_bat_stochastic <- function(mu_A, sigma_A, mu_B, sigma_B, sigma, t_seq, t_exp = 0, MFI_min = 0, MFI_max = 50000, shrink_factor) {
    
    # draw bat-level parameters from population
    #shrink_factor = 0 is fully shrunk to mean
    #shrink_factor = 1 is full sigma variance
    A <- rnorm(1, mean = mu_A, sd = sigma_A*shrink_factor)
    B <- rnorm(1, mean = mu_B, sd = sigma_B*shrink_factor)
    
    tau <- pmax(t_seq - t_exp, 0)
    
    # simulate raw MFI
    MFI <- A * tau * exp(-B * tau) + rnorm(length(tau), 0, sigma)
    
    # constrain MFI to [MFI_min, MFI_max]
    MFI <- pmin(pmax(MFI, MFI_min), MFI_max)
    
    return(MFI)
  }
  
  
  n_new_bats <- 300
  t_seq <- seq(0, 150, by = 1)   # time in days
  
  # sample 300 posterior rows for variation
  posterior_rows <- hyperpars[sample(nrow(hyperpars), n_new_bats, replace = TRUE), ]
  
  sim_list <- lapply(1:n_new_bats, function(i) {
    row <- posterior_rows[i, ]
    data.frame(
      t = t_seq,
      MFI_sim = simulate_new_bat_stochastic(row$mu_A, row$sigma_A, row$mu_B, row$sigma_B, row$sigma, t_seq, shrink_factor = 0),
      bat = i
    )
  })
  
  
  sim_df <- do.call(rbind, sim_list)
  
  #and the deterministic
  simulate_new_bat_deterministic <- function(A, B, t_seq, t_exp = 0) {
    
    A <- mean(posterior$mu_A)
    B <- mean(posterior$mu_B)
    t_exp <- 0   # or mean(posterior$t_exp[bat_id]) if you want average exposure
    
    
    tau <- pmax(t_seq - t_exp, 0)
    
    # simulate raw MFI
    MFI <- A * tau * exp(-B * tau) 
    
    
    
    return(MFI)
  }
  
  det.bat <- cbind.data.frame(t=t_seq,
                              MFI=simulate_new_bat_deterministic(A=row$mu_A, B=row$mu_B, t_seq))
  
  
  # #plot sims
  # p3 <- ggplot() +
  #   geom_line(data= sim_df, aes(x = t, y = MFI_sim, group = bat),alpha = 0.2, color = "steelblue") +
  #   labs(title = "Simulated antibody trajectories for 300 new bats",
  #        x = "Days since exposure",
  #        y = "MFI") +
  #   geom_line(data=det.bat, aes(x=t, y=MFI),size=1) +
  #   theme(panel.grid = element_blank()) +
  #   theme_bw() + coord_cartesian(ylim=c(0,40000)) + 
  #   geom_hline(aes(yintercept =cutoff.df$MFI[cutoff.df$bat_species=="Eidolon dupreanum"]), linetype=2); p3
  # 
  #get the peak of the deterministic model
  peak.det = det.bat$t[det.bat$MFI==max(det.bat$MFI)]
  
  # sub_df <- subset(sim_df, bat<=100 )
  # sub_df <- subset(sim_df, bat>100 & bat<=200)
  sub_df <- subset(sim_df, bat>200 & bat<=300)
  # 
  
  
  #align stan_df bats around the peak
  stan_df$t_exp_max[stan_df$bat_id==2] <- 91
  stan_df$new_ts <- stan_df$t_obs- stan_df$t_exp_max + peak.det
  
  #and save 
  write.csv(stan_df, file = filename_fitted_bats, row.names = F)
  write.csv(sub_df, file = filename_stochastic_sims, row.names = F)
  write.csv(det.bat, file = filename_deterministic_sim, row.names = F)
}


run.adultAB.deterministic(catch.dat <- read.csv(file = paste0(homewd, "/master_data_files/catching_data.csv"), header = T, stringsAsFactors = F),
                       lum.dat <- read.csv(file = paste0(homewd, "/serology/data/mix_sub_data_75_90_96_long.csv")),
                       homewd = homewd,
                       filename_fitted_bats = paste0(homewd, "/serology/data/stan_fitted_bats.csv"), 
                       filename_deterministic_df= paste0(homewd, "/serology/data/stan_fit_deterministic_df_adultAB.csv"))


run.adultAB.stochastic(catch.dat <- catch.dat <- read.csv(file = paste0(homewd, "/master_data_files/catching_data.csv"), header = T, stringsAsFactors = F),
            lum.dat <- read.csv(file = paste0(homewd, "/serology/data/mix_sub_data_75_90_96_long.csv")),
            homewd = homewd,
            filename_fitted_bats = paste0(homewd, "/serology/data/stan_fitted_bats.csv"), 
            filename_stochastic_sims= paste0(homewd, "/serology/data/stan_fit_stochastic_sim_adultAB.csv"), 
            filename_deterministic_sim= paste0(homewd, "/serology/data/stan_fit_deterministic_sim_adultAB.csv"))









#below is all exploratory and produces some temporary plots, but above is all that is needed

# first, collect all the adults that have a 3-prong peak


#load the recaptures
catch.dat <- read.csv(paste0(batwd,"/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)
head(catch.dat)
recap = subset(catch.dat, recap_new=="yes")
head(recap)
names(recap)
#slim to columns of interest
recap <- dplyr::select(recap, barcode_id, sampleid, roost_site, processing_date, bat_species, bat_age_class, bat_sex, mass_forearm_residual)
head(recap)


#just take rousettus and eidolon
recap=subset(recap, bat_species!="Pteropus rufus")
recap$processing_date <- as.Date(recap$processing_date)
min(recap$processing_date)#"2013-11-24"
max(recap$processing_date)#"2026-01-24"



#and now load the new luminex data 
lum.dat <- read.csv(file = paste0(homewd, "/serology/data/mix_sub_data_75_90_96_long.csv"))
head(lum.dat)

#slim
lum.dat <- dplyr::select(lum.dat, sampleid,bat_species, antigen, MFI, seropos, seropos_lci, seropos_uci)

head(lum.dat)
length(unique(lum.dat$sampleid)) #1055
head(recap)

lum.dat = subset(lum.dat, bat_species=="Eidolon dupreanum" & antigen=="GhV")
recap = subset(recap, bat_species=="Eidolon dupreanum")

recap <- merge(recap, lum.dat, by=c("bat_species", "sampleid"), all.x = T)
head(recap)
min(recap$processing_date)#"2013-11-24"
max(recap$processing_date)#"2026-01-22"


recap = subset(recap, antigen=="GhV")

# split by antigen and barcode and make sure there are multiple timepoints. 
# color as bats that stay neg, seroconvert, serorevert or do both
id.recap <- function(df){
  if(nrow(df)>1){
    df <- arrange(df, processing_date)
    #df.new <- dplyr::select(df, -(keep))
    df.new <- df
    df.new$class <- NA
    if(sum(as.numeric(as.character(df$seropos)))==0){
      df.new$class <- "always seronegative" # all seroneg
    }else if(sum(as.numeric(as.character(df$seropos)))==length(df$seropos)){
      df.new$class <- "always seropositive" #all seropos
    }else{
      vec.sero <- unique(diff(as.numeric(as.character(df$seropos))))
      if(sum(unique(vec.sero))==0){
        df.new$class <- "seroconversion-seroreversion" #both
      }else if((sum(unique(match(vec.sero,-1)),na.rm=T)==0) & (sum(unique(match(vec.sero,1)),na.rm=T)==1)){
        df.new$class <- "seroconversion" #both
      }else if((sum(unique(match(vec.sero,1)),na.rm=T)==0) & (sum(unique(match(vec.sero,-1)),na.rm=T)==1)){
        df.new$class <- "seroreversion" #both
      }
    }
    df.new$N_cap<- nrow(df.new)
    return(df.new)
  }
  
}
recap.list <- dlply(recap, .(antigen, barcode_id))

recap.sum <- data.table::rbindlist(lapply(recap.list, id.recap))

recap.sum$seropos <- as.factor(recap.sum$seropos)

#and plot
recap.sum$is_seroneg <- 1
recap.sum$is_seroneg[recap.sum$class!="always seronegative"] <- 0
recap.sum$is_seroneg <- as.factor(recap.sum$is_seroneg)

unique(recap.sum$class)
recap.sum$seropos <- as.character(recap.sum$seropos)
recap.sum$seropos[recap.sum$seropos==0] <- "seronegative"
recap.sum$seropos[recap.sum$seropos==1] <- "seropositive"
recap.sum$seropos <- factor(recap.sum$seropos, levels=c("seronegative", "seropositive"))

filz = c('seronegative' = "navy",'seropositive'= "violetred")

colorz1 = c('always seronegative' = "navy", 
            "always seropositive" = "violetred4",
            "seroconversion-seroreversion"="goldenrod2", 
            "seroreversion" = "lightseagreen",
            "seroconversion" = "orangered")

shapez = c('female' = 21, 'male' = 24)
alphaz <- c('0' = 1, '1'=.1)


recap.sum$class <- factor(recap.sum$class, levels = c("always seronegative", "always seropositive", "seroconversion","seroreversion", "seroconversion-seroreversion" ))

Eid.recap <- ggplot(subset(recap.sum, bat_species=="Eidolon dupreanum")) + theme_bw()+ 
  geom_line(aes(x=processing_date, y=barcode_id, group=barcode_id, color=class, alpha=is_seroneg), size=1) + 
  geom_point(aes(x=processing_date, y=barcode_id, group=barcode_id, shape=bat_sex, fill=seropos, alpha=is_seroneg), size=3, color="gray65") + 
  scale_shape_manual(name="bat sex", values = shapez) + 
  scale_fill_manual(name = "serostatus", values = filz) +
  scale_color_manual(name= "serological\ntrajectory", values = colorz1) + ylab("individual barcode id\n\n") +
  facet_nested(bat_sex~antigen, scales = "free", space = "free") + scale_alpha_manual(values = alphaz) +
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        legend.box = "vertical",
        legend.position = "left",
        panel.spacing.y = unit(.1, "lines"),
        panel.spacing.x = unit(1, "lines"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        plot.margin = unit(c(.3,.3,0,2),"lines"),
        axis.text.x = element_text(size=14),  strip.text.x = element_text(size=16),
        strip.text.y = element_text(size=16),#, face="italic"),
        panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
  scale_x_continuous(breaks = c(as.Date("2014-07-01"),as.Date("2015-07-01"), as.Date("2016-07-01"), as.Date("2017-07-01"), as.Date("2018-07-01"), as.Date("2019-07-01")),labels = c(2014,2015,2016,2017,2018,2019)) +
  guides(alpha="none",  shape=guide_legend(order=1, override.aes = list(size=3, color="black")), fill=guide_legend(order=2, override.aes = list(shape=21,size=3)), color=guide_legend(order=3, ncol=1));Eid.recap




head(recap.sum)

#get those with < 1 year bewteen sampling dates
recap.list <- dlply(recap.sum, .(barcode_id))

count_spacing <- function(df, time_days){
  df <- arrange(df, processing_date)
  diff= c(diff(df$processing_date))
  
  df$days_diff <- c(NA, diff)
  
  df$time_series <- cumsum(c(1,diff))
  if(sum(df$days_diff[!is.na(df$days_diff)]<time_days)){
    return(df)
  }
  
}

recap.short <-  data.table::rbindlist(lapply(recap.list, count_spacing, time_days=6*31))

#recap.short <- subset(recap.short, N_cap>2)


recap.short$simple_class <- recap.short$class
recap.short$simple_class <- as.character(recap.short$simple_class)
recap.short$simple_class[recap.short$simple_class!="always seronegative"] <- "dynamic"
recap.short = subset(recap.short, simple_class=="dynamic")
alphaz2 = c("always seronegative" = .3, "dynamic" = 1 )

recap.short$era <- "2014-2016"
recap.short$era[recap.short$processing_date>as.Date("2018-01-01")] <- "2018-2020"

# add <- recap.short[1,]
# add$era <- "2014-2016"
# add$processing_date <- as.Date("2017-01-01")
# add$MFI <- NA
# add$antigen <- "NiV"
# add2 <- add
# add2$antigen <- "HeV"
# add <- rbind(add, add2)
# add2 <- add
# add2$era <- "2018-2020"
# add2$processing_date <- as.Date("2017-01-02")
# 
# recap.short <- rbind(recap.short, add, add2)

shapez = c("female" = 21, "male" =24)

#and include plot of the cutoff line by species/antigen/era
cutoff.df <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 6), rep("Rousettus madagascariensis",3)), antigen=c("NiV", "HeV", "NiV", "HeV", "GhV", "AngV", "GhV", "BdbV", "RavV"), era=c(rep("2014-2016",2), rep("2018-2020",7)))
cutoff.df$MFI = c(757,215,936,2851,1141,1763,2642,5647,6617)
cutoff.df$antigen <- factor(cutoff.df$antigen, levels=c("NiV", "HeV", "GhV", "AngV", "BdbV", "RavV"))
cutoff.df = subset(cutoff.df, antigen=="GhV")

EidMFIall <- ggplot(data=subset(recap.short, bat_species=="Eidolon dupreanum")) + 
  facet_nested(bat_sex~., scales = "free", space = "free") + theme_bw() +
  scale_shape_manual(name="bat sex", values = shapez) + 
  scale_alpha_manual(values = alphaz2) + ylab("mean fluorescence intensity (MFI)")+
  theme(axis.title.y = element_text(size=16), axis.title.x = element_blank(),
        legend.box = "horizontal", strip.text.y = element_text(size=16),#, face = "italic"),
        strip.background.y = element_rect(fill="white"),
        strip.background.x = element_blank(),
        axis.text = element_text(size=14), 
        panel.spacing.y = unit(.1, "lines"),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,.3,.3,.3),"lines"),
        panel.spacing.x = unit(c(0,1,0,1,1),"lines"),
        panel.grid = element_blank()) +
  geom_hline(data = subset(cutoff.df, bat_species=="Eidolon dupreanum"), aes(yintercept = MFI), color="red", linetype=2) + 
  geom_line(aes(x=processing_date, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1)+
  geom_point(aes(x=processing_date, y=MFI, fill=barcode_id, shape = bat_sex, alpha = simple_class), show.legend = F, size=3) +
  scale_x_continuous(breaks = c(as.Date("2014-07-01"),as.Date("2015-07-01"), as.Date("2016-07-01"), as.Date("2017-07-01"), as.Date("2018-07-01"), as.Date("2019-07-01")),labels = c(2014,2015,2016,2017,2018,2019)) +
  guides(alpha="none",  shape=guide_legend(order=1));EidMFIall


#now try restricting to only triptics that go "low-high-low". center on the peak
#also need the first and second point to be within 6 months of each other
recap.list <- dlply(recap.short, .(barcode_id))
find.tryptic <- function(df){
  if(nrow(df)>=3){
    
    if(df$MFI[1]<df$MFI[2]){
      df.out = df[1:3,]
      
      
      df.out$new_ts <- NA
      df.out$new_ts[df.out$MFI==max(df.out$MFI)] <- 0
      
      df.out$new_ts[is.na(df.out$new_ts)] <- df.out$time_series[is.na(df.out$new_ts)] + -df.out$time_series[df.out$MFI==max(df.out$MFI)]
      
      return(df.out)
    }else if(df$MFI[2]<df$MFI[3] & nrow(df)>3){
      df.out = df[2:nrow(df),]
      
      #remove anything separated by >6 months and recompute diff
      df.out =subset(df.out, days_diff<(6*31))
      df.out$days_diff <- c(NA, diff(df.out$processing_date))
      
      if(sum(df.out$days_diff[!is.na(df.out$days_diff)]<(6*31))>1){ #otherwise you lose it
        
 
      
      df.out$new_ts <- NA
      df.out$new_ts[df.out$MFI==max(df.out$MFI)] <- 0
      
      df.out$new_ts[is.na(df.out$new_ts)] <- df.out$time_series[is.na(df.out$new_ts)] + -df.out$time_series[df.out$MFI==max(df.out$MFI)]
      
      return(df.out)
    }     }
    
    
  }
  
  
}


triptic_df <- data.table::rbindlist(lapply(recap.list, find.tryptic))

triptic_df$seropos <- as.character(triptic_df$seropos)

tryptic.plot <- ggplot(data=triptic_df) + 
  facet_nested(bat_sex~., scales = "free", space = "free") + theme_bw() +
  scale_shape_manual(name="bat sex", values = shapez) + 
  scale_alpha_manual(values = alphaz2) + ylab("mean fluorescence intensity (MFI)")+
  theme(axis.title.y = element_text(size=16), axis.title.x = element_blank(),
        legend.box = "horizontal", strip.text.y = element_text(size=16),#, face = "italic"),
        strip.background.y = element_rect(fill="white"),
        strip.background.x = element_blank(),
        axis.text = element_text(size=14), 
        panel.spacing.y = unit(.1, "lines"),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,.3,.3,.3),"lines"),
        panel.spacing.x = unit(c(0,1,0,1,1),"lines"),
        panel.grid = element_blank()) +
  geom_hline(data = subset(cutoff.df, bat_species=="Eidolon dupreanum"), aes(yintercept = MFI), color="red", linetype=2) + 
  geom_line(aes(x=processing_date, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1)+
  geom_point(aes(x=processing_date, y=MFI, fill=barcode_id, shape = bat_sex, alpha = simple_class), show.legend = F, size=3) +
  scale_x_continuous(breaks = c(as.Date("2014-07-01"),as.Date("2015-07-01"), as.Date("2016-07-01"), as.Date("2017-07-01"), as.Date("2018-07-01"), as.Date("2019-07-01")),labels = c(2014,2015,2016,2017,2018,2019)) +
  guides(alpha="none",  shape=guide_legend(order=1));tryptic.plot

#and plot with centered peak
tryptic.plot <- ggplot(data=triptic_df) + 
  facet_nested(bat_sex~., scales = "free", space = "free") + theme_bw() +
  scale_shape_manual(name="bat sex", values = shapez) + 
  xlab("Time relative to observed peak (t_rel)") +
  ylab("MFI") +
  scale_alpha_manual(values = alphaz2) + ylab("mean fluorescence intensity (MFI)")+
  theme(axis.title.y = element_text(size=16), 
        legend.box = "horizontal", strip.text.y = element_text(size=16),#, face = "italic"),
        strip.background.y = element_rect(fill="white"),
        strip.background.x = element_blank(),
        axis.text = element_text(size=14), 
        panel.spacing.y = unit(.1, "lines"),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,.3,.3,.3),"lines"),
        panel.spacing.x = unit(c(0,1,0,1,1),"lines"),
        panel.grid = element_blank()) +
  geom_hline(data = subset(cutoff.df, bat_species=="Eidolon dupreanum"), aes(yintercept = MFI), color="red", linetype=2) + 
  geom_line(aes(x=new_ts, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1)+
  geom_point(aes(x=new_ts, y=MFI, fill=barcode_id, shape = bat_sex, alpha = simple_class), show.legend = F, size=3) +
  guides(alpha="none",  shape=guide_legend(order=1));tryptic.plot

#interpolate that curve for the top peak

plot1 <- ggplot(data=subset(triptic_df, barcode_id=="956000006105427")) + 
  facet_nested(bat_sex~., scales = "free", space = "free") + theme_bw() +
  scale_shape_manual(name="bat sex", values = shapez) + 
  xlab("Time relative to observed peak (t_rel)") +
  ylab("MFI") +
  scale_alpha_manual(values = alphaz2) + ylab("mean fluorescence intensity (MFI)")+
  theme(axis.title.y = element_text(size=16), 
        legend.box = "horizontal", strip.text.y = element_text(size=16),#, face = "italic"),
        strip.background.y = element_rect(fill="white"),
        strip.background.x = element_blank(),
        axis.text = element_text(size=14), 
        panel.spacing.y = unit(.1, "lines"),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,.3,.3,.3),"lines"),
        panel.spacing.x = unit(c(0,1,0,1,1),"lines"),
        panel.grid = element_blank()) +
  geom_hline(data = subset(cutoff.df, bat_species=="Eidolon dupreanum"), aes(yintercept = MFI), color="red", linetype=2) + 
  geom_line(aes(x=new_ts, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1)+
  geom_point(aes(x=new_ts, y=MFI, fill=barcode_id, shape = bat_sex, alpha = simple_class), show.legend = F, size=3) +
  guides(alpha="none",  shape=guide_legend(order=1));plot1

pre.MFI <- lm(MFI~new_ts, data = subset(triptic_df, barcode_id=="956000006105427" & new_ts<1))
post.MFI <- lm(MFI~new_ts, data = subset(triptic_df, barcode_id=="956000006105427" & new_ts>-1))

#MFI = 1141 is the cutoff for seropositivity
add.pseudo <- triptic_df[triptic_df$barcode_id=="956000006105427",]
add.pseudo <- add.pseudo[1:2,]
add.pseudo$MFI <- c(predict(pre.MFI, newdata = data.frame(new_ts = -90)), predict(post.MFI, newdata = data.frame(new_ts = 100)))
add.pseudo$new_ts <- c(-90, 100)
add.pseudo$seropos <- "seronegative"

triptic_df <- rbind(triptic_df, add.pseudo)


# and plot with centered peak again 
tryptic.plot2 <- ggplot(data=triptic_df) + 
  facet_nested(bat_sex~., scales = "free", space = "free") + theme_bw() +
  scale_shape_manual(name="bat sex", values = shapez) + 
  xlab("Time relative to observed peak (t_rel)") +
  ylab("MFI") +
  scale_alpha_manual(values = alphaz2) + ylab("mean fluorescence intensity (MFI)")+
  theme(axis.title.y = element_text(size=16), 
        legend.box = "horizontal", strip.text.y = element_text(size=16),#, face = "italic"),
        strip.background.y = element_rect(fill="white"),
        strip.background.x = element_blank(),
        axis.text = element_text(size=14), 
        panel.spacing.y = unit(.1, "lines"),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,.3,.3,.3),"lines"),
        panel.spacing.x = unit(c(0,1,0,1,1),"lines"),
        panel.grid = element_blank()) +
  geom_hline(data = subset(cutoff.df, bat_species=="Eidolon dupreanum"), aes(yintercept = MFI), color="red", linetype=2) + 
  geom_line(aes(x=new_ts, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1)+
  geom_point(aes(x=new_ts, y=MFI, fill=barcode_id, shape = bat_sex, alpha = simple_class), show.legend = F, size=3) +
  guides(alpha="none",  shape=guide_legend(order=1)) +
  geom_vline(aes(xintercept = 28))+ #duration 2 seropositive timesteps
  geom_vline(aes(xintercept = 42))+ #duration 3 seropositive timesteps
  geom_vline(aes(xintercept = 135)) #duration 9 seropositive timesteps - this is the mean
tryptic.plot2


# and plot with centered peak again 
tryptic.plot3 <- ggplot(data=triptic_df) + 
  #facet_nested(bat_sex~., scales = "free", space = "free") + 
  theme_bw() +
  scale_shape_manual(name="bat sex", values = shapez) + 
  xlab("Time relative to observed peak (t_rel)") +
  ylab("MFI") +
  scale_alpha_manual(values = alphaz2) + ylab("mean fluorescence intensity (MFI)")+
  theme(axis.title = element_text(size=18), 
        legend.box = "horizontal", strip.text.y = element_text(size=16),#, face = "italic"),
        strip.background.y = element_rect(fill="white"),
        strip.background.x = element_blank(),
        axis.text = element_text(size=14), 
        panel.spacing.y = unit(.1, "lines"),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,.3,.3,.3),"lines"),
        panel.spacing.x = unit(c(0,1,0,1,1),"lines"),
        panel.grid = element_blank()) +
  geom_hline(data = subset(cutoff.df, bat_species=="Eidolon dupreanum"), aes(yintercept = MFI), linetype=2) + 
  geom_line(aes(x=new_ts, y=MFI, color=barcode_id, alpha = simple_class), show.legend = F, size=1)+
  geom_point(aes(x=new_ts, y=MFI, fill=barcode_id, shape = bat_sex, alpha = simple_class), show.legend = F, size=3) +
  guides(alpha="none",  shape=guide_legend(order=1)) #+
tryptic.plot3

#tthis would only be supplemental but would be good to keep 
#
ggsave(file = paste0(homewd, "/final_figures_tables/tmp/tryptic-plot.png"),
       plot = tryptic.plot3,
       units="mm",  
       width=65, 
       height=55, 
       scale=2.5, 
       dpi=300)

#save for stan and include marker of the last seroneg and first seropos for each bat

names(triptic_df)
triptic.list <- dlply(triptic_df,.(barcode_id))

prep.data.stan <- function(df){
  df$t_obs <- df$new_ts + abs(min(df$new_ts))+1
  df$last_neg <- max(df$t_obs[df$new_ts<0 & df$seropos=="seronegative"])
  df$first_pos <- min(df$t_obs[df$seropos=="seropositive"])
  
  return(df)
}

triptic_df <- data.table::rbindlist(lapply(triptic.list, prep.data.stan))

#now run stan
# Map barcode_id to consecutive integers for Stan
bat_levels <- unique(triptic_df$barcode_id)
triptic_df <- triptic_df %>%
  mutate(bat_idx = as.numeric(factor(barcode_id, levels = bat_levels)))

# Prepare Stan data
stan_data <- list(
  N = nrow(triptic_df),
  N_bats = length(bat_levels),
  bat_id = triptic_df$bat_idx,
  t_obs = triptic_df$t_obs,
  MFI = triptic_df$MFI,
  t_exp_min = triptic_df$last_neg,
  t_exp_max = triptic_df$first_pos
)

# Fit Stan model to tryptic data get trajectory of antibody dynamics
# in exposed adults
fit <- stan(file = paste0(homewd,"/model_fitting/GhV-fits/latent_exposure_model.stan"),
            data = stan_data,
            iter = 10000, chains = 4, cores = 4, warmup = 2000,
            control = list(adapt_delta = 0.99))
print(fit, pars = c("A", "B", "sigma", "t_exp"), probs = c(0.05, 0.5, 0.95))

posterior <- as.data.frame(fit)

#get the hyperparameters (population-level)
hyperpars <- posterior %>%
  select(mu_A, sigma_A, mu_B, sigma_B, sigma)

# quick summary
summary(hyperpars)

#get the bat-specific parameters 
#(A and B for each of the 4 bats included in the fitting)
bat_params <- posterior %>%
  select(matches("^A\\[|^B\\["))  # matches A[1], A[2], ... B[1], B[2], ...

# latent exposure times
t_exp_post <- posterior %>%
  select(matches("^t_exp\\["))     # t_exp[1], t_exp[2], ...

bat_summary <- bat_params %>%
  summarise(across(everything(), list(
    mean = ~mean(.),
    median = ~median(.),
    q05 = ~quantile(., 0.05),
    q95 = ~quantile(., 0.95)
  )))

head(bat_summary)

# to plot the individual bat fits:

sim_fun <- function(A, B, t_exp_val, sigma, t_seq) {
  tau <- pmax(t_seq - t_exp_val, 0)
  A * tau * exp(-B * tau) + rnorm(length(tau), 0, sigma)
}

# simulate trajectory for bat 1
sim_traj <- sim_fun(A=mean(posterior$`A[1]`), B=mean(posterior$`B[1]`), t_exp_val=mean(posterior$`t_exp[1]`), sigma=mean(posterior$sigma), t_seq=0:150)

plot(0:150, sim_traj) #this is the Ab trajectory for bat 1
t_seq = 0:150

sim_list <- lapply(1:4, function(i) {
  
  A_hat <- mean(posterior[[paste0("A[", i, "]")]])
  B_hat <- mean(posterior[[paste0("B[", i, "]")]])
  sigma_hat <- mean(posterior[["sigma"]])
  t_exp_val <- mean(posterior[[paste0("t_exp[", i, "]")]])
  data.frame(
    t = t_seq,
    MFI_sim = sim_fun(A_hat, B_hat, t_exp_val, sigma_hat, t_seq),
    bat = i
  )
})

sim_df <- do.call(rbind, sim_list)

stan_df <- as.data.frame(stan_data)

p.ind.bat <- ggplot(data=sim_df)+
  geom_line(aes(x=t, y=MFI_sim,group=bat, color=as.factor(bat) )) + theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "inside", legend.position.inside = c(.85,.85)) +
  geom_point(data=stan_df, (aes(x=t_obs, y=MFI, color=as.factor(bat_id)))); p.ind.bat
# here are the individual trajectories for our 4 tryptic bats

#now simulate many new bats from the population-level hyperparameters
hyperpars <- posterior %>%
  select(mu_A, sigma_A, mu_B, sigma_B, sigma)



simulate_new_bat_stochastic <- function(mu_A, sigma_A, mu_B, sigma_B, sigma, t_seq, t_exp = 0, MFI_min = 0, MFI_max = 50000, shrink_factor) {
  
  # draw bat-level parameters from population
  #shrink_factor = 0 is fully shrunk to mean
  #shrink_factor = 1 is full sigma variance
  A <- rnorm(1, mean = mu_A, sd = sigma_A*shrink_factor)
  B <- rnorm(1, mean = mu_B, sd = sigma_B*shrink_factor)
  
  tau <- pmax(t_seq - t_exp, 0)
  
  # simulate raw MFI
  MFI <- A * tau * exp(-B * tau) + rnorm(length(tau), 0, sigma)
  
  # constrain MFI to [MFI_min, MFI_max]
  MFI <- pmin(pmax(MFI, MFI_min), MFI_max)
  
  return(MFI)
}


n_new_bats <- 300
t_seq <- seq(0, 150, by = 1)   # time in days

# sample 300 posterior rows for variation
posterior_rows <- hyperpars[sample(nrow(hyperpars), n_new_bats, replace = TRUE), ]

sim_list <- lapply(1:n_new_bats, function(i) {
  row <- posterior_rows[i, ]
  data.frame(
    t = t_seq,
    MFI_sim = simulate_new_bat_stochastic(row$mu_A, row$sigma_A, row$mu_B, row$sigma_B, row$sigma, t_seq, shrink_factor = 0),
    bat = i
  )
})


sim_df <- do.call(rbind, sim_list)

#and the deterministic
simulate_new_bat_deterministic <- function(A, B, t_seq, t_exp = 0) {
  
  A <- mean(posterior$mu_A)
  B <- mean(posterior$mu_B)
  t_exp <- 0   # or mean(posterior$t_exp[bat_id]) if you want average exposure
  
  
  tau <- pmax(t_seq - t_exp, 0)
  
  # simulate raw MFI
  MFI <- A * tau * exp(-B * tau) 
  
  
  
  return(MFI)
}

det.bat <- cbind.data.frame(t=t_seq,
                            MFI=simulate_new_bat_deterministic(A=row$mu_A, B=row$mu_B, t_seq))


#plot sims
p3 <- ggplot() +
  geom_line(data= sim_df, aes(x = t, y = MFI_sim, group = bat),alpha = 0.2, color = "steelblue") +
  labs(title = "Simulated antibody trajectories for 300 new bats",
       x = "Days since exposure",
       y = "MFI") +
  geom_line(data=det.bat, aes(x=t, y=MFI),size=1) +
  theme(panel.grid = element_blank()) +
  theme_bw() + coord_cartesian(ylim=c(0,40000)) + 
  geom_hline(aes(yintercept =cutoff.df$MFI[cutoff.df$bat_species=="Eidolon dupreanum"]), linetype=2); p3

#get the peak of the deterministic model
peak.det = det.bat$t[det.bat$MFI==max(det.bat$MFI)]

sub_df <- subset(sim_df, bat<=100 )
sub_df <- subset(sim_df, bat>100 & bat<=200)
sub_df <- subset(sim_df, bat>200 & bat<=300)



#align stan_df bats around the peak
stan_df$t_exp_max[stan_df$bat_id==2] <- 91
stan_df$new_ts <- stan_df$t_obs- stan_df$t_exp_max + peak_det

#plot with data (100 sims)
p4 <- ggplot() +
  geom_line(data=sub_df , aes(x = t, y = MFI_sim, group = bat, alpha=as.factor(bat)), color = "gray", show.legend = F) +
  labs(title = "Simulated adult antibody trajectories for 100 bats",
       x = "Days since exposure",
       y = "MFI") +
  geom_line(data=det.bat, aes(x=t, y=MFI),size=1) +
  theme_bw() + coord_cartesian(ylim=c(0,40000), xlim=c(0,150)) + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=18),
        plot.title = element_text(size=18),
        axis.text = element_text(size=14)) +
  geom_hline(aes(yintercept =cutoff.df$MFI[cutoff.df$bat_species=="Eidolon dupreanum" & cutoff.df$type=="mean"]), linetype=2) +
  geom_point(data=stan_df, (aes(x=new_ts, y=MFI, color=as.factor(bat_id))), size=3, show.legend = F); p4
  #geom_line(data=stan_df, aes(x=t_obs, y=MFI, color=as.factor(bat_id)), linetype=2); p4

ggsave(file = paste0(homewd, "/final_figures_tables/tmp/within-host-adultAb-fit.png"),
       plot = p4,
       units="mm",  
       width=65, 
       height=55, 
       scale=2.5, 
       dpi=300)


min(det.bat$t[det.bat$MFI<cutoff.df$MFI[cutoff.df$bat_species=="Eidolon dupreanum"] & det.bat$t>25])
# 135 day minimum duration Ab immunity from the deterministic model fit
# but the stochastic simulations show a lot of variation
# Deterministic output is an Ab duration of 4.5 months - 
# but we start seeing many dropoff at ~75 days (2.5 months)

# for this reason, we decided to take the same approach for adult Ab immunity
# as for maternal: we fixed two timesteps of the R class (28 days total, 2x biweeks)
# as protected Ab-mediated immunity, then estimated the rate of exiting from 
# the third class


### now take a look at the active infection data to evaluate whether 
### there is enough shedding data to develop a similar process for the I class
### Verdict: THERE IS NOT - commented out
# rm(list=ls())
# https://files.osf.io/v1/resources/kqdsn_v1/providers/osfstorage/655c82f12ee3334478c4165b?action=download&direct&version=1
# ################################################################################
# ################################################################################
# homewd = "/Users/carabrook/Developer/luminex-serology-madabat"
# batwd = "/Users/carabrook/Developer/Madagascar-Bat-Data"
# ################################################################################
# ################################################################################
# 
# shed.dat <- read.csv(paste0(homewd,"/Gwen_AngV_GhV_viral_shedding/indiv_urine_qPCR_AngV_GhV_data.csv"), header = T, stringsAsFactors = F)
# head(shed.dat)
# 
# #and take the positives as any that are positive by both M and L (the two most sensitive targets for GhV)
# names(shed.dat)
# shed.dat <- dplyr::select(shed.dat, sampleid, qPCR_GhV_multi_L,  qPCR_GhV_multi_L_ct, qPCR_GhV_multi_M, qPCR_GhV_multi_M_ct)
# 
# catch.dat <- read.csv(paste0(batwd,"/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)
# head(catch.dat)
# 
# #and merge a few things
# catch.merge <- dplyr::select(catch.dat,barcode_id, sampleid, roost_site, processing_date, bat_species, bat_age_class, recap_new)
# head(catch.merge)
# 
# shed.merge <- merge(shed.dat, catch.merge, by="sampleid", all.x = T)
# head(shed.merge)
# #take eidolon only
# shed.merge = subset(shed.merge, bat_species == "Eidolon dupreanum")
# head(shed.merge)
# 
# #now ID those that are positive by the two targets
# shed.merge$shed_pos <- 0
# shed.merge$shed_pos[shed.merge$qPCR_GhV_multi_L==1 & shed.merge$qPCR_GhV_multi_M==1] <- 1
# 
# sum(shed.merge$shed_pos) # 70
# 
# #are there any recaps?
# shed.recap = subset(shed.merge, recap_new=="yes")
# #ANGB110 is a mess if I recall, so remove that
# shed.recap = subset(shed.recap, sampleid!="ANGB110")
# 
# shed.recap$ct <- 0
# shed.recap$ct[shed.recap$shed_pos==1] <- shed.recap$qPCR_GhV_multi_L_ct[shed.recap$shed_pos==1]
# 
# #also positive
# shed.recap.pos = subset(shed.recap, shed_pos==1)
# #a few... 
# #any repeat captures?
# duplicated(shed.recap.pos$barcode_id) #no... 
# # so any recaps we have do not have multiple Ct data corresponding
# shed.recap$processing_date <- as.Date(shed.recap$processing_date)
# 
# #nothing positive at multiple timepoints
# ggplot(shed.recap) + 
#   geom_point(aes(x=processing_date, y=ct, color=barcode_id), show.legend = F) +
#   geom_line(aes(x = processing_date,y=ct, color=barcode_id, group=barcode_id), show.legend = F)
# 
# 
# #and now load the new luminex data 
# lumdat <- read.csv(file = paste0(homewd, "/data/mix_sub_data_75_90_96.csv"))
# head(lumdat)
# 
# #slim
# lumdat <- dplyr::select(lumdat, sampleid,bat_species, antigen, MFI, seropos, seropos_lci, seropos_uci)
# 
# lumdat$log10MFI <- log10(lumdat$MFI)
# 
# #focus just on GhV
# lumdat = subset(lumdat, antigen=="GhV"  & bat_species=="Eidolon dupreanum")
# head(lumdat)
# lum.merge <- dplyr::select(lumdat, sampleid, MFI, seropos)
# 
# shed.recap <- merge(shed.recap, lum.merge, by="sampleid", all.x = T)
# head(shed.recap)
# shed.recap$seropos <- as.character(shed.recap$seropos)
# shed.recap$seropos[shed.recap$seropos=="NA"] <- "no data"
# shed.recap$seropos[shed.recap$seropos=="0"] <- "seroneg"
# shed.recap$seropos[shed.recap$seropos=="1"] <- "seropos"
# shed.recap$seropos[is.na(shed.recap$seropos)] <- "no data"
# shed.recap$seropos <- factor(shed.recap$seropos, levels=c("seropos", "seroneg", "no data"))
# 
# ggplot(shed.recap) + 
#   geom_point(aes(x=processing_date, y=ct, color=barcode_id, shape=seropos), size=3) +
#   geom_line(aes(x = processing_date,y=ct, color=barcode_id, group=barcode_id), show.legend = F) +
#   guides(color="none")
# 
# 
# #do any of these have values within 6 months of one another?
# head(shed.recap)
# shed.recap <- arrange(shed.recap, barcode_id, processing_date)
# 
# #and split
# shed.split <- dlply(shed.recap, .(barcode_id))
# 
# #find 6 month pairs that have 1 positive
# find.pairs <- function(df){
#   if(sum(df$shed_pos>0)){
#     ts.diff = c(as.numeric(diff(df$processing_date)))
#     if(length(ts.diff[ts.diff<180])>0){
#       return(df)
#     }
#   }
# }
# pos.pairs <-data.table::rbindlist(lapply(shed.split,find.pairs))
# pos.pairs # two pairs
# 
# # These are the only two recaps with one positive recaptured within 6 months
# # we can't really learn anything about the duration of infection from these...
# # The one pregnant female who is molecular positive and seropositive in early 
# # October and then molecular negative and seronegative when not pregnant in
# # December is fascinating
# 
# ggplot(pos.pairs) + #facet_grid(~barcode_id) +
#   geom_point(aes(x=processing_date, y=ct, color=barcode_id, shape=seropos), size=3) +
#   geom_line(aes(x = processing_date,y=ct, color=barcode_id, group=barcode_id), show.legend = F) +
#   guides(color="none")
# 
# # We'll keep the assumption of infection duration of two weeks
