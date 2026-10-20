rm(list=ls())
#set wd
#.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4")

library(plyr)
library(dplyr)   
library(epitools)
library(deSolve)
library(lubridate)
library(matrixcalc)
library(Matrix)
library(mvtnorm)
library(tidyr)
library(scales)
library(ggplot2)

#this is the version that builds a comparable matrix to fit the top five model structures
#load data
homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/"
setwd(homewd)
#dat <- read.csv(file =  paste0(homewd,"/data/age-seroprevalence-data.csv"), header = T, stringsAsFactors = F)
# dat <- read.csv(file =  "age-seroprevalence-data.csv", header = T, stringsAsFactors = F)
# sero.dat <- read.csv(file =  "biweekly-serology-data.csv", header = T, stringsAsFactors = F)
# shed.dat <- read.csv(file =  "biweekly-shedding-data.csv", header = T, stringsAsFactors = F)
# head(dat)
# head(sero.dat)
# head(shed.dat)




# w_1_star = proportion of bats in youngest age class at equilibrium - extracted from model simulations by taking the first element of 
# the dominant right eigenvector of your demographic transition matrix; w_1_star = 0.276
# w_2_star = proportion of bats in the youngest ADULT age class at equilibrium - extracted from model simulations by taking the
# second element of the dominant right eigenvector of your demographic transition matrix; w_2_star = 0.150
# age_rate = 1/26


calculate_R0_single <- function(par.df, w_1_star) {
  
  # Parameters
  beta    <- par.df$beta
  m_J     <- par.df$mort_juv     # Biweekly juvenile survival: (annual_s_J)^(1/26)
  m_A     <- par.df$mort      # Biweekly adult survival: (annual_s_A)^(1/26)
  #w_1_star <- hardcoded  # 1st element of normalized right eigenvector
  
  # Mortality multiplier (default to 1 if no elevated mortality)
  mu_I    <- par.df$mu_sick
  
  # State-specific biweekly survival probabilities on biweekly timestep
  s_J_I   <- (1 - (m_J*mu_I))
  if(is.na(s_J_I)|(s_J_I<0)){s_J_I=0}
  s_J_I_biwk   <- s_J_I^(1/26)
  s_A_I   <- (1 - (m_A*mu_I))
  if(is.na(s_A_I)|(s_A_I<0)){s_A_I=0}
  s_A_I_biwk   <- s_A_I^(1/26)
  
  # R0 Calculation
  R0 <- beta * (s_J_I_biwk * w_1_star +   s_A_I_biwk* (1 - w_1_star))
  return(R0)
}
calculate_R0_recrudescence <- function(par.df, w_1_star, age_rate) {
  
  # Base Parameters
  beta     <- par.df$beta
  m_A      <- par.df$mort
  m_J      <- par.df$mort_juv
  a        <- age_rate         # Aging rate: 1/26
  # w_2_star  hard-coded
  
  # Mortality multipliers (default to 1 if no elevated mortality)
  mu_I     <- par.df$mu_sick
  mu_R     <- par.df$mu_rec
  
  # Epidemiological transition rates
  gamma    <- par.df$recru
  gamma_2  <- par.df$recru2
  sigma_1  <- par.df$sigma1
  sigma_2  <- par.df$sigma2
  rho      <- par.df$loss
  epsilon  <- par.df$epsilon
  
  # State-specific biweekly survival probabilities
  s_J_I   <- (1 - (m_J*mu_I))
  if(is.na(s_J_I)|(s_J_I<0)){s_J_I=0}
  s_J_I_biwk   <- s_J_I^(1/26)
  
  s_A_I   <- (1 - (m_A*mu_I))
  if(is.na(s_A_I)|(s_A_I<0)){s_A_I=0}
  s_A_I_biwk   <- s_A_I^(1/26)
  
  s_J_R   <- (1 - (m_J*mu_R))
  if(is.na(s_J_R)|(s_J_R<0)){s_J_R=0}
  s_J_R_biwk   <- s_J_R^(1/26)
  
  s_A_R   <- (1 - (m_A*mu_R))
  if(is.na(s_A_R)|(s_A_R<0)){s_A_R=0}
  s_A_R_biwk   <- s_A_R^(1/26)
  
  s_J_N_biwk    <- (1 - m_J)^(1/26)  
  s_A_N_biwk    <- (1 - m_A)^(1/26)  # Baseline survival for carrier/waning class N
  
  # State-specific stay terms (survival * (1 - aging))
  s_stay_I <- s_A_I_biwk * (1 - a)
  s_stay_R <- s_A_R_biwk * (1 - a)
  s_stay_N <- s_A_N_biwk * (1 - a)
  
  # Loop denominators (expected duration in R3 and N1 states)
  denom_R3 <- 1 - (1 - sigma_1 - sigma_2 - gamma_2) * s_stay_R
  denom_N1 <- 1 - (1 - rho - gamma - epsilon) * s_stay_N
  
  # Recrudescence return pathways
  pathway_1 <- (gamma_2 * s_stay_R) / denom_R3
  pathway_2 <- (gamma * sigma_2 * s_stay_R * s_stay_N) / (denom_R3 * denom_N1)
  
  # Total return probability (includes 3-step delay through R1 -> R2 -> R3)
  R_return  <- (pathway_1 + pathway_2) * (s_stay_R^3)
  
  # Numerator (Primary infectious transmission)
  # Option B Numerator: Primary infection uses total biweekly survival
  numerator <- beta * (s_J_I_biwk * w_1_star + s_A_I_biwk * (1 - w_1_star))
  
  # Denominator: Continues using s_stay_R = s_A_R * (1 - a) and s_stay_N = s_A * (1 - a)
  denom_R3 <- 1 - (1 - sigma_1 - sigma_2 - gamma_2) * s_stay_R
  denom_N1 <- 1 - (1 - rho - gamma - epsilon) * s_stay_N
  
  pathway_1 <- (gamma_2 * s_stay_R) / denom_R3
  pathway_2 <- (gamma * sigma_2 * s_stay_R * s_stay_N) / (denom_R3 * denom_N1)
  
  R_return  <- (pathway_1 + pathway_2) * (s_stay_R^3)
  
  # Final R0 Calculation
  R0 <- numerator / (1 - R_return)
  
  
  
  return(R0)
}


wrap.R0 <- function(par.df, w_1_star,  age_rate){
  #adult and juvenile annual survival need to be expressed in biweeks
  
  #get the base model
  par.df$model_base <- sapply(strsplit(par.df$model, "-mu"),'[',1)
  
  par.df$R0_method <- NA
  for (i in 1:nrow(par.df)){
    if (par.df$model_base[i] == "MSIRN" | par.df$model_base[i] == "MSIRS" | par.df$model_base[i] == "MSIRNS" | par.df$model_base[i] == "MSIRNNS"){
      par.df$R0_method[i] <- "single-infected-class"
    }else{
      par.df$R0_method[i] <- "infection-recrudescence"
    }
  }
  
  par.df$R0 <- NA
  
  
  
  for (i in 1:length(par.df$model)){
    if (par.df$R0_method[i]=="single-infected-class"){
      
      par.df$R0[i] <- calculate_R0_single(par.df = par.df[i,], w_1_star = w_1_star)
      
    }else if(par.df$R0_method[i]=="infection-recrudescence"){
      
      par.df$R0[i] <- calculate_R0_recrudescence(par.df = par.df[i,], w_1_star = w_1_star, age_rate = age_rate)
    }
    
  }
  
  return(par.df)
  
  
  
}

#load the best fit parameters

#load the completed models with age-seroprevalence and time series
all.par <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/par_fits_GhV_CI_final.csv"), header = T, stringsAsFactors = F)


all.R0 <- wrap.R0(par.df = all.par, w_1_star = 0.276, age_rate = 1/26)


#now reorder columns and rearrange rows:
sort.fitted.pars <- function(df){
  par.list <- dlply(df,.(model_base))
  
  get.avg <- function(df){
    df$mean_AIC <- mean(df$AIC)
    df$mean_llik <- mean(df$llik)
  return(df)
  }

  
  df <- data.table::rbindlist(lapply(par.list, get.avg))
  df <- arrange(df, mean_AIC, llik)
  
  #and rearrange columns
  names(df)
  df <- dplyr::select(df, model_base, model, mort, mort_juv, adult_fec, recov, sigma_fix, wane_fix, wane,
                      beta, beta_lci, beta_uci, sigma1, sigma1_lci, sigma1_uci, sigma2, sigma2_lci, sigma2_uci,
                      recru, recru_lci, recru_uci, recru2, recru2_lci, recru2_uci, loss, loss_lci, loss_uci,
                      epsilon, epsilon_lci, epsilon_uci, mu_sick, mu_sick_lci, mu_sick_uci, mu_rec, mu_rec_lci, mu_rec_uci,
                      weightPCR, llik, AIC, mean_AIC, mean_llik, R0_method, R0)
      
return(df)
  }

save.R0 <- sort.fitted.pars(all.R0)

write.csv(save.R0, file = paste0(homewd,"/model_fitting/GhV-fits/final_par_fits_wR0.csv"), row.names = F)

#and include the additional
save.R0$mort_class <- gsub(pattern = "-mu", replacement = "mu", x= sapply(strsplit(save.R0$model, save.R0$model_base),'[',2))
save.R0$mort_class[is.na(save.R0$mort_class)] <- "base-model"
save.R0$mort_class <- gsub(pattern = "sickmu", replacement = "sick-mu", x= save.R0$mort_class)

save.R0 <- arrange(save.R0, mean_AIC, llik)
save.R0$model_base <- factor(save.R0$model_base, levels = unique(save.R0$model_base))

#and plot
p1 <- ggplot(save.R0) + geom_point(aes(x=llik, y=R0, color=model_base, shape=mort_class), size=3, show.legend = F) + scale_x_log10() +theme_bw() +
      xlab("negative log-likelihood") + ylab(bquote(R[0]))+ scale_color_discrete(name = "base model") + scale_shape_discrete("mortality class") +
      theme(panel.grid = element_blank(), axis.title = element_text(size=14), axis.text = element_text(size=12));p1

p2 <- ggplot(save.R0) + geom_point(aes(x=llik, y=R0, color=model_base, shape=mort_class), size=3) + scale_x_log10() +theme_bw() +
  xlab("negative log-likelihood") + ylab(bquote(R[0]))+ scale_color_discrete(name = "base model") + scale_shape_discrete("mortality class") +
  theme(panel.grid = element_blank(), axis.title = element_text(size=14), axis.text = element_text(size=12)) + coord_cartesian(xlim=c(628,646))



pR0 <- cowplot::plot_grid(p1,p2, labels = c("A", "B"), label_size = 22, rel_widths = c(1,1.3))

#and also include a plot of the fitted beta par by model for top ten model fits
all.pars <- read.csv(file = paste0(homewd,"/model_fitting/GhV-fits/top-10-matrix-round1.csv"), header = T, stringsAsFactors = F)
head(all.pars) #top 10 fits by llik for each model construction
unique(all.pars$model)

#now plot beta value for each model
all.pars$model_base <- sapply(strsplit(all.pars$model,"-"), '[',1) 
all.pars$mort_class <- gsub(pattern = "-mu", replacement = "mu", x= sapply(strsplit(all.pars$model, all.pars$model_base),'[',2))
all.pars$mort_class[is.na(all.pars$mort_class)] <- "base-model"
all.pars$mort_class <- gsub(pattern = "sickmu", replacement = "sick-mu", x= all.pars$mort_class)

#and rank them by mean AIC
list.pars <- dlply(all.pars,.(model_base))
get.mean.AIC.llik <- function(df){
  df$mean_AIC <- mean(df$AIC)
  df$mean_llik <- mean(df$llik)
  return(df)
}
all.pars <- data.table::rbindlist(lapply(list.pars,get.mean.AIC.llik))
all.pars <- arrange(all.pars, mean_AIC)
all.pars$model_base <- factor(all.pars$model_base, levels=unique(all.pars$model_base))

ggplot(all.pars) + geom_point(aes(x=llik, y=beta, color=model_base, shape=mort_class)) + coord_cartesian(xlim=c(628,700), ylim=c(0,100))


all.pars$log10beta <- log10(all.pars$beta)


beta.plot <- ggplot(all.pars) + #coord_cartesian(ylim=c(-1,10)) +# scale_x_log10() +
                  geom_boxplot(aes(x=model_base, y=log10beta), show.legend = F) +
                  geom_jitter(aes(x=model_base, y=log10beta, color=model_base, shape=mort_class), width = .15, show.legend = F) +
                  theme_bw() +coord_cartesian(ylim=c(-1,6)) +
                  ylab(bquote(log[10]~beta~', top 10 model fits'))+ scale_color_discrete(name = "base model") + 
                  scale_shape_discrete("mortality class") +
                  theme(panel.grid = element_blank(), 
                        axis.title.y = element_text(size=14),
                        axis.title.x = element_blank(),
                        axis.text.y = element_text(size=12),
                        axis.text.x = element_text(size=12, angle = 45, vjust=-.0001,margin = margin(t = -20)))


all.pars$any_sigma <- all.pars$sigma2
all.pars$any_sigma[all.pars$any_sigma==0] <- all.pars$sigma1[all.pars$any_sigma==0]
all.pars$any_sigma_log10 <- log10(all.pars$any_sigma) 
all.pars$any_sigma_log10[all.pars$any_sigma_log10 ==-Inf] <- NA

sigma.plot <- ggplot(all.pars) + #coord_cartesian(ylim=c(-1,10)) +# scale_x_log10() +
  geom_boxplot(aes(x=model_base, y=any_sigma_log10), show.legend = F) +
  geom_jitter(aes(x=model_base, y=any_sigma_log10, color=model_base, shape=mort_class), width = .15, show.legend = F) +
  theme_bw() + coord_cartesian(ylim=c(-1,23)) + coord_cartesian(ylim=c(-1,25)) +
  ylab(bquote(log[10]~sigma~', top 10 model fits'))+ scale_color_discrete(name = "base model") + 
  scale_shape_discrete("mortality class") +
  theme(panel.grid = element_blank(), 
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 45, vjust=-.0001,margin = margin(t = -20)))



all.pars$any_recru <- all.pars$recru2
all.pars$any_recru[all.pars$any_recru==0] <- all.pars$recru[all.pars$any_recru==0]
all.pars$any_recru_log10 <- log10(all.pars$any_recru) 
all.pars$any_recru_log10[all.pars$any_recru_log10==-Inf] <- NA

recru.plot <- ggplot(all.pars) + #coord_cartesian(ylim=c(-1,10)) +# scale_x_log10() +
  geom_boxplot(aes(x=model_base, y=any_recru_log10), show.legend = F) +
  geom_jitter(aes(x=model_base, y=any_recru_log10, color=model_base, shape=mort_class), width = .15, show.legend = F) +
  theme_bw() + coord_cartesian(ylim=c(-1,6)) +
  ylab(bquote(log[10]~gamma~', top 10 model fits'))+ scale_color_discrete(name = "base model") + 
  scale_shape_discrete("mortality class") +
  theme(panel.grid = element_blank(), 
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 45, vjust=-.0001, margin = margin(t = -20)))



all.pars$any_loss <- all.pars$loss
all.pars$any_loss_log10 <- log10(all.pars$any_loss) 
all.pars$any_loss_log10[all.pars$any_loss_log10==-Inf] <- NA

loss.plot <- ggplot(all.pars) + #coord_cartesian(ylim=c(-1,10)) +# scale_x_log10() +
  geom_boxplot(aes(x=model_base, y=any_loss_log10), show.legend = F) +
  geom_jitter(aes(x=model_base, y=any_loss_log10, color=model_base, shape=mort_class), width = .15, show.legend = F) +
  theme_bw() + coord_cartesian(ylim=c(-1,6)) +
  ylab(bquote(log[10]~rho~', top 10 model fits'))+ scale_color_discrete(name = "base model") + 
  scale_shape_discrete("mortality class") +
  theme(panel.grid = element_blank(), 
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 45, vjust=-.0001, margin = margin(t = -20)))

#and combined



pcomppar <- cowplot::plot_grid(beta.plot, sigma.plot, recru.plot, loss.plot, ncol = 2, nrow = 2, labels = c("C", "D", "E", "F"), label_size = 22, scale = .97) +
  theme(plot.margin = margin(t=.1,r=8,b=2,l=.1, unit= "lines"))
  
            

pall <- cowplot::plot_grid(pR0, pcomppar, nrow=2, ncol=1, rel_heights  = c(1,2))

ggsave(file = paste0(homewd,"/final_figures_tables/supplemental/FigSX-model-comparison.pdf"),
       plot = pall,
       units="mm",  
       width=100, 
       height=140, 
       scale=3, 
       dpi=300)
