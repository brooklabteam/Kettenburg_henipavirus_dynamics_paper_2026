rm(list=ls())


library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(mgcv)
library(reshape2)
library(Hmisc)
library(exactci)
library(lme4)
library(sjPlot)

homewd = "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/"
batwd = "/Users/carabrook/Developer/Madagascar-Bat-Data"

#homewd = "/Users/gwenddolenkettenburg/Desktop/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
#batwd = "/Users/gwenddolenkettenburg/Desktop/Developer/Madagascar-Bat-Data"


# load comprehensive catching data 
catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))


head(catch.dat)
names(catch.dat)

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, pup_or_mom, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, recap_new)

#remove duplicate bats from dataset (bats mis-identified)
catch.dat = subset(catch.dat, sampleid!="ANGB110")
catch.dat = subset(catch.dat, sampleid!="MIZ545")
catch.dat = subset(catch.dat, sampleid!="MIZ548")
catch.dat = subset(catch.dat, sampleid!="AND6")
catch.dat = subset(catch.dat, sampleid!="RAK32")
catch.dat = subset(catch.dat, sampleid!="MIZ305" &sampleid!="MIZ306" & sampleid!="MIZ307" & sampleid!="MIZ546" & sampleid!="MIZ547")

#sub to timeframe of interest
catch.dat$processing_date <- as.Date(catch.dat$processing_date)
catch.dat <- subset(catch.dat, processing_date < "2021-01-01")
#catch.dat <- subset(catch.dat, processing_date < "2023-08-01")
unique(catch.dat$young_of_year)
catch.dat$bat_age_class[catch.dat$young_of_year=="yes" & !is.na(catch.dat$bat_age_class)] <- "J"
unique(catch.dat$bat_species)


adult.dat = subset(catch.dat, bat_species=="Pteropus rufus" & bat_age_class!="J" & !is.na(bat_age_class) & young_of_year!="yes"  | bat_species=="Eidolon dupreanum" & bat_age_class!="J" & !is.na(bat_age_class) & young_of_year!="yes"| bat_species=="Rousettus madagascariensis" & bat_age_class!="J" & !is.na(bat_age_class) & young_of_year!="yes")
unique(adult.dat$bat_species)
unique(adult.dat$bat_age_class)
unique(adult.dat$young_of_year)

length(adult.dat$bat_weight_g[is.na(adult.dat$bat_weight_g)])
length(adult.dat$bat_forearm_mm[is.na(adult.dat$bat_forearm_mm)])



adult.dat$bat_species <- factor(adult.dat$bat_species, levels=c("Rousettus madagascariensis", "Eidolon dupreanum", "Pteropus rufus"))


#and redo the MFR with this
# To make sure it is important to allow for different slopes and y-intercepts,
# you can test all possible models.

# Here, just test a model that is weight by forearm but no species differentiation
# you ca

modBaseF <- lm(log10(bat_weight_g)~log10(bat_forearm_mm), data=subset(adult.dat, bat_sex=="female"))
summary(modBaseF)

# Now, allow different slope per species:
# You have to keep optimization to REML=FALSE to be able to compare with lm()
modSlopeF <- lm(log10(bat_weight_g)~log10(bat_forearm_mm):(bat_species),  data=subset(adult.dat, bat_sex=="female"))

# Now, different y-intercept
modYIntF <- lmer(log10(bat_weight_g)~log10(bat_forearm_mm) +(1|bat_species), REML=FALSE, data=subset(adult.dat, bat_sex=="female"))

# And, finally, both unique y-intercept and slope
modYIntSlopeF <-  lmer(log10(bat_weight_g)~log10(bat_forearm_mm):bat_species +(1|bat_species), REML=FALSE, data=subset(adult.dat, bat_sex=="female"))

# compare AIC
AIC.dat.F = AIC(modBaseF, modSlopeF, modYIntF, modYIntSlopeF)
print(AIC.dat.F)
AIC.dat.F[which(AIC.dat.F$AIC==min(AIC.dat.F$AIC)),] # Here is your best fit model
# This allows for unique slopes per species 

# You can visualize these slopes here:
plot_model(modSlopeF)

# You can also visualize how this differs from the mixed effect model:
#plot_model(modYIntSlopeF, type="int") #it won't run because it can't handle the logarithms
# But you can get around it by making a logged variable in a new dataset
dat.mod = adult.dat
dat.mod$log10forearm = log10(dat.mod$bat_forearm_mm)
dat.mod$log10mass = log10(dat.mod$bat_weight_g)

modYIntSlopeFPlot <-  lmer(log10mass~log10forearm:bat_species +(1|bat_species), REML=FALSE, data=subset(dat.mod, bat_sex=="female"))
plot_model(modYIntSlopeFPlot, type="int")


#### Now, repeat for the males
modBaseM <- lm(log10(bat_weight_g)~log10(bat_forearm_mm), data=subset(adult.dat, bat_sex=="male"))
summary(modBaseM)

# Now, allow different slope per species:
# You have to keep optimization to REML=FALSE to be able to compare with lm()
modSlopeM <- lm(log10(bat_weight_g)~log10(bat_forearm_mm):(bat_species),  data=subset(adult.dat, bat_sex=="male"))

# Now, different y-intercept
modYIntM <- lmer(log10(bat_weight_g)~log10(bat_forearm_mm) +(1|bat_species), REML=FALSE, data=subset(adult.dat, bat_sex=="male"))

# And, finally, both unique y-intercept and slope
modYIntSlopeM <-  lmer(log10(bat_weight_g)~log10(bat_forearm_mm):bat_species +(1|bat_species), REML=FALSE, data=subset(adult.dat, bat_sex=="male"))
#singular fit means the random effect is not doing anything--could fit as just a linear model

# compare AIC
AIC.dat.M = AIC(modBaseM, modSlopeM, modYIntM, modYIntSlopeM)
print(AIC.dat.M)
AIC.dat.M[which(AIC.dat.M$AIC==min(AIC.dat.M$AIC)),] 
# Males do best with just different slope by species, same as females

# plot the best two 
plot_model(modSlopeM)

# need to remake the male model with logged variables again
modYIntSlopeMPlot <-  lmer(log10mass~log10forearm:bat_species +(1|bat_species), REML=FALSE, data=subset(dat.mod, bat_sex=="male"))
plot_model(modYIntSlopeMPlot, type="int") #As in the females, the Pteropus curve is steeper

## Conclusion: Best fit model is just different slopes by species with no random effect for both
## males and females so go with that.

# Add a place for "predicted mass" in your dataset, and predict from this model:
adult.dat$predicted_mass_g <- NA
adult.dat$predicted_mass_g[adult.dat$bat_sex=="female" & !is.na(adult.dat$bat_weight_g) &!is.na(adult.dat$bat_forearm_mm)] <- 10^predict(modSlopeF) #you are predicting log10 mass, so you need to exponent it ere
adult.dat$predicted_mass_g[adult.dat$bat_sex=="male"& !is.na(adult.dat$bat_weight_g) &!is.na(adult.dat$bat_forearm_mm)] <- 10^predict(modSlopeM)

#and the confidence intervals around that
adult.dat$predicted_mass_g_lci <- adult.dat$predicted_mass_g_uci <- NA

adult.dat$predicted_mass_g_lci[adult.dat$bat_sex=="female" & !is.na(adult.dat$bat_weight_g) &!is.na(adult.dat$bat_forearm_mm)] <- 10^(predict(modSlopeF) -1.96*predict(modSlopeF,se.fit=T)$se)
adult.dat$predicted_mass_g_uci[adult.dat$bat_sex=="female" & !is.na(adult.dat$bat_weight_g) &!is.na(adult.dat$bat_forearm_mm)] <- 10^(predict(modSlopeF) +1.96*predict(modSlopeF,se.fit=T)$se)
adult.dat$predicted_mass_g_lci[adult.dat$bat_sex=="male" & !is.na(adult.dat$bat_weight_g) &!is.na(adult.dat$bat_forearm_mm)] <- 10^(predict(modSlopeM) -1.96*predict(modSlopeM,se.fit=T)$se)
adult.dat$predicted_mass_g_uci[adult.dat$bat_sex=="male" & !is.na(adult.dat$bat_weight_g) &!is.na(adult.dat$bat_forearm_mm)] <- 10^(predict(modSlopeM) +1.96*predict(modSlopeM,se.fit=T)$se)


#plot for reality check

#ColM<- c("Pteropus rufus"="blue", "Eidolon dupreanum"="lightgreen","Rousettus madagascariensis"="purple")

# plot: bat mass:forearm residuals
# blue line is your predicted mass from the model. 
# black dots are your actual data

p1 <- ggplot(data=subset(adult.dat, !is.na(predicted_mass_g))) +
  geom_point(aes(x=bat_forearm_mm, y=bat_weight_g),size=1) +
  geom_ribbon(aes(x=bat_forearm_mm, ymin=predicted_mass_g_lci, ymax= predicted_mass_g_uci, group=bat_species),  fill="blue", alpha=.3) +
  geom_line(aes(x=bat_forearm_mm, y=predicted_mass_g, group=bat_species),  color="blue", size=1) +
  facet_grid(bat_sex~bat_species, scales = "free") +theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        strip.text.x  = element_text(face="italic", size=14),
        strip.text.y  = element_text(size=14),
        axis.title = element_text(size=14), axis.text = element_text(size=12),
        strip.background = element_rect(fill="white"), 
        legend.title = element_blank()) +
  xlab("bat forearm (mm)") + ylab("bat mass (g)")
p1

ggsave(file = paste0(homewd,"/final_figures_tables/supplemental/SuppFig8.pdf"),
       plot = p1,
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=600)
# 
# # Now calculate mass residual and save data
# adult.dat$mass_forearm_residual <- adult.dat$bat_weight_g - adult.dat$predicted_mass_g
# hist(adult.dat$mass_forearm_residual) #histogram centers around 0. perfect.
# 
# # Now attach these mass:forearm residuals to the original dataset
# 
# merge.mass_resid <- dplyr::select(adult.dat, sampleid, mass_forearm_residual)
# 
# 
# catch.dat <- merge(catch.dat, merge.mass_resid, by="sampleid", all.x = T)
# head(catch.dat)
# 
# sort(unique(catch.dat$mass_forearm_residual))
# unique(catch.dat$mass_forearm_residual[catch.dat$young_of_year=="yes"]) 
# #these 5 with values are double-labelled
# subset(catch.dat, young_of_year=="yes" & !is.na(mass_forearm_residual))
# unique(catch.dat$mass_forearm_residual[catch.dat$bat_age_class=="J"]) #same
# # remove
# 
# catch.dat$mass_forearm_residual[catch.dat$young_of_year=="yes"] <- NA
# length(catch.dat$mass_forearm_residual[is.na(catch.dat$mass_forearm_residual)]) #427 - juveniles or YOY 
# 
# #and save to the master data files 
# write.csv(catch.dat, file = paste0(homewd, "/master_data_files/catching_data.csv"), row.names = F)
# 
# 
