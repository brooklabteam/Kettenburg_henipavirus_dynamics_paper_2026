#Script to make figure 3, which shows the mountain plots of henipavirus serology and the seasonal serology

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
library(heatmaply)
library(seqinr)
library(scales)
library(ggdendroplot)
library(ggimage)
library(ggtree)
library(ape)
library(ggnewscale)
library(gdata)
library(phylotools)
library(phylobase)
library(cowplot)

#use this script to calculate cutoffs for the luminex


#first record your wd for this repo (home) 
homewd= "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/"
#homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/"

#now load in the merged luminex data from the "data" folder

#make sure this is latest catching file
dat <- read.csv(paste0(homewd, "serology/data/lum_catch_merge.csv"), header = T, stringsAsFactors = F)
head(dat)

length(unique(dat$sampleid))
dat$sampleid[duplicated(dat$sampleid)]

# test <- ddply(dat, .(bat_species, bat_sex), summarise, Ntot=length(unique(sampleid)))
# test

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

# # #convert any negative values to 1
# lumdat.long$MFI[lumdat.long$MFI<0] <- 1
# lumdat.long$log10MFI <- log10(lumdat.long$MFI)

lum.split <- dlply(lumdat.long,.(bat_species, antigen))

correct.negs <- function(df){
  abs.min <- abs(min(df$MFI))
  df$MFI <- df$MFI + abs.min + 1
  
  #then, also take the log10
  df$log10MFI <- log10(df$MFI)
  return(df)
}

lumdat.long <- data.table::rbindlist(lapply(lum.split, correct.negs))


#add family
lumdat.long$virus_family <- "Paramyxovirus"
unique(lumdat.long$virus_family)

#and genus information
lumdat.long$virus_genus <- "Henipavirus"
lumdat.long$virus_genus[lumdat.long$antigen=="MojV"] <- "Parahenipavirus"
unique(lumdat.long$virus_genus)

#factor antigens by relatedness
lumdat.long$antigen <- factor(lumdat.long$antigen, levels=c( "NiV", "HeV","CedV", "GhV", "AngV", "MojV"))

lumdat.long$bat_species <- factor(lumdat.long$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))
 
lumdat.long.plot <- lumdat.long

#and heatmap MFI by species/antigen combination 
library(ggh4x)
p1 <- ggplot(lumdat.long.plot) + 
  geom_tile(aes(x=sampleid, y=antigen, fill=MFI)) + theme_bw()+
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(size=12), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        legend.key.width =  unit(c(1.5), "cm"),
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        legend.position = "bottom") +
  scale_fill_viridis_c(trans="log10") + facet_nested(virus_family+virus_genus~bat_species, scales = "free", space = "free");p1

# and also add a mountain plot for additional support
# here, with antigens arranged according to phylogenetic relatedness... 
# makes it a bit easier to see the relatedness between bats
# first pass is all antigens in order of true relatedness

# load aa alignment
henipa.align <- read.alignment(file = paste0(homewd,"serology/data/Henipavirus-Alignment.fasta"), format = "fasta")
# replace names
new.names <- strsplit(henipa.align$nam,"_")
length.names <- lapply(new.names, length)
new.names <- mapply(FUN='[', new.names, length.names)
henipa.align$nam <- new.names

# make distance matrix
dist.henipa <- dist.alignment(henipa.align)

colorz = c('AngV' = "#F8766D", 'HeV' ="#A58AFF", 
           'NiV' = "#DB72FB", 'GhV' = "#00C094", 
           'CedV'= "#53B400", "MojV" = "#619CFF")


#and put them together
image.df2 <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"), antigen=rep("NiV", 3), MFI=rep(21000,3), image=c(paste0(homewd,"/serology/images/eidolon_dupreanum_headshot.png"),paste0(homewd,"/serology/images/rousettus_madagascariensis_headshot.png"),paste0(homewd,"/serology/images/pteropus_rufus_headshot.png")))
image.df2$bat_species <- factor(image.df2$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))
image.df2$virus_family="Paramyxovirus"
image.df <- image.df2

lumdat.long$bat_species <- factor(lumdat.long$bat_species,  levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))

#some cutoffs were established by mixture modeling
cutoff.df <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 4), rep("Rousettus madagascariensis",3)), antigen=c("NiV", "HeV", "GhV", "AngV", "GhV", "BdbV", "RavV"))
cutoff.df$virus_family = "Paramyxovirus"
cutoff.df$MFI = c(936,2851,1141,1763,2642,5647,6617)
mean.cut <- ddply(cutoff.df,.(bat_species, virus_family), summarise, meanMFI = mean(MFI))

plot.cut <- cbind.data.frame(bat_species=rep(c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"),each=2),
                             virus_family=rep(c("Filovirus", "Paramyxovirus"), 3),
                             MFI= c(rep(1672,2), 2642, 6132, rep(1672,2))) #might come back and mixture model all of these

plot.cut$bat_species <- factor(plot.cut$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))

hc.henipa <- hclust(dist.henipa)

#edit these a bit for ease of viewing
hc.henipa$height[2]<-hc.henipa$height[2]-.3
hc.henipa$height[3]<-hc.henipa$height[3]-.2
hc.henipa$height[4]<-hc.henipa$height[4]-.1

p.h.1 <- ggplot() + geom_dendro(hc.henipa, size=1) + theme_bw() +
      scale_y_reverse() + scale_x_reverse()+ 
      theme(panel.grid = element_blank(), 
            panel.border = element_blank(),
            axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), plot.margin = unit(c(0,1,0,3.3), "cm")) 
p.h.1
# p.h.2 <- ggplot() + geom_dendro(hc.henipa, size=1) + theme_bw() +
#   scale_y_reverse() + scale_x_reverse()+
#   theme(panel.grid = element_blank(), panel.border = element_blank(),
#         axis.title = element_blank(), axis.text = element_blank(),
#         axis.ticks = element_blank(), plot.margin = unit(c(0,1,0,1), "cm")) 
# p.h.3 <- ggplot() + geom_dendro(hc.henipa, size=1) + theme_bw() +
#   scale_y_reverse() + scale_x_reverse()+
#   theme(panel.grid = element_blank(), panel.border = element_blank(),
#         axis.title = element_blank(), axis.text = element_blank(),
#         axis.ticks = element_blank(), plot.margin = unit(c(0,1,0,1), "cm")) 

# p.h.row <- cowplot::plot_grid(p.h.1,p.h.2,p.h.3, ncol = 3, nrow = 1, rel_widths = c(1.22,1,1)) + theme(plot.margin = unit(c(0,0,0,0), "cm")) 

p2bsave <- ggplot(subset(lumdat.long)) + 
  #coord_cartesian(ylim = c(0,27000)) +
  facet_nested(bat_species~., scales = "free") +
  geom_line(aes(x=antigen, y=MFI, group=sampleid), show.legend=F) + theme_bw()+
  geom_point(aes(x=antigen, y=MFI, color=antigen), show.legend=F) + 
  geom_hline(data = subset(plot.cut, virus_family=="Paramyxovirus"), aes(yintercept = MFI), color="red", linetype=2, size=.8) +
  scale_color_manual(values=colorz) + ylab("Mean fluorescence intensity (MFI)") +
  #geom_dendro(rowclus) +
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size=14), 
        axis.title.y = element_text(size=16), 
        axis.text.y = element_text(size=16), 
        strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=16, face="italic",color = "black"), 
        strip.text.x = element_text(size=16, face="italic", color="black"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")); p2bsave
  #geom_image(data=subset(image.df, virus_family=="Paramyxovirus"), aes(x=antigen, y=c(23000, 24000,23000), image = image), size=.3, nudge_x = .7); p2bsave
mountain <- cowplot::plot_grid(p2bsave, p.h.1, ncol = 1,nrow = 2, rel_heights = c(1,.2), rel_widths = c(1,0.2))+
  theme(plot.background = element_rect(fill = "white", colour = NA));mountain
mountain<-as.ggplot(mountain)

##Now make all the seasonal serology plots

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(mgcv)
library(reshape2)
library(Hmisc)
library(epitools)

homewd="/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
batwd="/Users/carabrook/Developer/Madagascar-Bat-Data"
#homewd = "/Users/gwenddolenkettenburg/Desktop/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
#batwd = "/Users/gwenddolenkettenburg/Desktop/Developer/Madagascar-Bat-Data"

## Section 1. Prep the data

# First, load the luminex data after determining who is seropositive/
# seronegative in the other scripts
lum.dat <- read.csv(file = paste0(homewd, "/serology/data/mix_sub_data_75_90_96_long.csv"))
head(lum.dat)
# Pteropus plots are all cross-reactivity and only two seropositives so no point examining them statistically
# For Eid and Rou, the only antigens with known viruses from sequence data are GhV and AngV
# Rou has no seropositives for AngV, so we focus exclusively on GhV in Rou and ghV/AngV in Eid

lum.dat = subset(lum.dat, bat_species=="Eidolon dupreanum" & antigen=="GhV"|  bat_species=="Eidolon dupreanum" & antigen=="AngV" | bat_species=="Rousettus madagascariensis" & antigen=="GhV")

lum.dat <- dplyr::select(lum.dat, sampleid, antigen, MFI, log10MFI, seropos, seropos_lci, seropos_uci)
## merge the luminex data with the catching data

# link this to the catching data to plot by sampling session
catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))
#catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/Historic-Catching-Data-Through2020/catching_data_historic.csv"))

head(catch.dat)
names(catch.dat)

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)

#and merge all three datasets
merge.dat <- merge(lum.dat, catch.dat, by="sampleid", all.x = T)
head(merge.dat)


#get the midpoint processing date per sampling session
merge.dat$processing_date <- as.Date(merge.dat$processing_date)
sampling.date <- ddply(merge.dat, .(bat_species, sampling_session), summarise, min_date = min(processing_date), max_date = max(processing_date))
sampling.date$mid_date <- (sampling.date$max_date-sampling.date$min_date)/2 + sampling.date$min_date

#and slim to merge
sampling.merge <- dplyr::select(sampling.date, bat_species, sampling_session, mid_date)

merge.dat <- merge(merge.dat, sampling.merge, by = c("bat_species", "sampling_session"))
merge.dat$year <- year(merge.dat$processing_date)
merge.dat$month <- as.character(month(merge.dat$processing_date))
merge.dat$month[merge.dat$month!="10"& merge.dat$month!="11" &merge.dat$month!="12"]<-paste0("0", merge.dat$month[merge.dat$month!="10"& merge.dat$month!="11" &merge.dat$month!="12"])
#merge.dat$mid_date <- as.Date(paste0(merge.dat$year, "-", merge.dat$month, "-15"))
head(merge.dat)

#focusing just on 2018-2020 window
#merge.dat <- subset(merge.dat, mid_date>"2018-01-01")
merge.dat <- subset(merge.dat, mid_date>"2018-01-15")

#change roost site to just Moramanga region or Ankarana
unique(merge.dat$roost_site)
merge.dat$roost_site[merge.dat$roost_site=="Angavobe_Edup" |merge.dat$roost_site=="Angavokely_Edup"| merge.dat$roost_site=="Maromizaha_Rmad"] <- "Moramanga"
merge.dat$roost_site[merge.dat$roost_site!="Moramanga"] <- "Ankarana"

#also add in quantitative age data where available:
#age.dat <- read.csv(file = paste0(batwd, "/Age-Data/aged_bats.csv"), header = T, stringsAsFactors = F)
#head(age.dat)

#slim to just the columns you want to add
#age.dat <- dplyr::select(age.dat, sampleid, age)

#and merge
#merge.dat <- merge(merge.dat, age.dat, by="sampleid", all.x = T)
head(merge.dat)
#names(merge.dat)[names(merge.dat)=="mid_date"] <- "mid_date_old"
#names(merge.dat)[names(merge.dat)=="processing_date"] <- "mid_date"

# #curate date once more
sampling.date <- ddply(merge.dat, .(bat_species, roost_site, mid_date), summarise)
sampling.date# examine
# #write over those that are super close
merge.dat$mid_date <- as.character(merge.dat$mid_date)
merge.dat$mid_date[merge.dat$bat_species=="Eidolon dupreanum" & merge.dat$roost_site=="Moramanga" & merge.dat$mid_date=="2018-02-15"] <- "2018-02-14"
merge.dat$mid_date[merge.dat$bat_species=="Rousettus madagascariensis" & merge.dat$roost_site=="Moramanga" & merge.dat$mid_date=="2019-01-09"] <- "2019-01-05"
merge.dat$mid_date <- as.Date(merge.dat$mid_date)

## Seasonality of seroprevalence 

# now look at seroprevalence by sampling session by species and antigen
# plot data first
seroprev.sum <- ddply(merge.dat, .(bat_species, roost_site, antigen, mid_date), summarise, N=length(seropos), seropos=sum(seropos))
#seroprev.sum <- ddply(merge.dat, .(bat_species, roost_site, antigen, quarter_year), summarise, N=length(seropos), seropos=sum(seropos))
#names(seroprev.sum)[names(seroprev.sum)=="quarter_year"] <- "mid_date"
#seroprev.sum <- ddply(merge.dat, .(bat_species, roost_site, antigen, mid_date), summarise, N=length(seropos), seropos=sum(seropos))
seroprev.sum$seroneg <- seroprev.sum$N - seroprev.sum$seropos
seroprev.sum$seroprev <- seroprev.sum$seropos/seroprev.sum$N
seroprev.sum$seroprev_lci <-binom.exact(x=seroprev.sum$seropos, n=seroprev.sum$N)$lower
seroprev.sum$seroprev_uci <-binom.exact(x=seroprev.sum$seropos, n=seroprev.sum$N)$upper
head(seroprev.sum)

#and plot

# it's pretty evident there are correlations in seasonality 
# between these antigens

colorz = c('AngV' = "#F8766D", 'HeV' ="#A58AFF", 
           'NiV' = "#DB72FB", 'GhV' = "#00C094", 
           'CedV'= "#53B400", "MojV" = "#619CFF")

p1 <- ggplot(data=seroprev.sum) + geom_point(aes(x=mid_date, y=seroprev, size=N, color=antigen)) +
  geom_line(aes(x=mid_date, y=seroprev, color=antigen)) + scale_color_manual(values = colorz) +
  facet_grid(bat_species~roost_site) + theme_bw() 
p1

# time series is super short for Ankarana and does not span a whole year,
# so let's remove it from the data for seasonal analysis
# we can come back to it for age and day of year

# add in the timing of the reproductive season in Moramanga for both species
# from angelo and santino's paper
# eidolon, gestation is from: Aug-3 - Nov 15; lactation Nov 16 - Feb 2
# rousettus, gestation is from: Sept 11 - Dec 11; lactation Dec 12 - Feb 19
# dry season = May 15 - October 10

gest.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis", "Rousettus madagascariensis"), start = c(as.Date("2018-08-03"), as.Date("2019-08-03"), as.Date("2018-09-11"), as.Date("2019-09-11")), stop=c(as.Date("2018-11-15"), as.Date("2019-11-15"), as.Date("2018-12-11"), as.Date("2019-12-11")))
gest.dat$season = "gestation"
lact.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis"), start = c(as.Date("2018-01-01"), as.Date("2018-11-16"), as.Date("2019-11-16"), as.Date("2018-01-01"), as.Date("2018-12-12"), as.Date("2019-12-12")), stop=c(as.Date("2018-02-02"), as.Date("2019-02-02"), as.Date("2020-02-02"), as.Date("2018-02-19"), as.Date("2019-02-19"), as.Date("2020-02-19")))
lact.dat$season = "lactation"

#and winter
cal.dat <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis", "Rousettus madagascariensis"), start = c(as.Date("2018-05-15"), as.Date("2019-05-15"), as.Date("2018-05-15"), as.Date("2019-05-15")), stop=c(as.Date("2018-10-10"), as.Date("2019-10-10"), as.Date("2018-10-10"), as.Date("2019-10-10")))
cal.dat$season <- "winter"

cal.dat <- rbind(cal.dat, gest.dat, lact.dat)
cal.dat.save=cal.dat

#now add to plot

p2 <- ggplot(data=subset(seroprev.sum, roost_site=="Moramanga")) + 
  geom_rect(data=cal.dat, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=mid_date, y=seroprev, size=N, color=antigen)) +
  geom_line(aes(x=mid_date, y=seroprev, color=antigen)) + scale_color_manual(values = colorz) +
  facet_grid(bat_species~.) + theme_bw()

p2

#what if you did males and females separately?
cal.dat$bat_sex <- "female"
cal.dat$bat_sex[cal.dat$season=="winter"] <- "male"

seroprev.sum.sex <- ddply(merge.dat, .(bat_species, roost_site, bat_sex, antigen,  mid_date), summarise, N=length(seropos), seropos=sum(seropos))
#seroprev.sum.sex <- ddply(merge.dat, .(bat_species, roost_site, bat_sex, antigen,  quarter_year), summarise, N=length(seropos), seropos=sum(seropos))
#names(seroprev.sum.sex)[names(seroprev.sum.sex)=="quarter_year"] <- "mid_date"
#seroprev.sum.sex <- ddply(merge.dat, .(bat_species, roost_site, bat_sex, antigen,  mid_date), summarise, N=length(seropos), seropos=sum(seropos))
seroprev.sum.sex$seroneg <- seroprev.sum.sex$N - seroprev.sum.sex$seropos
seroprev.sum.sex$seroprev <- seroprev.sum.sex$seropos/seroprev.sum.sex$N
seroprev.sum.sex$seroprev_lci <-binom.exact(x=seroprev.sum.sex$seropos, n=seroprev.sum.sex$N)$lower
seroprev.sum.sex$seroprev_uci <-binom.exact(x=seroprev.sum.sex$seropos, n=seroprev.sum.sex$N)$upper

head(seroprev.sum.sex)
seroprev.sum.sex = subset(seroprev.sum.sex, bat_sex=="female" | bat_sex=="male")

p3 <- ggplot(data=subset(seroprev.sum.sex, roost_site=="Moramanga")) + 
  geom_rect(data=cal.dat, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=mid_date, y=seroprev, size=N, color=antigen)) + scale_color_manual(values = colorz) +
  geom_line(aes(x=mid_date, y=seroprev, color=antigen)) +
  facet_grid(bat_species~bat_sex) + theme_bw()

p3
# interesting differences here

#check out Ankarana
p3b <- ggplot(data=subset(seroprev.sum.sex, roost_site=="Ankarana")) + 
  geom_rect(data=cal.dat, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=mid_date, y=seroprev, size=N, color=antigen)) +
  geom_line(aes(x=mid_date, y=seroprev, color=antigen)) + scale_color_manual(values = colorz) +
  facet_grid(bat_species~bat_sex) + theme_bw()

p3b

# time series is too short here, so we will allow the model to zero out data 
# from the Ankarana site if it does not find it useful. 

# we will need to make a bonferonni correction on p-values
# since we are doing multiple hypothesis tests on these data
# 4 correlated antigens for the eidolon
.001/length(unique(seroprev.sum.sex$antigen[seroprev.sum.sex$bat_species=="Eidolon dupreanum"])) #0.00025 threshold for significance for Eidolon
# 3 for the rousettus
.001/length(unique(seroprev.sum.sex$antigen[seroprev.sum.sex$bat_species=="Rousettus madagascariensis"]))#0.00025 threshold for significance for Rousettus


# you want to test whether including sex or antigen matters in driving seroprevalence
seroprev.sum.sex$antigen <- as.factor(seroprev.sum.sex$antigen)
seroprev.sum.sex$bat_sex <- as.factor(seroprev.sum.sex$bat_sex)
seroprev.sum.sex$roost_site <- as.factor(seroprev.sum.sex$roost_site)

# prior work suggests important sex-based differences, so we will 
# model sex interactions using the "by" term and as a fixed effect
# we will also model all antigens together but allow for antigen-specific 
# effects of sex and site

# research question to be answered: is there sex-specific seasonality in
# the serological data for each antigen?

# for each species and antigen, we first test: 
# is there seasonality in seroprevalence by antigen without accounting for sex?
# we can include both sites (Ankarana and Moramanga) but allow the model to drop 
# terms internally if one dataset (e.g. Ankarana) is not robust enough.

seroprev.sum$year <- year(seroprev.sum$mid_date)
seroprev.sum$doy <- yday(seroprev.sum$mid_date)
seroprev.sum$antigen <- as.factor(seroprev.sum$antigen)
seroprev.sum$roost_site <- as.factor(seroprev.sum$roost_site)
seroprev.sum$year_f <- as.factor(seroprev.sum$year)


##Now get the interannual data

#GhV
Eid.sum.GhV <- subset(seroprev.sum, bat_species=="Eidolon dupreanum" & antigen=="GhV")
Eid.sum.GhV$month <- month(Eid.sum.GhV$mid_date)
Eid.sum.GhV$roost_site <- droplevels(Eid.sum.GhV$roost_site)
Eid.sum.GhV$antigen <- droplevels(Eid.sum.GhV$antigen)
Eid.sum.GhV$year_f <- droplevels(Eid.sum.GhV$year_f)

# 
# gamEidNoSexGhV <- gam(cbind(seropos, seroneg) ~ 
#                         roost_site + 
#                         s(as.numeric(doy), bs = "cc", k =10, by = roost_site) +
#                         s(year_f,  bs="fs", k=3, by= roost_site),
#                       data= Eid.sum.GhV, 
#                       family = "binomial", 
#                       method = "REML", 
#                       gamma = 1,
#                       select = TRUE) #


gamEidNoSexGhV <- gam(cbind(seropos, seroneg) ~ 
                        #roost_site + 
                        s(as.numeric(doy), bs = "cc", k =10, by = roost_site) +
                        s(as.numeric(mid_date), bs = "tp", k =6, by = roost_site),
                      data= Eid.sum.GhV, 
                      family = "binomial", 
                      method = "REML", 
                      gamma = 1,
                      select = TRUE) #



summary(gamEidNoSexGhV) 
#plot(gamEidNoSexGhV)
# sig seasonality inter- and intra-annual for GhV Moramanga
# Ankarana interannual edf = 0

# 
# #now try with sex too
# seroprev.sum.sex$doy <- yday(seroprev.sum.sex$mid_date)
# seroprev.sum.sex$year <- year(seroprev.sum.sex$mid_date)
# 
# Eid.sex.GhV = subset(seroprev.sum.sex, bat_species=="Eidolon dupreanum" & antigen=="GhV")
# Eid.sex.GhV$year_f <- as.factor(Eid.sex.GhV$year)
# Eid.sex.GhV$antigen <- droplevels(Eid.sex.GhV$antigen)
# Eid.sex.GhV$roost_site <- droplevels(Eid.sex.GhV$roost_site)
# Eid.sex.GhV$bat_sex <- droplevels(Eid.sex.GhV$bat_sex)
# 
# group <- interaction( Eid.sex.GhV$bat_sex, Eid.sex.GhV$roost_site, drop = TRUE)
# 
# Eid.sex.GhV$group <- group
# Eid.sex.GhV$year <- as.numeric(as.character(Eid.sex.GhV$year))
# 
# gamEidSexGhV <- gam(cbind(seropos, seroneg)~ 
#                       roost_site + bat_sex +
#                       s(as.numeric(doy), bs = c("cc"), k = 10, by = group) +
#                       s(year_f, bs="fs", k=3, by=group),
#                     data=Eid.sex.GhV, 
#                     family = "binomial", 
#                     method = "REML", 
#                     gamma =1,
#                     select = TRUE) #this term allows non-impactful terms to essentially be removed from model - check edf in summary to see if this happened to any
# summary(gamEidSexGhV)
# significant interannual trends for female Moramanga and significant intraanual trends for male Moramanga 
# most Ankarana  data reduced to edf = 0, especial interannual


#AIC(gamEidNoSexGhV, gamEidSexGhV) # a lot better with no sex

#now try for AngV as well
Eid.sum.AngV <- subset(seroprev.sum, bat_species=="Eidolon dupreanum" & antigen=="AngV")
# 
Eid.sum.AngV$year <- year(Eid.sum.AngV$mid_date)
Eid.sum.AngV$year_f <- as.factor(Eid.sum.AngV$year)
Eid.sum.AngV$month <- month(Eid.sum.AngV$mid_date)
Eid.sum.AngV$roost_site <- droplevels(Eid.sum.AngV$roost_site)
Eid.sum.AngV$antigen <- droplevels(Eid.sum.AngV$antigen)
Eid.sum.AngV$year_f <- droplevels(Eid.sum.AngV$year_f)

gamEidNoSexAngV <- gam(cbind(seropos, seroneg) ~ 
                       #  roost_site + 
                         s(as.numeric(doy), bs = "cc", k =10, by = roost_site) +
                         s(as.numeric(mid_date),  bs="tp", k=6, by= roost_site),
                       data= Eid.sum.AngV, 
                       family = "binomial", 
                       method = "REML", 
                       gamma = 1,
                       select = TRUE) #



# gamEidNoSexAngVPlot <- gam(cbind(seropos, seroneg) ~ 
#                              roost_site + 
#                              s(as.numeric(doy), bs = "cc", k =10, by = roost_site) +
#                              s(as.numeric(mid_date),  bs="tp", k=6, by= roost_site),
#                            data= Eid.sum.AngV, 
#                            family = "binomial", 
#                            method = "REML", 
#                            gamma = .48,
#                            select = TRUE) #

summary(gamEidNoSexAngV) #Moramanga interannual sig
# #plot(gamEidNoSexAngV)
# # Ankarana edf = 0
# 
# 
# #now try with sex too
# 
# seroprev.sum.sex$doy <- yday(seroprev.sum.sex$mid_date)
# seroprev.sum.sex$year <- year(seroprev.sum.sex$mid_date)
# 
# Eid.sex.AngV = subset(seroprev.sum.sex, bat_species=="Eidolon dupreanum" & antigen=="AngV")
# Eid.sex.AngV$antigen <- droplevels(Eid.sex.AngV$antigen)
# Eid.sex.AngV$roost_site <- droplevels(Eid.sex.AngV$roost_site)
# Eid.sex.AngV$bat_sex <- droplevels(Eid.sex.AngV$bat_sex)
# 
# group <- interaction( Eid.sex.AngV$bat_sex, Eid.sex.AngV$roost_site, drop = TRUE)
# 
# Eid.sex.AngV$group <- group
# Eid.sex.AngV$year_f <- as.factor(Eid.sex.AngV$year)
# Eid.sex.AngV$year <- as.numeric(as.character(Eid.sex.AngV$year))
# 
# gamEidSexAngV <- gam(cbind(seropos, seroneg)~ 
#                        roost_site + bat_sex +
#                        s(as.numeric(doy), bs = c("cc"), k = 10, by = group) +
#                        s(as.numeric(mid_date), bs="tp", k=6, by=group),
#                      data=Eid.sex.AngV, 
#                      family = "binomial", 
#                      method = "REML", 
#                      gamma =1,
#                      select = TRUE) #this term allows non-impactful terms to essentially be removed from model - check edf in summary to see if this happened to any
# summary(gamEidSexAngV)
# significant interannual trends for female Moramanga and significant intraanual trends for male Moramanga 
# most Ankarana  data reduced to edf = 0, especial interannual


#AIC(gamEidNoSexAngV, gamEidSexAngV) # better with no sex





#use above to predict and overlay on the data
#make dataset to capture predictions
pred.df <- cbind.data.frame(mid_date = rep(seq(as.Date("2018-01-01"), as.Date("2020-04-01"), by = 1), (2*length(unique(subset(seroprev.sum.sex, bat_species=="Eidolon dupreanum")$antigen)))), bat_sex=rep(c("female", "male"), each=length(seq(as.Date("2018-01-01"), as.Date("2020-04-01"), by = 1))*length(unique(subset(seroprev.sum.sex, bat_species=="Eidolon dupreanum")$antigen))), antigen= rep(rep(unique(subset(seroprev.sum.sex, bat_species=="Eidolon dupreanum")$antigen), each = length(seq(as.Date("2018-01-01"), as.Date("2020-04-01"), by = 1))),2))
head(pred.df)
pred.df$antigen <- droplevels(pred.df$antigen)
pred.df$doy <- yday(pred.df$mid_date)
pred.df$year <- year(pred.df$mid_date)
pred.df$year_f <- as.factor(pred.df$year)
pred.df$roost_site <- "Moramanga"

pred.df$bat_sex <- as.factor(pred.df$bat_sex)
pred.df$roost_site <- as.factor(pred.df$roost_site)

group.pred <- interaction(pred.df$bat_sex,pred.df$roost_site,  drop=TRUE)

pred.df$group <- group.pred


levels(pred.df$group)
#levels(Eid.sex.GhV$group)


# #GhV
# pred.df$predicted_seroprevalence <- NA
# pred.df$predicted_seroprevalence_lci <- NA
# pred.df$predicted_seroprevalence_uci <- NA
# #pred.df$predicted_seroprevalence <- predict.gam(gamEidSexGhV, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_seroprevalence[pred.df$antigen=="GhV"] <- predict.gam(gamEidSexGhV, newdata=subset(pred.df,antigen=="GhV"), type = "response", se.fit = T)$fit
# pred.df$predicted_seroprevalence_lci[pred.df$antigen=="GhV"] <-pred.df$predicted_seroprevalence[pred.df$antigen=="GhV"] - 1.96*predict.gam(gamEidSexGhV, newdata=subset(pred.df,antigen=="GhV"), type = "response", se.fit = T)$se
# pred.df$predicted_seroprevalence_uci[pred.df$antigen=="GhV"] <-pred.df$predicted_seroprevalence[pred.df$antigen=="GhV"] + 1.96*predict.gam(gamEidSexGhV, newdata=subset(pred.df,antigen=="GhV"), type = "response", se.fit = T)$se
# 
# 
# #AngV
# pred.df$predicted_seroprevalence[pred.df$antigen=="AngV"] <- predict.gam(gamEidSexAngV, newdata=subset(pred.df,antigen=="AngV"), type = "response", se.fit = T)$fit
# pred.df$predicted_seroprevalence_lci[pred.df$antigen=="AngV"] <-pred.df$predicted_seroprevalence[pred.df$antigen=="AngV"] - 1.96*predict.gam(gamEidSexAngV, newdata=subset(pred.df,antigen=="AngV"), type = "response", se.fit = T)$se
# pred.df$predicted_seroprevalence_uci[pred.df$antigen=="AngV"] <-pred.df$predicted_seroprevalence[pred.df$antigen=="AngV"] + 1.96*predict.gam(gamEidSexAngV, newdata=subset(pred.df,antigen=="AngV"), type = "response", se.fit = T)$se
# 
# 
# 
# 
# pred.df$bat_species = "Eidolon dupreanum"
# pred.df$predicted_seroprevalence_lci[pred.df$predicted_seroprevalence_lci <0] <- 0
# pred.df$predicted_seroprevalence_uci[pred.df$predicted_seroprevalence_uci >1] <- 1
# 
# 

# #alphaz = c('0' = 0, '1' = 1, '2' = 0.5)
# 
# p5 <- ggplot(data=subset(seroprev.sum.sex, roost_site=="Moramanga")) +
#   geom_rect(data=cal.dat, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3) +
#   geom_point(aes(x=mid_date, y=seroprev, size=N, color=antigen)) + scale_color_manual(values=colorz) + scale_fill_manual(values=filz) +
#   geom_linerange(aes(x=mid_date, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen), linetype=2, alpha=.3, show.legend = F) +
#   geom_ribbon(data=pred.df, aes(x=mid_date, ymin= predicted_seroprevalence_lci, ymax=predicted_seroprevalence_uci, fill=antigen), alpha=.3) +
#   geom_line(data=pred.df, aes(x=mid_date, y= predicted_seroprevalence, color=antigen)) +
#   facet_grid(bat_sex~bat_species) + theme_bw() #+ scale_alpha_discrete(values=alphaz)
# p5
# 
# #look at lag between females and males
# head(pred.df)

#predict males and females all together here
pred.df <- cbind.data.frame(mid_date = rep(seq(as.Date("2018-01-01"), as.Date("2020-04-01"), by = 1), (length(unique(subset(seroprev.sum.sex, bat_species=="Eidolon dupreanum")$antigen)))), antigen= rep(unique(subset(seroprev.sum.sex, bat_species=="Eidolon dupreanum")$antigen), each = length(seq(as.Date("2018-01-01"), as.Date("2020-04-01"), by = 1))))
head(pred.df)
pred.df$antigen <- droplevels(pred.df$antigen)
pred.df$doy <- yday(pred.df$mid_date)
pred.df$year <- year(pred.df$mid_date)

pred.df$month <- month(pred.df$mid_date)
pred.df$year_f <- as.factor(pred.df$year)
pred.df$roost_site <- "Moramanga"


pred.df$roost_site <- as.factor(pred.df$roost_site)



#GhV
pred.df$predicted_seroprevalence <- NA
pred.df$predicted_seroprevalence_lci <- NA
pred.df$predicted_seroprevalence_uci <- NA
pred.df$predicted_seroprevalence[pred.df$antigen=="GhV"] <- predict.gam(gamEidNoSexGhV, newdata=subset(pred.df,antigen=="GhV"), type = "response", se.fit = T)$fit
pred.df$predicted_seroprevalence_lci[pred.df$antigen=="GhV"] <-pred.df$predicted_seroprevalence[pred.df$antigen=="GhV"] - 1.96*predict.gam(gamEidNoSexGhV, newdata=subset(pred.df,antigen=="GhV"), type = "response", se.fit = T)$se
pred.df$predicted_seroprevalence_uci[pred.df$antigen=="GhV"] <-pred.df$predicted_seroprevalence[pred.df$antigen=="GhV"] + 1.96*predict.gam(gamEidNoSexGhV, newdata=subset(pred.df,antigen=="GhV"), type = "response", se.fit = T)$se


#AngV
pred.df$predicted_seroprevalence[pred.df$antigen=="AngV"] <- predict.gam(gamEidNoSexAngV, newdata=subset(pred.df,antigen=="AngV"), type = "response", se.fit = T)$fit
pred.df$predicted_seroprevalence_lci[pred.df$antigen=="AngV"] <-pred.df$predicted_seroprevalence[pred.df$antigen=="AngV"] - 1.96*predict.gam(gamEidNoSexAngV, newdata=subset(pred.df,antigen=="AngV"), type = "response", se.fit = T)$se
pred.df$predicted_seroprevalence_uci[pred.df$antigen=="AngV"] <-pred.df$predicted_seroprevalence[pred.df$antigen=="AngV"] + 1.96*predict.gam(gamEidNoSexAngV, newdata=subset(pred.df,antigen=="AngV"), type = "response", se.fit = T)$se





pred.df$bat_species = "Eidolon dupreanum"
pred.df$predicted_seroprevalence_lci[pred.df$predicted_seroprevalence_lci <0] <- 0
pred.df$predicted_seroprevalence_uci[pred.df$predicted_seroprevalence_uci >1] <- 1


#and plot with data:
pred.df$group <- paste0(pred.df$bat_species, "-", pred.df$antigen, "-", pred.df$year)

 filz <- c(colorz, 'gestation' = "yellow", 'lactation' = "lightpink", 'winter' = "skyblue")

p6 <- ggplot(data=subset(seroprev.sum, roost_site=="Moramanga")) +
  geom_rect(data=cal.dat, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=mid_date, y=seroprev, size=N, color=antigen)) + scale_color_manual(values=colorz) + scale_fill_manual(values=filz) +
  geom_linerange(aes(x=mid_date, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen), linetype=2, alpha=.3, show.legend = F) +
  geom_ribbon(data=pred.df, aes(x=mid_date, ymin= predicted_seroprevalence_lci, ymax=predicted_seroprevalence_uci, fill=antigen), alpha=.3) +
  geom_line(data=pred.df, aes(x=mid_date, y= predicted_seroprevalence, color=antigen, group=group)) +
  facet_grid(~bat_species) + theme_bw()
p6

#plot just GhV

p6b <- ggplot(data=subset(seroprev.sum, roost_site=="Moramanga" & antigen=="GhV" & bat_species =="Eidolon dupreanum")) +
  geom_rect(data=subset(cal.dat, bat_species =="Eidolon dupreanum"), aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=mid_date, y=seroprev, size=N, color=antigen)) + scale_color_manual(values=colorz) + scale_fill_manual(values=filz) +
  geom_linerange(aes(x=mid_date, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen), linetype=2, alpha=.3, show.legend = F) +
  geom_ribbon(data=subset(pred.df,  antigen=="GhV"& bat_species =="Eidolon dupreanum"), aes(x=mid_date, ymin= predicted_seroprevalence_lci, ymax=predicted_seroprevalence_uci, fill=antigen), alpha=.3) +
  geom_line(data=subset(pred.df, antigen=="GhV"& bat_species =="Eidolon dupreanum"), aes(x=mid_date, y= predicted_seroprevalence, color=antigen, group=year)) +
  facet_grid(~bat_species) + theme_bw()
p6b




# And do the same for Rousettus
# "GhV" only

Rou.sum.GhV <- subset(seroprev.sum, bat_species=="Rousettus madagascariensis" & antigen=="GhV")
Rou.sum.GhV$antigen <- droplevels(Rou.sum.GhV$antigen)
Rou.sum.GhV$roost_site <- droplevels(Rou.sum.GhV$roost_site)
Rou.sum.GhV$year_f <- as.factor(Rou.sum.GhV$year)



gamRouNoSexGhV <- gam(cbind(seropos, seroneg) ~# roost_site +
                        s(as.numeric(doy), bs = "cc", k =10, by = roost_site) +
                        s(as.numeric(mid_date),  bs="tp", k=6, by= roost_site),
                      data=Rou.sum.GhV, 
                      family = "binomial", 
                      method = "REML", 
                      gamma = 1,
                      select = TRUE) #this term allows non-impactful terms to essentially be removed from model - check edf in summary to see if this happened to any
summary(gamRouNoSexGhV)  
# no sig

# #and with sex
# Rou.sex.GhV = subset(seroprev.sum.sex, bat_species=="Rousettus madagascariensis" & antigen=="GhV" )
# Rou.sex.GhV$antigen <- droplevels(Rou.sex.GhV$antigen)
# Rou.sex.GhV$roost_site <- droplevels(Rou.sex.GhV$roost_site)
# Rou.sex.GhV$bat_sex <- droplevels(Rou.sex.GhV$bat_sex)
# Rou.sex.GhV$year_f <- as.factor(Rou.sex.GhV$year)
# 
# group4 <- interaction(Rou.sex.GhV$roost_site, Rou.sex.GhV$bat_sex)
# 
# gamRouSexGhV <- gam(cbind(seropos, seroneg)~ roost_site + bat_sex +
#                       s(as.numeric(doy),  k=10, bs="cc", by = group4) +
#                       s(year_f, k=3, bs="fs", by = group4), 
#                     data=Rou.sex.GhV, 
#                     family = "binomial", 
#                     method = "REML", 
#                     gamma =1,
#                     select = TRUE) #this term allows non-impactful terms to essentially be removed from model - check edf in summary to see if this happened to any
# summary(gamRouSexGhV) 
# # no sig
# 
# # compare model with and without sex
# AIC(gamRouNoSexGhV, gamRouSexGhV) # simpler model without sex is still more parsimonious


#use above to predict and overlay on the data
# #make dataset to capture predictions
# pred.df2 <- cbind.data.frame(mid_date = rep(seq(as.Date("2018-01-01"), as.Date("2020-04-01"), by = 1), (2*length(unique(subset(seroprev.sum.sex, bat_species=="Rousettus madagascariensis")$antigen)))), bat_sex=rep(c("female", "male"), each=length(seq(as.Date("2018-01-01"), as.Date("2020-04-01"), by = 1))*length(unique(subset(seroprev.sum.sex, bat_species=="Rousettus madagascariensis")$antigen))), antigen= rep(rep(unique(subset(seroprev.sum.sex, bat_species=="Rousettus madagascariensis")$antigen), each = length(seq(as.Date("2018-01-01"), as.Date("2020-04-01"), by = 1))),2))
# head(pred.df2)
# pred.df2$doy <- yday(pred.df2$mid_date)
# pred.df2$year <- year(pred.df2$mid_date)
# pred.df2$year_f <- as.factor(pred.df2$year)
# pred.df2$antigen <- droplevels(pred.df2$antigen) 
# pred.df2$roost_site <- "Moramanga"
# pred.df2$group4 <- interaction(pred.df2$roost_site,  pred.df2$bat_sex)
# 
# 
# 
# #GhV
# pred.df2$predicted_seroprevalence <- NA
# pred.df2$predicted_seroprevalence_lci <- NA
# pred.df2$predicted_seroprevalence_uci <- NA
# pred.df2$predicted_seroprevalence[pred.df2$antigen=="GhV"] <- predict.gam(gamRouSexGhV, newdata=subset(pred.df2,antigen=="GhV"), type = "response", se.fit = T)$fit
# pred.df2$predicted_seroprevalence_lci[pred.df2$antigen=="GhV"] <-pred.df2$predicted_seroprevalence[pred.df2$antigen=="GhV"] - 1.96*predict.gam(gamRouSexGhV, newdata=subset(pred.df2,antigen=="GhV"), type = "response", se.fit = T)$se
# pred.df2$predicted_seroprevalence_uci[pred.df2$antigen=="GhV"] <-pred.df2$predicted_seroprevalence[pred.df2$antigen=="GhV"] + 1.96*predict.gam(gamRouSexGhV, newdata=subset(pred.df2,antigen=="GhV"), type = "response", se.fit = T)$se
# 
# 
# 
# 
# pred.df2$bat_species = "Rousettus madagascariensis"
# pred.df2$predicted_seroprevalence_lci[pred.df2$predicted_seroprevalence_lci <0] <- 0
# pred.df2$predicted_seroprevalence_uci[pred.df2$predicted_seroprevalence_uci >1] <- 1
# 
# #then, mark those projections that were from significant model fits only for plotting
# pred.df2$sig <- 0
# pred.df2$sig[pred.df2$antigen=="RavV"] <- 1
# pred.df2$sig[pred.df2$antigen=="BdbV"] <- 1 # not sig here but is for the intra-annual
# #color significant those for which seroprevalence is estimated in Marana's paper
# 
# 
# #correct any predictions < 0:
# pred.df2$predicted_seroprevalence_lci[pred.df2$predicted_seroprevalence_lci<0] <- 0
# 

# #and plot with data:
# library(scales)
# setdiff(pred.df2$antigen, seroprev.sum.sex$antigen)
# setdiff(seroprev.sum.sex$antigen, pred.df2$antigen)

filz1 <- filz[1:(length(filz)-3)]
filz2 <- filz[(length(filz)-2):length(filz)]



#or, as a collective
#predict males and females all together here
pred.df2 <- cbind.data.frame(mid_date = rep(seq(as.Date("2018-01-01"), as.Date("2020-04-01"), by = 1), (length(unique(subset(seroprev.sum.sex, bat_species=="Rousettus madagascariensis")$antigen)))), antigen= rep(unique(subset(seroprev.sum.sex, bat_species=="Rousettus madagascariensis")$antigen), each = length(seq(as.Date("2018-01-01"), as.Date("2020-04-01"), by = 1))))
head(pred.df2)
pred.df2$antigen <- droplevels(pred.df2$antigen)
pred.df2$doy <- yday(pred.df2$mid_date)
pred.df2$year <- year(pred.df2$mid_date)
pred.df2$month <- month(pred.df2$mid_date)
pred.df2$year_f <- as.factor(pred.df2$year)
pred.df2$roost_site <- "Moramanga"


pred.df2$roost_site <- as.factor(pred.df2$roost_site)


#GhV
pred.df2$predicted_seroprevalence <- NA
pred.df2$predicted_seroprevalence_lci <- NA
pred.df2$predicted_seroprevalence_uci <- NA
pred.df2$predicted_seroprevalence[pred.df2$antigen=="GhV"] <- predict.gam(gamRouNoSexGhV, newdata=subset(pred.df2,antigen=="GhV"), type = "response", se.fit = T)$fit
pred.df2$predicted_seroprevalence_lci[pred.df2$antigen=="GhV"] <-pred.df2$predicted_seroprevalence[pred.df2$antigen=="GhV"] - 1.96*predict.gam(gamRouNoSexGhV, newdata=subset(pred.df2,antigen=="GhV"), type = "response", se.fit = T)$se
pred.df2$predicted_seroprevalence_uci[pred.df2$antigen=="GhV"] <-pred.df2$predicted_seroprevalence[pred.df2$antigen=="GhV"] + 1.96*predict.gam(gamRouNoSexGhV, newdata=subset(pred.df2,antigen=="GhV"), type = "response", se.fit = T)$se




pred.df2$bat_species = "Rousettus madagascariensis"
pred.df2$predicted_seroprevalence_lci[pred.df2$predicted_seroprevalence_lci <0] <- 0
pred.df2$predicted_seroprevalence_uci[pred.df2$predicted_seroprevalence_uci >1] <- 1
pred.df2$group <- paste0(pred.df2$bat_species, "-", pred.df2$antigen, "-", pred.df2$year)


pred.df <- rbind(pred.df, pred.df2)
pred.df$viralcat <- NA
pred.df$viralcat[pred.df$antigen == "HeV"] <- "HeV and NiV"
pred.df$viralcat[pred.df$antigen == "NiV"] <- "HeV and NiV"
pred.df$viralcat[pred.df$antigen == "AngV"] <- "AngV"
pred.df$viralcat[pred.df$antigen == "GhV"] <- "GhV"

dat_text <- data.frame(
  label = c("p = NS"),
  x     = c(Inf),
  y     = c(1)
)

#and plot with data:
pInterannual_new_Rou <- ggplot(data=subset(seroprev.sum, roost_site=="Moramanga" & bat_species=="Rousettus madagascariensis")) +
  geom_rect(data=subset(cal.dat,bat_species=="Rousettus madagascariensis"), aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=mid_date, y=seroprev, size=N, color=antigen)) + scale_color_manual(values=colorz) + scale_fill_manual(values=filz) +
  geom_linerange(aes(x=mid_date, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen), linetype=2, alpha=.3, show.legend = F) +
  geom_ribbon(data=subset(pred.df, bat_species=="Rousettus madagascariensis"), aes(x=mid_date, ymin= predicted_seroprevalence_lci, ymax=predicted_seroprevalence_uci, fill=antigen), alpha=.3) +
  geom_line(data=subset(pred.df, bat_species=="Rousettus madagascariensis"), aes(x=mid_date, y= predicted_seroprevalence, color=antigen, group = group)) +
  facet_grid(~bat_species) + theme_bw() +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  theme(axis.text.x = element_text(size=14), axis.title.x = element_blank(), 
        axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),
        strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        legend.position = "none",
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  scale_x_continuous(breaks=c(as.Date("2018-01-01"),as.Date("2019-01-01"), as.Date("2020-01-01")), labels=c("2018", "2019", "2020")) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + guides(alpha="none")
pInterannual_new_Rou
pInterannual_new_Rou2<-pInterannual_new_Rou+geom_text(data= dat_text, aes(x = x, y = y, label = label, hjust=1.1, vjust=1.5),
                                              fontface="italic");pInterannual_new_Rou2

#save for supplement



ggsave(file = paste0(homewd, "/final_figures_tables/tmp/Fig3_supplement_mountainplot_seasonal_serology_Rousettus.png"),
       plot= pInterannual_new_Rou2,
       units="mm",  
       width=60, 
       height=45, 
       scale=2, 
       dpi=500)




seroprev.sum$sig <- 2
seroprev.sum$sig[seroprev.sum$antigen=="AngV" | seroprev.sum$antigen=="GhV"] <- 1

seroprev.sum$sig <- as.factor(seroprev.sum$sig)

alphaz = c('0' = 0, '1' = 1, '2' = 0.5)

filz.season <- filz[(length(filz)-2):length(filz)]
filz.antigen <- filz[1:(length(filz)-2)]

# dat_text <- data.frame(
#   label = c("p = 0.05", "p = 0.03", "p = NS HeV\n p = NS NiV"),
#   viralcat   = c("AngV", "GhV", "HeV and NiV"),
#   x     = c(Inf, Inf, Inf),
#   y     = c(1, 1, 1)
# )

dat_text <- data.frame(
  label = c("p = 0.03", "p = 0.03"),
  viralcat   = c("AngV", "GhV"),
  x     = c(Inf, Inf),
  y     = c(1, 1)
)

seroprev.sum$viralcat <- NA
#seroprev.sum$viralcat[seroprev.sum$antigen == "HeV"] <- "HeV and NiV"
#seroprev.sum$viralcat[seroprev.sum$antigen == "NiV"] <- "HeV and NiV"
seroprev.sum$viralcat[seroprev.sum$antigen == "AngV"] <- "AngV"
seroprev.sum$viralcat[seroprev.sum$antigen == "GhV"] <- "GhV"

pInterannual_new <- ggplot(data=subset(seroprev.sum, roost_site=="Moramanga" & bat_species=="Eidolon dupreanum")) +
  geom_rect(data=subset(cal.dat,bat_species=="Eidolon dupreanum"), aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_linerange(aes(x=mid_date, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen, alpha=sig), linetype=2, show.legend = F) +
  geom_point(aes(x=mid_date, y=seroprev, size=N, color=antigen, alpha = sig)) + scale_color_manual(values=colorz) + scale_fill_manual(values=filz.season) +scale_alpha_manual(values=alphaz) +
  ggnewscale::new_scale_fill() +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  geom_ribbon(data=subset(pred.df, bat_species=="Eidolon dupreanum"), aes(x=mid_date, ymin= predicted_seroprevalence_lci, ymax=predicted_seroprevalence_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df, bat_species=="Eidolon dupreanum"), aes(x=mid_date, y= predicted_seroprevalence, color=antigen, group = group)) +
  facet_grid(viralcat~bat_species) + theme_bw() +
  theme(axis.text.x = element_text(size=14), 
        axis.title.x = element_blank(), 
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        # axis.title.y = element_blank(),
        axis.text.y = element_text(size=14),
        legend.position="none",
        legend.direction = "vertical",
        axis.title.y = element_text(size=16),
        plot.title = element_text(size=16),
        strip.background = element_rect(fill="white"), 
        strip.text.x = element_text(size=16, face="italic"), 
        strip.text.y = element_text(size=16),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")
        ) +
  #ggtitle("Interannual seroprevalence (2018-2020 dataset) in" ~ italic("E. dupreanum"))+
  scale_x_continuous(breaks=c(as.Date("2018-01-01"),as.Date("2019-01-01"), as.Date("2020-01-01")), labels=c("2018", "2019", "2020")) +
  coord_cartesian(expand = F) + ylab("Seroprevalence") + guides(alpha="none") + scale_fill_manual(values=filz.antigen) 
pInterannual_new

#legend<-get_legend(pInterannual_new)

pInterannual_new2<-pInterannual_new+geom_text(data= dat_text, aes(x = x, y = y, label = label, hjust=1.1, vjust=1.5),
                                              fontface="italic");pInterannual_new2

 
# now test the predictors of seropositivity from the raw data
# reload the serology and catching

lum.dat <- read.csv(file = paste0(homewd, "/serology/data/mix_sub_data_75_90_96_long.csv"))
head(lum.dat)
# Pteropus plots are all cross-reactivity and only two seropositives so no point examining them statistically
# For Eid and Rou, the only antigens with known viruses from sequence data are GhV and AngV
# Rou has no seropositives for AngV, so we focus exclusively on GhV in Rou and ghV/AngV in Eid

lum.dat = subset(lum.dat, bat_species=="Eidolon dupreanum" & antigen=="GhV"|  bat_species=="Eidolon dupreanum" & antigen=="AngV" | bat_species=="Rousettus madagascariensis" & antigen=="GhV")

lum.dat <- dplyr::select(lum.dat, sampleid, antigen, MFI, log10MFI, seropos, seropos_lci, seropos_uci)
## merge the luminex data with the catching data

# link this to the catching data to plot by sampling session
catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))
#catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/Historic-Catching-Data-Through2020/catching_data_historic.csv"))

head(catch.dat)
names(catch.dat)

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)

#and merge all three datasets
merge.dat <- merge(lum.dat, catch.dat, by="sampleid", all.x = T)
head(merge.dat)


#get the midpoint processing date per sampling session
merge.dat$processing_date <- as.Date(merge.dat$processing_date)
sampling.date <- ddply(merge.dat, .(bat_species, sampling_session), summarise, min_date = min(processing_date), max_date = max(processing_date))
sampling.date$mid_date <- (sampling.date$max_date-sampling.date$min_date)/2 + sampling.date$min_date

#and slim to merge
sampling.merge <- dplyr::select(sampling.date, bat_species, sampling_session, mid_date)

merge.dat <- merge(merge.dat, sampling.merge, by = c("bat_species", "sampling_session"))
merge.dat$year <- year(merge.dat$processing_date)
merge.dat$month <- as.character(month(merge.dat$processing_date))
merge.dat$month[merge.dat$month!="10"& merge.dat$month!="11" &merge.dat$month!="12"]<-paste0("0", merge.dat$month[merge.dat$month!="10"& merge.dat$month!="11" &merge.dat$month!="12"])
#merge.dat$mid_date <- as.Date(paste0(merge.dat$year, "-", merge.dat$month, "-15"))
head(merge.dat)

#focusing just on 2018-2020 window
#merge.dat <- subset(merge.dat, mid_date>"2018-01-01")
merge.dat <- subset(merge.dat, mid_date>"2018-01-15")

#change roost site to just Moramanga region or Ankarana
unique(merge.dat$roost_site)
merge.dat$roost_site[merge.dat$roost_site=="Angavobe_Edup" |merge.dat$roost_site=="Angavokely_Edup"| merge.dat$roost_site=="Maromizaha_Rmad"] <- "Moramanga"
merge.dat$roost_site[merge.dat$roost_site!="Moramanga"] <- "Ankarana"
head(merge.dat)

# #curate date once more
sampling.date <- ddply(merge.dat, .(bat_species, roost_site, mid_date), summarise)
sampling.date# examine
# #write over those that are super close
merge.dat$mid_date <- as.character(merge.dat$mid_date)
merge.dat$mid_date[merge.dat$bat_species=="Eidolon dupreanum" & merge.dat$roost_site=="Moramanga" & merge.dat$mid_date=="2018-02-15"] <- "2018-02-14"
merge.dat$mid_date[merge.dat$bat_species=="Rousettus madagascariensis" & merge.dat$roost_site=="Moramanga" & merge.dat$mid_date=="2019-01-09"] <- "2019-01-05"
merge.dat$mid_date <- as.Date(merge.dat$mid_date)

head(merge.dat)

#as with shedding data, we count everything as a yoy as a juvenile
unique(merge.dat$young_of_year)
unique(merge.dat$bat_age_class[merge.dat$young_of_year=="yes" & !is.na(merge.dat$young_of_year)]) #lots of NL and A in young of year
merge.dat$bat_age_class[merge.dat$young_of_year=="yes"& !is.na(merge.dat$young_of_year)]<- "J" #all yoy are juveniles
# 
unique(merge.dat$bat_age_class)

# make a few extra variables
merge.dat$combo_class <- paste0(merge.dat$bat_sex, "-", merge.dat$bat_age_class)
#combine juveniles 
merge.dat$combo_class[merge.dat$bat_age_class=="J"] <- "juvenile"
merge.dat$combo_class[merge.dat$combo_class=="male-NA" | merge.dat$combo_class=="female-NA" | merge.dat$combo_class=="NA-NA"] <- NA
unique(merge.dat$combo_class) # don't use NAs in this analysis
merge.dat$MFR_class <- merge.dat$combo_class
merge.dat$MFR_class[merge.dat$MFR_class=="female-L" | merge.dat$MFR_class=="male-A" | merge.dat$MFR_class=="female-NL"] <- "adult"
merge.dat$MFR_class <- as.factor(merge.dat$MFR_class)
merge.dat$combo_class <- as.factor(merge.dat$combo_class)
unique(merge.dat$MFR_class) # don't use NAs in this analysis
merge.dat$mass_forearm_residual[merge.dat$MFR_class=="juvenile"] <- 0
merge.dat$doy <- yday(merge.dat$mid_date)
merge.dat$year <- as.factor(year(merge.dat$processing_date))
merge.dat$roost_site <- as.factor(merge.dat$roost_site)

#Eidolon
#GhV
gamGhV <- gam(seropos ~ combo_class + mass_forearm_residual:MFR_class +
                s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                s(year, bs="re"),
              family = "binomial",
              method = "REML", 
              gamma = 1,
              select = TRUE,
              data= subset(merge.dat, bat_species =="Eidolon dupreanum" & antigen=="GhV"))

summary(gamGhV) # moramanga and year

gamAngV <- gam(seropos ~ combo_class + mass_forearm_residual:MFR_class +
                s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                s(year, bs="re"),
              family = "binomial",
              method = "REML", 
              gamma = 1,
              select = TRUE,
              data= subset(merge.dat, bat_species =="Eidolon dupreanum" & antigen=="AngV"))

summary(gamAngV) # nothing

#Rousettus 
gamRGhV <- gam(seropos ~ combo_class + mass_forearm_residual:MFR_class +
                s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                s(year, bs="re"),
              family = "binomial",
              method = "REML", 
              gamma = 1,
              select = TRUE,
              data= subset(merge.dat, bat_species =="Rousettus madagascariensis" & antigen=="GhV"))

summary(gamRGhV) # nothing sig

#because there was clear seasonality, at least for GhV in Eidolon, now sum across intraannual

#then, look at predictors of MFI among seropositives
merge.sum <- ddply(merge.dat,.(bat_species, antigen, roost_site, mid_date, year, doy), summarise, Nseropos=sum(seropos), Ntot=length(seropos))
merge.sum$seroneg = merge.sum$Ntot- merge.sum$Nseropos



#and new gams
gamGhV <- gam(cbind(Nseropos, seroneg) ~ 
                s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                s(year, bs="re"),
              family = "binomial",
              method = "REML", 
              gamma = 1,
              select = TRUE,
              data= subset(merge.sum, bat_species =="Eidolon dupreanum" & antigen=="GhV"))

summary(gamGhV) # moramanga and year


gamAngV <- gam(cbind(Nseropos, seroneg) ~ 
                s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                s(year, bs="re"),
              family = "binomial",
              method = "REML", 
              gamma = 1,
              select = TRUE,
              data= subset(merge.sum, bat_species =="Eidolon dupreanum" & antigen=="AngV"))

summary(gamAngV) # nothing


#Rousettus madagascariensis
gamRGhV <- gam(cbind(Nseropos, seroneg) ~ 
                s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                s(year, bs="re"),
              family = "binomial",
              method = "REML", 
              gamma = 1,
              select = TRUE,
              data= subset(merge.sum, bat_species =="Rousettus madagascariensis" & antigen=="GhV"))

summary(gamRGhV) # nothing


#only significant seasonality for GhV in Eidolon in Moramanga

#plot that next


pred.df <- cbind.data.frame(mid_date = seq(from= as.Date("2019-01-01"), to = as.Date("2019-12-31"), by=1), doy=1:365)
pred.df$bat_species <- "Eidolon dupreanum"
pred.df$antigen <- "GhV"
# pred.df2 <- pred.df
# pred.df2$antigen <- "AngV"
# pred.df <- rbind(pred.df, pred.df2)
pred.df$roost_site <- "Moramanga"
pred.df$year <- 2019

pred.df$prevalence <- NA
pred.df$prev_lci <- NA
pred.df$prev_uci <- NA



pred.df$prevalence[pred.df$antigen=="GhV"] <- predict(gamGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)"))
pred.df$prev_lci[pred.df$antigen=="GhV"] <- predict(gamGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)")) - 1.96*predict(gamGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="GhV"] <- predict(gamGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)")) + 1.96*predict(gamGhV, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(year)"), se.fit = TRUE)$se
# 
# pred.df$prevalence[pred.df$antigen=="AngV"] <- predict(gamAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(year)"))
# pred.df$prev_lci[pred.df$antigen=="AngV"] <- predict(gamAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(year)")) - 1.96*predict(gamAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(year)"), se.fit = TRUE)$se
# pred.df$prev_uci[pred.df$antigen=="AngV"] <- predict(gamAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(year)")) + 1.96*predict(gamAngV, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(year)"), se.fit = TRUE)$se

pred.df$prev_lci[pred.df$prev_lci<0] <- 0
pred.df$prev_uci[pred.df$prev_uci>1] <- 1
pred.df$sig <- 1


#and make the gestation database
gest.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis"), start = c(as.Date("2019-08-03"), as.Date("2019-09-11")), stop=c(as.Date("2019-11-15"),  as.Date("2019-12-11")))
gest.dat$season = "gestation"
lact.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis",  "Rousettus madagascariensis"), start = c(as.Date("2019-11-16"), as.Date("2019-01-01"), as.Date("2019-12-12"),  as.Date("2019-01-01")), stop=c(as.Date("2019-12-31"), as.Date("2019-02-02"), as.Date("2019-12-31"), as.Date("2015-02-19")))
lact.dat$season = "lactation"

#and winter
cal.dat <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis"), start = c( as.Date("2019-05-15"),as.Date("2019-05-15")), stop=c(as.Date("2019-10-10"),as.Date("2019-10-10")))
cal.dat$season <- "winter"
cal.dat <- rbind(cal.dat, gest.dat, lact.dat)
cal.dat$year <- year(cal.dat$stop)

cal.dat$doy_start <- yday(cal.dat$start)
cal.dat$doy_stop <- yday(cal.dat$stop)




#collect data on prevalence just for plotting
head(merge.sum)

merge.sum$seroprevalence_lci <- merge.sum$seroprevalence_uci <- NA

for (i in 1:length(merge.sum$bat_species)){
  merge.sum$seroprevalence_lci[i] <-  binom.test(merge.sum$Nseropos[i],merge.sum$Ntot[i])$conf.int[1]
  merge.sum$seroprevalence_uci[i] <-  binom.test(merge.sum$Nseropos[i],merge.sum$Ntot[i])$conf.int[2]
  
}

merge.sum$seroprevalence <- merge.sum$Nseropos/merge.sum$Ntot


merge.sum$seroprevalence_lcimerge.sum$sig = 1




colorz = c('AngV' = "#F8766D", 'GhV' = "#00C094")

sitez = c('Ankarana' ="#A58AFF",  'Moramanga' = "#00C094", 'Mahabo' = "#619CFF")

filz <- c(colorz, 'gestation' = "yellow", 'lactation' = "lightpink", 'winter' = "skyblue")

filz1 <- filz[1:(length(filz)-3)]
filz2 <- filz[(length(filz)-2):length(filz)]

yearz <- c('2013'=3, '2014' = 7, '2015'=8, '2016'=18, '2018' = 16, '2019' = 17, '2020' = 15, '2021' = 4, '2023' = 6)


pIntra <- ggplot(pred.df) + 
  geom_point(data=subset(merge.sum, bat_species=="Eidolon dupreanum" & antigen=="GhV"), aes(x=doy, y=seroprevalence, size=Ntot, color=antigen, shape=year)) +
  geom_linerange(data=subset(merge.sum, bat_species=="Eidolon dupreanum" & antigen=="GhV"), aes(x=doy, ymin=seroprevalence_lci, ymax=seroprevalence_uci, color=antigen,  group=year), linetype=2, show.legend = F) +
  geom_rect(data=subset(cal.dat, bat_species=="Eidolon dupreanum"), aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3) +
  scale_fill_manual(values=filz2) +
  scale_shape_manual(values=yearz) +
  ggnewscale::new_scale_fill() +
  geom_ribbon(data = subset(pred.df), aes(x=doy, ymin=prev_lci, ymax= prev_uci, fill = antigen), alpha=.3, show.legend = F) +
  geom_line(data = subset(pred.df),aes(x=doy, y=prevalence, color=antigen)) +
  # geom_ribbon(aes(x=doy, ymin=prev_lci, ymax= prev_uci, fill = antigen), alpha=.3, show.legend = F) +
  # geom_line(aes(x=doy, y=prevalence, color=antigen)) + 
  facet_grid(bat_species~antigen) + theme_bw() +
  scale_color_manual(values=colorz) + scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  scale_fill_manual(values=filz1) +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  theme(axis.text = element_text(size=14), #axis.title.x = element_blank(), 
        axis.title = element_text(size=16),
        strip.background.x = element_rect(fill="white"), 
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        #strip.text.y = element_text(size=16),panel.spacing = unit(1.5, "lines"),
        strip.text.x = element_text(size=16), 
        panel.grid = element_blank()) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + xlab("month of year") + guides(color="none") ;pIntra


legend <- cowplot::get_legend(pIntra)

pIntra_noLeg <- ggplot(pred.df) + 
  geom_point(data=subset(merge.sum, bat_species=="Eidolon dupreanum" & antigen=="GhV"), aes(x=doy, y=seroprevalence, size=Ntot, color=antigen, shape=year), show.legend = F) +
  geom_linerange(data=subset(merge.sum, bat_species=="Eidolon dupreanum" & antigen=="GhV"), aes(x=doy, ymin=seroprevalence_lci, ymax=seroprevalence_uci, color=antigen,  group=year), linetype=2, show.legend = F) +
  geom_rect(data=subset(cal.dat, bat_species=="Eidolon dupreanum"), aes(xmin=doy_start, xmax=doy_stop, ymin=0, ymax=1, fill=season), alpha=.3, show.legend = F) +
  scale_fill_manual(values=filz2) +
  scale_shape_manual(values=yearz) +
  ggnewscale::new_scale_fill() +
  geom_ribbon(data = subset(pred.df), aes(x=doy, ymin=prev_lci, ymax= prev_uci, fill = antigen), alpha=.3, show.legend = F) +
  geom_line(data = subset(pred.df),aes(x=doy, y=prevalence, color=antigen)) +
  # geom_ribbon(aes(x=doy, ymin=prev_lci, ymax= prev_uci, fill = antigen), alpha=.3, show.legend = F) +
  # geom_line(aes(x=doy, y=prevalence, color=antigen)) + 
  facet_grid(bat_species~antigen) + theme_bw() +
  scale_color_manual(values=colorz) + scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  scale_fill_manual(values=filz1) +
  scale_size_continuous(name=c("N bats"), breaks = c(10,20,30)) +
  theme(axis.text = element_text(size=14), #axis.title.x = element_blank(), 
        axis.title = element_text(size=16),
        strip.background.x = element_rect(fill="white"), 
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        #strip.text.y = element_text(size=16),panel.spacing = unit(1.5, "lines"),
        strip.text.x = element_text(size=16), 
        panel.grid = element_blank()) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + xlab("month of year") + guides(color="none") ;pIntra_noLeg


#Next, test predictors of MFI in seropositives  and plot that for GhV in Eid
merge.dat$combo_class <- factor(merge.dat$combo_class, levels=c("female-NL", "female-P", "female-L", "male-A", "juvenile"))

#try the same structure as for Ct data
merge.dat$mass_cat <- NA
merge.dat$mass_cat[merge.dat$mass_forearm_residual<0 & !is.na(merge.dat$mass_forearm_residual)] <- "poor"
merge.dat$mass_cat[merge.dat$mass_forearm_residual>0 & !is.na(merge.dat$mass_forearm_residual)] <- "good"
merge.dat$mass_cat[merge.dat$combo_class=="juvenile"] <- "juvenile"
#merge.dat$mass_cat[merge.dat$combo_class=="subadult"] <- "subadult"
merge.dat$mass_cat[merge.dat$combo_class=="female-P"] <- "pregnant"
merge.dat$mass_cat <- as.factor(merge.dat$mass_cat)

library(lme4)
library(lmerTest)
#Eidolon
#GhV

lmGhV <- lmer(MFI ~ 
                mass_cat + (1| year),
              REML=T,
              data=  subset(merge.dat, bat_species =="Eidolon dupreanum" & antigen=="GhV" & seropos==1)) #singular fit, no need for random effect

lmGhV <- lm(MFI ~ mass_cat,
              data=  subset(merge.dat, bat_species =="Eidolon dupreanum" & antigen=="GhV" & seropos==1)) #singular fit, no need for random effect

summary(lmGhV) #juveniles and poor body condition adults associated with high MFI values as compared with healthy adults
# 

#angv
lmAngV <- lmer(MFI ~ 
                mass_cat + (1| year),
              REML=T,
              data=  subset(merge.dat, bat_species =="Eidolon dupreanum" & antigen=="AngV" & seropos==1)) #singular fit, no need for random effect


summary(lmAngV) #juveniles associated with high MFI values as compared with healthy adults

# Rousettus
#GhV

lmRGhV <- lmer(MFI ~ 
                mass_cat + (1| year),
              REML=T,
              data=  subset(merge.dat, bat_species =="Rousettus madagascariensis" & antigen=="GhV" & seropos==1)) #singular fit, no need for random effect

lmRGhV <- lm(MFI ~ mass_cat,
            data=  subset(merge.dat, bat_species =="Rousettus madagascariensis" & antigen=="GhV" & seropos==1)) #singular fit, no need for random effect

summary(lmRGhV) #not sig
# 

# 
# #Eidolon
# #GhV
# gamMFIGhV <- gam(MFI ~ combo_class + mass_forearm_residual:MFR_class +
#                 s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
#                 s(year, bs="re"),
#               method = "REML", 
#               gamma = 1,
#               select = TRUE,
#               data= subset(merge.dat, bat_species =="Eidolon dupreanum" & antigen=="GhV" & seropos==1 & combo_class!="female-P"))
# 
# summary(gamMFIGhV) # moramanga,  juveniles and lactating females associated with higher MFI as compared with female-NL
# 
# 
# gamMFIAngV <- gam(MFI ~ combo_class + mass_forearm_residual:MFR_class +
#                  s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
#                  s(year, bs="re"),
#                method = "REML", 
#                gamma = 1,
#                select = TRUE,
#                data= subset(merge.dat, bat_species =="Eidolon dupreanum" & antigen=="AngV" & combo_class!="female-P"))
# 
# summary(gamMFIAngV) # adult males associated with lower MFI as compared with female-NL

# #Rousettus 
# gamRMFIGhV <- gam(MFI ~ combo_class + mass_forearm_residual:MFR_class +
#                  s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
#                  s(year, bs="re"),
#                method = "REML", 
#                gamma = 1,
#                select = TRUE,
#                data= subset(merge.dat, bat_species =="Rousettus madagascariensis" & antigen=="GhV" & combo_class!="female-P"))
# 
# summary(gamRMFIGhV) # nothing sig
# 
# #and try the same analysis as with the shedding data
# merge.dat$mass_cat <- NA
# merge.dat$mass_cat[merge.dat$mass_forearm_residual<0 &!is.na(merge.dat$mass_forearm_residual) & merge.dat$combo_class!="juvenile"|merge.dat$mass_forearm_residual<0 &!is.na(merge.dat$mass_forearm_residual) & merge.dat$combo_class!="female-P"] <- "poor"
# merge.dat$mass_cat[merge.dat$mass_forearm_residual>0 &!is.na(merge.dat$mass_forearm_residual) & merge.dat$combo_class!="juvenile"|merge.dat$mass_forearm_residual>0 &!is.na(merge.dat$mass_forearm_residual) & merge.dat$combo_class!="female-P"] <- "good"
# merge.dat$mass_cat[ merge.dat$combo_class=="juvenile"] <- "juvenile"
# merge.dat$mass_cat[ merge.dat$combo_class=="female-P"] <- "female-P"
# merge.dat$mass_cat <- as.factor(merge.dat$mass_cat)
# 
# 
# #and redo these with simler predictors
# #GhV
# gamMFIGhV <- lmer(MFI ~ mass_cat + (1|year), REML = T,
#                  data= subset(merge.dat, bat_species =="Eidolon dupreanum" & antigen=="GhV" & seropos==1 & mass_cat!="female-P")) # singlular fit, remove random effect
# 
# gamMFIGhV <- lm(MFI ~ mass_cat,
#                   data= subset(merge.dat, bat_species =="Eidolon dupreanum" & antigen=="GhV" & seropos==1 & mass_cat!="female-P")) 
# 
# summary(gamMFIGhV) #juveniles and poor body condition associated with higher MFI vs good body condition
# 
# 
# gamMFIAngV <- lmer(MFI ~ mass_cat + (1|year), REML = T,
#                   data= subset(merge.dat, bat_species =="Eidolon dupreanum" & antigen=="AngV" & mass_cat!="female-P"))# singlular fit, remove random effect
# 
# gamMFIAngV <- lm(MFI ~ mass_cat, 
#                    data= subset(merge.dat, bat_species =="Eidolon dupreanum" & antigen=="AngV" & mass_cat!="female-P"))
# 
# summary(gamMFIAngV) # nothing
# 
# #Rousettus 
# gamRMFIGhV <- lmer(MFI ~ mass_cat + (1|year), REML = T,
#                   data= subset(merge.dat, bat_species =="Rousettus madagascariensis" & antigen=="GhV" & mass_cat!="female-P"))# singlular fit, remove random effect
# 
# gamRMFIGhV <- lm(MFI ~ mass_cat,
#                    data= subset(merge.dat, bat_species =="Rousettus madagascariensis" & antigen=="GhV" & mass_cat!="female-P"))
# 
# summary(gamRMFIGhV) # nothing sig


#and plot
merge.dat$cat_plot <- as.character(merge.dat$mass_cat)
merge.dat$cat_plot[merge.dat$cat_plot=="good"] <- "adult,\nMFR>0"
merge.dat$cat_plot[merge.dat$cat_plot=="poor"] <- "adult,\nMFR<0"
merge.dat$cat_plot[merge.dat$cat_plot=="pregnant"] <- "pregnant\nfemale"
merge.dat$cat_plot <- factor(merge.dat$cat_plot, levels = c("adult,\nMFR>0", "adult,\nMFR<0", "pregnant\nfemale", "juvenile"))

#and plot MFI by condition/juvenile status - leave out pregnant females since there are so few
pMFR_MFI_GhV <- ggplot(data= subset(merge.dat, bat_species =="Eidolon dupreanum" & antigen=="GhV" & seropos==1 & mass_cat!="female-P")) +
  geom_violin(aes(x=cat_plot, y=MFI), draw_quantiles = c(0,.25,.5,.75)) + 
  geom_jitter(aes(x=cat_plot, y=MFI, color=roost_site, shape=as.factor(year)), width = .2, size=3) +
  theme_bw() + ylab("MFI") + facet_grid(bat_species~antigen)+
  scale_shape_manual(name="year", values=yearz) +
  scale_color_manual(values = sitez, name="roost site") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=14), 
        strip.background.y = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.23,.85),
        strip.text.y = element_blank(),
        strip.background.x = element_rect(fill="white"),
        strip.text.x = element_text(size = 16),
        axis.title.x = element_blank(), 
        #axis.title.x = element_text(size=16, colour = "transparent"),
        axis.title.y = element_text(size=16)) + guides(shape="none");pMFR_MFI_GhV 

#now put these all together
pCD <- cowplot::plot_grid(pIntra_noLeg, pMFR_MFI_GhV, nrow = 1, ncol = 2, labels = c("C", "D"), label_size = 22);pCD


#and complete the rest of fig3 

# inter <- cowplot::plot_grid(pInterannual_new2, NULL, ncol = 1, nrow = 2, 
#                             rel_widths = c(1,1), rel_heights = c(1,0.2), labels = c("B"), label_size = 22);inter
# as.ggplot(inter)

inter <- cowplot::plot_grid(pInterannual_new2, 
                             labels = c("B"), label_size = 22);inter
as.ggplot(inter)

# inter <- cowplot::plot_grid(pInterannual_new2, legend, ncol = 1, nrow = 2, 
#                             rel_widths = c(1,1), rel_heights = c(1,0.2), labels = c("", ""));inter
# as.ggplot(inter)

#right column
pRight <- cowplot::plot_grid(inter, pCD, ncol = 1, nrow = 2, 
                             rel_widths = c(1,1), rel_heights = c(1.8,1), labels = c("", ""));pRight

#add overall legend
pRightwLeg <- cowplot::plot_grid(pRight, legend, ncol = 2, nrow = 1, 
                             rel_widths = c(1,.2), labels = c("", ""));pRightwLeg


#and the whole plot

fig3<-cowplot::plot_grid(mountain, pRightwLeg,
                         labels=c("A",""),
                         rel_widths=c(1,1.6), rel_heights = c(1,1),
                         ncol=2,align="hv", axis="l", label_size = 23)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        plot.background = element_rect(fill = "white", colour = NA)); fig3


ggsave(file = paste0(homewd, "/final_figures_tables/tmp/Fig3_mountainplot_seasonal_serology_new.png"),
       plot= fig3,
       units="mm",  
       width=180, 
       height=140, 
       scale=2, 
       dpi=500)


##Add the historical IPM time series for comparison
 
 # load the IPM data and the JAE data
 
IPM.dat <- read.csv(file = paste0(batwd,"/Serological-Data/IPM-data-2006-2007.csv"),header = T, stringsAsFactors = F)
head(IPM.dat) 

# and do the same male/female specific serology for both for the antigens of 
# interest: HeV, NiV 

unique(IPM.dat$prev)
IPM.dat = subset(IPM.dat, prev=="0" | prev=="1")
unique(IPM.dat$species)
unique(IPM.dat$type)
IPM.dat = subset(IPM.dat, type=="ELISA_HeV" | type =="ELISA_NiV")

sort(unique(IPM.dat$date)) #each of these is a discrete sampling event
unique(IPM.dat$site)#"AngavoKely" "AngavoBe"  
subset(IPM.dat, date=="2007-04-05" ) #angavokely
subset(IPM.dat, date=="2007-04-11" ) #angavobe

#summarise to seropositivity
IPM.dat$prev <- as.numeric(IPM.dat$prev)
IPM.sero <- ddply(IPM.dat,.(type, date, sex),summarise, N_seropos=sum(prev), N_tot=length(prev))
IPM.sero.nosex <- ddply(IPM.dat,.(type, date),summarise, N_seropos=sum(prev), N_tot=length(prev))
IPM.sero$N_seroneg = IPM.sero$N_tot-IPM.sero$N_seropos
IPM.sero$seroprev <- IPM.sero$N_seropos/IPM.sero$N_tot
IPM.sero$seroprev_lci <-binom.exact(x=IPM.sero$N_seropos, n=IPM.sero$N_tot)$lower
IPM.sero$seroprev_uci <-binom.exact(x=IPM.sero$N_seropos, n=IPM.sero$N_tot)$upper
IPM.sero$date <- as.Date(IPM.sero$date)
head(IPM.sero)
unique(IPM.sero$type)
IPM.sero$type[IPM.sero$type=="ELISA_HeV"] <- "HeV"
IPM.sero$type[IPM.sero$type=="ELISA_NiV"] <- "NiV"
names(IPM.sero)[names(IPM.sero)=="type"] <- "antigen"
IPM.sero$antigen <- as.factor(IPM.sero$antigen)
unique(IPM.sero$sex)

IPM.sero$sex[IPM.sero$sex=="F"] <- "female"
IPM.sero$sex[IPM.sero$sex=="M"] <- "male"

IPM.sero$sex <- as.factor(IPM.sero$sex)


IPM.sero.nosex$N_seroneg = IPM.sero.nosex$N_tot-IPM.sero.nosex$N_seropos
IPM.sero.nosex$seroprev <- IPM.sero.nosex$N_seropos/IPM.sero.nosex$N_tot
IPM.sero.nosex$seroprev_lci <-binom.exact(x=IPM.sero.nosex$N_seropos, n=IPM.sero.nosex$N_tot)$lower
IPM.sero.nosex$seroprev_uci <-binom.exact(x=IPM.sero.nosex$N_seropos, n=IPM.sero.nosex$N_tot)$upper
IPM.sero.nosex$date <- as.Date(IPM.sero.nosex$date)
head(IPM.sero.nosex)
unique(IPM.sero.nosex$type)
IPM.sero.nosex$type[IPM.sero.nosex$type=="ELISA_HeV"] <- "HeV"
IPM.sero.nosex$type[IPM.sero.nosex$type=="ELISA_NiV"] <- "NiV"
names(IPM.sero.nosex)[names(IPM.sero.nosex)=="type"] <- "antigen"
IPM.sero.nosex$antigen <- as.factor(IPM.sero.nosex$antigen)

IPM.sero.nosex$doy <- yday(IPM.sero.nosex$date)
IPM.sero.nosex$year <- year(IPM.sero.nosex$date)
IPM.sero$doy <- yday(IPM.sero$date)
IPM.sero$year <- year(IPM.sero$date)

#because there is only one full calendar year, we cannot plot "cc" for doy
# so we just do the "tp" on date


# first, is there seasoyear()# first, is there seasonality without accounting for sex?
gamIPMNiVNoSex <- gam(cbind(N_seropos, N_seroneg)~ 
                      s(as.numeric(date), k=-1, bs="tp"), 
                      data=subset(IPM.sero.nosex, antigen=="NiV"), 
                      family = "binomial", 
                      select = TRUE, 
                      method = "REML")
summary(gamIPMNiVNoSex) # yes, very significant for both antigens and both timescales

gamIPMHeVNoSex <- gam(cbind(N_seropos, N_seroneg)~ 
                        s(as.numeric(date), k=-1, bs="tp"), 
                      data=subset(IPM.sero.nosex, antigen=="HeV"), 
                      family = "binomial", 
                      select = TRUE, 
                      method = "REML")
summary(gamIPMHeVNoSex) # yes, very significant for both antigens and both timescales


# now, with sex
gamIPMNiVSex <- gam(cbind(N_seropos, N_seroneg)~  sex +
                    s(as.numeric(date), k=-1, bs="tp", by=interaction(antigen, sex)),
                   data=subset(IPM.sero, antigen=="NiV"), family = "binomial", select = TRUE,
                  method = "REML")
summary(gamIPMNiVSex) # yes, very significant for both antigens and both timescales
#plot(gamIPMSex)
#intraanual for HeV male, NiV male, borderline NiV female
# interannual for HeV male and NiV male

gamIPMHeVSex <- gam(cbind(N_seropos, N_seroneg)~  sex +
                      s(as.numeric(date), k=-1, bs="tp", by=interaction(antigen, sex)),
                    data=subset(IPM.sero, antigen=="HeV"), family = "binomial", select = TRUE,
                    method = "REML")
summary(gamIPMHeVSex) # yes, very significant for both antigens and both timescales
#plot(gamIPMSex)
#intraanual for HeV male, NiV male, borderline NiV female
# interannual for HeV male and NiV male

AIC(gamIPMNiVNoSex, gamIPMNiVSex) # simpler, no sex model prevails
AIC(gamIPMHeVNoSex, gamIPMHeVSex) # simpler, no sex model prevails


# #and plot the predictions
pred.df.IPM <- cbind.data.frame(date = rep(seq(as.Date("2005-11-01"), as.Date("2007-04-30"), by = 1), (2*length(unique(IPM.sero$antigen)))), sex=rep(c("female", "male"), each=length(seq(as.Date("2005-11-01"), as.Date("2007-04-30"), by = 1))*length(unique(IPM.sero$antigen))), antigen= rep(rep(unique(IPM.sero$antigen), each = length(seq(as.Date("2005-11-01"), as.Date("2007-04-30"), by = 1))),2))
pred.df.IPM$sex = "female"
pred.df.tmp <- pred.df.IPM
pred.df.tmp$sex <- "male"
pred.df.IPM <- rbind(pred.df.IPM, pred.df.tmp)
#pred.df.IPM$doy <- yday(pred.df.IPM$date)
#pred.df.IPM$year <- year(pred.df.IPM$date)

head(pred.df.IPM)

#pred.df.IPM = subset(pred.df.IPM, antigen=="NiV")

#predict without effect of year for plotting
#NiV
pred.df.IPM$predicted_seropos <- NA
pred.df.IPM$predicted_seropos_lci <- NA
pred.df.IPM$predicted_seropos_uci <- NA
pred.df.IPM$predicted_seropos[pred.df.IPM$antigen=="NiV"] <- predict.gam(gamIPMNiVNoSex, newdata = pred.df.IPM[pred.df.IPM$antigen=="NiV",], type="response")
pred.df.IPM$predicted_seropos_lci[pred.df.IPM$antigen=="NiV"] <- pred.df.IPM$predicted_seropos[pred.df.IPM$antigen=="NiV"] - 1.96*predict.gam(gamIPMNiVNoSex, newdata = pred.df.IPM[pred.df.IPM$antigen=="NiV",], type="response", se.fit = T)$se
pred.df.IPM$predicted_seropos_uci[pred.df.IPM$antigen=="NiV"] <- pred.df.IPM$predicted_seropos[pred.df.IPM$antigen=="NiV"] + 1.96*predict.gam(gamIPMNiVNoSex, newdata = pred.df.IPM[pred.df.IPM$antigen=="NiV",], type="response", se.fit = T)$se

#HeV
pred.df.IPM$predicted_seropos[pred.df.IPM$antigen=="HeV"] <- predict.gam(gamIPMHeVNoSex, newdata = pred.df.IPM[pred.df.IPM$antigen=="HeV",], type="response")
pred.df.IPM$predicted_seropos_lci[pred.df.IPM$antigen=="HeV"] <- pred.df.IPM$predicted_seropos[pred.df.IPM$antigen=="HeV"] - 1.96*predict.gam(gamIPMHeVNoSex, newdata = pred.df.IPM[pred.df.IPM$antigen=="HeV",], type="response", se.fit = T)$se
pred.df.IPM$predicted_seropos_uci[pred.df.IPM$antigen=="HeV"] <- pred.df.IPM$predicted_seropos[pred.df.IPM$antigen=="HeV"] + 1.96*predict.gam(gamIPMHeVNoSex, newdata = pred.df.IPM[pred.df.IPM$antigen=="HeV",], type="response", se.fit = T)$se

pred.df.IPM$predicted_seropos_lci[pred.df.IPM$predicted_seropos_lci<0]<-0
pred.df.IPM$predicted_seropos_uci[pred.df.IPM$predicted_seropos_uci>1]<-1
pred.df.IPM$bat_species="Eidolon dupreanum"

#plot


gest.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum"), start = c(as.Date("2005-11-01"), as.Date("2006-08-03")), stop=c(as.Date("2005-11-15"), as.Date("2006-11-15")))
gest.dat$season = "gestation"
gest.dat$sex="female"
lact.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum"), start = c(as.Date("2005-11-16"), as.Date("2006-11-16")), stop=c(as.Date("2006-02-02"), as.Date("2007-02-02")))
lact.dat$season = "lactation"
lact.dat$sex="female"

#and winter
cal.dat <-  cbind.data.frame(bat_species=c("Eidolon dupreanum"), start = c(as.Date("2006-05-15")), stop=c(as.Date("2006-10-10")))
cal.dat$season <- "winter"
cal.dat$sex="male"

cal.dat <- rbind(cal.dat, gest.dat, lact.dat)
cal.dat$bat_species=="Eidolon dupreanum"
IPM.sero$bat_species=="Eidolon dupreanum"
IPM.sero$sig <- 0
IPM.sero$sig[IPM.sero$antigen=="NiV"] <- 1
IPM.sero$sig[IPM.sero$antigen=="HeV"] <- 1
IPM.sero$sig <- as.factor(IPM.sero$sig)

alphaz = c('0' = 0, '1' = 1, '2' = 0.5)
#with legend
pIPMLeg <- ggplot(data=IPM.sero) +
  scale_color_manual(values=colorz) +
  scale_fill_manual(values=filz2) +
  scale_size_continuous(name="N bats", breaks = c(10,20,30)) +
  geom_rect(data=cal.dat, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=date, y=seroprev, size=N_tot, color=antigen, alpha = sig)) +
  #geom_line(aes(x=date, y=seroprev, color=antigen, alpha = sig), linetype=2) +
  geom_linerange(aes(x=date, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen, alpha=sig), linetype=2, show.legend = F) +
  scale_alpha_manual(values=alphaz) +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values=filz1) +
  geom_ribbon(data=subset(pred.df.IPM, antigen=="HeV"), aes(x=date, ymin= predicted_seropos_lci, ymax=predicted_seropos_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df.IPM, antigen=="HeV"), aes(x=date, y= predicted_seropos, color=antigen), alpha=.3) +
  geom_ribbon(data=subset(pred.df.IPM, antigen=="NiV"), aes(x=date, ymin= predicted_seropos_lci, ymax=predicted_seropos_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df.IPM, antigen=="NiV"), aes(x=date, y= predicted_seropos, color=antigen)) +
  facet_grid(sex~bat_species) + theme_bw()+
  theme(axis.text = element_text(size=14), axis.title.x = element_blank(),
        axis.title.y = element_text(size=16), strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines")) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + guides(alpha="none"); pIPMLeg

#strip
pIPM <- ggplot(data=IPM.sero) +
  scale_color_manual(values=colorz) +
  scale_fill_manual(values=filz2) +
  scale_size_continuous(name="N bats", breaks = c(10,20,30)) +
  geom_rect(data=cal.dat, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3, show.legend = F) +
  geom_point(aes(x=date, y=seroprev, size=N_tot, color=antigen, alpha = sig), show.legend = F) +
  #geom_line(aes(x=date, y=seroprev, color=antigen, alpha = sig), linetype=2, show.legend = F) +
  geom_linerange(aes(x=date, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen, alpha=sig), linetype=2, show.legend = F) +
  scale_alpha_manual(values=alphaz) +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values=filz1) +
  geom_ribbon(data=subset(pred.df.IPM, antigen=="HeV"), aes(x=date, ymin= predicted_seropos_lci, ymax=predicted_seropos_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df.IPM, antigen=="HeV"), aes(x=date, y= predicted_seropos, color=antigen), show.legend = F) +
  geom_ribbon(data=subset(pred.df.IPM, antigen=="NiV"), aes(x=date, ymin= predicted_seropos_lci, ymax=predicted_seropos_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df.IPM, antigen=="NiV"), aes(x=date, y= predicted_seropos, color=antigen), show.legend = F) +
  facet_grid(sex~bat_species) + theme_bw()+
  theme(axis.text = element_text(size=14), axis.title.x = element_blank(), legend.text = element_blank(), legend.title = element_blank(),
        strip.background.x = element_rect(fill="white"), strip.background.y = element_blank(), strip.text.y = element_blank(),
        strip.text.x = element_text(size=16, face="italic"),
        axis.title.y = element_text(size=16),
        plot.margin = unit(c(.5,1,.5,1), "lines"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines")) +
  scale_x_continuous(breaks=c(as.Date("2006-01-01"),as.Date("2007-01-01")), labels=c("2006", "2007")) +
  coord_cartesian(expand = F) + ylab("seroprevalence\n") + guides(alpha="none"); pIPM

#and replot with the no sex version


#and plot the predictions
pred.df.IPM <- cbind.data.frame(date = rep(seq(as.Date("2005-11-01"), as.Date("2007-04-30"), by = 1), (length(unique(IPM.sero$antigen)))), antigen= rep(rep(unique(IPM.sero$antigen), each = length(seq(as.Date("2005-11-01"), as.Date("2007-04-30"), by = 1)))))
head(pred.df.IPM)
#pred.df.IPM =subset(pred.df.IPM , antigen=="NiV")
#pred.df.IPM$doy <- yday(pred.df.IPM$date)
#pred.df.IPM$year_f <- as.factor(year(pred.df.IPM$date))
#predict without effect of year for plotting
#NiV
pred.df.IPM$predicted_seropos <- NA
pred.df.IPM$predicted_seropos_lci <- NA
pred.df.IPM$predicted_seropos_uci <- NA
pred.df.IPM$predicted_seropos[pred.df.IPM$antigen=="NiV"] <- predict.gam(gamIPMNiVNoSex, newdata = pred.df.IPM[pred.df.IPM$antigen=="NiV",], type="response")
pred.df.IPM$predicted_seropos_lci[pred.df.IPM$antigen=="NiV"] <- pred.df.IPM$predicted_seropos[pred.df.IPM$antigen=="NiV"] - 1.96*predict.gam(gamIPMNiVNoSex, newdata = pred.df.IPM[pred.df.IPM$antigen=="NiV",], type="response", se.fit = T)$se
pred.df.IPM$predicted_seropos_uci[pred.df.IPM$antigen=="NiV"] <- pred.df.IPM$predicted_seropos[pred.df.IPM$antigen=="NiV"] + 1.96*predict.gam(gamIPMNiVNoSex, newdata = pred.df.IPM[pred.df.IPM$antigen=="NiV",], type="response", se.fit = T)$se

#HeV
pred.df.IPM$predicted_seropos[pred.df.IPM$antigen=="HeV"] <- predict.gam(gamIPMHeVNoSex, newdata = pred.df.IPM[pred.df.IPM$antigen=="HeV",], type="response")
pred.df.IPM$predicted_seropos_lci[pred.df.IPM$antigen=="HeV"] <- pred.df.IPM$predicted_seropos[pred.df.IPM$antigen=="HeV"] - 1.96*predict.gam(gamIPMHeVNoSex, newdata = pred.df.IPM[pred.df.IPM$antigen=="HeV",], type="response", se.fit = T)$se
pred.df.IPM$predicted_seropos_uci[pred.df.IPM$antigen=="HeV"] <- pred.df.IPM$predicted_seropos[pred.df.IPM$antigen=="HeV"] + 1.96*predict.gam(gamIPMHeVNoSex, newdata = pred.df.IPM[pred.df.IPM$antigen=="HeV",], type="response", se.fit = T)$se


pred.df.IPM$predicted_seropos_lci[pred.df.IPM$predicted_seropos_lci<0]<-0
pred.df.IPM$predicted_seropos_uci[pred.df.IPM$predicted_seropos_uci>1]<-1
pred.df.IPM$bat_species="Eidolon dupreanum"

IPM.sero.nosex$sig = 1
#IPM.sero.nosex$sig[IPM.sero.nosex$antigen!="NiV"] <- 0
IPM.sero.nosex$sig <- as.factor(IPM.sero.nosex$sig)

dat_text <- data.frame(
  label = c("p = 0 HeV\n p = 0 NiV"),
  x     = c(Inf),
  y     = c(1))

IPM.sero.nosex$dataset <- "Eidolon dupreanum - IPM"

pIPM_new <- ggplot(data=IPM.sero.nosex) +
  scale_color_manual(values=colorz) +
  scale_fill_manual(values=filz2) +
  scale_size_continuous(name="N bats", breaks = c(10,20,30)) +
  geom_rect(data=cal.dat, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3, show.legend = F) +
  geom_point(aes(x=date, y=seroprev, size=N_tot, color=antigen, alpha = sig), show.legend = F) +
  #geom_line(aes(x=date, y=seroprev, color=antigen, alpha = sig), linetype=2, show.legend = F) +
  geom_linerange(aes(x=date, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen, alpha=sig), linetype=2, show.legend = F) +
  scale_alpha_manual(values=alphaz) +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values=filz1) +
  geom_ribbon(data=subset(pred.df.IPM, antigen=="HeV"), aes(x=date, ymin= predicted_seropos_lci, ymax=predicted_seropos_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df.IPM, antigen=="HeV"), aes(x=date, y= predicted_seropos, color=antigen), show.legend = F) +
  geom_ribbon(data=subset(pred.df.IPM, antigen=="NiV"), aes(x=date, ymin= predicted_seropos_lci, ymax=predicted_seropos_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df.IPM, antigen=="NiV"), aes(x=date, y= predicted_seropos, color=antigen), show.legend = F) +
  facet_grid(.~dataset) + theme_bw()+
  theme(axis.text.x = element_text(size=14), 
        axis.title.x = element_blank(), 
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        # axis.title.y = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title.y = element_text(size=16),
        plot.title = element_text(size=16),
        strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  #ggtitle("Interannual seroprevalence (historical IPM dataset) in" ~ italic("E. dupreanum"))+
  scale_x_continuous(breaks=c(as.Date("2006-01-01"),as.Date("2007-01-01")), labels=c("2006", "2007")) +
  coord_cartesian(expand = F) + ylab("Seroprevalence") + guides(alpha="none"); pIPM_new

pIPM_new2<-pIPM_new+geom_text(data= dat_text,
        mapping = aes(x = x, y = y, label = label, hjust=1.1, vjust=1.5, fontface="italic"));pIPM_new2

## Add the historical JAE time series for comparison

JAE.dat <-  read.csv(file = paste0(homewd,"/serology/data/mix_2019_lum_dat.csv"),header = T, stringsAsFactors = F)

# new.dat <-  read.csv(file = paste0(homewd,"/data/mix_2019_lum_dat.csv"),header = T, stringsAsFactors = F)
head(JAE.dat)

head(catch.dat)
catch.merge <- dplyr::select(catch.dat,sampleid, processing_date, bat_sex, roost_site, sampling_session)

JAE.long <- merge(JAE.dat, catch.merge, by="sampleid", all.x = T)
head(JAE.long)
JAE.long$processing_date <- as.Date(JAE.long$processing_date)
unique(JAE.long$antigen)
unique(JAE.long$roost_site)
JAE.long$roost_site[is.na(JAE.long$roost_site)] <- "Moramanga"
JAE.long$roost_site[JAE.long$roost_site=="Ankarana_Canyon_Edup"] <- "Ankarana"
JAE.long$roost_site[JAE.long$roost_site=="Mahabo_Edup"] <- "Mahabo"
JAE.long$roost_site[JAE.long$roost_site=="Ankarana_Canyon_Rmad"] <- "Ankarana"
JAE.long$roost_site[JAE.long$roost_site=="Makira_Rmad"] <- "Makira"
JAE.long$roost_site[JAE.long$roost_site!="Makira" &JAE.long$roost_site!="Ankarana" & JAE.long$roost_site!="Mahabo" ] <- "Moramanga"

#get sampling session
sample.df <- ddply(JAE.long,.(sampling_session), summarise, min_date = min(processing_date), max_date =max(processing_date))
sample.df <- sample.df[complete.cases(sample.df),]
sample.df$mid_date <- ((sample.df$max_date-sample.df$min_date)/2) + sample.df$min_date

JAE.long <- merge(JAE.long, sample.df, by="sampling_session")

#now, sub to just NiV in Eidolon
JAE.long$plot <- 0
JAE.long$plot[JAE.long$bat_species=="Eidolon dupreanum" & JAE.long$antigen=="NiV"]<- 1
JAE.long$plot[JAE.long$bat_species=="Eidolon dupreanum" & JAE.long$antigen=="HeV"]<- 1


#and do sex- and antigen-specific GAMS
JAE.merge = subset(JAE.long, plot==1)
head(JAE.merge)

###here better compile data

#if in the same month and year, make the same datapoint
JAE.merge$month <- as.character(month(JAE.merge$mid_date))
JAE.merge$month[JAE.merge$month!="10" &JAE.merge$month!="11" & JAE.merge$month!="12"] <- paste0("0", JAE.merge$month[JAE.merge$month!="10" &JAE.merge$month!="11" & JAE.merge$month!="12"])
JAE.merge$year <- year(JAE.merge$mid_date)
JAE.merge$plot_date <- as.Date(as.character(paste0(JAE.merge$year, "-", JAE.merge$month, "-15")))

#unique(JAE.merge$site)
unique(JAE.merge$roost_site)
# JAE.merge$roost_site[JAE.merge$roost_site=="Angavokely_Edup"| JAE.merge$roost_site=="Angavobe_Edup"|JAE.merge$roost_site=="Maromizaha_Rmad"] <- "Moramanga"
# JAE.merge = subset(JAE.merge, roost_site=="Moramanga" )
JAE.sum <- ddply(JAE.merge,.(bat_species, roost_site, bat_sex, antigen, plot_date), summarise, N_seropos=sum(seropos), N_tot=length(seropos))
head(JAE.sum)
names(JAE.sum)[names(JAE.sum)=="plot_date"] <- "mid_date"


JAE.sum.all <- ddply(JAE.merge,.(bat_species, roost_site, antigen, plot_date), summarise, N_seropos=sum(seropos), N_tot=length(seropos))
head(JAE.sum.all)
names(JAE.sum.all)[names(JAE.sum.all)=="plot_date"] <-"mid_date"

JAE.sum.all$N_seroneg = JAE.sum.all$N_tot-JAE.sum.all$N_seropos

JAE.sum$N_seroneg = JAE.sum$N_tot-JAE.sum$N_seropos
JAE.sum$seroprev <- JAE.sum$N_seropos/JAE.sum$N_tot
JAE.sum$seroprev_lci <- binom.exact(JAE.sum$N_seropos, JAE.sum$N_tot)$lower
JAE.sum$seroprev_uci <- binom.exact(JAE.sum$N_seropos, JAE.sum$N_tot)$upper
JAE.sum$bat_sex[JAE.sum$bat_sex=="M"]<- "male"
JAE.sum$bat_sex[JAE.sum$bat_sex=="F"]<- "female"
JAE.sum$bat_sex <- as.factor(JAE.sum$bat_sex)
JAE.sum$mid_date <- as.Date(JAE.sum$mid_date)
JAE.sum$antigen <- as.factor(JAE.sum$antigen)
JAE.sum$year <- year(JAE.sum$mid_date)
JAE.sum$doy <- yday(JAE.sum$mid_date)
JAE.sum$year_f <- as.factor(JAE.sum$year)

#is there seasonality without sex by antigen?
#largest possible k that runs (k=-1 is actually 10)
JAE.sum.all$antigen <- as.factor(JAE.sum.all$antigen)
JAE.sum.all$roost_site <- as.factor(JAE.sum.all$roost_site)
JAE.sum.all$doy <- yday(JAE.sum.all$mid_date)
JAE.sum.all$year <- year(JAE.sum.all$mid_date)
JAE.sum.all$year_f <- as.factor(JAE.sum.all$year)

# all antigens
JAE.Eid <- subset(JAE.sum.all, bat_species=="Eidolon dupreanum")
JAE.Eid$antigen <- droplevels(JAE.Eid$antigen)
JAE.Eid$roost_site <- droplevels(JAE.Eid$roost_site)


gamJEidNoSex <- gam(cbind(N_seropos, N_seroneg)~   roost_site + 
                      s(as.numeric(doy), k=9, bs="cc",  by = roost_site) +
                      s(year_f, k=3, bs="fs",  by = roost_site),
                    data=JAE.Eid, 
                    family = "binomial", 
                    gamma=.1,
                    select = TRUE,
                    method = "REML")


gamJEidNoSex <- gam(cbind(N_seropos, N_seroneg)~   roost_site + 
                      s(as.numeric(doy), k=9, bs="cc",  by = roost_site) +
                      s(as.numeric(mid_date), k=6, bs="tp",  by = roost_site),
                    data=JAE.Eid, 
                    family = "binomial", 
                    gamma=1,
                    select = TRUE,
                    method = "REML")


#plot this version
gamJEidNoSexPlot <- gam(cbind(N_seropos, N_seroneg)~   roost_site + 
                      s(as.numeric(doy), k=9, bs="cc",  by = roost_site) +
                      s(as.numeric(mid_date), k=6, bs="tp",  by = roost_site),
                    data=JAE.Eid, 
                    family = "binomial", 
                    gamma=.25,
                    select = TRUE,
                    method = "REML")




summary(gamJEidNoSex)
#plot(gamJEidNoSex)

# and with sex
# all antigens
JAE.Eid.sex <- subset(JAE.sum, bat_species=="Eidolon dupreanum")
JAE.Eid.sex$antigen <- droplevels(JAE.Eid.sex$antigen)
JAE.Eid.sex$roost_site <- as.factor(JAE.Eid.sex$roost_site)
JAE.Eid.sex$roost_site <- droplevels(JAE.Eid.sex$roost_site)
JAE.Eid.sex$bat_sex <- droplevels(JAE.Eid.sex$bat_sex)


gamJEidSex <- gam(cbind(N_seropos, N_seroneg)~   roost_site + bat_sex +
                    s(as.numeric(doy),  k=9, bs="cc", by = interaction(roost_site, bat_sex)) +
                    s(year_f, k=3, bs="fs",  by = interaction(roost_site, bat_sex)), 
                  #s(as.numeric(mid_date), k=7, bs="cr"), 
                  data=JAE.Eid.sex, 
                  family = "binomial", 
                  select = TRUE,
                  method = "REML")


summary(gamJEidSex)

AIC(gamJEidNoSex, gamJEidSex) #no sex is better


#and plot the predictions
pred.df.JAE <- cbind.data.frame(mid_date = rep(seq(as.Date("2013-11-01"), as.Date("2016-01-30"), by = 1), 2), sex=rep(c("female", "male"), each=length(seq(as.Date("2013-11-01"), as.Date("2016-01-30"), by = 1))))
pred.df.JAE$bat_sex = "female"
pred.df.tmp <- pred.df.JAE
pred.df.tmp$bat_sex <- "male"
pred.df.JAE <- rbind(pred.df.JAE, pred.df.tmp)

pred.df.JAE$doy <- yday(pred.df.JAE$mid_date)
pred.df.JAE$year <- year(pred.df.JAE$mid_date)
pred.df.JAE$year_f <- as.factor(pred.df.JAE$year)
pred.df.JAE$roost_site <- "Moramanga"
pred.df.JAE$antigen <- "NiV"
pred.df.JAE2 <- pred.df.JAE
pred.df.JAE2$antigen <- "HeV"
pred.df.JAE <- rbind(pred.df.JAE, pred.df.JAE2)
#predict 

pred.df.JAE$predicted_seropos <- NA
pred.df.JAE$predicted_seropos_lci <- NA
pred.df.JAE$predicted_seropos_uci <- NA
#NiV
pred.df.JAE$predicted_seropos[pred.df.JAE$antigen=="NiV"] <- predict.gam(gamJEidSex, newdata = pred.df.JAE[pred.df.JAE$antigen=="NiV",], type="response")
pred.df.JAE$predicted_seropos_lci[pred.df.JAE$antigen=="NiV"] <- pred.df.JAE$predicted_seropos[pred.df.JAE$antigen=="NiV"] - 1.96*predict.gam(gamJEidSex, newdata = pred.df.JAE[pred.df.JAE$antigen=="NiV",], type="response", se.fit = T)$se
pred.df.JAE$predicted_seropos_uci[pred.df.JAE$antigen=="NiV"] <- pred.df.JAE$predicted_seropos[pred.df.JAE$antigen=="NiV"] + 1.96*predict.gam(gamJEidSex, newdata = pred.df.JAE[pred.df.JAE$antigen=="NiV",], type="response", se.fit = T)$se

# # #HeV
pred.df.JAE$predicted_seropos[pred.df.JAE$antigen=="HeV"] <- predict.gam(gamJEidSex, newdata = pred.df.JAE[pred.df.JAE$antigen=="HeV",], type="response")
pred.df.JAE$predicted_seropos_lci[pred.df.JAE$antigen=="HeV"] <- pred.df.JAE$predicted_seropos[pred.df.JAE$antigen=="HeV"] - 1.96*predict.gam(gamJEidSex, newdata = pred.df.JAE[pred.df.JAE$antigen=="HeV",], type="response", se.fit = T)$se
pred.df.JAE$predicted_seropos_uci[pred.df.JAE$antigen=="HeV"] <- pred.df.JAE$predicted_seropos[pred.df.JAE$antigen=="HeV"] + 1.96*predict.gam(gamJEidSex, newdata = pred.df.JAE[pred.df.JAE$antigen=="HeV",], type="response", se.fit = T)$se

pred.df.JAE$predicted_seropos_lci[pred.df.JAE$predicted_seropos_lci<0]<-0
pred.df.JAE$predicted_seropos_uci[pred.df.JAE$predicted_seropos_uci>1]<-1
pred.df.JAE$bat_species="Eidolon dupreanum"

# 
# #and do Rousettus too
# # seasonality without sex?
# gamJRouEboV <- gam(cbind(N_seropos, N_seroneg)~s(as.numeric(mid_date), k=5, bs="cr"), data=subset(JAE.sum.all, bat_species=="Rousettus madagascariensis" ), family = "binomial", select = TRUE)
# plot(gamJRouEboV, pages = 1)
# summary(gamJRouEboV) #no seasonality detected
# 
# 
# #by sex
# gamJRouEboVSex <- gam(cbind(N_seropos, N_seroneg)~bat_sex+s(as.numeric(mid_date), k=5, bs="cr", by=bat_sex), data=subset(JAE.sum, bat_species=="Rousettus madagascariensis" ), family = "binomial", select = TRUE)
# plot(gamJRouEboVSex, pages = 1)
# summary(gamJRouEboVSex) #no seasonality detected in males but females borderline
# #no difference in parametric terms
# 
# AIC(gamJRouEboV, gamJRouEboVSex) #simpler no sex model better supported
# 
# #and plot the predictions
# pred.df.JAE2 <- cbind.data.frame(mid_date = rep(seq(as.Date("2014-08-01"), as.Date("2016-01-30"), by = 1),2), sex=rep(c("female", "male"), each=length(seq(as.Date("2014-08-01"), as.Date("2016-01-30"), by = 1))))
# pred.df.JAE2$antigen <- "EboV"
# pred.df.JAE2$bat_sex = "female"
# pred.df.tmp <- pred.df.JAE2
# pred.df.tmp$bat_sex <- "male"
# pred.df.JAE2 <- rbind(pred.df.JAE2, pred.df.tmp)
# 
# #predict 
# #females
# pred.df.JAE2$predicted_seropos <- NA
# pred.df.JAE2$predicted_seropos_lci <- NA
# pred.df.JAE2$predicted_seropos_uci <- NA
# pred.df.JAE2$predicted_seropos<- predict.gam(gamJRouEboVSex, newdata = pred.df.JAE2, type="response")
# pred.df.JAE2$predicted_seropos_lci<- pred.df.JAE2$predicted_seropos - 1.96*predict.gam(gamJRouEboVSex, newdata = pred.df.JAE2, type="response", se.fit = T)$se
# pred.df.JAE2$predicted_seropos_uci <- pred.df.JAE2$predicted_seropos + 1.96*predict.gam(gamJRouEboVSex, newdata = pred.df.JAE2, type="response", se.fit = T)$se
# 
# pred.df.JAE2$predicted_seropos_lci[pred.df.JAE2$predicted_seropos_lci<0]<-0
# pred.df.JAE2$predicted_seropos_uci[pred.df.JAE2$predicted_seropos_uci>1]<-1
# pred.df.JAE2$bat_species="Rousettus madagascariensis"
# 
# 
# pred.df.JAE <- rbind(pred.df.JAE, pred.df.JAE2)



gest.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis", "Rousettus madagascariensis"), start = c(as.Date("2018-08-03"), as.Date("2019-08-03"), as.Date("2018-09-11"), as.Date("2019-09-11")), stop=c(as.Date("2018-11-15"), as.Date("2019-11-15"), as.Date("2018-12-11"), as.Date("2019-12-11")))
gest.dat$season = "gestation"
lact.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis"), start = c(as.Date("2018-01-01"), as.Date("2018-11-16"), as.Date("2019-11-16"), as.Date("2018-01-01"), as.Date("2018-12-12"), as.Date("2019-12-12")), stop=c(as.Date("2018-02-02"), as.Date("2019-02-02"), as.Date("2020-02-02"), as.Date("2018-02-19"), as.Date("2019-02-19"), as.Date("2020-02-19")))
lact.dat$season = "lactation"


#plot
gest.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis", "Rousettus madagascariensis"), start = c(as.Date("2013-11-01"), as.Date("2014-08-03"), as.Date("2015-08-03"), as.Date("2014-09-11"), as.Date("2015-09-11")), stop=c(as.Date("2013-11-15"), as.Date("2014-11-15"), as.Date("2015-11-15"), as.Date("2014-12-11"), as.Date("2015-12-11")))
gest.dat$season = "gestation"
gest.dat$bat_sex="female"

lact.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum","Eidolon dupreanum", "Rousettus madagascariensis", "Rousettus madagascariensis"), start = c(as.Date("2013-11-16"), as.Date("2014-11-16"), as.Date("2015-11-16"), as.Date("2014-12-12"), as.Date("2015-12-12")), stop=c(as.Date("2014-02-02"), as.Date("2015-02-02"),  as.Date("2016-01-30"),as.Date("2015-02-19"), as.Date("2016-01-30")))
lact.dat$season = "lactation"
lact.dat$bat_sex="female"

#and winter
cal.dat <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis",  "Rousettus madagascariensis"), start = c(as.Date("2014-05-15"), as.Date("2015-05-15"), as.Date("2014-08-01"), as.Date("2015-05-15")), stop=c(as.Date("2014-10-10"), as.Date("2015-10-10"), as.Date("2014-10-10"), as.Date("2015-10-10")))
cal.dat$season <- "winter"
cal.dat$bat_sex="male"

cal.dat <- rbind(cal.dat, gest.dat, lact.dat)


JAE.sum$sig <- 0
JAE.sum$sig[JAE.sum$antigen=="NiV"] <- 1
JAE.sum$sig[JAE.sum$antigen=="HeV"] <- 1

JAE.sum$sig<- as.factor(JAE.sum$sig)
#only plot the Eidolon data since the rousettus did not test the same antigens
#and came out negative for any seasonal signature

#with legend
pJAELeg <- ggplot(data=subset(JAE.sum, bat_species=="Eidolon dupreanum")) +
  scale_color_manual(values=colorz) +
  scale_fill_manual(values=filz2) +
  scale_size_continuous(name="N bats", breaks = c(10,20,30)) +
  geom_rect(data=subset(cal.dat, bat_species=="Eidolon dupreanum"), aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=mid_date, y=seroprev, size=N_tot, color=antigen, alpha=sig)) +
  #geom_line(aes(x=mid_date, y=seroprev, color=antigen, alpha=sig), linetype=2) +
  geom_linerange(aes(x=mid_date, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen, alpha=sig), linetype=2, show.legend = F) +
  scale_alpha_manual(values=alphaz) +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values=filz1) +
  geom_ribbon(data=subset(pred.df.JAE, bat_species=="Eidolon dupreanum" & antigen=="HeV") , aes(x=mid_date, ymin= predicted_seropos_lci, ymax=predicted_seropos_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df.JAE, bat_species=="Eidolon dupreanum"& antigen=="HeV"), aes(x=mid_date, y= predicted_seropos, color=antigen), show.legend = F) +
  geom_ribbon(data=subset(pred.df.JAE, bat_species=="Eidolon dupreanum" & antigen=="NiV") , aes(x=mid_date, ymin= predicted_seropos_lci, ymax=predicted_seropos_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df.JAE, bat_species=="Eidolon dupreanum"& antigen=="NiV"), aes(x=mid_date, y= predicted_seropos, color=antigen, group=year)) +
  facet_grid(bat_sex~bat_species, space = "free", scales = "free") + theme_bw() +
  theme(axis.text = element_text(size=14), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16), strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines")) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + guides(alpha="none"); pJAELeg

leg<-get_legend(pJAELeg)

#strip
pJAE <- ggplot(data=subset(JAE.sum, bat_species=="Eidolon dupreanum")) +
  scale_color_manual(values=colorz) +
  scale_fill_manual(values=filz2) +
  scale_size_continuous(name="N bats", breaks = c(10,20,30)) +
  geom_rect(data=subset(cal.dat, bat_species=="Eidolon dupreanum"), aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3, show.legend = F) +
  geom_point(aes(x=mid_date, y=seroprev, size=N_tot, color=antigen, alpha=sig), show.legend = F) +
  #geom_line(aes(x=mid_date, y=seroprev, color=antigen, alpha=sig), linetype=2, show.legend = F) +
  geom_linerange(aes(x=mid_date, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen, alpha=sig), linetype=2, show.legend = F) +
  scale_alpha_manual(values=alphaz) +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values=filz1) +
  geom_ribbon(data=subset(pred.df.JAE, bat_species=="Eidolon dupreanum" & antigen=="HeV") , aes(x=mid_date, ymin= predicted_seropos_lci, ymax=predicted_seropos_uci, fill=antigen), alpha=.1, show.legend = F) +
  geom_line(data=subset(pred.df.JAE, bat_species=="Eidolon dupreanum"& antigen=="HeV"), aes(x=mid_date, y= predicted_seropos, color=antigen), alpha=.3, show.legend = F) +
  geom_ribbon(data=subset(pred.df.JAE, bat_species=="Eidolon dupreanum" & antigen=="NiV") , aes(x=mid_date, ymin= predicted_seropos_lci, ymax=predicted_seropos_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df.JAE, bat_species=="Eidolon dupreanum"& antigen=="NiV"), aes(x=mid_date, y= predicted_seropos, color=antigen, group=year), show.legend = F) +
  facet_grid(bat_sex~bat_species, space = "free", scales = "free") + theme_bw() +
  theme(axis.text.x = element_text(size=14), axis.title.x = element_blank(), 
        strip.background.x = element_rect(fill="white"), 
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_blank(), strip.background.y = element_blank(), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,1,.5,1), "lines")) +
  scale_x_continuous(breaks=c(as.Date("2014-01-01"),as.Date("2015-01-01"),as.Date("2016-01-01")), labels=c("2014", "2015", "2016")) +
  coord_cartesian(expand = F) + ylab("seroprevalence") + guides(alpha="none"); pJAE




# and repeat as one

#and plot the predictions
pred.df.JAE <- cbind.data.frame(mid_date = seq(as.Date("2013-11-01"), as.Date("2016-01-30"), by = 1) )

pred.df.JAE$doy <- yday(pred.df.JAE$mid_date)
pred.df.JAE$year <- year(pred.df.JAE$mid_date)
pred.df.JAE$year_f <- as.factor(pred.df.JAE$year)
pred.df.JAE$roost_site <- "Moramanga"
pred.df.JAE$antigen="NiV"
pred.df.JAE2 <- pred.df.JAE
pred.df.JAE2$antigen <- "HeV"
pred.df.JAE <- rbind(pred.df.JAE, pred.df.JAE2)
#predict 

pred.df.JAE$predicted_seropos <- NA
pred.df.JAE$predicted_seropos_lci <- NA
pred.df.JAE$predicted_seropos_uci <- NA
#NiV
pred.df.JAE$predicted_seropos[pred.df.JAE$antigen=="NiV"] <- predict.gam(gamJEidNoSexPlot, newdata = pred.df.JAE[pred.df.JAE$antigen=="NiV",], type="response")
pred.df.JAE$predicted_seropos_lci[pred.df.JAE$antigen=="NiV"] <- pred.df.JAE$predicted_seropos[pred.df.JAE$antigen=="NiV"] - 1.96*predict.gam(gamJEidNoSexPlot, newdata = pred.df.JAE[pred.df.JAE$antigen=="NiV",], type="response", se.fit = T)$se
pred.df.JAE$predicted_seropos_uci[pred.df.JAE$antigen=="NiV"] <- pred.df.JAE$predicted_seropos[pred.df.JAE$antigen=="NiV"] + 1.96*predict.gam(gamJEidNoSexPlot, newdata = pred.df.JAE[pred.df.JAE$antigen=="NiV",], type="response", se.fit = T)$se

#HeV
pred.df.JAE$predicted_seropos[pred.df.JAE$antigen=="HeV"] <- predict.gam(gamJEidNoSex, newdata = pred.df.JAE[pred.df.JAE$antigen=="HeV",], type="response")
pred.df.JAE$predicted_seropos_lci[pred.df.JAE$antigen=="HeV"] <- pred.df.JAE$predicted_seropos[pred.df.JAE$antigen=="HeV"] - 1.96*predict.gam(gamJEidNoSex, newdata = pred.df.JAE[pred.df.JAE$antigen=="HeV",], type="response", se.fit = T)$se
pred.df.JAE$predicted_seropos_uci[pred.df.JAE$antigen=="HeV"] <- pred.df.JAE$predicted_seropos[pred.df.JAE$antigen=="HeV"] + 1.96*predict.gam(gamJEidNoSex, newdata = pred.df.JAE[pred.df.JAE$antigen=="HeV",], type="response", se.fit = T)$se
# # 
pred.df.JAE$predicted_seropos_lci[pred.df.JAE$predicted_seropos_lci<0]<-0
pred.df.JAE$predicted_seropos_uci[pred.df.JAE$predicted_seropos_uci>1]<-1
pred.df.JAE$bat_species="Eidolon dupreanum"


JAE.sum.all$seroprevalence <- JAE.sum.all$N_seropos/JAE.sum.all$N_tot
JAE.sum.all$seroprev_lci <- binom.exact(JAE.sum.all$N_seropos, JAE.sum.all$N_tot)$lower
JAE.sum.all$seroprev_uci <- binom.exact(JAE.sum.all$N_seropos, JAE.sum.all$N_tot)$upper

dat_text <- data.frame(
  label = c("p = NS HeV\n p = NS NiV"),
  x     = c(Inf),
  y     = c(1))

JAE.sum.all$dataset <- "Eidolon dupreanum - JAE"

#strip
pJAE_inter <- ggplot(data=subset(JAE.sum.all, bat_species=="Eidolon dupreanum")) +
  scale_color_manual(values=colorz) +
  scale_fill_manual(values=filz2) +
  scale_size_continuous(name="N bats", breaks = c(10,20,30)) +
  geom_rect(data=subset(cal.dat, bat_species=="Eidolon dupreanum"), aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=season), alpha=.3, show.legend = F) +
  geom_point(aes(x=mid_date, y=seroprevalence, size=N_tot, color=antigen, alpha=.3,), show.legend = F) +
  #geom_line(aes(x=mid_date, y=seroprev, color=antigen, alpha=sig), linetype=2, show.legend = F) +
  geom_linerange(aes(x=mid_date, ymin=seroprev_lci, ymax=seroprev_uci, color=antigen, alpha=.3), linetype=2, show.legend = F) +
  #scale_alpha_manual(values=alphaz) +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values=filz1) +
  geom_ribbon(data=subset(pred.df.JAE, bat_species=="Eidolon dupreanum" & antigen=="HeV") , aes(x=mid_date, ymin= predicted_seropos_lci, ymax=predicted_seropos_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df.JAE, bat_species=="Eidolon dupreanum"& antigen=="HeV"), aes(x=mid_date, y= predicted_seropos, color=antigen), show.legend = F) +
  geom_ribbon(data=subset(pred.df.JAE, bat_species=="Eidolon dupreanum" & antigen=="NiV") , aes(x=mid_date, ymin= predicted_seropos_lci, ymax=predicted_seropos_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.df.JAE, bat_species=="Eidolon dupreanum"& antigen=="NiV"), aes(x=mid_date, y= predicted_seropos, color=antigen, group=year), show.legend = F) +
  facet_grid(.~dataset, space = "free", scales = "free") + theme_bw() +
  theme(axis.text.x = element_text(size=14), 
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        #axis.text.y = element_text(size=14),
        #axis.title.y = element_text(size=16),
        plot.title = element_text(size=16),
        strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  scale_x_continuous(breaks=c(as.Date("2014-01-01"),as.Date("2015-01-01"),as.Date("2016-01-01")), labels=c("2014", "2015", "2016")) +
  coord_cartesian(expand = F) + ylab("Seroprevalence") + guides(alpha="none"); pJAE_inter

pJAE_inter2<-pJAE_inter+geom_text(data= dat_text,
                              mapping = aes(x = x, y = y, label = label, hjust=1.1, vjust=1.5, fontface="italic"));pJAE_inter2

####put figures together
supp<-cowplot::plot_grid(pIPM_new2, pJAE_inter2, leg, 
                          labels=c("A","B","",""),
                          rel_widths=c(0.6,1,  0.2), rel_heights = c(1,1,1),
                          ncol=3,align="hv", axis="l", label_size = 23)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        plot.background = element_rect(fill = "white", colour = NA)); supp


#save figures
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/final_figures_tables/gwen_tmp/"
ggsave(file = paste0(homewd, "Fig3_mountainplot_seasonal_serology.png"),
       plot= fig3,
       units="mm",  
       width=150, 
       height=140, 
       scale=2, 
       dpi=500)

ggsave(file = paste0(homewd, "FigX_historical_and_rousettus_serology.png"),
       plot= supp,
       units="mm",  
       width=150, 
       height=50, 
       scale=2, 
       dpi=500)

