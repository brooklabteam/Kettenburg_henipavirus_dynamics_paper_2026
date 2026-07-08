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

#use this script to shwo qPCR cutoffs

homewd = "/Users/gwenddolenkettenburg/Desktop/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
batwd = "/Users/gwenddolenkettenburg/Desktop/Developer/Madagascar-Bat-Data"

#now load pcr data
dat <- read.csv(file = paste0(homewd, "/virus_shedding/data/henipavirus_shedding_data.csv"))
head(dat)

#subset to just the data we need for further analysis
dat.ct <- dplyr::select(dat, sampleid, 
                        qPCR_GhV_multi_M_ct_ur, qPCR_GhV_multi_N_ct_ur, qPCR_GhV_multi_L_ct_ur,
                        qPCR_AndV_multi_M_ct_ur, qPCR_AndV_multi_N_ct_ur, qPCR_AndV_multi_L_ct_ur,
                        qPCR_AngV_multi_M_ct_ur, qPCR_AngV_multi_N_ct_ur, qPCR_AngV_multi_L_ct_ur)
dat.pos <- dplyr::select(dat, sampleid, 
                         qPCR_GhV_multi_M_ur, qPCR_GhV_multi_N_ur, qPCR_GhV_multi_L_ur,
                         qPCR_AndV_multi_M_ur, qPCR_AndV_multi_N_ur, qPCR_AndV_multi_L_ur,
                         qPCR_AngV_multi_M_ur, qPCR_AngV_multi_N_ur, qPCR_AngV_multi_L_ur)

#make long versions
shed.long.pos <- melt(dat.pos, id.vars = c("sampleid"), value.name = "pos", variable.name = "virus")
head(shed.long.pos)
shed.long.pos$virus <- as.character(shed.long.pos$virus) 
shed.long.pos$virus <- gsub(pattern="qPCR_", replacement = "", x=shed.long.pos$virus, fixed = T)
names(shed.long.pos)[1] <- "sampleid"
#complete cases
shed.long.pos <- shed.long.pos[complete.cases(shed.long.pos),]

shed.long.ct <- melt(dat.ct, id.vars = c("sampleid"), value.name = "ct", variable.name = "virus")
head(shed.long.ct)

shed.long.ct$virus <- as.character(shed.long.ct$virus) 
shed.long.ct$virus <- gsub(pattern="qPCR_", replacement = "", x=shed.long.ct$virus, fixed = T)
names(shed.long.ct)[1] <- "sampleid"
#complete cases
shed.long.ct <- shed.long.ct[complete.cases(shed.long.ct),]

#restructure 
shed.long.pos.ghvm<-subset(shed.long.pos, virus=="GhV_multi_M_ur")
shed.long.pos.ghvn<-subset(shed.long.pos, virus=="GhV_multi_N_ur")
shed.long.pos.ghvl<-subset(shed.long.pos, virus=="GhV_multi_L_ur")
shed.long.pos.andvm<-subset(shed.long.pos, virus=="AndV_multi_M_ur")
shed.long.pos.andvn<-subset(shed.long.pos, virus=="AndV_multi_N_ur")
shed.long.pos.andvl<-subset(shed.long.pos, virus=="AndV_multi_L_ur")
shed.long.pos.angvm<-subset(shed.long.pos, virus=="AngV_multi_M_ur")
shed.long.pos.angvn<-subset(shed.long.pos, virus=="AngV_multi_N_ur")
shed.long.pos.angvl<-subset(shed.long.pos, virus=="AngV_multi_L_ur")

shed.long.pos.m<-rbind( shed.long.pos.ghvm, shed.long.pos.andvm, shed.long.pos.angvm)
shed.long.pos.m$gene= "M"
shed.long.pos.n<-rbind( shed.long.pos.ghvn, shed.long.pos.andvn, shed.long.pos.angvn)
shed.long.pos.n$gene= "N"
shed.long.pos.l<-rbind( shed.long.pos.ghvl, shed.long.pos.andvl, shed.long.pos.angvl)
shed.long.pos.l$gene= "L"

shed.long.pos<-rbind(shed.long.pos.m, shed.long.pos.n, shed.long.pos.l)
shed.long.pos$gene <- as.factor(shed.long.pos$gene)

shed.long.ct.ghvm<-subset(shed.long.ct, virus=="GhV_multi_M_ct_ur")
shed.long.ct.ghvn<-subset(shed.long.ct, virus=="GhV_multi_N_ct_ur")
shed.long.ct.ghvl<-subset(shed.long.ct, virus=="GhV_multi_L_ct_ur")
shed.long.ct.andvm<-subset(shed.long.ct, virus=="AndV_multi_M_ct_ur")
shed.long.ct.andvn<-subset(shed.long.ct, virus=="AndV_multi_N_ct_ur")
shed.long.ct.andvl<-subset(shed.long.ct, virus=="AndV_multi_L_ct_ur")
shed.long.ct.angvm<-subset(shed.long.ct, virus=="AngV_multi_M_ct_ur")
shed.long.ct.angvn<-subset(shed.long.ct, virus=="AngV_multi_N_ct_ur")
shed.long.ct.angvl<-subset(shed.long.ct, virus=="AngV_multi_L_ct_ur")

shed.long.ct.m<-rbind(shed.long.ct.ghvm, shed.long.ct.andvm, shed.long.ct.angvm)
shed.long.ct.m$gene= "M"
shed.long.ct.n<-rbind(shed.long.ct.ghvn, shed.long.ct.andvn, shed.long.ct.angvn)
shed.long.ct.n$gene= "N"
shed.long.ct.l<-rbind( shed.long.ct.ghvl, shed.long.ct.andvl, shed.long.ct.angvl)
shed.long.ct.l$gene= "L"

shed.long.ct<-rbind(shed.long.ct.m, shed.long.ct.n, shed.long.ct.l)
shed.long.ct$gene <- as.factor(shed.long.ct$gene)

#simplify names
shed.long.pos$virus[shed.long.pos$virus=="GhV_multi_M_ur"] <- "GhV_cov"
shed.long.pos$virus[shed.long.pos$virus=="GhV_multi_N_ur"] <- "GhV_cov"
shed.long.pos$virus[shed.long.pos$virus=="GhV_multi_L_ur"] <- "GhV_cov"
shed.long.pos$virus[shed.long.pos$virus=="AndV_multi_M_ur"] <- "AndV_cov"
shed.long.pos$virus[shed.long.pos$virus=="AndV_multi_N_ur"] <- "AndV_cov"
shed.long.pos$virus[shed.long.pos$virus=="AndV_multi_L_ur"] <- "AndV_cov"
shed.long.pos$virus[shed.long.pos$virus=="AngV_multi_M_ur"] <- "AngV_cov"
shed.long.pos$virus[shed.long.pos$virus=="AngV_multi_N_ur"] <- "AngV_cov"
shed.long.pos$virus[shed.long.pos$virus=="AngV_multi_L_ur"] <- "AngV_cov"

shed.long.pos$virus <- as.factor(shed.long.pos$virus)

shed.long.ct$virus[shed.long.ct$virus=="GhV_multi_M_ct_ur"] <- "GhV_cov"
shed.long.ct$virus[shed.long.ct$virus=="GhV_multi_N_ct_ur"] <- "GhV_cov"
shed.long.ct$virus[shed.long.ct$virus=="GhV_multi_L_ct_ur"] <- "GhV_cov"
shed.long.ct$virus[shed.long.ct$virus=="AndV_multi_M_ct_ur"] <- "AndV_cov"
shed.long.ct$virus[shed.long.ct$virus=="AndV_multi_N_ct_ur"] <- "AndV_cov"
shed.long.ct$virus[shed.long.ct$virus=="AndV_multi_L_ct_ur"] <- "AndV_cov"
shed.long.ct$virus[shed.long.ct$virus=="AngV_multi_M_ct_ur"] <- "AngV_cov"
shed.long.ct$virus[shed.long.ct$virus=="AngV_multi_N_ct_ur"] <- "AngV_cov"
shed.long.ct$virus[shed.long.ct$virus=="AngV_multi_L_ct_ur"] <- "AngV_cov"

shed.long.ct$virus <- as.factor(shed.long.ct$virus)

#combine
shed.long<-merge(shed.long.pos, shed.long.ct)

#merge with the catching data
catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))
#catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/Historic-Catching-Data-Through2020/catching_data_historic.csv"))

head(catch.dat)
names(catch.dat)

# drop most coshedns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)

#and merge all three datasets
merge.dat <- merge(shed.long, catch.dat, by="sampleid", all.x = T)

#get the midpoint processing date per sampling session
merge.dat$processing_date <- as.Date(merge.dat$processing_date)
sampling.date <- ddply(merge.dat, .(bat_species, sampling_session), summarise, min_date = min(processing_date), max_date = max(processing_date))
sampling.date$mid_date <- (sampling.date$max_date-sampling.date$min_date)/2 + sampling.date$min_date
sampling.merge <- dplyr::select(sampling.date, bat_species, sampling_session, mid_date)
merge.dat <- merge(merge.dat, sampling.merge, by = c("bat_species", "sampling_session"))
merge.dat$year <- year(merge.dat$processing_date)
merge.dat$month <- as.character(month(merge.dat$processing_date))
merge.dat$month[merge.dat$month!="10"& merge.dat$month!="11" &merge.dat$month!="12"]<-paste0("0", merge.dat$month[merge.dat$month!="10"& merge.dat$month!="11" &merge.dat$month!="12"])
#merge.dat$mid_date <- as.Date(paste0(merge.dat$year, "-", merge.dat$month, "-15"))
head(merge.dat)

#also add in quantitative age data where available:
age.dat <- read.csv(file = paste0(batwd, "/Age-Data/aged_bats.csv"), header = T, stringsAsFactors = F)
age.dat <- dplyr::select(age.dat, sampleid, age)

#and merge
merge.dat <- merge(merge.dat, age.dat, by="sampleid", all.x = T)

#curate date once more
merge.dat$mid_date <- as.Date(merge.dat$mid_date)

colorz = c(
           'Eidolon dupreanum'="goldenrod1",
           'Rousettus madagascariensis'="aquamarine3",
           'Pteropus rufus'="#EE5C42")

#subset to just those that are positive
merge.dat2<-subset(merge.dat, pos=="1")
#and remove hipposideros samples
merge.dat2<-subset(merge.dat2, bat_species=="Eidolon dupreanum" | bat_species=="Pteropus rufus" |bat_species=="Rousettus madagascariensis")

#plot just those that are positive in each gene
#title<-expression(paste("ED_cov sensitivity: L>E>N\nPR_cov sensitivity: L>E>N\nRM1_cov sensitivity: L>N>E\nRM2_cov sensitivity: L>E>N"))
dat_text <- data.frame(
  label = c(
            "M: R^2 = 0.998\n L: R^2 = 0.999\n N: R^2 = 0.990", 
            "M: R^2 = 0.992\n L: R^2 = 0.998\n N: R^2 = 0.993", 
            "M: R^2 = 0.902\n L: R^2 = 0.885\n N: R^2 = 0.867"),
  virus   = c("GhV", "AngV", "AndV"),
  x     = c(Inf, Inf, Inf, Inf),
    y     = c(12, 12, 12, 12)
)

p1 <- ggplot(subset(merge.dat2)) + 
  facet_nested(.~virus, scales = "free", space = "free") +
  geom_line(aes(x=gene, y=ct, group=sampleid, color=bat_species), show.legend=F) + theme_bw()+
  geom_point(aes(x=gene, y=ct, color=bat_species), show.legend=T) + 
  scale_color_manual(values=colorz) + 
  ylab("Ct") +
  xlab("Gene") +
  #scale_y_continuous(limits = c(10, 40))+
  scale_y_reverse(limits=c(12,40))+
  theme(axis.title.x = element_text(size=14) , panel.grid = element_blank(),
        axis.text.x = element_text(size=14), 
        axis.title.y = element_text(size=14), 
        legend.position = c(0.91,0.7),
        #legend.position = "right",
        legend.text = element_text(face="italic"),
        axis.text.y = element_text(size=14), 
        strip.background = element_rect(fill="white"),
        strip.text.x = element_text(size=14, color = "black"), 
        strip.text.y = element_text(size=14, face="italic", color="black"),
        plot.margin = unit(c(0.1,0.1,0,0), "cm")); p1

pcr<-p1+geom_text(data= dat_text, aes(x = x, y = y, label = label, hjust=1.2, vjust=1),
                 fontface="italic");pcr

##put them together
library(cowplot)

pcr2<-cowplot::plot_grid(NULL,pcr,NULL, labels=c("","", ""),
                          rel_widths=c(0.02,1,0), rel_heights = c(1,1),
                          ncol=3,align="hv",  label_size = 23)+
  theme(plot.background = element_rect(fill = "white", colour = NA)); pcr2

plot2<-cowplot::plot_grid(a, pcr2, labels=c("", "C"),
                rel_widths=c(1,1), rel_heights = c(1.3,1),
                ncol=1,align="h",axis="lr",  label_size = 23)+
  theme(plot.background = element_rect(fill = "white", colour = NA)); plot2

homewd = "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/final_figures_tables/gwen_tmp/"
ggsave(file = paste0(homewd, "/FigX_qPCR_cutoffs.png"),
       plot= plot2,
       units="mm",
       width=180,
       height=200,
       scale=2,
       dpi=500)



