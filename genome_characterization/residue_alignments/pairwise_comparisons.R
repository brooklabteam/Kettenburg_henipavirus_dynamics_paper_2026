rm(list=ls())

#packages
library(plyr) 
library(dplyr) 
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(ggrepel)
library(hrbrthemes)
library(reshape)
library(paletteer)
library(ggrepel)
library(hrbrthemes)
library(reshape)
library(paletteer)
library(ggh4x)
#####################################################################
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_circulating_clades/genome_characterization/"
setwd(paste0(homewd, "residue_alignments")) 

#load data
dat <- read.csv(file = "V_W_pairwise_comparisons.csv", header = T, stringsAsFactors = F)
head(dat)
names(dat)

dat$var2[dat$var2=="ROU125_NODE_1_czid2023_add2"] <- "PX127342 AndV_add2"
dat$var2[dat$var2=="ROU148_manual_assembly_add2"] <- "PX127343 AndV_add2"
dat$var2[dat$var2=="ROU173_manual_assembly_add2"] <- "PX127344 AndV_add2"
dat$var2[dat$var2=="ROU125_NODE_1_czid2023_add1"] <- "PX127342 AndV_add1"
dat$var2[dat$var2=="ROU148_manual_assembly_add1"] <- "PX127343 AndV_add1"
dat$var2[dat$var2=="ROU173_manual_assembly_add1"] <- "PX127344 AndV_add1"

dat$gene<-factor(dat$gene, levels=c("V protein", "W protein"))

#subset the data
V<-subset(dat, gene=="V protein")
W<-subset(dat, gene=="W protein")

#plot comparison to known sequences
compare_V <- ggplot(V, aes(x=var1_alt, y=var2, fill=id)) +
  geom_tile() +
  geom_text(aes(label = id), color = "black", size = 3) +
 # scale_fill_viridis_c(option = "C", direction = -1) +
  scale_fill_gradient(low="white", high="tomato1")+
  facet_nested(gene~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Identity",
    title="")+
  theme_linedraw()+
  scale_y_discrete(position="right", limits=rev)+
  scale_x_discrete(position="top")+
  theme(plot.margin = margin(0, 0, 2, 2, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour=c("black","black","black","black","black","black","black","black","black","black","darkcyan"), face="bold"),
        axis.text.y = element_text(colour=c("darkcyan","darkcyan","darkcyan"), face="bold"),
        plot.title = element_text(size=12, face="italic"), 
        legend.position = "left",
        strip.text = element_text(size=12))
compare_V

compare_W <- ggplot(W, aes(x=var1_alt, y=var2, fill=id)) +
  geom_tile() +
  geom_text(aes(label = id), color = "black", size = 3) +
  # scale_fill_viridis_c(option = "C", direction = -1) +
  scale_fill_gradient(low="white", high="tomato1")+
  facet_nested(gene~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Identity",
    title="")+
  theme_linedraw()+
  scale_y_discrete(position="right", limits=rev)+
  scale_x_discrete(position="top")+
  theme(plot.margin = margin(0, 0, 2, 2, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour=c("black","black","black","black","black","black","black","black","black","black","darkcyan"), face="bold"),
        axis.text.y = element_text(colour=c("darkcyan","darkcyan","darkcyan"), face="bold"),
        plot.title = element_text(size=12, face="italic"), 
        legend.position = "left",
        strip.text = element_text(size=12))
compare_W

#add together: 
pair<-plot_grid(compare_V, compare_W, labels=c("",""),
                rel_widths=c(1,1), rel_heights = c(1,1),
                ncol=1,align="hv", axis="l", label_size = 23)
pair



