#Script to make bat phylogeny for supplemental

rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
#install.packages('gdata')
library(gdata)
library(phylotools)
library(phylobase)
library(cowplot)

##First plot the phylogenetic tree

##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/phylo_trees/"
#homewd= "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/phylo_trees/"

#henipa tree
setwd(paste0(homewd,"/polymerase"))
tree2 <-  read.tree("batLmid.treefile") 
rooted.tree2 <- root(tree2, which(tree2$tip.label == "YP_009094051.1"))
#drop root and other reference sequences for aesthetics
rooted.tree2<-drop.tip(rooted.tree2, "YP_009094051.1")
# rooted.tree2<-drop.tip(rooted.tree2, "YP_010797981.1")
# rooted.tree2<-drop.tip(rooted.tree2, "YP_010790495.1")
# rooted.tree2<-drop.tip(rooted.tree2, "YP_010799241.1")
# rooted.tree2<-drop.tip(rooted.tree2, "YP_010798268.1")
# rooted.tree2<-drop.tip(rooted.tree2, "NP_056924.1")
# rooted.tree2<-drop.tip(rooted.tree2, "NP_054697.1")
# rooted.tree2<-drop.tip(rooted.tree2, "YP_138518.1")
# rooted.tree2<-drop.tip(rooted.tree2, "YP_010790514.1")
# rooted.tree2<-drop.tip(rooted.tree2, "YP_009094339.1")
# rooted.tree2<-drop.tip(rooted.tree2, "YP_010790479.1")
# rooted.tree2<-drop.tip(rooted.tree2, "YP_010790486.1")
# rooted.tree2<-drop.tip(rooted.tree2, "YP_010774627.1")
# rooted.tree2<-drop.tip(rooted.tree2, "YP_010790504.1")
# rooted.tree2<-drop.tip(rooted.tree2, "AET43339.1")
#load tree data prepared from elsewhere
dat2 <- read.csv(("L_bat_metadata.csv"), header = T, stringsAsFactors = F)

#check subgroup names
unique(dat2$Genus)
dat2$novel <- as.factor(dat2$novel)

# colz = c("Orthoavulavirus"="darkorange","Paraavulavirus"="darkorange3","Metaavulavirus"="darkgoldenrod1","Jeilongvirus"="coral",
#          "Morbillivirus"="coral3","Respirovirus"="deeppink1", "Parahenipavirus"="deeppink4",
#          "Henipavirus"="firebrick1","Cynoglossusvirus"="grey50","Pararubulavirus"="firebrick4","Salemvirus"="indianred1",
#          "Ferlavirus"="indianred3","Orthorubulavirus"="indianred4",
#          "Hoplichthysvirus"="lightsalmon1","Narmovirus"="lightsalmon3","Aquaparamyxovirus"="lightpink2","Scoliodonvirus"="lightpink3",
#          "Synodonvirus"="pink1","Unclassified"="black",
#          "Sunshinevirus"="chartreuse1")

colz = c("Orthoavulavirus"="black","Paraavulavirus"="black","Metaavulavirus"="black","Jeilongvirus"="black","Parajeilongvirus"="black",
         "Morbillivirus"="black","Respirovirus"="black", "Parahenipavirus"="black",
         "Henipavirus"="deepskyblue2","Cynoglossusvirus"="black","Pararubulavirus"="goldenrod","Salemvirus"="black",
         "Ferlavirus"="black","Orthorubulavirus"="black",
         "Hoplichthysvirus"="black","Narmovirus"="black","Aquaparamyxovirus"="black","Scoliodonvirus"="black",
         "Synodonvirus"="black","Unclassified"="black", "Unclassified henipavirus"="black", "Tupaivirus"="black", "Narmovirus"="black",
         "Hippocavirus"="black", "Bovinavirus"="black","Sunshinevirus"="black")

#now get new tip labels
dat2$old_tip_label <- dat2$tip_label
dat2$new_label <- NA

#now, select the name based on what components are present for each sample
#all components with values:
dat2$new_label[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)] <- paste(dat2$Accession[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                  dat2$Organism_Name[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                  dat2$Isolate[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                  dat2$source[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                  dat2$Country[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                  dat2$Collection_Date[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)])
# 
# #and if there is an NA just drop it
# 
# #here NA in Isolate only:
dat2$new_label[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)] <- paste(dat2$Accession[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                 dat2$Organism_Name[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                 #dat2$Isolate[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                 dat2$source[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                 dat2$Country[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                 dat2$Collection_Date[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)])


# dat2$new_label[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)] <- paste(dat2$Accession[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
#                                                                                                                   dat2$Organism_Name[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)])
# 
# #and if there is an NA just drop it
# 
# #here NA in Isolate only:
# dat2$new_label[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)] <- paste(dat2$Accession[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|", 
#                                                                                                                  dat2$Organism_Name[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)])



#make sure to sort in order
colnames(dat2)[colnames(dat2)=="Accession"]="old_tip_label"
tree.dat2 <- data.frame(old_tip_label=rooted.tree2$tip.label, num =1:length(rooted.tree2$tip.label))
head(tree.dat2)
head(dat2)
tree.dat2 <- merge(tree.dat2, dat2, by = "old_tip_label", all.x = T, sort = F)
names(tree.dat2)
tree.dat2$tip_label <- tree.dat2$new_label
tree.dat2 <- dplyr::select(tree.dat2, tip_label, Isolate, Host,Host2, source, Country, Collection_Date, Genus, Organism_Name, novel, old_tip_label)
rooted.tree2$tip.label <- tree.dat2$tip_label
head(tree.dat2)
cbind(tree.dat2$old_tip_label, rooted.tree2$tip.label)
tree.dat2$tip_label#all look good

#shaped and highlighting
tree.dat2$Host[tree.dat2$Host==0] <- "Non-bat host"
tree.dat2$Host[tree.dat2$Host==1] <- "Eidolon dupreanum"
tree.dat2$Host[tree.dat2$Host==2] <- "Rousettus madagascariensis"
tree.dat2$Host[tree.dat2$Host==3] <- "Other bat host"
tree.dat2$Host[tree.dat2$Host==4] <- "Collapsed clade - multihost"
tree.dat2$Host <- factor(tree.dat2$Host, levels = c("Eidolon dupreanum", "Rousettus madagascariensis", 
                                                    "Other bat host", "Non-bat host", "Collapsed clade - multihost"))
#tree.dat2$Host <- as.factor(tree.dat2$Host)
shapez = c("Non-bat host" =  18, "Eidolon dupreanum" = 17, "Rousettus madagascariensis" = 19, "Other bat host" =  15, "Collapsed clade - multihost"=25)

##Get the clade numbers so we can collapse unnecessary clades
ggtree(rooted.tree2) + geom_text(aes(label=node), hjust=-.3)

#henipavirus focused tree
title<-expression(paste("                                Polymerase"))
p2 <- ggtree(rooted.tree2) %<+% tree.dat2 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0.5,show.legend = T) +
  guides(color="none")+
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz, guide =
                       guide_legend(label.theme = element_text(angle = 0, face = "italic")))+
  scale_shape_manual(values=shapez, guide =
                       guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
  new_scale_color() +
  #ggtitle(title)+
  geom_tiplab(aes(color = novel, show.legend=F), geom = "text", Genus="Helvetica", size=4, nudge_x=0.05) +
  guides(color="none")+#
  scale_color_manual(values=c("0"="black","1"="red")) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm")) +
  geom_cladelabel(node= 534, label='HeV/NiV sister clade', offset=8.2, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=520, label='AndV clade', offset=8.2, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=1028, label='AndV-2 clade', offset=8.2, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=1032, label='AndV-3.1 clade', offset=8.2, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=511, label='AndV-3.2 clade', offset=8.2, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=528, label='GhV clade', offset=8.2, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=549, label='AngV clade', offset=8.2, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=30, label='AngV-2.1 clade', offset=8.2, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=552, label='AngV-2.2 clade', offset=8.2, fontsize=4, color="black", align=TRUE)+
  #geom_cladelabel(node=553, label='AngV-2.3 clade', offset=8.2, fontsize=4, color="black", align=TRUE)+

  geom_cladelabel(node = 693, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 973, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 999, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 1005, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 678, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 691, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 674, label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  geom_cladelabel(node = 651, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 662, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 672, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 535, label = "Hendra/Nipah henipavirus clade (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 530, label = "Ghanaian bat henipavirus clade | Eidolon helvum | Zambia (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 1013, label = "Unclassified bat paramyxovirus clade | Pteropus sp. | Asia (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  geom_cladelabel(node = 1031, label = "Unclassified bat paramyxovirus clade | Eidolon helvum | Tanzania (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 568, label = "Rubulavirus clade (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  xlim(c(0,14))
p2
p2<-flip(p2, 692, 561)
p2<-flip(p2, 526, 524)

p2.2<-ggtree::collapse(p2, 693)+geom_point2(aes(subset=(node==693)), size=4, shape=25, fill="black")
p2.3<-ggtree::collapse(p2.2, 973)+geom_point2(aes(subset=(node==973)), size=4, shape=25, fill="black")
p2.4<-ggtree::collapse(p2.3, 999)+geom_point2(aes(subset=(node==999)), size=4, shape=25, fill="black")
p2.5<-ggtree::collapse(p2.4, 1005)+geom_point2(aes(subset=(node==1005)), size=4, shape=25, fill="black")
p2.6<-ggtree::collapse(p2.5, 678)+geom_point2(aes(subset=(node==678)), size=4, shape=25, fill="black")
p2.7<-ggtree::collapse(p2.6, 691)+geom_point2(aes(subset=(node==691)), size=4, shape=25, fill="black")
p2.8<-ggtree::collapse(p2.7, 674)+geom_point2(aes(subset=(node==674)), size=4, shape=25, fill="black")
p2.9<-ggtree::collapse(p2.8, 651)+geom_point2(aes(subset=(node==651)), size=4, shape=25, fill="black")
p2.10<-ggtree::collapse(p2.9, 662)+geom_point2(aes(subset=(node==662)), size=4, shape=25, fill="black")
p2.11<-ggtree::collapse(p2.10, 672)+geom_point2(aes(subset=(node==672)), size=4, shape=25, fill="black")
p2.12<-ggtree::collapse(p2.11, 535)+geom_point2(aes(subset=(node==535)), size=4, shape=25, fill="deepskyblue2")
p2.13<-ggtree::collapse(p2.12, 530)+geom_point2(aes(subset=(node==530)), size=4, shape=25, fill="deepskyblue2")
p2.14<-ggtree::collapse(p2.13, 1013)+geom_point2(aes(subset=(node==1013)), size=4, shape=25, fill="deepskyblue2")
p2.15<-ggtree::collapse(p2.14, 1031)+geom_point2(aes(subset=(node==1031)), size=4, shape=25, fill="deepskyblue2")
p2.1<-ggtree::collapse(p2.15, 568)+geom_point2(aes(subset=(node==568)), size=4, shape=25, fill="black")
p2.1

p2.1<- p2.1 +  
  geom_hilight(node=534, fill="#A58AFF", alpha=0.3, extendto=12)+ #HeV/NiV sister
  geom_hilight(node=520, fill="darkslategray1", alpha=0.3, extendto=12)+ #AndV
  geom_hilight(node=1028, fill="dodgerblue1", alpha=0.3, extendto=12)+ #AndV2
  geom_hilight(node=1032, fill="dodgerblue4", alpha=0.3, extendto=12)+ #AndV3.1
  geom_hilight(node=511, fill="dodgerblue4", alpha=0.3, extendto=12)+ #AndV3.2
  geom_hilight(node=528, fill="#00C094", alpha=0.3, extendto=12)+ #GhVlike
  geom_hilight(node=549, fill="tomato1", alpha=0.3, extendto=12)+ #AngV
  geom_hilight(node=30, fill="darkorange1", alpha=0.3, extendto=12) + #AngV2.1
  geom_hilight(node=552, fill="darkorange1", alpha=0.3, extendto=12) #AngV2.2
  #geom_hilight(node=553, fill="darkorange1", alpha=0.3, extendto=12) #AngV2.3
p2.1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree2)
p0.dat1 <- p0$data
p0.dat1$Bootstrap <- NA
Bootstrap<-p0.dat1$Bootstrap[(length(tree.dat2$tip_label)+1):length(p0.dat1$label)] <- as.numeric(p0.dat1$label[(length(tree.dat2$tip_label)+1):length(p0.dat1$label)])#fill with label

#add bootstrap values to original plot
tree2 <- p2.1  %<+% p0.dat1 + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0, size=2)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1))
tree2

homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/final_figures_tables/extended-data/"
ggsave(file = paste0(homewd, "ExtendedDataFig2.png"),
       plot= tree2,
       units="mm",
       width=250,
       height=110,
       scale=2,
       dpi=500)





