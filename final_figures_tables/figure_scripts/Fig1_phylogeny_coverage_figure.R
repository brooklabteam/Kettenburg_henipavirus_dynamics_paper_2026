#Script to make figure 1, which shows the full genome phylogeny, polymerase phylogeny with the coverage plots

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
library(reshape2)
library(gggenes)
library(patchwork)
library(cowplot)
library(ggh4x)

##First plot the phylogenetic tree

##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/phylo_trees/"
#homewd= "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/phylo_trees/"

#full genome tree first
setwd(paste0(homewd,"/full_genome/full_genome_tree_GhVlike"))
tree1 <-  read.tree("fullgenome_align.fasta.treefile")
#drop root for aesthetics
rooted.tree1 <- root(tree1, which(tree1$tip.label == "NC_025345.1"))
rooted.tree1<-drop.tip(rooted.tree1, "NC_025345.1")
#load tree data prepared from elsewhere
dat1 <- read.csv(("fullgenome_metadata.csv"), header = T, stringsAsFactors = F)

#check subgroup names
unique(dat1$Genus)
dat1$novel <- as.factor(dat1$novel)

# colz = c("Orthoavulavirus"="darkorange","Paraavulavirus"="darkorange3","Metaavulavirus"="darkgoldenrod1","Jeilongvirus"="coral",
#          "Morbillivirus"="coral3","Respirovirus"="deeppink1", "Parahenipavirus"="deeppink4",
#          "Henipavirus"="firebrick1","Cynoglossusvirus"="grey50","Pararubulavirus"="firebrick4","Salemvirus"="indianred1",
#          "Ferlavirus"="indianred3","Orthorubulavirus"="indianred4",
#          "Hoplichthysvirus"="lightsalmon1","Narmovirus"="lightsalmon3","Aquaparamyxovirus"="lightpink2","Scoliodonvirus"="lightpink3",
#          "Synodonvirus"="pink1","Unclassified"="black",
#          "Sunshinevirus"="chartreuse1")

colz = c("Orthoavulavirus"="black","Paraavulavirus"="black","Metaavulavirus"="black","Jeilongvirus"="black",
         "Morbillivirus"="black","Respirovirus"="black", "Parahenipavirus"="black",
         "Henipavirus"="deepskyblue2","Cynoglossusvirus"="black","Pararubulavirus"="goldenrod","Salemvirus"="black",
         "Ferlavirus"="black","Orthorubulavirus"="black",
         "Hoplichthysvirus"="black","Narmovirus"="black","Aquaparamyxovirus"="black","Scoliodonvirus"="black",
         "Synodonvirus"="black","Unclassified"="black",
         "Sunshinevirus"="black")


#now get new tip labels
dat1$old_tip_label <- dat1$tip_label
dat1$new_label <- NA

#now, select the name based on what components are present for each sample
#all components with values:
# dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                   dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                   dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
#                                                                                                                   dat1$source[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                   dat1$Country[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                   dat1$Collection_Date[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
# 
# #and if there is an NA just drop it
# 
# #here NA in Isolate only:
# dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
#                                                                                                                  dat1$Organism_Name[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
#                                                                                                                  #dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
#                                                                                                                  dat1$source[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                  dat1$Country[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                  dat1$Collection_Date[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                                  dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                                 dat1$Organism_Name[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])



#make sure to sort in order
colnames(dat1)[colnames(dat1)=="Accession"]="old_tip_label"
tree.dat1 <- data.frame(old_tip_label=rooted.tree1$tip.label, num =1:length(rooted.tree1$tip.label))
head(tree.dat1)
head(dat1)
tree.dat1 <- merge(tree.dat1, dat1, by = "old_tip_label", all.x = T, sort = F)
names(tree.dat1)
tree.dat1$tip_label <- tree.dat1$new_label
tree.dat1 <- dplyr::select(tree.dat1, tip_label, Isolate, Host,Host2, source, Country, Collection_Date, Genus, Organism_Name, novel, old_tip_label)
rooted.tree1$tip.label <- tree.dat1$tip_label
head(tree.dat1)
cbind(tree.dat1$old_tip_label, rooted.tree1$tip.label)
tree.dat1$tip_label#all look good

#shaped and highlighting
tree.dat1$Host2[tree.dat1$Host2==0] <- "Non-bat host"
tree.dat1$Host2[tree.dat1$Host2==1] <- "Bat host"
tree.dat1$Host2 <- as.factor(tree.dat1$Host2)
#shapez = c("Bat host" =  17, "Non-bat host" = 19)
tree.dat1$Host[tree.dat1$Host==0] <- "Non-bat host"
tree.dat1$Host[tree.dat1$Host==1] <- "Eidolon dupreanum"
tree.dat1$Host[tree.dat1$Host==2] <- "Rousettus madagascariensis"
tree.dat1$Host[tree.dat1$Host==3] <- "Other bat host"
tree.dat1$Host[tree.dat1$Host==4] <- "Collapsed clade - multihost"
tree.dat1$Host <- factor(tree.dat1$Host, levels = c("Eidolon dupreanum", "Rousettus madagascariensis", 
                                                    "Other bat host", "Non-bat host", "Collapsed clade - multihost"))
#tree.dat1$Host <- as.factor(tree.dat1$Host)
shapez = c("Non-bat host" =  18, "Eidolon dupreanum" = 17, "Rousettus madagascariensis" = 19, "Other bat host" =  15, "Collapsed clade - multihost"=25)

colz2 = c('1' =  "yellow", '0' = "white")

##Get the clade numbers so we can collapse unnecessary clades
ggtree(rooted.tree1) + geom_text(aes(label=node), hjust=-.3)

#full genome alignment
title<-expression(paste("                    Full genome"))
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=4,stroke=0.5,show.legend = T) +
  guides(color="none")+
  #scale_fill_manual(values=colz) +
  scale_color_manual(values=colz, guide =
                       guide_legend(label.theme = element_text(angle = 0, face = "italic")))+
  scale_shape_manual(values=shapez, guide =
                       guide_legend(label.theme = element_text(angle = 0, face = "italic"), nrow=3)) +
  new_scale_color() +
  geom_tiplab(aes(color = novel, show.legend=F), geom = "text", Genus="Helvetica", size=4, nudge_x=0.05) +
  guides(color="none")+#
  #ggtitle(title)+
  scale_color_manual(values=c("0"="black","1"="red")) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  # geom_cladelabel(node=357, label='HeV/NiV sister clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=425, label='AndV clade', offset=4.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=213, label='GhV-like', offset=4.5, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=113, label='AndV-2 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=350, label='AndV-3 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=250, label='AngV clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=247, label='AngV-2 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=31, label='AchPV2-like virus', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=268, label='ThkPV1-like viruses', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # 
  geom_cladelabel(node = 275, label = "Hendra/Nipah virus clade (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 362, label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 351, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 341, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 338, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 333, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 324, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 394, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 391, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 372, label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 385, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 404, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  geom_cladelabel(node = 380, label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  xlim(c(0,12))

p1

p1.2<-ggtree::collapse(p1, 275)+geom_point2(aes(subset=(node==275)), size=4, shape=25, fill="deepskyblue2")
p1.3<-ggtree::collapse(p1.2, 362)+geom_point2(aes(subset=(node==362)), size=4, shape=25, fill="black")
p1.4<-ggtree::collapse(p1.3, 351)+geom_point2(aes(subset=(node==351)), size=4, shape=25, fill="black")
p1.5<-ggtree::collapse(p1.4, 341)+geom_point2(aes(subset=(node==341)), size=4, shape=25, fill="black")
p1.6<-ggtree::collapse(p1.5, 338)+geom_point2(aes(subset=(node==338)), size=4, shape=25, fill="black")
p1.7<-ggtree::collapse(p1.6, 333)+geom_point2(aes(subset=(node==333)), size=4, shape=25, fill="black")
p1.8<-ggtree::collapse(p1.7, 324)+geom_point2(aes(subset=(node==324)), size=4, shape=25, fill="black")
p1.9<-ggtree::collapse(p1.8, 394)+geom_point2(aes(subset=(node==394)), size=4, shape=25, fill="black")
p1.10<-ggtree::collapse(p1.9, 391)+geom_point2(aes(subset=(node==391)), size=4, shape=25, fill="black")
p1.11<-ggtree::collapse(p1.10, 372)+geom_point2(aes(subset=(node==372)), size=4, shape=25, fill="black")
p1.17<-ggtree::collapse(p1.11, 385)+geom_point2(aes(subset=(node==385)), size=4, shape=25, fill="black")
p1.18<-ggtree::collapse(p1.17, 404)+geom_point2(aes(subset=(node==404)), size=4, shape=25, fill="black")
p1.1<-ggtree::collapse(p1.18, 380)+geom_point2(aes(subset=(node==380)), size=4, shape=25, fill="black")
p1.1

p1.1<- p1.1 +  
  #geom_hilight(node=357, fill="deeppink3", alpha=0.3, extendto=12.8)+ #HeV/NiV sister
  geom_hilight(node=425, fill="darkslategray1", alpha=0.3, extendto=8)+ #AndV
  # geom_hilight(node=113, fill="darkturquoise", alpha=0.3, extendto=12.8)+ #AndV2
  # geom_hilight(node=350, fill="dodgerblue", alpha=0.3, extendto=12.8)+ #AndV3
  # geom_hilight(node=250, fill="springgreen3", alpha=0.3, extendto=12.8)+ #AngV
  # geom_hilight(node=247, fill="chartreuse", alpha=0.3, extendto=12.8)+ #AngV2
  geom_hilight(node=213, fill="#00C094", alpha=0.3, extendto=8) #GhVlike
p1.1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree1)
p0.dat1 <- p0$data
p0.dat1$Bootstrap <- NA
Bootstrap<-p0.dat1$Bootstrap[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)] <- as.numeric(p0.dat1$label[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)])#fill with label

#add bootstrap values to original plot
fulltree <- p1.1  %<+% p0.dat1 + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0, size=2)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  theme(legend.position = "none",
        legend.key.spacing.y = unit(0.001, "cm"),
        legend.spacing.y = unit(0.001, "cm"),
        legend.direction = "horizontal",
        legend.text = element_text(size=7),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1),
        plot.margin = unit(c(0,0,0,0), "cm"))
fulltree

leg<-get_legend(fulltree)

fulltree2<-cowplot::plot_grid(leg,NULL, fulltree, 
                         labels=c("","",""),
                         rel_widths=c(1,1,1), rel_heights = c(0.4,0.001,1),
                         ncol=1,nrow=3,align="hv", axis="l", label_size = 23)+
  theme(plot.background = element_rect(fill = "white", colour = NA)); fulltree2
fulltree2<-as.ggplot(fulltree2)

#then polymerase tree
setwd(paste0(homewd,"/polymerase"))
tree1 <-  read.tree("pmxLmid.treefile")
#drop root for aesthetics
rooted.tree1 <- root(tree1, which(tree1$tip.label == "YP_009094051.1"))
rooted.tree1<-drop.tip(rooted.tree1, "YP_009094051.1")
#load tree data prepared from elsewhere
dat1 <- read.csv(("L_refseq_metadata.csv"), header = T, stringsAsFactors = F)

#check subgroup names
unique(dat1$Genus)
dat1$novel <- as.factor(dat1$novel)

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
dat1$old_tip_label <- dat1$tip_label
dat1$new_label <- NA

#now, select the name based on what components are present for each sample
#all components with values:
# dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                   dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                   dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
#                                                                                                                   dat1$source[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                   dat1$Country[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                   dat1$Collection_Date[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
# 
# #and if there is an NA just drop it
# 
# #here NA in Isolate only:
# dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
#                                                                                                                  dat1$Organism_Name[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
#                                                                                                                  #dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
#                                                                                                                  dat1$source[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                  dat1$Country[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                  dat1$Collection_Date[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                                  dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                                 dat1$Organism_Name[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])



#make sure to sort in order
colnames(dat1)[colnames(dat1)=="Accession"]="old_tip_label"
tree.dat1 <- data.frame(old_tip_label=rooted.tree1$tip.label, num =1:length(rooted.tree1$tip.label))
head(tree.dat1)
head(dat1)
tree.dat1 <- merge(tree.dat1, dat1, by = "old_tip_label", all.x = T, sort = F)
names(tree.dat1)
tree.dat1$tip_label <- tree.dat1$new_label
tree.dat1 <- dplyr::select(tree.dat1, tip_label, Isolate, Host,Host2, source, Country, Collection_Date, Genus, Organism_Name, novel, old_tip_label)
rooted.tree1$tip.label <- tree.dat1$tip_label
head(tree.dat1)
cbind(tree.dat1$old_tip_label, rooted.tree1$tip.label)
tree.dat1$tip_label#all look good

#shaped and highlighting
tree.dat1$Host2[tree.dat1$Host2==0] <- "Non-bat host"
tree.dat1$Host2[tree.dat1$Host2==1] <- "Bat host"
tree.dat1$Host2 <- as.factor(tree.dat1$Host2)
#shapez = c("Bat host" =  17, "Non-bat host" = 19)
tree.dat1$Host[tree.dat1$Host==0] <- "Non-bat host"
tree.dat1$Host[tree.dat1$Host==1] <- "Eidolon dupreanum"
tree.dat1$Host[tree.dat1$Host==2] <- "Rousettus madagascariensis"
tree.dat1$Host[tree.dat1$Host==3] <- "Other bat host"
tree.dat1$Host[tree.dat1$Host==4] <- "Collapsed clade - multihost"
tree.dat1$Host <- factor(tree.dat1$Host, levels = c("Eidolon dupreanum", "Rousettus madagascariensis", 
                                                    "Other bat host", "Non-bat host", "Collapsed clade - multihost"))
#tree.dat1$Host <- as.factor(tree.dat1$Host)
shapez = c("Non-bat host" =  18, "Eidolon dupreanum" = 17, "Rousettus madagascariensis" = 19, "Other bat host" =  15, "Collapsed clade - multihost"=25)

colz2 = c('1' =  "yellow", '0' = "white")

##Get the clade numbers so we can collapse unnecessary clades
ggtree(rooted.tree1) + geom_text(aes(label=node), hjust=-.3)

#Refseq alignment
title<-expression(paste("                    Polymerase"))
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0.5,show.legend = T) +
  #scale_fill_manual(values=colz) +
  guides(color="none")+
  scale_color_manual(values=colz, guide =
                       guide_legend(label.theme = element_text(angle = 0, face = "italic")))+
  scale_shape_manual(values=shapez, guide =
                       guide_legend(label.theme = element_text(angle = 0, face = "italic"), nrow=3)) +
  new_scale_color() +
  geom_tiplab(aes(color = novel, show.legend=F), geom = "text", Genus="Helvetica", size=4, nudge_x=0.05) +
  guides(color="none")+#
  #ggtitle(title)+
  scale_color_manual(values=c("0"="black","1"="red")) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  geom_cladelabel(node=185, label='HeV/NiV sister clade', offset=4.6, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=176, label='AndV clade', offset=4.6, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=343, label='AndV-2 clade', offset=4.6, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=344, label='AndV-3 clade', offset=4.6, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=190, label='AngV clade', offset=4.6, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=191, label='AngV-2 clade', offset=4.6, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=182, label='GhV clade', offset=4.6, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=26, label='AchPV2-like virus', offset=4.6, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=205, label='ThkPV1-like viruses', offset=4.6, fontsize=4, color="black", align=TRUE)+
  # 
  #geom_cladelabel(node = , label = "Hendra/Nipah virus clade (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  #geom_cladelabel(node = , label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 235, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 225, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 245, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 260, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 251, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 284, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 264, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 214, label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 270, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 334, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  geom_cladelabel(node = 247, label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  xlim(c(0,12))

p1

p1<-flip(p1, 194, 261)
p1<-flip(p1, 178, 180)

#p1.2<-ggtree::collapse(p1, 236)+geom_point2(aes(subset=(node==236)), size=4, shape=25, fill="deepskyblue2")
#p1.3<-ggtree::collapse(p1.2, 344)+geom_point2(aes(subset=(node==344)), size=4, shape=25, fill="firebrick4")
p1.4<-ggtree::collapse(p1, 235)+geom_point2(aes(subset=(node==235)), size=4, shape=25, fill="black")
p1.5<-ggtree::collapse(p1.4, 225)+geom_point2(aes(subset=(node==225)), size=4, shape=25, fill="black")
p1.6<-ggtree::collapse(p1.5, 245)+geom_point2(aes(subset=(node==245)), size=4, shape=25, fill="black")
p1.7<-ggtree::collapse(p1.6, 260)+geom_point2(aes(subset=(node==260)), size=4, shape=25, fill="black")
p1.8<-ggtree::collapse(p1.7, 251)+geom_point2(aes(subset=(node==251)), size=4, shape=25, fill="black")
p1.9<-ggtree::collapse(p1.8, 284)+geom_point2(aes(subset=(node==284)), size=4, shape=25, fill="black")
p1.10<-ggtree::collapse(p1.9, 264)+geom_point2(aes(subset=(node==264)), size=4, shape=25, fill="black")
p1.11<-ggtree::collapse(p1.10, 214)+geom_point2(aes(subset=(node==214)), size=4, shape=25, fill="black")
p1.17<-ggtree::collapse(p1.11, 270)+geom_point2(aes(subset=(node==270)), size=4, shape=25, fill="black")
p1.18<-ggtree::collapse(p1.17, 334)+geom_point2(aes(subset=(node==334)), size=4, shape=25, fill="black")
p1.1<-ggtree::collapse(p1.18, 247)+geom_point2(aes(subset=(node==247)), size=4, shape=25, fill="black")
p1.1

p1.1<- p1.1 +  
  geom_hilight(node=185, fill="#A58AFF", alpha=0.3, extendto=9)+ #HeV/NiV sister
  geom_hilight(node=176, fill="darkslategray1", alpha=0.3, extendto=9)+ #AndV
  geom_hilight(node=343, fill="dodgerblue1", alpha=0.3, extendto=9)+ #AndV2
  geom_hilight(node=344, fill="dodgerblue4", alpha=0.3, extendto=9)+ #AndV3
  geom_hilight(node=190, fill="#F8766D", alpha=0.3, extendto=9)+ #AngV
  geom_hilight(node=191, fill="darkorange1", alpha=0.3, extendto=9)+ #AngV2
  geom_hilight(node=182, fill="#00C094", alpha=0.3, extendto=9)+ #GhVlike
  geom_hilight(node=26, fill="grey40", alpha=0.3, extendto=9)+ #AchPV2
  geom_hilight(node=205, fill="grey40", alpha=0.3, extendto=9) #ThkPV1
p1.1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree1)
p0.dat1 <- p0$data
p0.dat1$Bootstrap <- NA
Bootstrap<-p0.dat1$Bootstrap[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)] <- as.numeric(p0.dat1$label[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)])#fill with label

#add bootstrap values to original plot
polytree <- p1.1  %<+% p0.dat1 + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0, size=2)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1),
        plot.margin = unit(c(0.5,0,0,0.5), "cm"))
polytree


##Now do the coverage plots
#set directory
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/genome_characterization/full_genome_coverage"
setwd(paste0(homewd))

#GhV-like
dat <- read.csv(file= "KEL86_KEL94_coverage_noUTR.csv", header = T, stringsAsFactors = F) #load data

dat_text <- data.frame(
    label = c("GhV-like"),
  x     = c(Inf),
  y     = c(Inf))

coverage_ghvlike<-  ggplot()+ 
  annotate("rect", xmin=50, xmax=4179, ymin=0, ymax=550, alpha=0.6,  fill="#D29292")+
  annotate("rect", xmin=4603, xmax=7636, ymin=0, ymax=550, alpha=0.6,  fill="#525832")+
  annotate("rect", xmin=8096, xmax=10940, ymin=0, ymax=550, alpha=0.6,  fill="#F98F5B")+
  annotate("rect", xmin=11664, xmax=16846, ymin=0, ymax=550, alpha=0.6,  fill="#BD9A5A")+
  annotate("rect", xmin=17373, xmax=23508, ymin=0, ymax=550, alpha=0.6,  fill="#F17983")+
  annotate("rect", xmin=24156, xmax=48178, ymin=0, ymax=550, alpha=0.6,  fill="#6B8C81")+
  geom_area(data=dat, aes(x=Position, y=rPM), fill="black") +
  #geom_line(data=dat, aes(x=Position, y=rPM), size=1) +
  theme_linedraw() +
  theme(legend.position = "none")+
  #ggtitle(title)+
  scale_y_continuous(position="left", expand=c(0,0))+
  geom_hline(yintercept = 104.2684, linetype = "dashed", color="red")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  scale_x_continuous(breaks=c(0,5578,10396,17208,24521, 33205, 44590),
                     labels = c(0,3000,6000,9000,12000, 15000, 17000),expand=c(0,0))+
  #scale_y_continuous(expand=c(0,0))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.title.y = element_text(margin = margin(0,0,0,1)),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        axis.title.x = element_blank(),
        plot.title = element_text(color = "darkorange1"),
        axis.line.x = element_blank())
coverage_ghvlike

coverage_ghvlike<-coverage_ghvlike+geom_text(data= dat_text,
                           mapping = aes(x = x, y = y, label = label, hjust=1.1, vjust=1.5, fontface="italic"));coverage_ghvlike


#ROU148
dat <- read.csv(file= "ROU148_coverage.csv", header = T, stringsAsFactors = F) #load data

dat_text <- data.frame(
  label = c("AndV PX127343"),
  x     = c(Inf),
  y     = c(Inf))

title<-expression(paste(italic("PX127343 AndV genome coverage")))
coverage_rou148<-  ggplot()+ 
  annotate("rect", xmin=156, xmax=1763, ymin=0, ymax=174, alpha=0.6,  fill="#D29292")+
  annotate("rect", xmin=2159, xmax=3868, ymin=0, ymax=174, alpha=0.6,  fill="#525832")+
  annotate("rect", xmin=4329, xmax=5360, ymin=0, ymax=174, alpha=0.6,  fill="#F98F5B")+
  annotate("rect", xmin=5595, xmax=7280, ymin=0, ymax=174, alpha=0.6,  fill="#BD9A5A")+
  annotate("rect", xmin=7782, xmax=10361, ymin=0, ymax=174, alpha=0.6,  fill="#F17983")+
  annotate("rect", xmin=10794, xmax=17546, ymin=0, ymax=174, alpha=0.6,  fill="#6B8C81")+
  geom_area(data=dat, aes(x=Position, y=rPM), fill="black") +
  theme_linedraw() +
  theme(legend.position = "none")+
  #ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  geom_hline(yintercept = 49.40666, linetype = "dashed", color="red")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        plot.title = element_text(color = "darkcyan"),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())
coverage_rou148
coverage_rou148<-coverage_rou148+geom_text(data= dat_text,
                                             mapping = aes(x = x, y = y, label = label, hjust=1.1, vjust=1.5, fontface="italic"));coverage_rou148


#ROU173 and 125 together
dat <- read.csv(file= "ROU125_coverage.csv", header = T, stringsAsFactors = F) #load data
dat2 <- read.csv(file= "ROU173_coverage.csv", header = T, stringsAsFactors = F) #load data

dat_text <- data.frame(
  label = c("AndV PX127342 (black) \nAndV PX127344 (grey)"),
  x     = c(Inf),
  y     = c(Inf))

title<-expression(paste(italic("PX127342 (black) and PX127344 (grey) AndV genome coverage")))
coverage_andv<-  ggplot()+ 
  annotate("rect", xmin=111, xmax=1718, ymin=0, ymax=1.6, alpha=0.6,  fill="#D29292")+
  annotate("rect", xmin=2114, xmax=3823, ymin=0, ymax=1.6, alpha=0.6,  fill="#525832")+
  annotate("rect", xmin=4284, xmax=5315, ymin=0, ymax=1.6, alpha=0.6,  fill="#F98F5B")+
  annotate("rect", xmin=5550, xmax=7235, ymin=0, ymax=1.6, alpha=0.6,  fill="#BD9A5A")+
  annotate("rect", xmin=7737, xmax=10316, ymin=0, ymax=1.6, alpha=0.6,  fill="#F17983")+
  annotate("rect", xmin=10749, xmax=17501, ymin=0, ymax=1.6, alpha=0.6,  fill="#6B8C81")+
  geom_area(data=dat, aes(x=Position, y=rPM), fill="black") +
  geom_area(data=dat2, aes(x=Position, y=rPM), fill="grey50") +
  theme_linedraw() +
  theme(legend.position = "none")+
  #ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  geom_hline(yintercept = 0.5368673, linetype = "dashed", color="red")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        #axis.title.y = element_text(margin = margin(0,0,0,1)),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(color = "darkcyan"),
        axis.line.x = element_blank())
coverage_andv

coverage_andv<-coverage_andv+geom_text(data= dat_text,
                                           mapping = aes(x = x, y = y, label = label, hjust=1.1, vjust=1.5, fontface="italic"));coverage_andv


#and a version with lines instead of area plots
# dat <- read.csv(file= "AndV_merged_coverage.csv", header = T, stringsAsFactors = F) #load data
# 
# 
# 
# colzpalette<-c("black","darkslategray4","darkslategray")
# 
# title<-expression(paste(italic("AndV genome coverage")))
# coverage_andv2<-  ggplot()+ 
#   annotate("rect", xmin=111, xmax=1718, ymin=0, ymax=180, alpha=0.6,  fill="#D29292")+
#   annotate("rect", xmin=2114, xmax=3823, ymin=0, ymax=180, alpha=0.6,  fill="#525832")+
#   annotate("rect", xmin=4284, xmax=5315, ymin=0, ymax=180, alpha=0.6,  fill="#F98F5B")+
#   annotate("rect", xmin=5550, xmax=7235, ymin=0, ymax=180, alpha=0.6,  fill="#BD9A5A")+
#   annotate("rect", xmin=7737, xmax=10316, ymin=0, ymax=180, alpha=0.6,  fill="#F17983")+
#   annotate("rect", xmin=10749, xmax=17501, ymin=0, ymax=180, alpha=0.6,  fill="#6B8C81")+
#   #geom_area(data=dat, aes(x=Position, y=rPM, color=virus)) +
#   geom_line(data=dat, aes(x=Position, y=rPM, color=virus), size=1) +
#   theme_linedraw() +
#   theme(legend.position = "none")+
#   ggtitle(title)+
#   xlab("Genome position") + ylab("rPM (reads per million)")+
#   scale_x_continuous(expand=c(0,0))+
#   scale_y_break(c(1.6, 10), scales = 1.5, expand=c(0,0))+
#   geom_hline(yintercept = 0.5368673, linetype = "dashed", color="red")+
#   scale_color_manual(values=colzpalette) + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         legend.direction = "horizontal",
#         legend.margin=margin(),
#         legend.justification = "left",
#         legend.text = element_text(face="italic", size = 10),
#         legend.title = element_blank(),
#         legend.position="top", 
#         legend.key.height= unit(3.5, 'mm'),
#         legend.key.width= unit(3.5, 'mm'),
#         legend.background =element_rect(fill = alpha("white", 0)),
#         axis.title.y = element_text(margin = margin(0,0,0,1)),
#         axis.ticks.x=element_blank(),
#         plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank(),
#         plot.title = element_text(color = "darkcyan"),
#         axis.line.x = element_blank())
# coverage_andv2


##Now make genome maps to go with the coverage plots

#GhV-like first
map <- read.csv("genome_genes.csv", header = T, stringsAsFactors = F)
map$gene<-factor(map$gene, levels = c("N","P","M","F","G","L","gap"))
#Load the peptide files
map_pep <- read.csv("genome_peptides.csv", header = T, stringsAsFactors = F)
map_pep$gene<-factor(map_pep$gene, levels = c("N","P","M","F","G","L","gap"))
map_pep$string<-factor(map_pep$string, levels = c("first","second"))
#Load the feature file in case its needed
map_feat <- read.csv("genome_features.csv", header = T, stringsAsFactors = F)
map_feat$gene<-factor(map_feat$gene, levels = c("N","P","M","F","G","L","gap"))

colz=c("N"="#D29292","P"="#525832","M"="#F98F5B","F"="#BD9A5A","G"="#F17983","L"="#6B8C81","gap"="#0E171E","gap2"="#0E171E")
colz2=c("N"="white","P"="white","M"="white","F"="white","G"="white","L"="white","gap"="#0E171E")

#subset by the sequence you want
map_combo<-subset(map,type=="KEL86_KEL94")
map_combo_pep<-subset(map_pep,type=="KEL86_KEL94")
map_combo_pep1<-subset(map_combo_pep,string=="first")
map_combo_feat<-subset(map_feat,type=="KEL86_KEL94")
map_combo_pep2<-subset(map_combo_pep,string=="second")

#map the genomes
title<-expression(paste(italic("GhV-like")))
genome1ghvlike<-ggplot(map_combo_pep1, aes(xmin = from, xmax = to, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(0, "mm"),
                  arrowhead_height = grid::unit(5, "mm"),
                  arrow_body_height = grid::unit(5, "mm")
  ) +
  scale_fill_manual(values=colz2)+
  #ggtitle(title)+
  theme_genes()+
  scale_x_continuous(breaks=c(0,5578,10396,17208,24521, 33205, 44590),
                     labels = c(0,3000,6000,9000,12000, 15000, 17000),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
genome1ghvlike

genome2ghvlike<-ggplot(map_combo, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(5, "mm"),
                  arrow_body_height = grid::unit(5, "mm")) +
  geom_subgene_arrow(data = map_combo_pep2, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                        xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(5, "mm"),
                     arrow_body_height = grid::unit(5, "mm"))+
  geom_text(data=map_combo_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(breaks=c(0,5578,10396,17208,24521, 33205, 44590),
                     labels = c(0,3000,6000,9000,12000, 15000, 17000),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank()
        #axis.ticks.x=element_blank(),
        #axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        #axis.line.x = element_blank()
  )+
  xlab("Genome position") + ylab("")
genome2ghvlike

#put coverage and genome maps together
coverageplot_ghvlike<-coverage_ghvlike/genome1ghvlike/genome2ghvlike+plot_layout(nrow=3,  heights = c(1, 0.2, 0.4))
coverageplot_ghvlike

coverage_ghvlike2<-genome1ghvlike/genome2ghvlike+plot_layout(nrow=2,  heights = c(1,1.2))
coverage_ghvlike2


coverage_ghvlike2<-as.ggplot(coverage_ghvlike2)
coverage_ghvlike2

#ROU148
map <- read.csv("genome_genes.csv", header = T, stringsAsFactors = F)
map$gene<-factor(map$gene, levels = c("N","P","M","F","G","L","gap"))
#Load the peptide files
map_pep <- read.csv("genome_peptides.csv", header = T, stringsAsFactors = F)
map_pep$gene<-factor(map_pep$gene, levels = c("N","P","M","F","G","L","gap"))
map_pep$string<-factor(map_pep$string, levels = c("first","second"))
#Load the feature file in case its needed
map_feat <- read.csv("genome_features.csv", header = T, stringsAsFactors = F)
map_feat$gene<-factor(map_feat$gene, levels = c("N","P","M","F","G","L","gap"))

colz=c("N"="#D29292","P"="#525832","M"="#F98F5B","F"="#BD9A5A","G"="#F17983","L"="#6B8C81","gap"="#0E171E","gap2"="#0E171E")
colz2=c("N"="white","P"="white","M"="white","F"="white","G"="white","L"="white","gap"="#0E171E")

#subset by the sequence you want
map_combo<-subset(map,type=="ROU148")
map_combo_pep<-subset(map_pep,type=="ROU148")
map_combo_pep1<-subset(map_combo_pep,string=="first")
map_combo_feat<-subset(map_feat,type=="ROU148")
map_combo_pep2<-subset(map_combo_pep,string=="second")

#map the genomes
title<-expression(paste(italic("AndV PX127342 and AndV PX127343")))
genome1ROU148<-ggplot(map_combo_pep1, aes(xmin = from, xmax = to, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(0, "mm"),
                  arrowhead_height = grid::unit(5, "mm"),
                  arrow_body_height = grid::unit(5, "mm")
  ) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
genome1ROU148

genome2ROU148<-ggplot(map_combo, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(5, "mm"),
                  arrow_body_height = grid::unit(5, "mm")) +
  geom_subgene_arrow(data = map_combo_pep2, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                        xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(5, "mm"),
                     arrow_body_height = grid::unit(5, "mm"))+
  geom_text(data=map_combo_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank()
  )+
  xlab("Genome position") + ylab("")
genome2ROU148

#put coverage and genome maps together
genomeplot_ROU148<-coverage_rou148/genome1ROU148/genome2ROU148+plot_layout(nrow=3,  heights = c(1, 0.2, 0.4))
genomeplot_ROU148

genomeplot_ROU148.2<-genome1ROU148/genome2ROU148+plot_layout(nrow=2,  heights = c(1,1.2))
genomeplot_ROU148.2

genomeplot_ROU148.2<-as.ggplot(genomeplot_ROU148.2)
genomeplot_ROU148.2


#rou125 and rou148 ref
map <- read.csv("genome_genes.csv", header = T, stringsAsFactors = F)
map$gene<-factor(map$gene, levels = c("N","P","M","F","G","L","gap"))
#Load the peptide files
map_pep <- read.csv("genome_peptides.csv", header = T, stringsAsFactors = F)
map_pep$gene<-factor(map_pep$gene, levels = c("N","P","M","F","G","L","gap"))
map_pep$string<-factor(map_pep$string, levels = c("first","second"))
#Load the feature file in case its needed
map_feat <- read.csv("genome_features.csv", header = T, stringsAsFactors = F)
map_feat$gene<-factor(map_feat$gene, levels = c("N","P","M","F","G","L","gap"))

colz=c("N"="#D29292","P"="#525832","M"="#F98F5B","F"="#BD9A5A","G"="#F17983","L"="#6B8C81","gap"="#0E171E","gap2"="#0E171E")
colz2=c("N"="white","P"="white","M"="white","F"="white","G"="white","L"="white","gap"="#0E171E")

#subset by the sequence you want
map_combo<-subset(map,type=="ROU125")
map_combo_pep<-subset(map_pep,type=="ROU125")
map_combo_pep1<-subset(map_combo_pep,string=="first")
map_combo_feat<-subset(map_feat,type=="ROU125")
map_combo_pep2<-subset(map_combo_pep,string=="second")

#map the genomes
#title<-expression(paste(italic("PX127342 & PX127343")))
genome1andv<-ggplot(map_combo_pep1, aes(xmin = from, xmax = to, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(0, "mm"),
                  arrowhead_height = grid::unit(5, "mm"),
                  arrow_body_height = grid::unit(5, "mm")
  ) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  #ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
genome1andv

genome2andv<-ggplot(map_combo, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(5, "mm"),
                  arrow_body_height = grid::unit(5, "mm")) +
  geom_subgene_arrow(data = map_combo_pep2, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                        xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(5, "mm"),
                     arrow_body_height = grid::unit(5, "mm"))+
  geom_text(data=map_combo_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank()
  )+
  xlab("") + ylab("")
genome2andv


#ROU173
map <- read.csv("genome_genes.csv", header = T, stringsAsFactors = F)
map$gene<-factor(map$gene, levels = c("N","P","M","F","G","L","gap"))
#Load the peptide files
map_pep <- read.csv("genome_peptides.csv", header = T, stringsAsFactors = F)
map_pep$gene<-factor(map_pep$gene, levels = c("N","P","M","F","G","L","gap"))
map_pep$string<-factor(map_pep$string, levels = c("first","second"))
#Load the feature file in case its needed
map_feat <- read.csv("genome_features.csv", header = T, stringsAsFactors = F)
map_feat$gene<-factor(map_feat$gene, levels = c("N","P","M","F","G","L","gap"))

colz=c("N"="#D29292","P"="#525832","M"="#F98F5B","F"="#BD9A5A","G"="#F17983","L"="#6B8C81","gap"="#0E171E","gap2"="#0E171E")
colz2=c("N"="white","P"="white","M"="white","F"="white","G"="white","L"="white","gap"="#0E171E")

#subset by the sequence you want
map_combo<-subset(map,type=="ROU173")
map_combo_pep<-subset(map_pep,type=="ROU173")
map_combo_pep1<-subset(map_combo_pep,string=="first")
map_combo_feat<-subset(map_feat,type=="ROU173")
map_combo_pep2<-subset(map_combo_pep,string=="second")

#map the genomes
title<-expression(paste(italic("AndV PX127344")))
genome1ROU173<-ggplot(map_combo_pep1, aes(xmin = from, xmax = to, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(0, "mm"),
                  arrowhead_height = grid::unit(5, "mm"),
                  arrow_body_height = grid::unit(5, "mm")
  ) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
genome1ROU173

genome2ROU173<-ggplot(map_combo, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(5, "mm"),
                  arrow_body_height = grid::unit(5, "mm")) +
  geom_subgene_arrow(data = map_combo_pep2, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                        xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(5, "mm"),
                     arrow_body_height = grid::unit(5, "mm"))+
  geom_text(data=map_combo_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        # axis.ticks.x=element_blank(),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        # axis.line.x = element_blank()
  )+
  xlab("Genome position") + ylab("")
genome2ROU173

#put coverage and genome maps together
genomeplot_ROU173<-coverage_andv/genome1ROU173/genome2ROU173+plot_layout(nrow=3,  heights = c(1, 0.2, 0.2))
genomeplot_ROU173

genomeplot_ROU173.2<-genome1ROU173/genome2ROU173+plot_layout(nrow=2,  heights = c(1, 1.2))
genomeplot_ROU173.2

genomeplot_ROU173.2<-as.ggplot(genomeplot_ROU173.2)
genomeplot_ROU173.2

genomeplot_andv<-coverage_rou148/genome1ROU148/genome2ROU148/genome1ROU173/genome2ROU173+plot_layout(nrow=5,  heights = c(1, 0.2, 0.4, 0.2, 0.4))
genomeplot_andv



#simplot for GhVlike looks meh, maybe make a pairwise comparison plot instead
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/genome_characterization"
setwd(paste0(homewd))

dat <- read.csv(file= "pairwise_comparisons_heatmap_aa.csv", header = T, stringsAsFactors = F) #load data
head(dat)
names(dat)

dat$gene<-factor(dat$gene, levels=c("N", "P","M","F","G","L"))
dat$var1_alt<-factor(dat$var1_alt, levels=c("GhV", "YuV2"))
dat$var2<-factor(dat$var2, levels=c("GhV-like identity to GhV", "AndV identity to YuV2"))

#plot comparison to known sequences
compare_plot <- ggplot(dat, aes(x=var1_alt, y=gene, fill=id)) +
  geom_tile() +
  geom_text(aes(label = id), color = "black", size = 4) +
  # scale_fill_viridis_c(option = "C", direction = -1) +
  scale_fill_gradient(low="yellow", high="red")+
  facet_nested(.~var2, scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "Gene",
    fill="Amino acid \nidentity",
    title="")+
  theme_linedraw()+
  scale_y_discrete(position="left", limits=rev)+
  scale_x_discrete(position="bottom")+
  theme( 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=7),
        legend.title = element_text(size=10),
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        plot.margin = margin(0, 5, 5, 5, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12), 
        legend.position = "right",
        strip.text = element_text(size=12))
compare_plot


#full comparison

part1<-cowplot::plot_grid(compare_plot,NULL,fulltree2, 
                          labels=c("B","","C"),
                          rel_widths=c(1,1,1), rel_heights = c(1,0.01,2.3),
                          ncol=1,nrow=3,align="hv", axis="l", label_size = 23)+
  theme(plot.background = element_rect(fill = "white", colour = NA)); part1
part1<-as.ggplot(part1)

part2<-cowplot::plot_grid(part1, polytree, 
                            labels=c("","D"),
                            rel_widths=c(1,1), rel_heights = c(1.2,1),
                            ncol=2,nrow=1,align="hv", axis="l", label_size = 23)+
  theme(plot.background = element_rect(fill = "white", colour = NA)); part2


part3<-cowplot::plot_grid(coverageplot_ghvlike, genomeplot_andv,
                        labels=c("A",""),
                        rel_widths=c(1,1), rel_heights = c(1,1.4),
                        ncol=1,nrow=2,align="hv", axis="l", label_size = 23)+
  theme(plot.background = element_rect(fill = "white", colour = NA)); part3

# part4<-cowplot::plot_grid(coverage_ghvlike2, genomeplot_ROU148.2, genomeplot_ROU173.2,
#                           labels=c("","",""),
#                           rel_widths=c(1,1,1), rel_heights = c(1,1,1),
#                           ncol=1,nrow=3,align="hv", axis="r", label_size = 23)+
#   theme(plot.background = element_rect(fill = "white", colour = NA)); part4
# 
# part5<-cowplot::plot_grid(NULL, part4, NULL,
#                           labels=c("","",""),
#                           rel_widths=c(0.05,1,0.03), rel_heights = c(1,1,1),
#                           ncol=3,nrow=1,align="hv", axis="lr", label_size = 23)+
#   theme(plot.background = element_rect(fill = "white", colour = NA)); part5
# 
# part6<-cowplot::plot_grid(part3, part5,
#                           labels=c("",""),
#                           rel_widths=c(1,1), rel_heights = c(1,1.5),
#                           ncol=1,nrow=2,align="hv", axis="lr", label_size = 23)+
#   theme(plot.background = element_rect(fill = "white", colour = NA)); part6
# part6<-as.ggplot(part6)


fig1<-cowplot::plot_grid(part3, part2, 
                labels=c("",""),
                rel_widths=c(1,2), rel_heights = c(1,1),
                ncol=2,nrow=1,align="hv", axis="lr", label_size = 23)+
  theme(plot.background = element_rect(fill = "white", colour = NA)); fig1

#save
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/final_figures_tables/tmp/"
ggsave(file = paste0(homewd, "Fig1_phylogeny_coverage_pairwise.png"),
       plot= fig1,
       units="mm",
       width=290,
       height=110,
       scale=2,
       dpi=500)



