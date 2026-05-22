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

##First plot the phylogenetic tree

##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/phylo_trees/with_GhVlike/"
#homewd= "/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/phylo_trees/"

#full genome tree first
setwd(paste0(homewd,"/full_genome_tree"))
tree1 <-  read.tree("fullgenome_align.fasta.treefile")
#drop root for aesthetics
rooted.tree1 <- root(tree1, which(tree1$tip.label == "NC_025345.1"))
rooted.tree1<-drop.tip(rooted.tree1, "NC_025345.1")
#load tree data prepared from elsewhere
dat1 <- read.csv(("fullgenome_metadata.csv"), header = T, stringsAsFactors = F)

#check subgroup names
unique(dat1$Genus)
dat1$novel <- as.factor(dat1$novel)

colz = c("Orthoavulavirus"="darkorange","Paraavulavirus"="darkorange3","Metaavulavirus"="darkgoldenrod1","Jeilongvirus"="coral",
         "Morbillivirus"="coral3","Respirovirus"="deeppink1", "Parahenipavirus"="deeppink4",
         "Henipavirus"="firebrick1","Cynoglossusvirus"="grey50","Pararubulavirus"="firebrick4","Salemvirus"="indianred1",
         "Ferlavirus"="indianred3","Orthorubulavirus"="indianred4",
         "Hoplichthysvirus"="lightsalmon1","Narmovirus"="lightsalmon3","Aquaparamyxovirus"="lightpink2","Scoliodonvirus"="lightpink3",
         "Synodonvirus"="pink1","Unclassified"="black",
         "Sunshinevirus"="chartreuse1")

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
tree.dat1$Host[tree.dat1$Host==0] <- "Other host"
tree.dat1$Host[tree.dat1$Host==1] <- "Eidolon dupreanum"
tree.dat1$Host[tree.dat1$Host==2] <- "Rousettus madagascariensis"
tree.dat1$Host[tree.dat1$Host==3] <- "Bat host"
tree.dat1$Host <- as.factor(tree.dat1$Host)
shapez = c("Other host" =  15, "Eidolon dupreanum" = 17, "Rousettus madagascariensis" = 19, "Bat host" =  18)

colz2 = c('1' =  "yellow", '0' = "white")

##Get the clade numbers so we can collapse unnecessary clades
ggtree(rooted.tree1) + geom_text(aes(label=node), hjust=-.3)

#full genome alignment
title<-expression(paste("          Full genome"))
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=4,stroke=0,show.legend = T) +
  #scale_fill_manual(values=colz) +
  scale_color_manual(values=colz, guide =
                       guide_legend(label.theme = element_text(angle = 0, face = "italic")))+
  scale_shape_manual(values=shapez, guide =
                       guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
  new_scale_color() +
  geom_tiplab(aes(color = novel, show.legend=F), geom = "text", Genus="Helvetica", size=4, nudge_x=0.05) +
  guides(color="none")+#
  ggtitle(title)+
  scale_color_manual(values=c("0"="black","1"="red")) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm")) +
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
  geom_cladelabel(node = 275, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 362, label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 351, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 341, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 338, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 333, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 324, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 394, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 391, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 372, label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 385, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 404, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 380, label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  xlim(c(0,12))

p1

p1.2<-ggtree::collapse(p1, 275)+geom_point2(aes(subset=(node==275)), size=4, shape=22, fill="firebrick1")
p1.3<-ggtree::collapse(p1.2, 362)+geom_point2(aes(subset=(node==362)), size=4, shape=22, fill="firebrick4")
p1.4<-ggtree::collapse(p1.3, 351)+geom_point2(aes(subset=(node==351)), size=4, shape=22, fill="darkgoldenrod1")
p1.5<-ggtree::collapse(p1.4, 341)+geom_point2(aes(subset=(node==341)), size=4, shape=22, fill="darkorange")
p1.6<-ggtree::collapse(p1.5, 338)+geom_point2(aes(subset=(node==338)), size=4, shape=22, fill="darkorange3")
p1.7<-ggtree::collapse(p1.6, 333)+geom_point2(aes(subset=(node==333)), size=4, shape=22, fill="lightpink2")
p1.8<-ggtree::collapse(p1.7, 324)+geom_point2(aes(subset=(node==324)), size=4, shape=22, fill="deeppink1")
p1.9<-ggtree::collapse(p1.8, 394)+geom_point2(aes(subset=(node==394)), size=4, shape=22, fill="coral")
p1.10<-ggtree::collapse(p1.9, 391)+geom_point2(aes(subset=(node==391)), size=4, shape=22, fill="lightsalmon3")
p1.11<-ggtree::collapse(p1.10, 372)+geom_point2(aes(subset=(node==372)), size=4, shape=22, fill="indianred3")
p1.17<-ggtree::collapse(p1.11, 385)+geom_point2(aes(subset=(node==385)), size=4, shape=22, fill="coral3")
p1.18<-ggtree::collapse(p1.17, 404)+geom_point2(aes(subset=(node==404)), size=4, shape=22, fill="deeppink4")
p1.1<-ggtree::collapse(p1.18, 380)+geom_point2(aes(subset=(node==380)), size=4, shape=22, fill="lightpink3")
p1.1

p1.1<- p1.1 +  
  #geom_hilight(node=357, fill="deeppink3", alpha=0.3, extendto=12.8)+ #HeV/NiV sister
  geom_hilight(node=425, fill="darkcyan", alpha=0.3, extendto=8)+ #AndV
  # geom_hilight(node=113, fill="darkturquoise", alpha=0.3, extendto=12.8)+ #AndV2
  # geom_hilight(node=350, fill="dodgerblue", alpha=0.3, extendto=12.8)+ #AndV3
  # geom_hilight(node=250, fill="springgreen3", alpha=0.3, extendto=12.8)+ #AngV
  # geom_hilight(node=247, fill="chartreuse", alpha=0.3, extendto=12.8)+ #AngV2
  geom_hilight(node=213, fill="orange", alpha=0.3, extendto=8) #GhVlike
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
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1),
        plot.margin = unit(c(0,0,0,2), "cm"))
fulltree


#then polymerase tree
setwd(paste0(homewd,"/pmx_gene_trees/polymerase"))
tree1 <-  read.tree("polymerase_mid_align.fasta.treefile")
#drop root for aesthetics
rooted.tree1 <- root(tree1, which(tree1$tip.label == "YP_009094051.1"))
rooted.tree1<-drop.tip(rooted.tree1, "YP_009094051.1")
#load tree data prepared from elsewhere
dat1 <- read.csv(("L_metadata_mid.csv"), header = T, stringsAsFactors = F)

#check subgroup names
unique(dat1$Genus)
dat1$novel <- as.factor(dat1$novel)

colz = c("Orthoavulavirus"="darkorange","Paraavulavirus"="darkorange3","Metaavulavirus"="darkgoldenrod1","Jeilongvirus"="coral",
         "Morbillivirus"="coral3","Respirovirus"="deeppink1", "Parahenipavirus"="deeppink4",
         "Henipavirus"="firebrick1","Cynoglossusvirus"="grey50","Pararubulavirus"="firebrick4","Salemvirus"="indianred1",
         "Ferlavirus"="indianred3","Orthorubulavirus"="indianred4",
         "Hoplichthysvirus"="lightsalmon1","Narmovirus"="lightsalmon3","Aquaparamyxovirus"="lightpink2","Scoliodonvirus"="lightpink3",
         "Synodonvirus"="pink1","Unclassified"="black",
         "Sunshinevirus"="chartreuse1")

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
tree.dat1$Host[tree.dat1$Host==0] <- "Other host"
tree.dat1$Host[tree.dat1$Host==1] <- "Eidolon dupreanum"
tree.dat1$Host[tree.dat1$Host==2] <- "Rousettus madagascariensis"
tree.dat1$Host[tree.dat1$Host==3] <- "Bat host"
tree.dat1$Host <- as.factor(tree.dat1$Host)
shapez = c("Other host" =  15, "Eidolon dupreanum" = 17, "Rousettus madagascariensis" = 19, "Bat host" =  18)

colz2 = c('1' =  "yellow", '0' = "white")

##Get the clade numbers so we can collapse unnecessary clades
ggtree(rooted.tree1) + geom_text(aes(label=node), hjust=-.3)

#Refseq alignment
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +
  #scale_fill_manual(values=colz) +
  scale_color_manual(values=colz, guide =
                       guide_legend(label.theme = element_text(angle = 0, face = "italic")))+
  scale_shape_manual(values=shapez, guide =
                       guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
  new_scale_color() +
  geom_tiplab(aes(color = novel, show.legend=F), geom = "text", Genus="Helvetica", size=4, nudge_x=0.05) +
  guides(color="none")+#
  scale_color_manual(values=c("0"="black","1"="red")) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm")) +
  geom_cladelabel(node=235, label='HeV/NiV sister clade', offset=8.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=221, label='AndV clade', offset=8.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=10, label='AndV-2 clade', offset=8.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=225, label='AndV-3 clade', offset=8.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=339, label='AngV clade', offset=8.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=336, label='AngV-2 clade', offset=8.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=110, label='GhV-like clade', offset=8.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=131, label='AchPV2-like virus', offset=8.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=350, label='ThkPV1-like viruses', offset=8.5, fontsize=4, color="black", align=TRUE)+
  # 
  geom_cladelabel(node = 236, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 344, label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 370, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 378, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 383, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 392, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 387, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 404, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 395, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 362, label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 399, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 414, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 384, label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  xlim(c(0,12))

p1

p1.2<-ggtree::collapse(p1, 236)+geom_point2(aes(subset=(node==236)), size=4, shape=22, fill="firebrick1")
#p1.3<-ggtree::collapse(p1.2, 344)+geom_point2(aes(subset=(node==344)), size=4, shape=22, fill="firebrick4")
p1.4<-ggtree::collapse(p1.2, 370)+geom_point2(aes(subset=(node==370)), size=4, shape=22, fill="darkgoldenrod1")
p1.5<-ggtree::collapse(p1.4, 378)+geom_point2(aes(subset=(node==378)), size=4, shape=22, fill="darkorange")
p1.6<-ggtree::collapse(p1.5, 383)+geom_point2(aes(subset=(node==383)), size=4, shape=22, fill="darkorange3")
p1.7<-ggtree::collapse(p1.6, 392)+geom_point2(aes(subset=(node==392)), size=4, shape=22, fill="lightpink2")
p1.8<-ggtree::collapse(p1.7, 387)+geom_point2(aes(subset=(node==387)), size=4, shape=22, fill="deeppink1")
p1.9<-ggtree::collapse(p1.8, 404)+geom_point2(aes(subset=(node==404)), size=4, shape=22, fill="coral")
p1.10<-ggtree::collapse(p1.9, 395)+geom_point2(aes(subset=(node==395)), size=4, shape=22, fill="lightsalmon3")
p1.11<-ggtree::collapse(p1.10, 362)+geom_point2(aes(subset=(node==362)), size=4, shape=22, fill="indianred3")
p1.17<-ggtree::collapse(p1.11, 399)+geom_point2(aes(subset=(node==399)), size=4, shape=22, fill="coral3")
p1.18<-ggtree::collapse(p1.17, 414)+geom_point2(aes(subset=(node==414)), size=4, shape=22, fill="deeppink4")
p1.1<-ggtree::collapse(p1.18, 384)+geom_point2(aes(subset=(node==384)), size=4, shape=22, fill="lightpink3")
p1.1

p1.1<- p1.1 +  
  geom_hilight(node=235, fill="deeppink3", alpha=0.3, extendto=11.7)+ #HeV/NiV sister
  geom_hilight(node=221, fill="darkcyan", alpha=0.3, extendto=11.7)+ #AndV
  geom_hilight(node=10, fill="darkturquoise", alpha=0.3, extendto=11.7)+ #AndV2
  geom_hilight(node=225, fill="dodgerblue", alpha=0.3, extendto=11.7)+ #AndV3
  geom_hilight(node=339, fill="springgreen3", alpha=0.3, extendto=11.7)+ #AngV
  geom_hilight(node=336, fill="chartreuse", alpha=0.3, extendto=11.7)+ #AngV2
  geom_hilight(node=110, fill="orange", alpha=0.3, extendto=11.7)+ #GhVlike
  geom_hilight(node=131, fill="grey40", alpha=0.3, extendto=11.7)+ #AchPV2
  geom_hilight(node=350, fill="grey40", alpha=0.3, extendto=11.7) #ThkPV1
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
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1),
        plot.margin = unit(c(0,0,0,2), "cm"))
polytree











