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
title<-expression(paste("                    Full genome"))
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
  geom_cladelabel(node = 275, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4,fontface=3, color="black") +
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
        plot.margin = unit(c(0,0,0,0), "cm"))
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
title<-expression(paste("                    Polymerase"))
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +
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
        legend.key.size = unit(0.2, "cm"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  geom_cladelabel(node=235, label='HeV/NiV sister clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=221, label='AndV clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=10, label='AndV-2 clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=225, label='AndV-3 clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=339, label='AngV clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=336, label='AngV-2 clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=110, label='GhV-like clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=131, label='AchPV2-like virus', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=350, label='ThkPV1-like viruses', offset=4.8, fontsize=4, color="black", align=TRUE)+
  # 
  geom_cladelabel(node = 236, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4,fontface=3, color="black") +
  #geom_cladelabel(node = 344, label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 370, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 378, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 383, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 392, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 387, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 404, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 395, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 362, label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 399, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 414, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  geom_cladelabel(node = 384, label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
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
  geom_hilight(node=235, fill="deeppink3", alpha=0.3, extendto=8)+ #HeV/NiV sister
  geom_hilight(node=221, fill="darkcyan", alpha=0.3, extendto=8)+ #AndV
  geom_hilight(node=10, fill="darkturquoise", alpha=0.3, extendto=8)+ #AndV2
  geom_hilight(node=225, fill="dodgerblue", alpha=0.3, extendto=8)+ #AndV3
  geom_hilight(node=339, fill="springgreen3", alpha=0.3, extendto=8)+ #AngV
  geom_hilight(node=336, fill="chartreuse", alpha=0.3, extendto=8)+ #AngV2
  geom_hilight(node=110, fill="orange", alpha=0.3, extendto=8)+ #GhVlike
  geom_hilight(node=131, fill="grey40", alpha=0.3, extendto=8)+ #AchPV2
  geom_hilight(node=350, fill="grey40", alpha=0.3, extendto=8) #ThkPV1
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
        plot.margin = unit(c(0,0,0,0), "cm"))
polytree


##Now do the coverage plots
#set directory
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/genome_characterization/full_genome_coverage"
setwd(paste0(homewd))

#GhV-like
dat <- read.csv(file= "KEL86_KEL94_coverage_noUTR.csv", header = T, stringsAsFactors = F) #load data

title<-expression(paste(italic("GhV-like genome coverage")))
coverage_ghvlike<-  ggplot()+ 
  annotate("rect", xmin=50, xmax=4179, ymin=0, ymax=550, alpha=0.6,  fill="#D29292")+
  annotate("rect", xmin=4603, xmax=7636, ymin=0, ymax=550, alpha=0.6,  fill="#525832")+
  annotate("rect", xmin=8096, xmax=10940, ymin=0, ymax=550, alpha=0.6,  fill="#F98F5B")+
  annotate("rect", xmin=11664, xmax=16846, ymin=0, ymax=550, alpha=0.6,  fill="#BD9A5A")+
  annotate("rect", xmin=17373, xmax=23508, ymin=0, ymax=550, alpha=0.6,  fill="#F17983")+
  annotate("rect", xmin=24156, xmax=48178, ymin=0, ymax=550, alpha=0.6,  fill="#6B8C81")+
  geom_area(data=dat, aes(x=Position, y=rPM), fill="black") +
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  scale_x_continuous(breaks=c(0,5578,10396,17208,24521, 33205, 44590),
                     labels = c(0,3000,6000,9000,12000, 15000, 17000),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
        axis.title.x = element_blank(),
        plot.title = element_text(color = "darkorange1"),
        axis.line.x = element_blank())
coverage_ghvlike

#ROU148
dat <- read.csv(file= "ROU148_coverage.csv", header = T, stringsAsFactors = F) #load data

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
  ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
        plot.title = element_text(color = "darkcyan"),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())
coverage_rou148

#ROU173 and 125 together
dat <- read.csv(file= "ROU125_coverage.csv", header = T, stringsAsFactors = F) #load data
dat2 <- read.csv(file= "ROU173_coverage.csv", header = T, stringsAsFactors = F) #load data

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
  ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(color = "darkcyan"),
        axis.line.x = element_blank())
coverage_andv

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
genome1ghvlike<-ggplot(map_combo_pep1, aes(xmin = from, xmax = to, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(0, "mm"),
                  arrowhead_height = grid::unit(5, "mm"),
                  arrow_body_height = grid::unit(5, "mm")
  ) +
  scale_fill_manual(values=colz2)+
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
        axis.text.y = element_blank(),
        # axis.ticks.x=element_blank(),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        # axis.line.x = element_blank()
  )+
  xlab("Genome position") + ylab("")
genome2ghvlike

#put coverage and genome maps together
coverage_ghvlike<-coverage_ghvlike/genome1ghvlike/genome2ghvlike+plot_layout(nrow=3,  heights = c(1, 0.2, 0.2))
coverage_ghvlike

coverage_ghvlike<-as.ggplot(coverage_ghvlike)
coverage_ghvlike

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
genome1ROU148<-ggplot(map_combo_pep1, aes(xmin = from, xmax = to, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(0, "mm"),
                  arrowhead_height = grid::unit(5, "mm"),
                  arrow_body_height = grid::unit(5, "mm")
  ) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
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
        # axis.ticks.x=element_blank(),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        # axis.line.x = element_blank()
  )+
  xlab("Genome position") + ylab("")
genome2ROU148

#put coverage and genome maps together
genomeplot_ROU148<-coverage_rou148/genome1ROU148/genome2ROU148+plot_layout(nrow=3,  heights = c(1, 0.2, 0.2))
genomeplot_ROU148

genomeplot_ROU148<-as.ggplot(genomeplot_ROU148)
genomeplot_ROU148


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
genome1ROU173<-ggplot(map_combo_pep1, aes(xmin = from, xmax = to, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(0, "mm"),
                  arrowhead_height = grid::unit(5, "mm"),
                  arrow_body_height = grid::unit(5, "mm")
  ) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
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

genomeplot_ROU173<-as.ggplot(genomeplot_ROU173)
genomeplot_ROU173


#Figure with all the AndV figures
genomeplot_andv<-coverage_andv/genome1ROU173/genome2ROU173+plot_layout(nrow=3,  heights = c(1, 0.2, 0.2))
genomeplot_andv

genomeplot_andv<-as.ggplot(genomeplot_andv)


#AndV all together
andv_coverage<-plot_grid(genomeplot_ROU148, genomeplot_andv, labels=c("",""),
                         rel_widths=c(1,1), rel_heights = c(1,1),
                         ncol=1,align="hv", axis="l", label_size = 23)
andv_coverage

andv_coverage_long<-plot_grid(genomeplot_ROU148, genomeplot_andv, labels=c("",""),
                         rel_widths=c(1,1), rel_heights = c(1,1),
                         ncol=2,align="hv", axis="l", label_size = 23)
andv_coverage_long

#all coverage
all_coverage<-plot_grid(coverage_ghvlike, genomeplot_ROU148, genomeplot_andv, labels=c("","",""),
                         rel_widths=c(1,1,1), rel_heights = c(1,1,1),
                         ncol=3,align="hv", axis="l", label_size = 23)
all_coverage


#Now add the simplots in case they end up going with the figure
#set directory
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/genome_characterization/simplot"
setwd(paste0(homewd))

#aa plots
andv_map <- read.csv(file = "aa_align_AndV_noUTRtest.csv", header = T, stringsAsFactors = F)
head(andv_map)

#move to long
long.sim_aa <- melt(andv_map, id.vars = c("pointer"), measure.vars = c("AndV_ROU148","AndV_ROU173","AndV_ROU125"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "AndV_ROU148"] <- "PX127343 AndV"
long.sim_aa$accession[long.sim_aa$accession == "AndV_ROU173"] <- "PX127344 AndV"
long.sim_aa$accession[long.sim_aa$accession == "AndV_ROU125"] <- "PX127342 AndV"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c( 
  "PX127342 AndV",
  "PX127343 AndV", 
  "PX127344 AndV"))
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

colzpalette<-c("grey","grey30","black")

title<-expression(paste(italic("Reference: YuV2")))
## amino acid
andv_map_aa <- ggplot(long.sim_aa) + 
  annotate("rect", xmin=57, xmax=592, ymin=0, ymax=1, alpha=0.6,  fill="#D29292")+
  annotate("rect", xmin=780, xmax=1527, ymin=0, ymax=1, alpha=0.6,  fill="#525832")+
  annotate("rect", xmin=1702, xmax=2035, ymin=0, ymax=1, alpha=0.6,  fill="#F98F5B")+
  annotate("rect", xmin=2143, xmax=2742, ymin=0, ymax=1, alpha=0.6,  fill="#BD9A5A")+
  annotate("rect", xmin=2934, xmax=3847, ymin=0, ymax=1, alpha=0.6,  fill="#F17983")+
  annotate("rect", xmin=4041, xmax=6340, ymin=0, ymax=1, alpha=0.6,  fill="#6B8C81")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Amino acid identity")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.direction = "horizontal",
        legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 10),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position=c(0.8,0.1), 
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  ggtitle(title)+
  scale_color_manual(values=colzpalette) + 
  coord_cartesian(ylim=c(0,1.0))+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055, 10000/3.055, 12000/3.055, 14000/3.055, 16000/3.055, 18000/3.055), 
                     labels = c(0,2000, 4000,6000,8000, 10000, 12000, 14000, 16000, 18000),
                     expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0)); andv_map_aa

#ghvlike
ghv_map <- read.csv(file = "aa_align_GhVlike_noUTRtest.csv", header = T, stringsAsFactors = F)
head(ghv_map)

#move to long
long.sim_aa <- melt(ghv_map, id.vars = c("pointer"), measure.vars = c("GhVlike"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "GhVlike"] <- "GhVlike"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c( 
  "GhVlike"))
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

colzpalette<-c("black","grey30","grey")

title<-expression(paste(italic("Reference: GhV")))
## amino acid
ghv_map_aa <- ggplot(long.sim_aa) + 
  annotate("rect", xmin=1, xmax=554, ymin=0, ymax=1, alpha=0.6,  fill="#D29292")+
  annotate("rect", xmin=730, xmax=1628, ymin=0, ymax=1, alpha=0.6,  fill="#525832")+
  annotate("rect", xmin=1722, xmax=2197, ymin=0, ymax=1, alpha=0.6,  fill="#F98F5B")+
  annotate("rect", xmin=2448, xmax=3001, ymin=0, ymax=1, alpha=0.6,  fill="#BD9A5A")+
  annotate("rect", xmin=3117, xmax=3889, ymin=0, ymax=1, alpha=0.6,  fill="#F17983")+
  annotate("rect", xmin=4061, xmax=6369, ymin=0, ymax=1, alpha=0.6,  fill="#6B8C81")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  #annotate("rect", xmin=3519, xmax=4168, ymin=0, ymax=1, alpha=1,  fill="grey30")+
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Amino acid identity")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.direction = "horizontal",
        legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 10),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position=c(0.8,0.05), 
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  ggtitle(title)+
  scale_color_manual(values=colzpalette) + 
  coord_cartesian(ylim=c(0,1.0))+
  scale_x_continuous(
    breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055, 10000/3.055, 12000/3.055, 14000/3.055, 16000/3.055, 18000/3.055), 
    labels = c(0,2000, 4000,6000,8000, 10000, 12000, 14000, 16000, 18000),
    expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0)); ghv_map_aa



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
map_combo<-subset(map,type=="KEL86_KEL94" & alignment=="nt")
map_combo_pep<-subset(map_pep,type=="KEL86_KEL94"& alignment=="nt")
map_combo_pep1<-subset(map_combo_pep,string=="first"& alignment=="nt")
map_combo_feat<-subset(map_feat,type=="KEL86_KEL94"& alignment=="nt")
map_combo_pep2<-subset(map_combo_pep,string=="second"& alignment=="nt")

#map the genomes
genome1ghvlike<-ggplot(map_combo_pep1, aes(xmin = from, xmax = to, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(0, "mm"),
                  arrowhead_height = grid::unit(5, "mm"),
                  arrow_body_height = grid::unit(5, "mm")
  ) +
  scale_fill_manual(values=colz2)+
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
  scale_x_continuous(
    #breaks=c(0,5578,10396,17208,24521, 33205, 44590),
    #labels = c(0,3000,6000,9000,12000, 15000, 17000),
    expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        # axis.ticks.x=element_blank(),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        # axis.line.x = element_blank()
  )+
  xlab("Genome position") + ylab("")
genome2ghvlike

#put coverage and genome maps together
genomeplot_ghvlike<-ghv_map_aa/genome1ghvlike/genome2ghvlike+plot_layout(nrow=3,  heights = c(1, 0.2,0.2))
genomeplot_ghvlike

genomeplot_ghvlike<-as.ggplot(genomeplot_ghvlike)
genomeplot_ghvlike


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
map_combo<-subset(map,type=="ROU173" & alignment=="nt")
map_combo_pep<-subset(map_pep,type=="ROU173"& alignment=="nt")
map_combo_pep1<-subset(map_combo_pep,string=="first"& alignment=="nt")
map_combo_feat<-subset(map_feat,type=="ROU173"& alignment=="nt")
map_combo_pep2<-subset(map_combo_pep,string=="second"& alignment=="nt")

#map the genomes
title<-expression(paste(italic("PX127344")))
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
map_combo<-subset(map,type=="AndV"& alignment=="nt")
map_combo_pep<-subset(map_pep,type=="AndV"& alignment=="nt")
map_combo_pep1<-subset(map_combo_pep,string=="first"& alignment=="nt")
map_combo_feat<-subset(map_feat,type=="AndV"& alignment=="nt")
map_combo_pep2<-subset(map_combo_pep,string=="second"& alignment=="nt")

#map the genomes
title<-expression(paste(italic("PX127342 & PX127343")))
genome1andv<-ggplot(map_combo_pep1, aes(xmin = from, xmax = to, y = molecule, fill=gene)) +
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

#Figure with all the AndV figures
genomeplot_andv<-andv_map_aa/genome1andv/genome2andv/genome1ROU173/genome2ROU173+plot_layout(nrow=5,  heights = c(1, 0.2,0.2,0.2,0.2))
genomeplot_andv

genomeplot_andv<-as.ggplot(genomeplot_andv)

genomeplot_andv_reduced<-andv_map_aa/genome1ROU173/genome2ROU173+plot_layout(nrow=3,  heights = c(1, 0.2,0.2))
genomeplot_andv_reduced

#both genome plots together
simplots<-plot_grid(genomeplot_ghvlike, genomeplot_andv_reduced,
                 rel_widths=c(1,1), rel_heights = c(1,1),
                 ncol=2,align="hv", axis="r", label_size = 23); simplots
simplots<-as.ggplot(simplots)

#all together in figure 1

part1<-plot_grid(all_coverage, simplots,
                labels=c("B","C"),
                rel_widths=c(1,1), rel_heights = c(1,1),
                ncol=1,align="hv", axis="l", label_size = 23)+
  theme(plot.margin = unit(c(0,0,0,0), "cm")); part1
part1<-as.ggplot(part1)


fig1<-plot_grid(polytree, part1, 
                labels=c("A",""),
                rel_widths=c(1,1.5), rel_heights = c(1,1),
                ncol=2,nrow=1,align="hv", axis="l", label_size = 23)+
  theme(plot.background = element_rect(fill = "white", colour = NA)); fig1

fig1.2<-plot_grid(polytree, part1, 
                labels=c("A",""),
                rel_widths=c(1,1), rel_heights = c(1,1.2),
                ncol=1,nrow=2,align="hv", axis="r", label_size = 23)+
  theme(plot.background = element_rect(fill = "white", colour = NA)); fig1.2


#save
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/final_figures_tables/gwen_tmp/"
ggsave(file = paste0(homewd, "Fig1_phylogeny_coverage_simplot_v2.png"),
       plot= fig1.2,
       units="mm",
       width=210,
       height=230,
       scale=2,
       dpi=500)



