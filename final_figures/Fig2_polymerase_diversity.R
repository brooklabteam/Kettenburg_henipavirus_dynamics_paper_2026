#Script to make figure 2

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

##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_circulating_clades/phylo_trees/"
#refseq tree first
setwd(paste0(homewd,"/pmx_gene_trees/polymerase"))
tree1 <-  read.tree("polymerase_mid_align.fasta.treefile")
#drop root for aesthetics
rooted.tree1 <- root(tree1, which(tree1$tip.label == "YP_009094051.1"))
rooted.tree1<-drop.tip(rooted.tree1, "YP_009094051.1")
#load tree data prepared from elsewhere
dat1 <- read.csv(("L_metadata_mid.csv"), header = T, stringsAsFactors = F)

#henipa tree next
setwd(paste0(homewd,"/henipa_focus_L_tree"))
tree2 <-  read.tree("ref_right.fasta.treefile") 
rooted.tree2 <- root(tree2, which(tree2$tip.label == "YP_009094051.1"))
#drop root and other reference sequences for aesthetics
rooted.tree2<-drop.tip(rooted.tree2, "YP_009094051.1")
rooted.tree2<-drop.tip(rooted.tree2, "YP_010797981.1")
rooted.tree2<-drop.tip(rooted.tree2, "YP_010790495.1")
rooted.tree2<-drop.tip(rooted.tree2, "YP_010799241.1")
rooted.tree2<-drop.tip(rooted.tree2, "YP_010798268.1")
rooted.tree2<-drop.tip(rooted.tree2, "NP_056924.1")
rooted.tree2<-drop.tip(rooted.tree2, "NP_054697.1")
rooted.tree2<-drop.tip(rooted.tree2, "YP_138518.1")
rooted.tree2<-drop.tip(rooted.tree2, "YP_010790514.1")
rooted.tree2<-drop.tip(rooted.tree2, "YP_009094339.1")
rooted.tree2<-drop.tip(rooted.tree2, "YP_010790479.1")
rooted.tree2<-drop.tip(rooted.tree2, "YP_010790486.1")
rooted.tree2<-drop.tip(rooted.tree2, "YP_010774627.1")
rooted.tree2<-drop.tip(rooted.tree2, "YP_010790504.1")
#load tree data prepared from elsewhere
dat2 <- read.csv(("ref_right_metadata.csv"), header = T, stringsAsFactors = F)

#check subgroup names
unique(dat1$Genus)
dat1$novel <- as.factor(dat1$novel)
dat2$novel <- as.factor(dat2$novel)

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
dat2$old_tip_label <- dat2$tip_label
dat2$new_label <- NA

#now, select the name based on what components are present for each sample
#all components with values:
dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                                  dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                                  dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                                  dat1$source[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                                  dat1$Country[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                                  dat1$Collection_Date[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                                 dat1$Organism_Name[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                                 #dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                                 dat1$source[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                                 dat1$Country[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                                 dat1$Collection_Date[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
dat2$new_label[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)] <- paste(dat2$Accession[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                  dat2$Organism_Name[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                  dat2$Isolate[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|", 
                                                                                                                  dat2$source[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                  dat2$Country[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                  dat2$Collection_Date[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat2$new_label[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)] <- paste(dat2$Accession[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|", 
                                                                                                                 dat2$Organism_Name[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|", 
                                                                                                                 #dat2$Isolate[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|", 
                                                                                                                 dat2$source[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                 dat2$Country[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                 dat2$Collection_Date[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)])
#make sure to sort in order
colnames(dat1)[colnames(dat1)=="Accession"]="old_tip_label"
tree.dat1 <- data.frame(old_tip_label=rooted.tree1$tip.label, num =1:length(rooted.tree1$tip.label))
head(tree.dat1)
head(dat1)
tree.dat1 <- merge(tree.dat1, dat1, by = "old_tip_label", all.x = T, sort = F)
names(tree.dat1)
tree.dat1$tip_label <- tree.dat1$new_label
tree.dat1 <- dplyr::select(tree.dat1, tip_label, Isolate, Host, source, Country, Collection_Date, Genus, Organism_Name, novel, old_tip_label)
rooted.tree1$tip.label <- tree.dat1$tip_label
head(tree.dat1)
cbind(tree.dat1$old_tip_label, rooted.tree1$tip.label)
tree.dat1$tip_label#all look good

colnames(dat2)[colnames(dat2)=="Accession"]="old_tip_label"
tree.dat2 <- data.frame(old_tip_label=rooted.tree2$tip.label, num =1:length(rooted.tree2$tip.label))
head(tree.dat2)
head(dat2)
tree.dat2 <- merge(tree.dat2, dat2, by = "old_tip_label", all.x = T, sort = F)
names(tree.dat2)
tree.dat2$tip_label <- tree.dat2$new_label
tree.dat2 <- dplyr::select(tree.dat2, tip_label, Isolate, Host, source, Country, Collection_Date, Genus, Organism_Name, novel, old_tip_label)
rooted.tree2$tip.label <- tree.dat2$tip_label
head(tree.dat2)
cbind(tree.dat2$old_tip_label, rooted.tree2$tip.label)
tree.dat2$tip_label#all look good

#shaped and highlighting
tree.dat1$Host[tree.dat1$Host==0] <- "Non-bat host"
tree.dat1$Host[tree.dat1$Host==1] <- "Bat host"
tree.dat1$Host <- as.factor(tree.dat1$Host)
shapez = c("Bat host" =  17, "Non-bat host" = 19)
colz2 = c('1' =  "yellow", '0' = "white")

tree.dat2$Host[tree.dat2$Host==0] <- "Non-bat host"
tree.dat2$Host[tree.dat2$Host==1] <- "Bat host"
tree.dat2$Host <- as.factor(tree.dat2$Host)
shapez = c("Bat host" =  17, "Non-bat host" = 19)
colz2 = c('1' =  "yellow", '0' = "white")

##Get the clade numbers so we can collapse unnecessary clades
ggtree(rooted.tree1) + geom_text(aes(label=node), hjust=-.3)
ggtree(rooted.tree2) + geom_text(aes(label=node), hjust=-.3)

#Refseq alignment
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +
  #scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_color() +
  geom_tiplab(aes(color = novel, show.legend=F), geom = "text", Genus="Helvetica", size=4, nudge_x=0.05) +
  guides(color="none")+#
  scale_color_manual(values=c("0"="black","1"="red")) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm")) +
  geom_cladelabel(node=390, label='HeV/NiV sister clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=158, label='GhV-like clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=386, label='AndV clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=152, label='AndV-2 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=382, label='AndV-3 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=282, label='AngV clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=279, label='AngV-2 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=279, label='AngV-2 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=70, label='AchPV2-like virus', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=300, label='ThkPV1-like viruses', offset=8.55, fontsize=4, color="black", align=TRUE)+
  
  geom_cladelabel(node = 271, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 344, label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 325, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 320, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 333, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 287, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 288, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 347, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 344, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 312, label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 339, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 357, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 334, label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  xlim(c(0,14))

p1

p1.2<-collapse(p1, 271)+geom_point2(aes(subset=(node==271)), size=4, shape=22, fill="firebrick1")
#p1.3<-collapse(p1.2, 344)+geom_point2(aes(subset=(node==344)), size=4, shape=22, fill="firebrick4")
p1.4<-collapse(p1.2, 325)+geom_point2(aes(subset=(node==325)), size=4, shape=22, fill="darkgoldenrod1")
p1.5<-collapse(p1.4, 320)+geom_point2(aes(subset=(node==320)), size=4, shape=22, fill="darkorange")
p1.6<-collapse(p1.5, 333)+geom_point2(aes(subset=(node==333)), size=4, shape=22, fill="darkorange3")
p1.7<-collapse(p1.6, 287)+geom_point2(aes(subset=(node==287)), size=4, shape=22, fill="lightpink2")
p1.8<-collapse(p1.7, 288)+geom_point2(aes(subset=(node==288)), size=4, shape=22, fill="deeppink1")
p1.9<-collapse(p1.8, 347)+geom_point2(aes(subset=(node==347)), size=4, shape=22, fill="coral")
p1.10<-collapse(p1.9, 344)+geom_point2(aes(subset=(node==344)), size=4, shape=22, fill="lightsalmon3")
p1.11<-collapse(p1.10, 312)+geom_point2(aes(subset=(node==312)), size=4, shape=22, fill="indianred3")
p1.17<-collapse(p1.11, 339)+geom_point2(aes(subset=(node==339)), size=4, shape=22, fill="coral3")
p1.18<-collapse(p1.17, 357)+geom_point2(aes(subset=(node==357)), size=4, shape=22, fill="deeppink4")
p1.1<-collapse(p1.18, 334)+geom_point2(aes(subset=(node==334)), size=4, shape=22, fill="lightpink3")
p1.1

p1.1<- p1.1 +  
  geom_hilight(node=390, fill="deeppink3", alpha=0.3, extendto=12.8)+ #HeV/NiV sister
  geom_hilight(node=158, fill="darkorange1", alpha=0.3, extendto=12.8)+ #GhVlike
  geom_hilight(node=386, fill="darkcyan", alpha=0.3, extendto=12.8)+ #AndV
  geom_hilight(node=152, fill="darkturquoise", alpha=0.3, extendto=12.8)+ #AndV2
  geom_hilight(node=382, fill="dodgerblue", alpha=0.3, extendto=12.8)+ #AndV3
  geom_hilight(node=282, fill="springgreen3", alpha=0.3, extendto=12.8)+ #AngV
  geom_hilight(node=279, fill="chartreuse", alpha=0.3, extendto=12.8)+ #AngV2
  geom_hilight(node=70, fill="grey40", alpha=0.3, extendto=12.8)+ #AngV
  geom_hilight(node=300, fill="grey40", alpha=0.3, extendto=12.8) #AngV2
p1.1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree1)
p0.dat1 <- p0$data
p0.dat1$Bootstrap <- NA
Bootstrap<-p0.dat1$Bootstrap[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)] <- as.numeric(p0.dat1$label[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)])#fill with label

#add bootstrap values to original plot
tree1 <- p1.1  %<+% p0.dat1 + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0, size=2)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1))
tree1

#henipavirus focused tree
p2 <- ggtree(rooted.tree2) %<+% tree.dat2 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_color() +
  geom_tiplab(aes(color = novel, show.legend=F), geom = "text", Genus="Helvetica", size=4, nudge_x=0.05) +
  guides(color="none")+#
  scale_color_manual(values=c("0"="black","1"="red")) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm")) +
  geom_cladelabel(node=467, label='HeV/NiV sister clade', offset=9.15, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=145, label='GhV-like clade', offset=9.15, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=512, label='AndV clade', offset=9.15, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=169, label='AndV-2 clade', offset=9.15, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=502, label='AndV-3.1 clade', offset=9.15, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=495, label='AndV-3.2 clade', offset=9.15, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=354, label='AngV clade', offset=9.15, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=363, label='AngV-2.1 clade', offset=9.15, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=370, label='AngV-2.2 clade', offset=9.15, fontsize=4, color="black", align=TRUE)+
  
  geom_cladelabel(node = 381, label = "Hendra/Nipah henipavirus clade (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 561, label = "Hendra/Nipa sister clade | Eidolon helvum | Zambia (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 476, label = "GhV-like clade | Eidolon helvum | Zambia (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 514, label = "Unclassified henipavirus clade | Eidolon helvum | Ghana and DRC (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 521, label = "Parahenipavirus clade (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 636, label = "Rubulavirus clade (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 531, label = "Unclassified bat paramyxovirus clade (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 344, label = "Unclassified bat paramyxovirus clade (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 372, label = "Unclassified bat paramyxovirus clade (collapsed)",offset=0.1, fontsize = 4, color="black")+
  xlim(c(0,14))
p2

p2.2<-collapse(p2, 381)+geom_point2(aes(subset=(node==381)), size=4, shape=22, fill="firebrick1")
#p2.3<-collapse(p2.2, 561)+geom_point2(aes(subset=(node==561)), size=4, shape=22, fill="firebrick1")
p2.4<-collapse(p2.2, 476)+geom_point2(aes(subset=(node==476)), size=4, shape=22, fill="firebrick1")
p2.5<-collapse(p2.4, 514)+geom_point2(aes(subset=(node==514)), size=4, shape=22, fill="firebrick1")
p2.6<-collapse(p2.5, 521)+geom_point2(aes(subset=(node==521)), size=4, shape=22, fill="deeppink4")
p2.7<-collapse(p2.6, 636)+geom_point2(aes(subset=(node==636)), size=4, shape=22, fill="firebrick4")
p2.8<-collapse(p2.7, 531)+geom_point2(aes(subset=(node==531)), size=4, shape=22, fill="black")
p2.9<-collapse(p2.8, 344)+geom_point2(aes(subset=(node==344)), size=4, shape=22, fill="black")
p2.1<-collapse(p2.9, 372)+geom_point2(aes(subset=(node==372)), size=4, shape=22, fill="black")
p2.1

p2.1<- p2.1 +  
  geom_hilight(node=467, fill="deeppink3", alpha=0.3, extendto=12.8)+ #HeV/NiV sister
  geom_hilight(node=145, fill="darkorange1", alpha=0.3, extendto=12.8)+ #GhVlike
  geom_hilight(node=512, fill="darkcyan", alpha=0.3, extendto=12.8)+ #AndV
  geom_hilight(node=169, fill="darkturquoise", alpha=0.3, extendto=12.8)+ #AndV2
  geom_hilight(node=502, fill="dodgerblue", alpha=0.3, extendto=12.8)+ #AndV3.1
  geom_hilight(node=495, fill="dodgerblue", alpha=0.3, extendto=12.8)+ #AndV3.2
  geom_hilight(node=354, fill="springgreen3", alpha=0.3, extendto=12.8)+ #AngV
  geom_hilight(node=363, fill="chartreuse", alpha=0.3, extendto=12.8) + #AngV2.1
  geom_hilight(node=370, fill="chartreuse", alpha=0.3, extendto=12.8) #AngV2.2

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
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1))
tree2

###Put together final figure
Phylo<-plot_grid(tree1, tree2, labels=c("A","B"),
                 rel_widths=c(1,1.22), rel_heights = c(1,1),
                 ncol=2,align="hv", axis="l", label_size = 23)
Phylo

#save trees
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_circulating_clades/final_figures/"
ggsave(file = paste0(homewd, "Fig2_polymerase_diversity.pdf"),
       plot= Phylo,
       units="mm",  
       width=380, 
       height=120, 
       scale=2, 
       dpi=500)














