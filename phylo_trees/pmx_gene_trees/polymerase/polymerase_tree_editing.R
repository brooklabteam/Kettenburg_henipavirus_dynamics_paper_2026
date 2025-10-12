#plotting of gene trees for project

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
homewd= "/Users/gwenddolenkettenburg/Desktop/brook_lab/henipavirus_shedding_cross_correlation_project/phylo_trees/"

#polymerase files first
setwd(paste0(homewd,"/pmx_gene_trees/polymerase"))

#load the trees and root them
tree1 <-  read.tree("polymerase_full_align.fasta.treefile")
tree2 <-  read.tree("polymerase_left_align.fasta.treefile") 
tree4 <-  read.tree("polymerase_mid_align.fasta.treefile") 
tree5 <-  read.tree("polymerase_right_align.fasta.treefile") 
rooted.tree1 <- root(tree1, which(tree1$tip.label == "YP_009094051.1"))
rooted.tree2 <- root(tree2, which(tree2$tip.label == "YP_009094051.1"))
rooted.tree4 <- root(tree4, which(tree4$tip.label == "YP_009094051.1"))
rooted.tree5 <- root(tree5, which(tree5$tip.label == "YP_009094051.1"))

#Remove root from displaying, still calculates changes correctly without it
rooted.tree1<-drop.tip(rooted.tree1, "YP_009094051.1")
rooted.tree2<-drop.tip(rooted.tree2, "YP_009094051.1")
rooted.tree2<-drop.tip(rooted.tree2, "MIZ193_NODE_81")
rooted.tree2<-drop.tip(rooted.tree2, "KEL277_NODE_6397")
rooted.tree4<-drop.tip(rooted.tree4, "YP_009094051.1")
rooted.tree4<-drop.tip(rooted.tree4, "MIZ193_NODE_166")
rooted.tree5<-drop.tip(rooted.tree5, "YP_009094051.1")
rooted.tree5<-drop.tip(rooted.tree5, "KEL277_NODE_812")

#load tree data prepared from elsewhere
dat1 <- read.csv(("L_metadata_full.csv"), header = T, stringsAsFactors = F)
dat2 <- read.csv(("L_metadata_left.csv"), header = T, stringsAsFactors = F)
dat4 <- read.csv(("L_metadata_mid.csv"), header = T, stringsAsFactors = F)
dat5 <- read.csv(("L_metadata_right.csv"), header = T, stringsAsFactors = F)

#check subgroup names
unique(dat1$Genus)
dat1$novel <- as.factor(dat1$novel)
dat2$novel <- as.factor(dat2$novel)
dat4$novel <- as.factor(dat4$novel)
dat5$novel <- as.factor(dat5$novel)

colz = c("Orthoavulavirus"="darkorange","Paraavulavirus"="darkorange3","Metaavulavirus"="darkgoldenrod1","Jeilongvirus"="coral",
         "Morbillivirus"="coral3","Respirovirus"="deeppink1", "Parahenipavirus"="deeppink4",
         "Henipavirus"="firebrick1","Cynoglossusvirus"="grey50","Pararubulavirus"="firebrick4","Salemvirus"="indianred1",
         "Ferlavirus"="indianred3","Orthorubulavirus"="indianred4",
         "Hoplichthysvirus"="lightsalmon1","Narmovirus"="lightsalmon3","Aquaparamyxovirus"="lightpink2","Scoliodonvirus"="lightpink3",
         "Synodonvirus"="pink1","Unclassified"="black",
         "Sunshinevirus"="chartreuse1")


#take a glance
p <- ggtree(rooted.tree2) %<+% dat2 + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=2) + geom_nodelab(size=2) +  #scale_color_manual(values=colz) + 
  theme(legend.position = "none", legend.title = element_blank())
p #looks great

#now get new tip labels
dat1$old_tip_label <- dat1$tip_label
dat1$new_label <- NA
dat2$old_tip_label <- dat2$tip_label
dat2$new_label <- NA
dat4$old_tip_label <- dat4$tip_label
dat4$new_label <- NA
dat5$old_tip_label <- dat5$tip_label
dat5$new_label <- NA

#now, select the name based on what components are present for each sample
#all components with values:
dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                             dat1$Species[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                             dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                             dat1$source[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                             dat1$Country[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                             dat1$Collection_Date[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                            dat1$Species[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                            #dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                            dat1$source[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                            dat1$Country[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                            dat1$Collection_Date[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
dat2$new_label[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)] <- paste(dat2$Accession[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                             dat2$Species[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                             dat2$Isolate[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|", 
                                                                                                             dat2$source[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                             dat2$Country[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                             dat2$Collection_Date[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat2$new_label[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)] <- paste(dat2$Accession[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|", 
                                                                                                            dat2$Species[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|", 
                                                                                                            #dat2$Isolate[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|", 
                                                                                                            dat2$source[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                            dat2$Country[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                            dat2$Collection_Date[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)])
dat4$new_label[!is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)] <- paste(dat4$Accession[!is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)], "|",
                                                                                                             dat4$Species[!is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)], "|",
                                                                                                             dat4$Isolate[!is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)], "|", 
                                                                                                             dat4$source[!is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)], "|",
                                                                                                             dat4$Country[!is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)], "|",
                                                                                                             dat4$Collection_Date[!is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat4$new_label[is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)] <- paste(dat4$Accession[is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)], "|", 
                                                                                                            dat4$Species[is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)], "|", 
                                                                                                            #dat4$Isolate[!is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)], "|", 
                                                                                                            dat4$source[is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)], "|",
                                                                                                            dat4$Country[is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)], "|",
                                                                                                            dat4$Collection_Date[is.na(dat4$Isolate) & !is.na(dat4$Accession) &!is.na(dat4$source) &!is.na(dat4$Country)])
dat5$new_label[!is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)] <- paste(dat5$Accession[!is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)], "|",
                                                                                                             dat5$Species[!is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)], "|",
                                                                                                             dat5$Isolate[!is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)], "|", 
                                                                                                             dat5$source[!is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)], "|",
                                                                                                             dat5$Country[!is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)], "|",
                                                                                                             dat5$Collection_Date[!is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat5$new_label[is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)] <- paste(dat5$Accession[is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)], "|", 
                                                                                                            dat5$Species[is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)], "|", 
                                                                                                            #dat5$Isolate[!is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)], "|", 
                                                                                                            dat5$source[is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)], "|",
                                                                                                            dat5$Country[is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)], "|",
                                                                                                            dat5$Collection_Date[is.na(dat5$Isolate) & !is.na(dat5$Accession) &!is.na(dat5$source) &!is.na(dat5$Country)])

#make sure to sort in order
colnames(dat1)[colnames(dat1)=="Accession"]="old_tip_label"
tree.dat1 <- data.frame(old_tip_label=rooted.tree1$tip.label, num =1:length(rooted.tree1$tip.label))
head(tree.dat1)
head(dat1)
tree.dat1 <- merge(tree.dat1, dat1, by = "old_tip_label", all.x = T, sort = F)
names(tree.dat1)
tree.dat1$tip_label <- tree.dat1$new_label
tree.dat1 <- dplyr::select(tree.dat1, tip_label, Isolate, Host, source, Country, Collection_Date, Genus, Species, novel, old_tip_label)
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
tree.dat2 <- dplyr::select(tree.dat2, tip_label, Isolate, Host, source, Country, Collection_Date, Genus, Species, novel, old_tip_label)
rooted.tree2$tip.label <- tree.dat2$tip_label
head(tree.dat2)
cbind(tree.dat2$old_tip_label, rooted.tree2$tip.label)
tree.dat2$tip_label#all look good

colnames(dat4)[colnames(dat4)=="Accession"]="old_tip_label"
tree.dat4 <- data.frame(old_tip_label=rooted.tree4$tip.label, num =1:length(rooted.tree4$tip.label))
head(tree.dat4)
head(dat4)
tree.dat4 <- merge(tree.dat4, dat4, by = "old_tip_label", all.x = T, sort = F)
names(tree.dat4)
tree.dat4$tip_label <- tree.dat4$new_label
tree.dat4 <- dplyr::select(tree.dat4, tip_label, Isolate, Host, source, Country, Collection_Date, Genus, Species, novel, old_tip_label)
rooted.tree4$tip.label <- tree.dat4$tip_label
head(tree.dat4)
cbind(tree.dat4$old_tip_label, rooted.tree4$tip.label)
tree.dat4$tip_label#all look good

colnames(dat5)[colnames(dat5)=="Accession"]="old_tip_label"
tree.dat5 <- data.frame(old_tip_label=rooted.tree5$tip.label, num =1:length(rooted.tree5$tip.label))
head(tree.dat5)
head(dat5)
tree.dat5 <- merge(tree.dat5, dat5, by = "old_tip_label", all.x = T, sort = F)
names(tree.dat5)
tree.dat5$tip_label <- tree.dat5$new_label
tree.dat5 <- dplyr::select(tree.dat5, tip_label, Isolate, Host, source, Country, Collection_Date, Genus, Species, novel, old_tip_label)
rooted.tree5$tip.label <- tree.dat5$tip_label
head(tree.dat5)
cbind(tree.dat5$old_tip_label, rooted.tree5$tip.label)
tree.dat5$tip_label#all look good


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

tree.dat4$Host[tree.dat4$Host==0] <- "Non-bat host"
tree.dat4$Host[tree.dat4$Host==1] <- "Bat host"
tree.dat4$Host <- as.factor(tree.dat4$Host)
shapez = c("Bat host" =  17, "Non-bat host" = 19)
colz2 = c('1' =  "yellow", '0' = "white")

tree.dat5$Host[tree.dat5$Host==0] <- "Non-bat host"
tree.dat5$Host[tree.dat5$Host==1] <- "Bat host"
tree.dat5$Host <- as.factor(tree.dat5$Host)
shapez = c("Bat host" =  17, "Non-bat host" = 19)
colz2 = c('1' =  "yellow", '0' = "white")

##Get the clade numbers so we can collapse unnecessary clades
ggtree(rooted.tree1) + geom_text(aes(label=node), hjust=-.3)
ggtree(rooted.tree2) + geom_text(aes(label=node), hjust=-.3)
ggtree(rooted.tree4) + geom_text(aes(label=node), hjust=-.3)
ggtree(rooted.tree5) + geom_text(aes(label=node), hjust=-.3)


##uncollapsed tree

#Full alignment
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Genus="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=2, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm")) +
  geom_cladelabel(node = 220, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 3, color="black") +
  xlim(c(0,14))

p1

p1.1<-collapse(p1, 220)+geom_point2(aes(subset=(node==220)), size=3, shape=22, fill="firebrick1")
p1.1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree1)
p0.dat1 <- p0$data
p0.dat1$Bootstrap <- NA
Bootstrap<-p0.dat1$Bootstrap[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)] <- as.numeric(p0.dat1$label[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)])#fill with label

#add bootstrap values to original plot
tree1 <- p1.1  %<+% p0.dat1 + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  ggtitle("polymerase tree full")+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1))
tree1


#Left alignment
p2 <- ggtree(rooted.tree2) %<+% tree.dat2 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Genus="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=2, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm")) +
  geom_cladelabel(node = 229, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 3, color="black") +
  xlim(c(0,14))

p2

p2.1<-collapse(p2, 229)+geom_point2(aes(subset=(node==229)), size=3, shape=22, fill="firebrick1")
p2.1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree2)
p0.dat1 <- p0$data
p0.dat1$Bootstrap <- NA
Bootstrap<-p0.dat1$Bootstrap[(length(tree.dat2$tip_label)+1):length(p0.dat1$label)] <- as.numeric(p0.dat1$label[(length(tree.dat2$tip_label)+1):length(p0.dat1$label)])#fill with label

#add bootstrap values to original plot
tree2 <- p2.1  %<+% p0.dat1 + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  ggtitle("polymerase tree left")+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1))
tree2



#Mid alignment
p4 <- ggtree(rooted.tree4) %<+% tree.dat4 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Genus="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=2, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm")) +
  geom_cladelabel(node = 249, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 3, color="black") +
  xlim(c(0,14))

p4

p4.1<-collapse(p4, 249)+geom_point2(aes(subset=(node==249)), size=3, shape=22, fill="firebrick1")
p4.1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree4)
p0.dat1 <- p0$data
p0.dat1$Bootstrap <- NA
Bootstrap<-p0.dat1$Bootstrap[(length(tree.dat4$tip_label)+1):length(p0.dat1$label)] <- as.numeric(p0.dat1$label[(length(tree.dat4$tip_label)+1):length(p0.dat1$label)])#fill with label

#add bootstrap values to original plot
tree4 <- p4.1  %<+% p0.dat1 + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  ggtitle("polymerase tree mid")+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1))
tree4

#Right alignment
p5 <- ggtree(rooted.tree5) %<+% tree.dat5 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Genus="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=2, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm")) +
  geom_cladelabel(node = 218, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 3, color="black") +
  xlim(c(0,14))

p5

p5.1<-collapse(p5, 218)+geom_point2(aes(subset=(node==218)), size=3, shape=22, fill="firebrick1")
p5.1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree5)
p0.dat1 <- p0$data
p0.dat1$Bootstrap <- NA
Bootstrap<-p0.dat1$Bootstrap[(length(tree.dat5$tip_label)+1):length(p0.dat1$label)] <- as.numeric(p0.dat1$label[(length(tree.dat5$tip_label)+1):length(p0.dat1$label)])#fill with label

#add bootstrap values to original plot
tree5 <- p5.1  %<+% p0.dat1 + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  ggtitle("polymerase tree right")+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1))
tree5


###Put together final figure
Phylo<-plot_grid(tree1, tree2, tree4, tree5, labels=c("","", "","",""),
                 rel_widths=c(1,1,1,1), rel_heights = c(1,1,1,1),
                 ncol=2,align="hv", axis="l", label_size = 23)
Phylo

#save trees
ggsave(file = paste0(homewd, "/figures/polymerase_phylogeny.pdf"),
       plot= Phylo,
       units="mm",  
       width=220, 
       height=250, 
       scale=2, 
       dpi=500)

