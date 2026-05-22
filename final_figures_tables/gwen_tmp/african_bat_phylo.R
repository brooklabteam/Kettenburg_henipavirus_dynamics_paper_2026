#Script to make african bat phylogeny for supplemental

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

#henipa tree
setwd(paste0(homewd,"/pmx_gene_trees/polymerase"))
setwd(paste0(homewd,"/henipa_focus_L_tree"))
tree2 <-  read.tree("ref_right_trim.fasta.treefile") 
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
rooted.tree2<-drop.tip(rooted.tree2, "AET43339.1")
#load tree data prepared from elsewhere
dat2 <- read.csv(("ref_right_metadata.csv"), header = T, stringsAsFactors = F)

#check subgroup names
unique(dat2$Genus)
dat2$novel <- as.factor(dat2$novel)

colz = c("Orthoavulavirus"="darkorange","Paraavulavirus"="darkorange3","Metaavulavirus"="darkgoldenrod1","Jeilongvirus"="coral",
         "Morbillivirus"="coral3","Respirovirus"="deeppink1", "Parahenipavirus"="deeppink4",
         "Henipavirus"="firebrick1","Cynoglossusvirus"="grey50","Pararubulavirus"="firebrick4","Salemvirus"="indianred1",
         "Ferlavirus"="indianred3","Orthorubulavirus"="indianred4",
         "Hoplichthysvirus"="lightsalmon1","Narmovirus"="lightsalmon3","Aquaparamyxovirus"="lightpink2","Scoliodonvirus"="lightpink3",
         "Synodonvirus"="pink1","Unclassified"="black",
         "Sunshinevirus"="chartreuse1")

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

#and if there is an NA just drop it

#here NA in Isolate only:
dat2$new_label[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)] <- paste(dat2$Accession[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|", 
                                                                                                                 dat2$Organism_Name[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|", 
                                                                                                                 #dat2$Isolate[!is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|", 
                                                                                                                 dat2$source[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                 dat2$Country[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)], "|",
                                                                                                                 dat2$Collection_Date[is.na(dat2$Isolate) & !is.na(dat2$Accession) &!is.na(dat2$source) &!is.na(dat2$Country)])
#make sure to sort in order
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
tree.dat2$Host[tree.dat2$Host==0] <- "Non-bat host"
tree.dat2$Host[tree.dat2$Host==1] <- "Bat host"
tree.dat2$Host <- as.factor(tree.dat2$Host)
shapez = c("Bat host" =  17, "Non-bat host" = 19)
colz2 = c('1' =  "yellow", '0' = "white")

##Get the clade numbers so we can collapse unnecessary clades
ggtree(rooted.tree2) + geom_text(aes(label=node), hjust=-.3)

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
  geom_cladelabel(node= 465, label='HeV/NiV sister clade', offset=9.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=508, label='AndV clade', offset=9.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=167, label='AndV-2 clade', offset=9.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=498, label='AndV-3.1 clade', offset=9.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=491, label='AndV-3.2 clade', offset=9.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=352, label='AngV clade', offset=9.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=361, label='AngV-2.1 clade', offset=9.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=368, label='AngV-2.2 clade', offset=9.8, fontsize=4, color="black", align=TRUE)+
  
  geom_cladelabel(node = 379, label = "Hendra/Nipah henipavirus clade (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 474, label = "Ghanaian bat henipavirus clade (collapsed) | Eidolon helvum | Zambia (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = , label = "Unclassified henipavirus clade | Eidolon helvum | Ghana and DRC (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 517, label = "Parahenipavirus clade (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 632, label = "Rubulavirus clade (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 342, label = "Unclassified bat paramyxovirus clade (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 510, label = "Unclassified bat paramyxovirus clade (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 527, label = "Unclassified bat paramyxovirus clade (collapsed)",offset=0.1, fontsize = 4, color="black")+
  xlim(c(0,14))
p2

p2.2<-ggtree::collapse(p2, 379)+geom_point2(aes(subset=(node==379)), size=4, shape=22, fill="firebrick1")
p2.3<-ggtree::collapse(p2.2, 474)+geom_point2(aes(subset=(node==474)), size=4, shape=22, fill="firebrick1")
#p2.4<-collapse(p2.3, )+geom_point2(aes(subset=(node==)), size=4, shape=22, fill="firebrick1")
#p2.5<-collapse(p2.4, )+geom_point2(aes(subset=(node==)), size=4, shape=22, fill="firebrick1")
p2.6<-ggtree::collapse(p2.3, 517)+geom_point2(aes(subset=(node==517)), size=4, shape=22, fill="deeppink4")
p2.7<-ggtree::collapse(p2.6, 632)+geom_point2(aes(subset=(node==632)), size=4, shape=22, fill="firebrick4")
p2.8<-ggtree::collapse(p2.7, 342)+geom_point2(aes(subset=(node==342)), size=4, shape=22, fill="black")
p2.9<-ggtree::collapse(p2.8, 510)+geom_point2(aes(subset=(node==510)), size=4, shape=22, fill="black")
p2.1<-ggtree::collapse(p2.9, 527)+geom_point2(aes(subset=(node==527)), size=4, shape=22, fill="black")
p2.1

p2.1<- p2.1 +  
  geom_hilight(node=465, fill="deeppink3", alpha=0.3, extendto=12.8)+ #HeV/NiV sister
  geom_hilight(node=508, fill="darkcyan", alpha=0.3, extendto=12.8)+ #AndV
  geom_hilight(node=167, fill="darkturquoise", alpha=0.3, extendto=12.8)+ #AndV2
  geom_hilight(node=498, fill="dodgerblue", alpha=0.3, extendto=12.8)+ #AndV3.1
  geom_hilight(node=491, fill="dodgerblue", alpha=0.3, extendto=12.8)+ #AndV3.2
  geom_hilight(node=352, fill="springgreen3", alpha=0.3, extendto=12.8)+ #AngV
  geom_hilight(node=361, fill="chartreuse", alpha=0.3, extendto=12.8) + #AngV2.1
  geom_hilight(node=368, fill="chartreuse", alpha=0.3, extendto=12.8) #AngV2.2

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





