#Script to make gene tree phylogenies

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

#fusion
setwd(paste0(homewd,"/pmx_gene_trees/fusion"))
tree1 <-  read.tree("fusion_left_align.fasta.treefile")
#drop root for aesthetics
rooted.tree1 <- root(tree1, which(tree1$tip.label == "YP_009094049.1"))
rooted.tree1<-drop.tip(rooted.tree1, "YP_009094049.1")
rooted.tree1<-drop.tip(rooted.tree1, "YP_010790476.1")
#load tree data prepared from elsewhere
dat1 <- read.csv(("F_metadata_left.csv"), header = T, stringsAsFactors = F)

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
# dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                   dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)],"|",
#                                                                                                                   dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
# 
# #and if there is an NA just drop it
# 
# #here NA in Isolate only:
# dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
#                                                                                                                  dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)],"|",
#                                                                                                                  dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
# 
# 
# 
dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                                  dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                                 dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])








#make sure to sort in order
colnames(dat1)[colnames(dat1)=="Accession"]="old_tip_label"
tree.dat1 <- data.frame(old_tip_label=rooted.tree1$tip.label, num =1:length(rooted.tree1$tip.label))
head(tree.dat1)
head(dat1)
tree.dat1 <- merge(tree.dat1, dat1, by = "old_tip_label", all.x = T, sort = F)
names(tree.dat1)
tree.dat1$tip_label <- tree.dat1$new_label
tree.dat1 <- dplyr::select(tree.dat1, tip_label, Isolate, Host, Host2, source, Country, Collection_Date, Genus, Organism_Name, novel, old_tip_label)
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

#plot
title<-expression(paste("          Fusion"))
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0.5,show.legend = T) +
  guides(color="none")+#
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
  geom_cladelabel(node=334, label='AndV clade', offset=4.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=170, label='GhV-like', offset=4.5, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=113, label='AndV-2 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=350, label='AndV-3 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=312, label='AngV clade', offset=4.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=313, label='AngV-2 clade', offset=4.5, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=31, label='AchPV2-like virus', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=268, label='ThkPV1-like viruses', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # 
  geom_cladelabel(node = 215, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 372, label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 398, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 392, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 389, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 358, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 361, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 336, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 355, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 382, label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 349, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 314, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  geom_cladelabel(node = 407, label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  xlim(c(0,12))

p1

p1.2<-ggtree::collapse(p1, 215)+geom_point2(aes(subset=(node==215)), size=4, shape=25, fill="black")
p1.3<-ggtree::collapse(p1.2, 372)+geom_point2(aes(subset=(node==372)), size=4, shape=25, fill="black")
p1.4<-ggtree::collapse(p1.3, 398)+geom_point2(aes(subset=(node==398)), size=4, shape=25, fill="black")
p1.5<-ggtree::collapse(p1.4, 392)+geom_point2(aes(subset=(node==392)), size=4, shape=25, fill="black")
p1.6<-ggtree::collapse(p1.5, 389)+geom_point2(aes(subset=(node==389)), size=4, shape=25, fill="black")
p1.7<-ggtree::collapse(p1.6, 358)+geom_point2(aes(subset=(node==358)), size=4, shape=25, fill="black")
p1.8<-ggtree::collapse(p1.7, 361)+geom_point2(aes(subset=(node==361)), size=4, shape=25, fill="black")
p1.9<-ggtree::collapse(p1.8, 336)+geom_point2(aes(subset=(node==336)), size=4, shape=25, fill="black")
p1.10<-ggtree::collapse(p1.9, 355)+geom_point2(aes(subset=(node==355)), size=4, shape=25, fill="black")
p1.11<-ggtree::collapse(p1.10, 382)+geom_point2(aes(subset=(node==382)), size=4, shape=25, fill="black")
p1.17<-ggtree::collapse(p1.11, 349)+geom_point2(aes(subset=(node==349)), size=4, shape=25, fill="black")
p1.18<-ggtree::collapse(p1.17, 314)+geom_point2(aes(subset=(node==314)), size=4, shape=25, fill="black")
p1.1<-ggtree::collapse(p1.18, 407)+geom_point2(aes(subset=(node==407)), size=4, shape=25, fill="black")
p1.1

p1.1<- p1.1 +  
  #geom_hilight(node=357, fill="deeppink3", alpha=0.3, extendto=12.8)+ #HeV/NiV sister
  geom_hilight(node=334, fill="darkslategray1", alpha=0.3, extendto=8.75)+ #AndV
  # geom_hilight(node=113, fill="darkturquoise", alpha=0.3, extendto=12.8)+ #AndV2
  # geom_hilight(node=350, fill="dodgerblue", alpha=0.3, extendto=12.8)+ #AndV3
  geom_hilight(node=312, fill="tomato1", alpha=0.3, extendto=8.75)+ #AngV
  geom_hilight(node=313, fill="darkorange1", alpha=0.3, extendto=8.75)+ #AngV2
  geom_hilight(node=170, fill="#00C094", alpha=0.3, extendto=8.75) #GhVlike
p1.1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree1)
p0.dat1 <- p0$data
p0.dat1$Bootstrap <- NA
Bootstrap<-p0.dat1$Bootstrap[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)] <- as.numeric(p0.dat1$label[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)])#fill with label

#add bootstrap values to original plot
fusiontree <- p1.1  %<+% p0.dat1 + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0, size=2)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1),
        plot.margin = unit(c(0,0,0,2), "cm"))
fusiontree

leg<-get_legend(fusiontree)


#glycoprotein
setwd(paste0(homewd,"/pmx_gene_trees/glycoprotein"))
tree1 <-  read.tree("glycoprotein_align.fasta.treefile")
#drop root for aesthetics
rooted.tree1 <- root(tree1, which(tree1$tip.label == "YP_009094050.1"))
rooted.tree1<-drop.tip(rooted.tree1, "YP_009094050.1")
#load tree data prepared from elsewhere
dat1 <- read.csv(("G_metadata_all.csv"), header = T, stringsAsFactors = F)

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
# dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                   dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)],"|",
#                                                                                                                   dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
# 
# #and if there is an NA just drop it
# 
# #here NA in Isolate only:
# dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
#                                                                                                                  dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)],"|",
#                                                                                                                  dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
# 
# 
# 
dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                                  dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                                 dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])








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

#plot
title<-expression(paste("          Glycoprotein"))
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0.5,show.legend = T) +
  guides(color="none")+
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
  geom_cladelabel(node=364, label='AndV clade', offset=5.5, fontsize=4, color="black", align=TRUE)+
  #geom_cladelabel(node=170, label='GhV-like', offset=4.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=166, label='AndV-2 clade', offset=5.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=165, label='AndV-3 clade', offset=5.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=360, label='AngV clade', offset=5.5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=359, label='AngV-2 clade', offset=5.5, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=31, label='AchPV2-like virus', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=268, label='ThkPV1-like viruses', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # 
  geom_cladelabel(node = 235, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 340, label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 325, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 332, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 338, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 296, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 297, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 303, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 353, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 315, label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 355, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 273, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  geom_cladelabel(node = 357, label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  xlim(c(0,12))

p1

p1.2<-ggtree::collapse(p1, 235)+geom_point2(aes(subset=(node==235)), size=4, shape=25, fill="black")
p1.3<-ggtree::collapse(p1.2, 340)+geom_point2(aes(subset=(node==340)), size=4, shape=25, fill="black")
p1.4<-ggtree::collapse(p1.3, 325)+geom_point2(aes(subset=(node==325)), size=4, shape=25, fill="black")
p1.5<-ggtree::collapse(p1.4, 332)+geom_point2(aes(subset=(node==332)), size=4, shape=25, fill="black")
p1.6<-ggtree::collapse(p1.5, 338)+geom_point2(aes(subset=(node==338)), size=4, shape=25, fill="black")
p1.7<-ggtree::collapse(p1.6, 296)+geom_point2(aes(subset=(node==296)), size=4, shape=25, fill="black")
p1.8<-ggtree::collapse(p1.7, 297)+geom_point2(aes(subset=(node==297)), size=4, shape=25, fill="black")
p1.9<-ggtree::collapse(p1.8, 303)+geom_point2(aes(subset=(node==303)), size=4, shape=25, fill="black")
p1.10<-ggtree::collapse(p1.9, 353)+geom_point2(aes(subset=(node==353)), size=4, shape=25, fill="black")
p1.11<-ggtree::collapse(p1.10, 315)+geom_point2(aes(subset=(node==315)), size=4, shape=25, fill="black")
p1.17<-ggtree::collapse(p1.11, 355)+geom_point2(aes(subset=(node==355)), size=4, shape=25, fill="black")
p1.18<-ggtree::collapse(p1.17, 273)+geom_point2(aes(subset=(node==273)), size=4, shape=25, fill="black")
p1.1<-ggtree::collapse(p1.18, 357)+geom_point2(aes(subset=(node==357)), size=4, shape=25, fill="black")
p1.1

p1.1<- p1.1 +  
  #geom_hilight(node=357, fill="deeppink3", alpha=0.3, extendto=12.8)+ #HeV/NiV sister
  geom_hilight(node=364, fill="darkslategray1", alpha=0.3, extendto=8.75)+ #AndV
  geom_hilight(node=166, fill="dodgerblue1", alpha=0.3, extendto=8.75)+ #AndV2
  geom_hilight(node=165, fill="dodgerblue4", alpha=0.3, extendto=8.75)+ #AndV3
  geom_hilight(node=360, fill="darkorange", alpha=0.3, extendto=8.75)+ #AngV
  geom_hilight(node=359, fill="tomato1", alpha=0.3, extendto=8.75) #AngV2
  #geom_hilight(node=170, fill="orange", alpha=0.3, extendto=8.75) #GhVlike
p1.1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree1)
p0.dat1 <- p0$data
p0.dat1$Bootstrap <- NA
Bootstrap<-p0.dat1$Bootstrap[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)] <- as.numeric(p0.dat1$label[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)])#fill with label

#add bootstrap values to original plot
glycotree <- p1.1  %<+% p0.dat1 + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0, size=2)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1),
        plot.margin = unit(c(0,0,0,2), "cm"))
glycotree




#matrix
setwd(paste0(homewd,"/pmx_gene_trees/matrix"))
tree1 <-  read.tree("matrix_right_align.fasta.treefile")
#drop root for aesthetics
rooted.tree1 <- root(tree1, which(tree1$tip.label == "YP_009094048.1"))
rooted.tree1<-drop.tip(rooted.tree1, "YP_009094048.1")
rooted.tree1<-drop.tip(rooted.tree1, "YP_010790482.1")
#load tree data prepared from elsewhere
dat1 <- read.csv(("M_metadata_right.csv"), header = T, stringsAsFactors = F)

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
# dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                   dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)],"|",
#                                                                                                                   dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
# 
# #and if there is an NA just drop it
# 
# #here NA in Isolate only:
# dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
#                                                                                                                  dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)],"|",
#                                                                                                                  dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
# 
# 
# 
dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                                  dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                                 dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])








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

#plot
title<-expression(paste("          Matrix"))
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0.5,show.legend = T) +
  guides(color="none")+
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
  geom_cladelabel(node=353, label='AndV clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=35, label='GhV-like', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=141, label='AndV-2 clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=355, label='AndV-3 clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=348, label='AngV clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=133, label='AngV-2 clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=31, label='AchPV2-like virus', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=268, label='ThkPV1-like viruses', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # 
  geom_cladelabel(node = 253, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 377, label = "Orthorubulavirus/Pararubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 399, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 408, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 397, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 375, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 366, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 356, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 225, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  #geom_cladelabel(node = , label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 220, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 228, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  geom_cladelabel(node = 415, label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  xlim(c(0,12))

p1

p1.2<-ggtree::collapse(p1, 253)+geom_point2(aes(subset=(node==253)), size=4, shape=25, fill="black")
p1.3<-ggtree::collapse(p1.2, 377)+geom_point2(aes(subset=(node==377)), size=4, shape=25, fill="black")
p1.4<-ggtree::collapse(p1.3, 399)+geom_point2(aes(subset=(node==399)), size=4, shape=25, fill="black")
p1.5<-ggtree::collapse(p1.4, 408)+geom_point2(aes(subset=(node==408)), size=4, shape=25, fill="black")
p1.6<-ggtree::collapse(p1.5, 397)+geom_point2(aes(subset=(node==397)), size=4, shape=25, fill="black")
p1.7<-ggtree::collapse(p1.6, 375)+geom_point2(aes(subset=(node==375)), size=4, shape=25, fill="black")
p1.8<-ggtree::collapse(p1.7, 366)+geom_point2(aes(subset=(node==366)), size=4, shape=25, fill="black")
p1.9<-ggtree::collapse(p1.8, 356)+geom_point2(aes(subset=(node==356)), size=4, shape=25, fill="black")
p1.10<-ggtree::collapse(p1.9, 225)+geom_point2(aes(subset=(node==225)), size=4, shape=25, fill="black")
#p1.11<-ggtree::collapse(p1.10, 382)+geom_point2(aes(subset=(node==382)), size=4, shape=25, fill="black")
p1.17<-ggtree::collapse(p1.10, 220)+geom_point2(aes(subset=(node==220)), size=4, shape=25, fill="black")
p1.18<-ggtree::collapse(p1.17, 228)+geom_point2(aes(subset=(node==228)), size=4, shape=25, fill="black")
p1.1<-ggtree::collapse(p1.18, 415)+geom_point2(aes(subset=(node==415)), size=4, shape=25, fill="black")
p1.1

p1.1<- p1.1 +  
  #geom_hilight(node=, fill="deeppink3", alpha=0.3, extendto=8.75)+ #HeV/NiV sister
  geom_hilight(node=353, fill="darkslategray1", alpha=0.3, extendto=8.75)+ #AndV
  geom_hilight(node=141, fill="dodgerblue1", alpha=0.3, extendto=8.75)+ #AndV2
  geom_hilight(node=355, fill="dodgerblue4", alpha=0.3, extendto=8.75)+ #AndV3
  geom_hilight(node=348, fill="tomato1", alpha=0.3, extendto=8.75)+ #AngV
  geom_hilight(node=133, fill="darkorange1", alpha=0.3, extendto=8.75)+ #AngV2
  geom_hilight(node=35, fill="#00C094", alpha=0.3, extendto=8.75) #GhVlike
p1.1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree1)
p0.dat1 <- p0$data
p0.dat1$Bootstrap <- NA
Bootstrap<-p0.dat1$Bootstrap[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)] <- as.numeric(p0.dat1$label[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)])#fill with label

#add bootstrap values to original plot
matrixtree <- p1.1  %<+% p0.dat1 + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0, size=2)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1),
        plot.margin = unit(c(0,0,0,2), "cm"))
matrixtree



#nucleoprotein
setwd(paste0(homewd,"/pmx_gene_trees/nucleoprotein"))
tree1 <-  read.tree("nucleoprotein_left_align.fasta.treefile")
#drop root for aesthetics
rooted.tree1 <- root(tree1, which(tree1$tip.label == "YP_009094046.1"))
rooted.tree1<-drop.tip(rooted.tree1, "YP_009094046.1")
rooted.tree1<-drop.tip(rooted.tree1, "YP_010790480.1")
rooted.tree1<-drop.tip(rooted.tree1, "YP_010790470.1")
rooted.tree1<-drop.tip(rooted.tree1, "KEL277_NODE_2164")
#load tree data prepared from elsewhere
dat1 <- read.csv(("N_metadata_left.csv"), header = T, stringsAsFactors = F)

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
# dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                   dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)],"|",
#                                                                                                                   dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
# 
# #and if there is an NA just drop it
# 
# #here NA in Isolate only:
# dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
#                                                                                                                  dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)],"|",
#                                                                                                                  dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
# 
# 
# 
dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                                  dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                                 dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])








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

#plot
title<-expression(paste("          Nucleoprotein"))
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0.5,show.legend = T) +
  guides(color="none")+
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
  geom_cladelabel(node=315, label='AndV clade', offset=1.7, fontsize=4, color="black", align=TRUE)+
  #geom_cladelabel(node=, label='GhV-like', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=109, label='AndV-2 clade', offset=1.7, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=312, label='AndV-3 clade', offset=1.7, fontsize=4, color="black", align=TRUE)+
  #geom_cladelabel(node=, label='AngV clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=320, label='AngV-2 clade', offset=1.7, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=31, label='AchPV2-like virus', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=268, label='ThkPV1-like viruses', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # 
  geom_cladelabel(node = 216, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 370, label = "Orthorubulavirus/Pararubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 399, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 392, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 407, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 367, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 358, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 349, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 347, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  #geom_cladelabel(node = , label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 342, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 321, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  geom_cladelabel(node = 409, label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  xlim(c(0,12))
p1

p1.2<-ggtree::collapse(p1, 216)+geom_point2(aes(subset=(node==216)), size=4, shape=25, fill="black")
p1.3<-ggtree::collapse(p1.2, 370)+geom_point2(aes(subset=(node==370)), size=4, shape=25, fill="black")
p1.4<-ggtree::collapse(p1.3, 399)+geom_point2(aes(subset=(node==399)), size=4, shape=25, fill="black")
p1.5<-ggtree::collapse(p1.4, 392)+geom_point2(aes(subset=(node==392)), size=4, shape=25, fill="black")
p1.6<-ggtree::collapse(p1.5, 407)+geom_point2(aes(subset=(node==407)), size=4, shape=25, fill="black")
p1.7<-ggtree::collapse(p1.6, 367)+geom_point2(aes(subset=(node==367)), size=4, shape=25, fill="black")
p1.8<-ggtree::collapse(p1.7, 358)+geom_point2(aes(subset=(node==358)), size=4, shape=25, fill="black")
p1.9<-ggtree::collapse(p1.8, 349)+geom_point2(aes(subset=(node==349)), size=4, shape=25, fill="black")
p1.10<-ggtree::collapse(p1.9, 347)+geom_point2(aes(subset=(node==347)), size=4, shape=25, fill="black")
#p1.11<-ggtree::collapse(p1.10, 382)+geom_point2(aes(subset=(node==382)), size=4, shape=25, fill="black")
p1.17<-ggtree::collapse(p1.10, 342)+geom_point2(aes(subset=(node==342)), size=4, shape=25, fill="black")
p1.18<-ggtree::collapse(p1.17, 321)+geom_point2(aes(subset=(node==321)), size=4, shape=25, fill="black")
p1.1<-ggtree::collapse(p1.18, 409)+geom_point2(aes(subset=(node==409)), size=4, shape=25, fill="black")
p1.1

p1.1<- p1.1 +  
  #geom_hilight(node=, fill="deeppink3", alpha=0.3, extendto=8.75)+ #HeV/NiV sister
  geom_hilight(node=315, fill="darkslategray1", alpha=0.3, extendto=8.75)+ #AndV
  geom_hilight(node=109, fill="dodgerblue1", alpha=0.3, extendto=8.75)+ #AndV2
  geom_hilight(node=312, fill="dodgerblue4", alpha=0.3, extendto=8.75)+ #AndV3
  #geom_hilight(node=348, fill="springgreen3", alpha=0.3, extendto=8.75)+ #AngV
  geom_hilight(node=320, fill="darkorange1", alpha=0.3, extendto=8.75) #AngV2
  #geom_hilight(node=35, fill="orange", alpha=0.3, extendto=8.75) #GhVlike
p1.1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree1)
p0.dat1 <- p0$data
p0.dat1$Bootstrap <- NA
Bootstrap<-p0.dat1$Bootstrap[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)] <- as.numeric(p0.dat1$label[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)])#fill with label

#add bootstrap values to original plot
nucleotree <- p1.1  %<+% p0.dat1 + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0, size=2)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1),
        plot.margin = unit(c(0,0,0,2), "cm"))
nucleotree



#phosphoprotein
setwd(paste0(homewd,"/pmx_gene_trees/phosphoprotein"))
tree1 <-  read.tree("phosphoprotein_long_align.fasta.treefile")
#drop root for aesthetics
rooted.tree1 <- root(tree1, which(tree1$tip.label == "YP_009094047.1"))
rooted.tree1<-drop.tip(rooted.tree1, "YP_009094047.1")
#rooted.tree1<-drop.tip(rooted.tree1, "YP_010790482.1")
#load tree data prepared from elsewhere
dat1 <- read.csv(("P_metadata_long.csv"), header = T, stringsAsFactors = F)

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
# dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
#                                                                                                                   dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)],"|",
#                                                                                                                   dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
# 
# #and if there is an NA just drop it
# 
# #here NA in Isolate only:
# dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
#                                                                                                                  dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)],"|",
#                                                                                                                  dat1$Isolate[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])
# 
# 
# 
dat1$new_label[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|",
                                                                                                                  dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat1$new_label[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)] <- paste(dat1$Accession[is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)], "|", 
                                                                                                                 dat1$Organism_Name[!is.na(dat1$Isolate) & !is.na(dat1$Accession) &!is.na(dat1$source) &!is.na(dat1$Country)])








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

#plot
title<-expression(paste("          Phosphoprotein"))
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0.5,show.legend = T) +
  guides(color="none")+
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
  geom_cladelabel(node=398, label='AndV clade', offset=5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=191, label='GhV-like', offset=5, fontsize=4, color="black", align=TRUE)+
  #geom_cladelabel(node=, label='AndV-2 clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  #geom_cladelabel(node=, label='AndV-3 clade', offset=4.8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=403, label='AngV clade', offset=5, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=202, label='AngV-2 clade', offset=5, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=31, label='AchPV2-like virus', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=268, label='ThkPV1-like viruses', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # 
  geom_cladelabel(node = 263, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 324, label = "Orthorubulavirus/Pararubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 344, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 355, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 353, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 322, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 362, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 383, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 373, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = , label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 377, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 301, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  #geom_cladelabel(node = , label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  xlim(c(0,12))

p1

p1.2<-ggtree::collapse(p1, 263)+geom_point2(aes(subset=(node==263)), size=4, shape=25, fill="black")
p1.3<-ggtree::collapse(p1.2, 324)+geom_point2(aes(subset=(node==324)), size=4, shape=25, fill="black")
p1.4<-ggtree::collapse(p1.3, 344)+geom_point2(aes(subset=(node==344)), size=4, shape=25, fill="black")
p1.5<-ggtree::collapse(p1.4, 355)+geom_point2(aes(subset=(node==355)), size=4, shape=25, fill="black")
p1.6<-ggtree::collapse(p1.5, 353)+geom_point2(aes(subset=(node==353)), size=4, shape=25, fill="black")
p1.7<-ggtree::collapse(p1.6, 322)+geom_point2(aes(subset=(node==322)), size=4, shape=25, fill="black")
p1.8<-ggtree::collapse(p1.7, 362)+geom_point2(aes(subset=(node==362)), size=4, shape=25, fill="black")
p1.9<-ggtree::collapse(p1.8, 383)+geom_point2(aes(subset=(node==383)), size=4, shape=25, fill="black")
p1.10<-ggtree::collapse(p1.9, 373)+geom_point2(aes(subset=(node==373)), size=4, shape=25, fill="black")
#p1.11<-ggtree::collapse(p1.10, 382)+geom_point2(aes(subset=(node==382)), size=4, shape=25, fill="black")
p1.17<-ggtree::collapse(p1.10, 377)+geom_point2(aes(subset=(node==377)), size=4, shape=25, fill="black")
p1.1<-ggtree::collapse(p1.17, 301)+geom_point2(aes(subset=(node==301)), size=4, shape=25, fill="black")
#p1.1<-ggtree::collapse(p1.18, 415)+geom_point2(aes(subset=(node==415)), size=4, shape=25, fill="black")
p1.1

p1.1<- p1.1 +  
  #geom_hilight(node=, fill="deeppink3", alpha=0.3, extendto=8.75)+ #HeV/NiV sister
  geom_hilight(node=398, fill="darkslategray1", alpha=0.3, extendto=8.75)+ #AndV
  #geom_hilight(node=141, fill="darkturquoise", alpha=0.3, extendto=8.75)+ #AndV2
  #geom_hilight(node=355, fill="dodgerblue", alpha=0.3, extendto=8.75)+ #AndV3
  geom_hilight(node=403, fill="tomato1", alpha=0.3, extendto=8.75)+ #AngV
  geom_hilight(node=202, fill="darkorange1", alpha=0.3, extendto=8.75)+ #AngV2
  geom_hilight(node=191, fill="#00C094", alpha=0.3, extendto=8.75) #GhVlike
p1.1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree1)
p0.dat1 <- p0$data
p0.dat1$Bootstrap <- NA
Bootstrap<-p0.dat1$Bootstrap[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)] <- as.numeric(p0.dat1$label[(length(tree.dat1$tip_label)+1):length(p0.dat1$label)])#fill with label

#add bootstrap values to original plot
phostree <- p1.1  %<+% p0.dat1 + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0, size=2)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1),
        plot.margin = unit(c(0,0,0,2), "cm"))
phostree


#put together
phylos <- cowplot::plot_grid(nucleotree, phostree, matrixtree,fusiontree, glycotree,leg,
                             rel_widths = c(1,1,1,1,1,1), align = "h", axis = "tb",
                             ncol=3, nrow = 2, labels = c("A", "B","C","D","E"), label_size = 22)+
                            theme(plot.background = element_rect(fill = "white", colour = NA)); phylos

#genetrees <- cowplot::plot_grid(phylos, legtree, ncol=2, nrow = 1, rel_widths = c(1,.2), align = "h", axis = "tb"); genetrees

homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/final_figures_tables/gwen_tmp/"
ggsave(file = paste0(homewd, "FigX_genetrees.png"),
       plot= phylos,
       units="mm",
       width=350,
       height=150,
       scale=2,
       dpi=500)



