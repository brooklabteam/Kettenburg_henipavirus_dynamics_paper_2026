#Script to make supplemental AndV figure

##Load packages
rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gggenes)
library(devtools)
library(gridExtra)
library(cowplot)
library(patchwork)
library(ggh4x)
library(grid)
library(ggplotify)
library(ggtree)
library(ape)
library(ggnewscale)
library(gdata)
library(phylotools)
library(phylobase)
library(seqinr)
library(ggmsa)
library(gdata)
library(cowplot)
library(lubridate)
library(tidyr)
library(reshape2)
library(Hmisc)
library(epitools)


#set directory
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/genome_characterization/full_genome_coverage"
setwd(paste0(homewd))

##Coverage plots first
dat2 <- read.csv(file= "AndV_merged_coverage.csv", header = T, stringsAsFactors = F) #load data
unique(dat2$virus)

library(ggbreak)

colz<-c("PX127342 AndV"="rosybrown4","PX127343 AndV"="tomato4","PX127344 AndV"="black")

title<-expression(paste("AndV genome coverage"))
coverage_andv<-  ggplot(dat2, aes(x=Position, y=rPM, color=virus))+ 
  annotate("rect", xmin=111, xmax=1718, ymin=0, ymax=Inf, alpha=0.6,  fill="#D29292")+
  annotate("rect", xmin=2114, xmax=3823, ymin=0, ymax=Inf, alpha=0.6,  fill="#525832")+
  annotate("rect", xmin=4384, xmax=5415, ymin=0, ymax=Inf, alpha=0.6,  fill="#F98F5B")+
  annotate("rect", xmin=5750, xmax=7435, ymin=0, ymax=Inf, alpha=0.6,  fill="#BD9A5A")+
  annotate("rect", xmin=8000, xmax=10716, ymin=0, ymax=Inf, alpha=0.6,  fill="#F17983")+
  annotate("rect", xmin=11300, xmax=18270, ymin=0, ymax=Inf, alpha=0.6,  fill="#6B8C81")+
  geom_line(linewidth=1) +
  theme_linedraw() +
  scale_color_manual(values=colz)+
  #facet_grid(virus~., scales="free")+
  theme(legend.position = "top")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million) log-scaled")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  scale_y_log10()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.title=element_blank(),
        axis.text = element_text(size=12), #axis.title.x = element_blank(), 
        axis.title = element_text(size=14), strip.background = element_rect(fill="white"), 
        plot.title = element_text(size=16),
        strip.text = element_text(size=16), panel.grid = element_blank())
coverage_andv

##Now make genome maps to go with the coverage plots
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
  scale_x_continuous(breaks=c(0,3000,6000,9000,12000, 15000, 17000),
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
map_combo<-subset(map,type=="ROU125")
map_combo_pep<-subset(map_pep,type=="ROU125")
map_combo_pep1<-subset(map_combo_pep,string=="first")
map_combo_feat<-subset(map_feat,type=="ROU125")
map_combo_pep2<-subset(map_combo_pep,string=="second")

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
  scale_x_continuous(breaks=c(0,3000,6000,9000,12000, 15000, 17000),
                     labels = c(0,3000,6000,9000,12000, 15000, 17000),expand=c(0,0))+
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
genomeplot_andv<-coverage_andv/genome1andv/genome2andv/genome1ROU173/genome2ROU173+plot_layout(nrow=5,  heights = c(1, 0.2, 0.2, 0.2, 0.2))
genomeplot_andv

genomeplot_andv<-as.ggplot(genomeplot_andv)



#Now do the phylogeny part
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/phylo_trees"

#fusion files first
setwd(paste0(homewd,"/full_genome_tree"))

#load the trees and root them
tree1 <-  read.tree("fullgenome_align.fasta.treefile")
rooted.tree1 <- root(tree1, which(tree1$tip.label == "NC_025345.1"))

#Remove root from displaying, still calculates changes correctly without it
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


#take a glance
p <- ggtree(rooted.tree1) %<+% dat1 + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=2) + geom_nodelab(size=2) +  #scale_color_manual(values=colz) + 
  theme(legend.position = "none", legend.title = element_blank())
p #looks great

#now get new tip labels
dat1$old_tip_label <- dat1$tip_label
dat1$new_label <- NA

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

#shaped and highlighting
tree.dat1$Host[tree.dat1$Host==0] <- "Non-bat host"
tree.dat1$Host[tree.dat1$Host==1] <- "Bat host"
tree.dat1$Host <- as.factor(tree.dat1$Host)
shapez = c("Bat host" =  17, "Non-bat host" = 19)
colz2 = c('1' =  "yellow", '0' = "white")

##Get the clade numbers so we can collapse unnecessary clades
ggtree(rooted.tree1) + geom_text(aes(label=node), hjust=-.3)

##uncollapsed tree

#Full alignment
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +
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
  geom_cladelabel(node=424, label='AndV clade', offset=9.32, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node = 274, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 360, label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 350, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 340, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 337, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 332, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 323, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 393, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 390, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 100, label = "Ferlavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 161, label = "Synodonvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 159, label = "Cynoglossusvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 158, label = "Hoplichthysvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 160, label = "Scoliodonvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 162, label = "Salemvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 384, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 403, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 379, label = "Cynoglossusvirus/Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  xlim(c(0,14))
p1

p1.2<-ggtree::collapse(p1, 274)+geom_point2(aes(subset=(node==274)), size=4, shape=22, fill="firebrick1")
p1.3<-ggtree::collapse(p1.2, 360)+geom_point2(aes(subset=(node==360)), size=4, shape=22, fill="firebrick4")
p1.4<-ggtree::collapse(p1.3, 350)+geom_point2(aes(subset=(node==350)), size=4, shape=22, fill="darkgoldenrod1")
p1.5<-ggtree::collapse(p1.4, 340)+geom_point2(aes(subset=(node==340)), size=4, shape=22, fill="darkorange")
p1.6<-ggtree::collapse(p1.5, 337)+geom_point2(aes(subset=(node==337)), size=4, shape=22, fill="darkorange3")
p1.7<-ggtree::collapse(p1.6, 332)+geom_point2(aes(subset=(node==332)), size=4, shape=22, fill="lightpink2")
p1.8<-ggtree::collapse(p1.7, 323)+geom_point2(aes(subset=(node==323)), size=4, shape=22, fill="deeppink1")
p1.9<-ggtree::collapse(p1.8, 393)+geom_point2(aes(subset=(node==393)), size=4, shape=22, fill="coral")
p1.10<-ggtree::collapse(p1.9, 390)+geom_point2(aes(subset=(node==390)), size=4, shape=22, fill="lightsalmon3")
#p1.11<-collapse(p1.10, 100)+geom_point2(aes(subset=(node==100)), size=4, shape=22, fill="indianred3")
#p1.12<-collapse(p1.10, 161)+geom_point2(aes(subset=(node==161)), size=4, shape=22, fill="pink1")
#p1.13<-collapse(p1.10, 159)+geom_point2(aes(subset=(node==159)), size=4, shape=22, fill="grey50")
#p1.14<-collapse(p1.10, 158)+geom_point2(aes(subset=(node==158)), size=4, shape=22, fill="lightsalmon1")
#p1.15<-collapse(p1.10, 160)+geom_point2(aes(subset=(node==160)), size=4, shape=22, fill="lightpink3")
#p1.16<-collapse(p1.10, 162)+geom_point2(aes(subset=(node==162)), size=4, shape=22, fill="indianred1")
p1.17<-ggtree::collapse(p1.10, 384)+geom_point2(aes(subset=(node==384)), size=4, shape=22, fill="coral3")
p1.18<-ggtree::collapse(p1.17, 403)+geom_point2(aes(subset=(node==403)), size=4, shape=22, fill="deeppink4")
p1.1<-ggtree::collapse(p1.18, 379)+geom_point2(aes(subset=(node==379)), size=4, shape=22, fill="grey40")
p1.1

p1.1<- p1.1 + 
  geom_hilight(node=424, fill="darkcyan", alpha=0.3, extendto=12.8)
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
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1),
        axis.text = element_text(size=12), #axis.title.x = element_blank(), 
        axis.title = element_text(size=14), strip.background = element_rect(fill="white"), 
        strip.text = element_text(size=16), panel.grid = element_blank())
tree1


##Next get the phosphoprotein data
#start with the residue alignments
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/genome_characterization/"

#load the msa of mRNA editing sites
setwd(paste0(homewd, "residue_alignments")) 
mRNAmsa = tidy_msa(msa="mRNA_alignment_AndVonly.fasta")
head(mRNAmsa)

#assign new names to each
mRNAmsa$name <- as.character(mRNAmsa$name)
unique(mRNAmsa$name)
mRNAmsa$name[mRNAmsa$name=="ROU125_NODE_1_czid2023"] <- "PX127342 AndV"
mRNAmsa$name[mRNAmsa$name=="ROU148_manual_assembly"] <- "PX127343 AndV"
mRNAmsa$name[mRNAmsa$name=="ROU173_manual_assembly"] <- "PX127344 AndV"
mRNAmsa$name[mRNAmsa$name=="CedV"] <- "*CedV"

unique(mRNAmsa$name)
names(mRNAmsa)[names(mRNAmsa)=="character"] <- "aa"

colz = c("A" = "darkorange1", "G"="deepskyblue", "T"="white")

#and control the order

#subset without partial sequence
#mRNAmsa = subset(mRNAmsa, name !="KEL153-virus")

mRNAmsa$name <- factor(mRNAmsa$name, levels = rev(c("HeV","NiV", "GhV","YuV1", "PX127342 AndV","PX127343 AndV","PX127344 AndV",
                                                    "YuV2" ,"AngV","SGV","MojV",        
                                                    "GAKV","DARV","LayV","*CedV")))
mRNAmsa.plot <- ggplot(data=mRNAmsa) + 
  geom_tile(aes(x=position, y=name, fill=aa, color=aa),width=0.9, height=0.9, show.legend = F) +
  geom_text(aes(x=position, y=name, label=aa), show.legend = F)+
  theme_bw() + 
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz) +
  theme(axis.title.y = element_blank(), 
        #panel.background = element_rect(fill="black"),
        #panel.border = element_line(color="black", size=1),
        panel.grid = element_blank(),
        axis.text = element_text(size=12),
        title = element_text(size=16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour=c("red","black","black","black","black","black","black", "cyan4", "cyan4",
                                            "cyan4","cyan4","black","black","black","black"), face="bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.text = element_text(size=10),
        #plot.background = element_rect(fill="transparent"),
        plot.margin = margin(0, 0, 0, 0, unit = "cm"))+
  #legend.position = "bottom", legend.direction = "horizontal")+
  #ggtitle("mRNA editing site")+
  #xlab("amino acid position\n(relative to AndV)")+ 
  scale_x_continuous(breaks = c(2,8), labels = c(1214, 1218),
                     sec.axis = dup_axis(breaks=c(8), labels=c("mRNA editing\n site")))
mRNAmsa.plot


#load the msa of c-terminal region of the v protein
Vmsa = tidy_msa(msa="v_align_trim.fasta")
head(Vmsa)

#assign new names to each
Vmsa$name <- as.character(Vmsa$name)
unique(Vmsa$name)
Vmsa$name[Vmsa$name=="ROU125_NODE_1_czid2023_add1"] <- "PX127342 AndV_add1"
Vmsa$name[Vmsa$name=="ROU148_manual_assembly_add1"] <- "PX127343 AndV_add1"
Vmsa$name[Vmsa$name=="ROU173_manual_assembly_add1"] <- "PX127344 AndV_add1"

unique(Vmsa$name)
names(Vmsa)[names(Vmsa)=="character"] <- "aa"

# colz = scales::hue_pal()(length(unique(Vmsa$aa)))
# names(colz) <- unique(Vmsa$aa)
# names(colz)[names(colz)=="-"] <- "white"

colz = c("A" = "paleturquoise1", "C"="palevioletred1","D"="plum1","E"="purple1","F"="orange",
         "G"="goldenrod1","H"="orchid1","I"="slateblue1","K"="sienna3","L"="royalblue2","M"="grey",
         "N"="yellow1","P"="khaki1","Q"="tomato1","R"="turquoise1","S"="thistle1","T"="chocolate",
         "V"="deeppink4","W"="firebrick3","Y"="lightsteelblue1","-"="white")

colz2 = c("A" = "white", "C"="yellow1","D"="white","E"="white","F"="white",
          "G"="green1","H"="orange1","I"="white","K"="dodgerblue","L"="white","M"="white",
          "N"="white","P"="white","Q"="white","R"="white","S"="white","T"="white",
          "V"="white","W"="white","Y"="white","-"="white", "*"="red")


#and control the order

Vmsa$name <- factor(Vmsa$name, levels = rev(c( "HeV_add1","NiV_add1", "GhV_add1","YuV1_add1", "PX127342 AndV_add1","PX127343 AndV_add1","PX127344 AndV_add1",
                                               "YuV2_add1","AngV_add1","SGV_add1","MojV_add1",        
                                               "GAKV_add1","DARV_add1","LayV_add1")))
Vmsa.plot <- ggplot(data=Vmsa) + 
  geom_tile(aes(x=position, y=name, fill=aa, color=aa),width=0.9, height=0.9, show.legend = F) +
  geom_text(aes(x=position, y=name, label=aa), show.legend = F)+
  theme_bw() + 
  scale_fill_manual(values=colz2) +
  scale_color_manual(values=colz2) +
  theme(axis.title.y = element_blank(), 
        #panel.background = element_rect(fill="black"),
        #panel.border = element_line(color="black", size=1),
        panel.grid = element_blank(),
        axis.text = element_text(size=12),
        title = element_text(size=16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour=c("black","black","black","black","black","black","cyan4", "cyan4",
                                            "cyan4","cyan4","black","black","black","black","black","black"), face="bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.text = element_text(size=10),
        plot.margin = margin(0, 0, 0, 2, unit = "cm"))+
  #legend.position = "bottom", legend.direction = "horizontal")+
  #ggtitle("V protein alignment")+
  scale_x_continuous(breaks = c(1,20,40,60), labels = c(409,429,449,469), 
                     sec.axis = dup_axis(breaks=c(4,23,27,39,41,44,48,51), labels=c("*","*","*","*","*","*","*","*")))

Vmsa.plot

#Now add the V and W pairwise similarity plots
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
        legend.position = "none",
        strip.text = element_text(size=10))
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
        legend.position = "bottom",
        strip.text = element_text(size=10))
compare_W

#add together: 
pair<-plot_grid(compare_V, compare_W, labels=c("",""),
                rel_widths=c(1,1), rel_heights = c(1,1),
                ncol=1,align="hv", axis="l", label_size = 23)
pair



##Now add all the virus shedding data
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/virus_shedding/data"
datawd = "/Users/gwenddolenkettenburg/Desktop/brook_lab/malagasy_fruit_bats_population_level_work/molecular_results"

# plot of seasonality of viruss for analysis
##  Prep the data

# First, load the viral shedding data - there is also a file for all molecular pos and all mNGS but we won't use them here
shed.dat <- read.csv(file = paste0(homewd, "/henipavirus_shedding_data.csv"))
head(shed.dat)

#slim to just AndV, AngV, GhV, and GhV-like data: one for just being positive and one for the ct values
shed.dat.pos <- dplyr::select(shed.dat, sampleid, AndV_shedding_ur)

## merge the shedding data with the catching data

#and make into long format for prevalence
molecular.long.pos <- melt(shed.dat.pos, id.vars = c("sampleid"), value.name = "pos", variable.name = "virus")
head(molecular.long.pos)
molecular.long.pos$virus <- as.character(molecular.long.pos$virus)
molecular.long.pos$virus <- gsub(pattern="qPCR_", replacement = "", x=molecular.long.pos$virus, fixed = T)
names(molecular.long.pos)[1] <- "sampleid"
#complete cases
molecular.long.pos <- molecular.long.pos[complete.cases(molecular.long.pos),]

# link this to the catching data to plot by sampling session, add the aging data as well
catch.dat <- read.csv(file = paste0(datawd, "/catching_data_clean_asofApril2025_edited.csv"), header = T, stringsAsFactors = F)
barcode.dat<- read.csv(file = paste0(datawd, "/barcode_ids.csv"), header = T, stringsAsFactors = F)
age.dat<-read.csv(file = paste0(datawd, "/aged_bats.csv"), header = T, stringsAsFactors = F)

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, region, roost_site, bat_species, processing_date, sampleid, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new,
                          year, quarter, month, doy, dataset, combo_class)
age.dat <-dplyr::select(age.dat, sampleid, age)

#and merge the catching data with the barcode and age data first
merge.dat <- merge(catch.dat, barcode.dat, by="sampleid", all.x = T)
head(merge.dat)

merge.dat <- merge(merge.dat, age.dat, by="sampleid", all.x = T)
head(merge.dat)

#now add the viral shedding data to the merged catching data
merge.dat.pos <- merge(molecular.long.pos, merge.dat, by="sampleid", all.x = T)
head(merge.dat.pos)

#get the midpoint processing date per sampling session
merge.dat.pos$processing_date <- as.Date(merge.dat.pos$processing_date)
sampling.date <- ddply(merge.dat.pos, .(bat_species, sampling_session), summarise, min_date = min(processing_date), max_date = max(processing_date))
sampling.date$mid_date <- (sampling.date$max_date-sampling.date$min_date)/2 + sampling.date$min_date

#and slim to merge
sampling.merge <- dplyr::select(sampling.date, bat_species, sampling_session, mid_date)

merge.dat.pos <- merge(merge.dat.pos, sampling.merge, by = c("bat_species", "sampling_session"))
merge.dat.pos$year <- year(merge.dat.pos$processing_date)
merge.dat.pos$month <- as.character(month(merge.dat.pos$processing_date))
head(merge.dat.pos)

#what species?
unique(merge.dat.pos$bat_species) #NA is the owl, get rid of that

merge.dat.pos<-subset(merge.dat.pos, bat_species=="Eidolon dupreanum" | bat_species=="Rousettus madagascariensis" | bat_species=="Pteropus rufus")
unique(merge.dat.pos$bat_species)

####################################################################################################
####################################################################################################
## INTERANNUAL PLOTTING

# now look at prevalence by sampling session by species and virus
# plot data first

#slim to those that are positive with any target for now, can come back and change

#rename virus for simplicity
merge.dat.pos$virus[merge.dat.pos$virus=="AndV_shedding_ur"] <- "AndV"
#merge.dat.pos$virus[merge.dat.pos$virus=="AndV_shedding_ct_ur"] <- "AndV_ct"

# merge.dat.pos <- rename(merge.dat.pos, pos = AndV_shedding_ur)
# merge.dat.pos <- rename(merge.dat.pos, ct = AndV_shedding_ct_ur)

#subset into the three time series datasets
merge.dat.pos_old <- merge.dat.pos[merge.dat.pos$mid_date > "2014-10-13" & merge.dat.pos$mid_date < "2016-01-19", ]
merge.dat.pos_mid <- merge.dat.pos[merge.dat.pos$mid_date > "2018-02-15" & merge.dat.pos$mid_date < "2019-09-21", ]
merge.dat.pos_new <- merge.dat.pos[merge.dat.pos$mid_date > "2022-07-19" & merge.dat.pos$mid_date < "2023-05-12", ]

prev.sum.old <- ddply(merge.dat.pos_old, .(bat_species, region, virus, mid_date), summarise, N=length(pos), pos=sum(pos))
prev.sum.old$neg <- prev.sum.old$N - prev.sum.old$pos
prev.sum.old$prev <- prev.sum.old$pos/prev.sum.old$N
prev.sum.old$prev_lci <-binom.exact(x=prev.sum.old$pos, n=prev.sum.old$N)$lower
prev.sum.old$prev_uci <-binom.exact(x=prev.sum.old$pos, n=prev.sum.old$N)$upper
head(prev.sum.old)

prev.sum.mid <- ddply(merge.dat.pos_mid, .(bat_species, region,virus, mid_date), summarise, N=length(pos), pos=sum(pos))
prev.sum.mid$neg <- prev.sum.mid$N - prev.sum.mid$pos
prev.sum.mid$prev <- prev.sum.mid$pos/prev.sum.mid$N
prev.sum.mid$prev_lci <-binom.exact(x=prev.sum.mid$pos, n=prev.sum.mid$N)$lower
prev.sum.mid$prev_uci <-binom.exact(x=prev.sum.mid$pos, n=prev.sum.mid$N)$upper
head(prev.sum.mid)

prev.sum.new <- ddply(merge.dat.pos_new, .(bat_species, region,virus, mid_date), summarise, N=length(pos), pos=sum(pos))
prev.sum.new$neg <- prev.sum.new$N - prev.sum.new$pos
prev.sum.new$prev <- prev.sum.new$pos/prev.sum.new$N
prev.sum.new$prev_lci <-binom.exact(x=prev.sum.new$pos, n=prev.sum.new$N)$lower
prev.sum.new$prev_uci <-binom.exact(x=prev.sum.new$pos, n=prev.sum.new$N)$upper
head(prev.sum.new)

#define colors per virus

colz1 = c('gestation' = "yellow", 'lactation' = "greenyellow", 'winter' = "skyblue",
          
          'AndV' = "black")

filz = c('gestation' = "yellow", 'lactation' = "greenyellow", 'winter' = "skyblue",
         
         'AndV' = "black")

p1 <- ggplot(data=prev.sum.mid) + geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) +
  geom_line(aes(x=mid_date, y=prev, color=virus)) + scale_color_manual(values = colz1) +
  facet_grid(bat_species~.) + theme_bw() 
p1

# add in the timing of the reproductive season for both species
# from angelo and santino's paper
# eidolon, gestation is from: Aug-3 - Nov 15; lactation Nov 16- 2 feb
# rousettus, gestation is from: sept 11 - Dec 11; lactation Dec 12- 19 feb
# pteropus, gestation is from: July 7 - Aug 2; lactation Aug 3- Jan 21
# dry season = may 15 - october 10

#old dataset 2013-2016
lact.dat.old <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum",
                                               "Rousettus madagascariensis", "Rousettus madagascariensis",
                                               "Pteropus rufus", "Pteropus rufus"), 
                                 start = c(as.Date("2014-11-16"), as.Date("2015-11-16"), 
                                           as.Date("2014-12-12"), as.Date("2015-12-12"),
                                           as.Date("2014-08-13"), as.Date("2015-08-03")), 
                                 stop=c(as.Date("2015-02-02"), as.Date("2016-02-02"),
                                        as.Date("2015-02-19"), as.Date("2016-02-19"),
                                        as.Date("2015-01-21"), as.Date("2016-01-21")))
lact.dat.old$type = "lactation"
lact.dat.old$bat_sex = "female"

gest.dat.old <- cbind.data.frame(bat_species=c("Eidolon dupreanum",
                                               "Rousettus madagascariensis", 
                                               "Pteropus rufus"), 
                                 start = c(as.Date("2015-08-03"), 
                                           as.Date("2015-09-11"),
                                           as.Date("2015-07-07")), 
                                 stop=c(as.Date("2015-11-15"),
                                        as.Date("2015-12-11"),
                                        as.Date("2015-08-02")))
gest.dat.old$type = "gestation"
gest.dat.old$bat_sex="female"

#and winter
cal.dat.old <-  cbind.data.frame(bat_species=c("Eidolon dupreanum",
                                               "Rousettus madagascariensis",
                                               "Pteropus rufus"), 
                                 start = c(as.Date("2015-05-15"), 
                                           as.Date("2015-05-15"),
                                           as.Date("2015-05-15")), 
                                 stop=c(as.Date("2015-10-10"), 
                                        as.Date("2015-10-10"),
                                        as.Date("2015-10-10")))
cal.dat.old$type <- "winter"
cal.dat.old$bat_sex <-"male"

library(reshape2)

repro.dat.old <- rbind(lact.dat.old, gest.dat.old)
demo.dat.old <- rbind(gest.dat.old, lact.dat.old, cal.dat.old)




#mid dataset 2018-2020
lact.dat.mid <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum",
                                               "Rousettus madagascariensis","Rousettus madagascariensis", 
                                               "Pteropus rufus", "Pteropus rufus"), 
                                 start = c(as.Date("2018-11-16"), as.Date("2019-11-16"),
                                           as.Date("2018-12-12"), as.Date("2019-12-12"),
                                           as.Date("2018-08-03"), as.Date("2019-08-03")), 
                                 stop=c(as.Date("2019-02-02"), as.Date("2020-02-02") ,
                                        as.Date("2019-02-19"), as.Date("2020-02-19") ,
                                        as.Date("2019-01-21"), as.Date("2020-01-21")))
lact.dat.mid$type = "lactation"
lact.dat.mid$bat_sex ="female"

gest.dat.mid <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", 
                                               "Rousettus madagascariensis", "Rousettus madagascariensis", 
                                               "Pteropus rufus", "Pteropus rufus"), 
                                 start = c(as.Date("2018-08-03"), as.Date("2019-08-03"),
                                           as.Date("2018-09-11"), as.Date("2019-09-11"),
                                           as.Date("2018-07-07"), as.Date("2019-07-07")), 
                                 stop=c(as.Date("2018-11-15"), as.Date("2019-11-15"),
                                        as.Date("2018-12-11"), as.Date("2019-12-11"),
                                        as.Date("2018-08-02"), as.Date("2019-08-02")))
gest.dat.mid$type = "gestation"
gest.dat.mid$bat_sex ="female"

#and winter
cal.dat.mid <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", 
                                               "Rousettus madagascariensis", "Rousettus madagascariensis",
                                               "Pteropus rufus","Pteropus rufus"), 
                                 start = c(as.Date("2018-05-15"), as.Date("2019-05-15"),
                                           as.Date("2018-05-15"), as.Date("2019-05-15"),
                                           as.Date("2018-05-15"), as.Date("2019-05-15")), 
                                 stop=c(as.Date("2018-10-10"), as.Date("2019-10-10"),
                                        as.Date("2018-10-10"), as.Date("2019-10-10"),
                                        as.Date("2018-10-10"), as.Date("2019-10-10")))
cal.dat.mid$type <- "winter"
cal.dat.mid$bat_sex="male"

library(reshape2)

repro.dat.mid <- rbind(lact.dat.mid, gest.dat.mid)
demo.dat.mid <- rbind(gest.dat.mid, lact.dat.mid, cal.dat.mid)

#now add to plot

p2 <- ggplot(data=subset(prev.sum.mid)) + 
  geom_rect(data=demo.dat.mid, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) +
  geom_line(aes(x=mid_date, y=prev, color=virus)) + scale_color_manual(values = colz1) +
  facet_grid(bat_species~.) + theme_bw()

p2


#what if you did males and females separately?
prev.sum.sex.old <- ddply(merge.dat.pos_old, .(bat_species, region, bat_sex, virus, sampling_session,  mid_date), summarise, N=length(pos), pos=sum(pos))
prev.sum.sex.old$neg <- prev.sum.sex.old$N - prev.sum.sex.old$pos
prev.sum.sex.old$prev <- prev.sum.sex.old$pos/prev.sum.sex.old$N
prev.sum.sex.old$prev_lci <-binom.exact(x=prev.sum.sex.old$pos, n=prev.sum.sex.old$N)$lower
prev.sum.sex.old$prev_uci <-binom.exact(x=prev.sum.sex.old$pos, n=prev.sum.sex.old$N)$upper
head(prev.sum.sex.old)
prev.sum.sex.old = subset(prev.sum.sex.old, bat_sex=="female" | bat_sex=="male")

prev.sum.sex.mid <- ddply(merge.dat.pos_mid, .(bat_species, region, bat_sex, virus, sampling_session, mid_date), summarise, N=length(pos), pos=sum(pos))
prev.sum.sex.mid$neg <- prev.sum.sex.mid$N - prev.sum.sex.mid$pos
prev.sum.sex.mid$prev <- prev.sum.sex.mid$pos/prev.sum.sex.mid$N
prev.sum.sex.mid$prev_lci <-binom.exact(x=prev.sum.sex.mid$pos, n=prev.sum.sex.mid$N)$lower
prev.sum.sex.mid$prev_uci <-binom.exact(x=prev.sum.sex.mid$pos, n=prev.sum.sex.mid$N)$upper
head(prev.sum.sex.mid)
prev.sum.sex.mid = subset(prev.sum.sex.mid, bat_sex=="female" | bat_sex=="male")

prev.sum.sex.new <- ddply(merge.dat.pos_new, .(bat_species, region, bat_sex, virus, sampling_session, mid_date), summarise, N=length(pos), pos=sum(pos))
prev.sum.sex.new$neg <- prev.sum.sex.new$N - prev.sum.sex.new$pos
prev.sum.sex.new$prev <- prev.sum.sex.new$pos/prev.sum.sex.new$N
prev.sum.sex.new$prev_lci <-binom.exact(x=prev.sum.sex.new$pos, n=prev.sum.sex.new$N)$lower
prev.sum.sex.new$prev_uci <-binom.exact(x=prev.sum.sex.new$pos, n=prev.sum.sex.new$N)$upper
head(prev.sum.sex.new)
prev.sum.sex.new = subset(prev.sum.sex.new, bat_sex=="female" | bat_sex=="male")

p3 <- ggplot(data=subset(prev.sum.sex.mid)) + 
  geom_rect(data=demo.dat.mid, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + scale_color_manual(values = colz1) +
  geom_line(aes(x=mid_date, y=prev, color=virus)) +
  facet_grid(bat_species~bat_sex) + theme_bw()

p3


#and plot with data:
#add category for the time period
prev.sum.sex.old$dataset <- "2014-2016"
prev.sum.sex.mid$dataset <- "2018-2020"
prev.sum.old$dataset <- "2014-2016"
prev.sum.mid$dataset <- "2018-2020"

demo.dat.old$dataset<-"2014-2016"
demo.dat.mid$dataset<-"2018-2020"

#bind the time periods together for plotting
prev.sum.sex.plot <- rbind(prev.sum.sex.old, prev.sum.sex.mid)
prev.sum.plot <- rbind(prev.sum.old, prev.sum.mid)

#bind together the demographic events
demo.dat.plot <- rbind(demo.dat.old, demo.dat.mid)

inter <- ggplot(data=subset(prev.sum.plot)) +
  geom_rect(data=demo.dat.plot, aes(xmin=as.Date(start), xmax=as.Date(stop), ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(aes(x=mid_date, y=prev, size=N, color=virus)) + 
  scale_color_manual(values=colz1) + 
  scale_fill_manual(values=filz, name = "season") +
  guides(size = guide_legend(title = "N bats"))+
  geom_line(aes(x=mid_date, y=prev, color=virus), linetype=2) +
  facet_grid(bat_species~dataset, scales="free") + 
  #scale_x_date(limits = as.Date(c("2018-02-18","2019-09-20")))+
  theme_bw()+
  ylab("Prevalence") + xlab("Sampling date")+
  guides(color="none")+
  ggtitle("AndV interannual shedding")+
  theme(axis.text = element_text(size=12), #axis.title.x = element_blank(), 
        axis.title = element_text(size=14), strip.background = element_rect(fill="white"), 
        plot.title = element_text(size=16),
        strip.text.y=element_text(face="italic"),
        strip.text = element_text(size=16), panel.grid = element_blank())
inter


####################################################################################################
####################################################################################################
## INTRAANNUAL PLOTTING

# # remove bats from the dataset with no age class listed
molecular.dat = subset(merge.dat.pos, !is.na(bat_age_class))
molecular.dat = subset(molecular.dat, !is.na(bat_sex))
molecular.dat = subset(molecular.dat, !is.na(bat_species))

#factor the categorical data
molecular.dat$year <- as.factor(molecular.dat$year)
molecular.dat$bat_sex <- as.factor(molecular.dat$bat_sex)
molecular.dat$roost_site <- as.factor(molecular.dat$roost_site)
molecular.dat$region <- as.factor(molecular.dat$region)
molecular.dat$virus <- as.factor(molecular.dat$virus)
molecular.dat$bat_species <- as.factor(molecular.dat$bat_species)
molecular.dat$bat_age_class <- as.factor(molecular.dat$bat_age_class)
molecular.dat$young_of_year <- as.factor(molecular.dat$young_of_year)
molecular.dat$pupmom_pair <- as.factor(molecular.dat$pupmom_pair)
molecular.dat$region <- as.factor(molecular.dat$region)
molecular.dat$combo_class<-as.factor(molecular.dat$combo_class)
molecular.dat$month<-as.numeric(molecular.dat$month)

molecular.dat <- molecular.dat[molecular.dat$mid_date > "2014-10-13" & molecular.dat$mid_date < "2019-09-21", ]

##New seasonality timing
#month
lact.dat.month <- cbind.data.frame(bat_species=rep(c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"), each=2), 
                                   start = c(11.5,1, 11.4,1, 8.1,1), 
                                   stop=c(12,2.05,12,2.7,12.9,1))
lact.dat.month$type = "lactation"
lact.dat.month$bat_sex="female"
gest.dat.month <- cbind.data.frame(bat_species=c("Eidolon dupreanum","Rousettus madagascariensis", "Pteropus rufus"), 
                                   start = c(8.15, 9.45, 7.2),
                                   stop=c(11.5, 11.39, 8.1))
gest.dat.month$type = "gestation"
gest.dat.month$bat_sex = "female"
#and winter
cal.dat.month <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"), 
                                   start = c(5.5,5.5,5.5),
                                   stop=c(10.3,10.3,10.3))
cal.dat.month$type <- "winter"
cal.dat.month$bat_sex="male"

library(reshape2)

repro.dat.month <- rbind(lact.dat.month, gest.dat.month)
demo.dat.month <- rbind(gest.dat.month, lact.dat.month, cal.dat.month)

#for plotting
merge.plot <- ddply(molecular.dat,.(bat_species, year, month, virus), summarise, N=length(pos), pos=sum(pos))
merge.plot$molecular <-   merge.plot$pos/merge.plot$N

merge.plot.noyear <- ddply(molecular.dat,.(bat_species, month, virus), summarise, N=length(pos), pos=sum(pos))
merge.plot.noyear$molecular <-   merge.plot.noyear$pos/merge.plot.noyear$N

merge.plot.sex <- ddply(molecular.dat,.(bat_species, bat_sex, year, month, virus), summarise, N=length(pos), pos=sum(pos))
merge.plot.sex$molecular <-   merge.plot.sex$pos/merge.plot.sex$N

merge.plot.sex.noyear <- ddply(molecular.dat,.(bat_species, bat_sex, month, virus), summarise, N=length(pos), pos=sum(pos))
merge.plot.sex.noyear$molecular <-   merge.plot.sex.noyear$pos/merge.plot.sex.noyear$N

shapez=c("2013"=16,"2014"=17,"2015"=15, "2016"=3, "2018"=7, "2019"=8, "2020"=18, "2022"=9, "2023"=10)

###plotting by day of year
#all datasets
intra <- 
  ggplot(merge.plot) + 
  #ggplot(merge.plot.noyear) + 
  geom_rect(data=demo.dat.month, aes(xmin=start, xmax=stop, ymin=0, ymax=1, fill=type), alpha=.3) +
  geom_point(data=merge.plot, aes(x=month, y=molecular, size=N, color=virus, shape=year)) +
  #geom_point(data=merge.plot.noyear, aes(x=month, y=molecular, size=N, color=virus)) +
  facet_grid(.~bat_species) + 
  theme_bw() +
  scale_color_manual(values=colz1) +
  scale_fill_manual(values=filz, name = "season") +
  guides(size = guide_legend(title = "N bats"))+  
  scale_shape_manual(values=shapez) +
  scale_x_continuous(breaks=c(2, 4, 6,8, 10, 12)) +
  theme(axis.text = element_text(size=12), #axis.title.x = element_blank(), 
        axis.title = element_text(size=14), 
        plot.title = element_text(size=16),
        strip.background = element_rect(fill="white"), 
        strip.text.x=element_text(face="italic"),
        strip.text = element_text(size=16), 
        panel.grid = element_blank()) +
  guides(color="none")+
  guides(fill="none")+
  guides(size="none")+
  ggtitle("AndV intraannual shedding")+
  coord_cartesian(expand = F) + ylab("Prevalence") + xlab("Month")
intra


##Put figure together
genome.coverage<-plot_grid(genomeplot_andv, tree1, labels=c("A","B"),rel_widths=c(1,1.3),
                ncol=2,  label_size = 23)
genome.coverage
genome.coverage<-as.ggplot(genome.coverage)

library(patchwork)

phospho.1<-mRNAmsa.plot+Vmsa.plot + plot_layout(widths = c(1, 6))
phospho.1

phospho<-phospho.1/pair+ plot_annotation(tag_levels = list(c('E','F','G'))) &
  theme(plot.tag = element_text(size = 23, face="bold"))
phospho
phospho<-as.ggplot(phospho)

seasonal<-plot_grid(inter, intra, labels=c("C","D"),rel_widths=c(1,1), rel_heights = c(3,1),
                           ncol=1,  label_size = 23)
seasonal
seasonal<-as.ggplot(seasonal)

part1<-plot_grid(genome.coverage, seasonal, labels=c("",""),rel_widths=c(2,1),
                    ncol=2,  label_size = 23)
part1
part1<-as.ggplot(part1)

Fig.AndV<-plot_grid(part1, phospho, labels=c("",""),rel_widths=c(1,1), rel_heights = c(1,1.4),
                ncol=1,  label_size = 23)
Fig.AndV







part1<-plot_grid(seasonal, phospho, labels=c("",""),rel_widths=c(1,1.2),
                 ncol=2,  label_size = 23)
part1
part1<-as.ggplot(part1)

Fig.AndV<-plot_grid(genome.coverage, part1, labels=c("",""),rel_widths=c(1,1), rel_heights = c(1,1.7),
                    ncol=1,  label_size = 23)
Fig.AndV


#save trees
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/final_figures_tables/supplemental/"
ggsave(file = paste0(homewd, "FigX_AndV_figure.png"),
       plot= Fig.AndV,
       units="mm",
       width=350,
       height=250,
       scale=2,
       dpi=500)







