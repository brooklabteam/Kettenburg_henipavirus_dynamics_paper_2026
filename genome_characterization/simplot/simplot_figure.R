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
library(gggenes)
library(grid)
library(ggplotify)
library(ggtree)
library(ape)
library(ggnewscale)
library(gdata)
library(phylotools)
library(phylobase)

#set directory
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/genome_characterization/simplot"
setwd(paste0(homewd))


#nucleotide first
andv_map <- read.csv(file = "nt_align_AndV_noUTRtest.csv", header = T, stringsAsFactors = F)
head(andv_map)

#move to long
long.sim_nt <- melt(andv_map, id.vars = c("pointer"), measure.vars = c("AndV_ROU148","AndV_ROU173","AndV_ROU125"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "AndV_ROU148"] <- "AndV_ROU148"
long.sim_nt$accession[long.sim_nt$accession == "AndV_ROU173"] <- "AndV_ROU173"
long.sim_nt$accession[long.sim_nt$accession == "AndV_ROU125"] <- "AndV_ROU125"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c(
  "AndV_ROU125",
  "AndV_ROU148", 
  "AndV_ROU173"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

title<-expression(paste(italic("Reference: YuV2")))

colzpalette<-c("grey","grey30","black")

#Plot nucleotide
andv_map_nt <- ggplot(long.sim_nt) + 
  annotate("rect", xmin=163, xmax=1844, ymin=0, ymax=1, alpha=0.6,  fill="#D29292")+
  annotate("rect", xmin=2299, xmax=4522, ymin=0, ymax=1, alpha=0.6,  fill="#525832")+
  annotate("rect", xmin=5047, xmax=6102, ymin=0, ymax=1, alpha=0.6,  fill="#F98F5B")+
  annotate("rect", xmin=6368, xmax=8136, ymin=0, ymax=1, alpha=0.6,  fill="#BD9A5A")+
  annotate("rect", xmin=8679, xmax=11408, ymin=0, ymax=1, alpha=0.6,  fill="#F17983")+
  annotate("rect", xmin=11862, xmax=18859, ymin=0, ymax=1, alpha=0.6,  fill="#6B8C81")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Nucleotide identity")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 10),
        legend.title = element_blank(),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 1))+
  scale_color_manual(values=colzpalette) + 
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0));andv_map_nt

#GhVlike
ghv_map <- read.csv(file = "nt_align_GhVlike_noUTRtest.csv", header = T, stringsAsFactors = F)
head(ghv_map)

#move to long
long.sim_nt <- melt(ghv_map, id.vars = c("pointer"), measure.vars = c("GhVlike"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "GhVlike"] <- "GhVlike"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c(
  "GhVlike"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

title<-expression(paste(italic("Reference: GhV")))

colzpalette<-c("black","grey30","grey")
ghvlike_map_nt <- ggplot(long.sim_nt) + 
  annotate("rect", xmin=1, xmax=1745, ymin=0, ymax=1, alpha=0.6,  fill="#D29292")+
  annotate("rect", xmin=2213, xmax=4887, ymin=0, ymax=1, alpha=0.6,  fill="#525832")+
  annotate("rect", xmin=5347, xmax=6403, ymin=0, ymax=1, alpha=0.6,  fill="#F98F5B")+
  annotate("rect", xmin=7112, xmax=9123, ymin=0, ymax=1, alpha=0.6,  fill="#BD9A5A")+
  annotate("rect", xmin=9390, xmax=11633, ymin=0, ymax=1, alpha=0.6,  fill="#F17983")+
  annotate("rect", xmin=12345, xmax=19379, ymin=0, ymax=1, alpha=0.6,  fill="#6B8C81")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
    theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Nucleotide identity")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 10),
        legend.title = element_blank(),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 1))+
  scale_color_manual(values=colzpalette) + 
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0));ghvlike_map_nt

##then the aa plots
#then aa
andv_map <- read.csv(file = "aa_align_AndV_noUTRtest.csv", header = T, stringsAsFactors = F)
head(andv_map)

#move to long
long.sim_aa <- melt(andv_map, id.vars = c("pointer"), measure.vars = c("AndV_ROU148","AndV_ROU173","AndV_ROU125"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "AndV_ROU148"] <- "AndV_ROU148"
long.sim_aa$accession[long.sim_aa$accession == "AndV_ROU173"] <- "AndV_ROU173"
long.sim_aa$accession[long.sim_aa$accession == "AndV_ROU125"] <- "AndV_ROU125"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c( 
  "AndV_ROU125",
  "AndV_ROU148", 
  "AndV_ROU177"))
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

colzpalette<-c("black","grey30","grey")

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
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 10),
        legend.title = element_blank(),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 1))+
  ggtitle(title)+
  scale_color_manual(values=colzpalette) + 
  coord_cartesian(ylim=c(0,1.02))+
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
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 10),
        legend.title = element_blank(),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 1))+
  ggtitle(title)+
  scale_color_manual(values=colzpalette) + 
  coord_cartesian(ylim=c(0,1.02))+
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
genomeplot_ghvlike<-genome1ghvlike/genome2ghvlike+plot_layout(nrow=2,  heights = c(1, 1))
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

#put coverage and genome maps together
genomeplot_ROU173<-genome1ROU173/genome2ROU173+plot_layout(nrow=2,  heights = c(1, 1))
genomeplot_ROU173

genomeplot_ROU173<-as.ggplot(genomeplot_ROU173)
genomeplot_ROU173


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
genomeplot_andv<-genome1andv/genome2andv/genome1ROU173/genome2ROU173+plot_layout(nrow=4,  heights = c(1, 1,1,1))
genomeplot_andv

genomeplot_andv<-as.ggplot(genomeplot_andv)


#put simplot and genome plots together

simGhV<-plot_grid(ghv_map_aa, genomeplot_ghvlike, labels=c("",""),
               rel_widths=c(1,1), rel_heights = c(1,1),
               ncol=1,align="v", axis="lr", label_size = 23)
simGhV
















