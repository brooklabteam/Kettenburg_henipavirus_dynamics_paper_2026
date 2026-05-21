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
andv_map <- read.csv(file = "nt_align_AndV.csv", header = T, stringsAsFactors = F)
head(andv_map)

#move to long
long.sim_nt <- melt(andv_map, id.vars = c("pointer"), measure.vars = c("AndV_ROU148","AndV_ROU173","AngV","YuV2"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "AndV_ROU148"] <- "AndV_ROU148"
long.sim_nt$accession[long.sim_nt$accession == "AndV_ROU173"] <- "AndV_ROU173"
long.sim_nt$accession[long.sim_nt$accession == "AngV"] <- "AngV"
long.sim_nt$accession[long.sim_nt$accession == "YuV2"] <- "YuV2"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c(
  "AndV_ROU148",
  "AndV_ROU173", 
  "AngV", 
  "YuV2"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

title<-expression(paste(italic("Reference: AndV ROU125")))

colzpalette<-c("#3B9AB2","#EBCC2A","#F21A00","darkslategray2","darkorange1","firebrick4", "black")

#Plot nucleotide
andv_map_nt <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  #annotate("rect", xmin=10416, xmax=12260, ymin=0, ymax=1, alpha=1,  fill="grey30")+
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
ghv_map <- read.csv(file = "nt_align_GhVlike.csv", header = T, stringsAsFactors = F)
head(ghv_map)

#move to long
long.sim_nt <- melt(ghv_map, id.vars = c("pointer"), measure.vars = c("AngV","GhV"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "AngV"] <- "AngV"
long.sim_nt$accession[long.sim_nt$accession == "GhV"] <- "GhV"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c(
  "AngV",
  "GhV"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

title<-expression(paste(italic("Reference: GhVlike")))

ghvlike_map_nt <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  annotate("rect", xmin=563, xmax=659, ymin=0, ymax=1, alpha=1,  fill="grey30")+
  annotate("rect", xmin=1171, xmax=3955, ymin=0, ymax=1, alpha=1,  fill="grey30")+
  annotate("rect", xmin=4696, xmax=5900, ymin=0, ymax=1, alpha=1,  fill="grey30")+
  annotate("rect", xmin=6317, xmax=7776, ymin=0, ymax=1, alpha=1,  fill="grey30")+
  annotate("rect", xmin=9138, xmax=9314, ymin=0, ymax=1, alpha=1,  fill="grey30")+
  annotate("rect", xmin=10209, xmax=10442, ymin=0, ymax=1, alpha=1,  fill="grey30")+
  annotate("rect", xmin=10989, xmax=12860, ymin=0, ymax=1, alpha=1,  fill="grey30")+
  annotate("rect", xmin=15133, xmax=15634, ymin=0, ymax=1, alpha=1,  fill="grey30")+
  annotate("rect", xmin=17902, xmax=18063, ymin=0, ymax=1, alpha=1,  fill="grey30")+
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
andv_map <- read.csv(file = "aa_align_AndV.csv", header = T, stringsAsFactors = F)
head(andv_map)

#move to long
long.sim_aa <- melt(andv_map, id.vars = c("pointer"), measure.vars = c("AndV_ROU148","AndV_ROU173","AngV","YuV2"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "AndV_ROU148"] <- "AndV_ROU148"
long.sim_aa$accession[long.sim_aa$accession == "AndV_ROU173"] <- "AndV_ROU173"
long.sim_aa$accession[long.sim_aa$accession == "AngV"] <- "AngV"
long.sim_aa$accession[long.sim_aa$accession == "YuV2"] <- "YuV2"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c( 
  "AndV_ROU148",
  "AndV_ROU173", 
  "AngV", 
  "YuV2"))
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

title<-expression(paste(italic("Reference: AndV ROU125")))
## amino acid
andv_map_aa <- ggplot(long.sim_aa) + 
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
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055, 10000/3.055, 12000/3.055, 14000/3.055, 16000/3.055, 18000/3.055), 
                     labels = c(0,2000, 4000,6000,8000, 10000, 12000, 14000, 16000, 18000),
                     expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0)); andv_map_aa

























dat <- read.csv(file= "KEL86_KEL94_coverage.csv", header = T, stringsAsFactors = F) #load data

title<-expression(paste(italic("GhV-like genome coverage")))
coverage_ghvlike<-  ggplot(dat, aes(x=Position, y=rPM))+ 
  annotate("rect", xmin=50, xmax=4179, ymin=0, ymax=550, alpha=0.6,  fill="#D29292")+
  annotate("rect", xmin=4603, xmax=7636, ymin=0, ymax=550, alpha=0.6,  fill="#525832")+
  annotate("rect", xmin=8096, xmax=10940, ymin=0, ymax=550, alpha=0.6,  fill="#F98F5B")+
  annotate("rect", xmin=11664, xmax=16846, ymin=0, ymax=550, alpha=0.6,  fill="#BD9A5A")+
  annotate("rect", xmin=17373, xmax=23508, ymin=0, ymax=550, alpha=0.6,  fill="#F17983")+
  annotate("rect", xmin=24156, xmax=48178, ymin=0, ymax=550, alpha=0.6,  fill="#6B8C81")+
  geom_line(linewidth=1) +
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
        axis.title.x = element_blank(),
        axis.line.x = element_blank())
coverage_ghvlike

dat2 <- read.csv(file= "AndV_merged_coverage.csv", header = T, stringsAsFactors = F) #load data
unique(dat2$virus)

library(ggbreak)

colz<-c("PX127342 AndV"="rosybrown4","PX127343 AndV"="tomato4","PX127344 AndV"="black")

title<-expression(paste(italic("AndV genome coverage")))
coverage_andv<-  ggplot(dat2, aes(x=Position, y=rPM, color=virus))+ 
  annotate("rect", xmin=111, xmax=1718, ymin=0, ymax=Inf, alpha=0.4,  fill="#D29292")+
  annotate("rect", xmin=2114, xmax=3823, ymin=0, ymax=Inf, alpha=0.4,  fill="#525832")+
  annotate("rect", xmin=4284, xmax=5315, ymin=0, ymax=Inf, alpha=0.4,  fill="#F98F5B")+
  annotate("rect", xmin=5550, xmax=7235, ymin=0, ymax=Inf, alpha=0.4,  fill="#BD9A5A")+
  annotate("rect", xmin=7737, xmax=10316, ymin=0, ymax=Inf, alpha=0.4,  fill="#F17983")+
  annotate("rect", xmin=10749, xmax=17501, ymin=0, ymax=Inf, alpha=0.4,  fill="#6B8C81")+
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
        legend.title=element_blank())
coverage_andv

#ROU125
dat <- read.csv(file= "ROU125_coverage.csv", header = T, stringsAsFactors = F) #load data

title<-expression(paste(italic("PX127342 AndV genome coverage")))
coverage_rou125<-  ggplot(dat, aes(x=Position, y=rPM))+ 
  annotate("rect", xmin=111, xmax=1718, ymin=0, ymax=1.6, alpha=0.6,  fill="#D29292")+
  annotate("rect", xmin=2114, xmax=3823, ymin=0, ymax=1.6, alpha=0.6,  fill="#525832")+
  annotate("rect", xmin=4284, xmax=5315, ymin=0, ymax=1.6, alpha=0.6,  fill="#F98F5B")+
  annotate("rect", xmin=5550, xmax=7235, ymin=0, ymax=1.6, alpha=0.6,  fill="#BD9A5A")+
  annotate("rect", xmin=7737, xmax=10316, ymin=0, ymax=1.6, alpha=0.6,  fill="#F17983")+
  annotate("rect", xmin=10749, xmax=17501, ymin=0, ymax=1.6, alpha=0.6,  fill="#6B8C81")+
  geom_area(linewidth=0.5) +
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
        axis.title.x = element_blank(),
        axis.line.x = element_blank())
coverage_rou125

#ROU148
dat <- read.csv(file= "ROU148_coverage.csv", header = T, stringsAsFactors = F) #load data

title<-expression(paste(italic("PX127343 AndV genome coverage")))
coverage_rou148<-  ggplot(dat, aes(x=Position, y=rPM))+ 
  annotate("rect", xmin=156, xmax=1763, ymin=0, ymax=174, alpha=0.6,  fill="#D29292")+
  annotate("rect", xmin=2159, xmax=3868, ymin=0, ymax=174, alpha=0.6,  fill="#525832")+
  annotate("rect", xmin=4329, xmax=5360, ymin=0, ymax=174, alpha=0.6,  fill="#F98F5B")+
  annotate("rect", xmin=5595, xmax=7280, ymin=0, ymax=174, alpha=0.6,  fill="#BD9A5A")+
  annotate("rect", xmin=7782, xmax=10361, ymin=0, ymax=174, alpha=0.6,  fill="#F17983")+
  annotate("rect", xmin=10794, xmax=17546, ymin=0, ymax=174, alpha=0.6,  fill="#6B8C81")+
  geom_area(linewidth=0.5) +
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
        axis.title.x = element_blank(),
        axis.line.x = element_blank())
coverage_rou148

#ROU173
dat <- read.csv(file= "ROU173_coverage.csv", header = T, stringsAsFactors = F) #load data

title<-expression(paste(italic("PX127344 AndV genome coverage")))
coverage_rou173<-  ggplot(dat, aes(x=Position, y=rPM))+ 
  annotate("rect", xmin=102, xmax=1709, ymin=0, ymax=0.85, alpha=0.6,  fill="#D29292")+
  annotate("rect", xmin=2105, xmax=3814, ymin=0, ymax=0.85, alpha=0.6,  fill="#525832")+
  annotate("rect", xmin=4275, xmax=5306, ymin=0, ymax=0.85, alpha=0.6,  fill="#F98F5B")+
  annotate("rect", xmin=5546, xmax=7231, ymin=0, ymax=0.85, alpha=0.6,  fill="#BD9A5A")+
  annotate("rect", xmin=7733, xmax=10198, ymin=0, ymax=0.85, alpha=0.6,  fill="#F17983")+
  annotate("rect", xmin=10631, xmax=17383, ymin=0, ymax=0.85, alpha=0.6,  fill="#6B8C81")+
  geom_area(linewidth=0.5) +
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
        axis.title.x = element_blank(),
        axis.line.x = element_blank())
coverage_rou173

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
genomeplot_ghvlike<-coverage_ghvlike/genome1ghvlike/genome2ghvlike+plot_layout(nrow=3,  heights = c(1, 0.2, 0.2))
genomeplot_ghvlike

genomeplot_ghvlike<-as.ggplot(genomeplot_ghvlike)
genomeplot_ghvlike


#ROU125
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
genome1ROU125<-ggplot(map_combo_pep1, aes(xmin = from, xmax = to, y = molecule, fill=gene)) +
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
genome1ROU125

genome2ROU125<-ggplot(map_combo, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
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
genome2ROU125

#put coverage and genome maps together
genomeplot_ROU125<-coverage_rou125/genome1ROU125/genome2ROU125+plot_layout(nrow=3,  heights = c(1, 0.2, 0.2))
genomeplot_ROU125

genomeplot_ROU125<-as.ggplot(genomeplot_ROU125)
genomeplot_ROU125


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
genomeplot_ROU173<-coverage_rou173/genome1ROU173/genome2ROU173+plot_layout(nrow=3,  heights = c(1, 0.2, 0.2))
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
genomeplot_andv<-coverage_andv/genome1andv/genome2andv/genome1ROU173/genome2ROU173+plot_layout(nrow=5,  heights = c(1, 0.2, 0.2, 0.2, 0.2))
genomeplot_andv

genomeplot_andv<-as.ggplot(genomeplot_andv)

#Put the whole figure together 
Fig<-plot_grid(genomeplot_ghvlike, genomeplot_ROU125, genomeplot_ROU148, genomeplot_ROU173, labels=c("A","B","C","D"),
                 rel_widths=c(1,1), rel_heights = c(1,1),
                 ncol=2,align="hv", axis="l", label_size = 23)
Fig


Fig<-plot_grid(genomeplot_ghvlike, genomeplot_andv, labels=c("A","B"),
               rel_widths=c(1,1), rel_heights = c(1,1.5),
               ncol=1,align="hv", axis="l", label_size = 23)
Fig



#homewd= "/Users/gwenddolenkettenburg/Desktop/GhV_like_genome_figure/"
ggsave(file = paste0(homewd, "/genome_coverage.pdf"),
       plot= Fig,
       units="mm",  
       width=200, 
       height=140, 
       scale=2, 
       dpi=500)

