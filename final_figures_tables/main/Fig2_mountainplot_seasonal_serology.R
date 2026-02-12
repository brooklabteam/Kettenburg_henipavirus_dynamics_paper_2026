#Script to make figure 1, which shows the mountain plots of henipavirus serology and the polymerase phylogeny side by side

rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(mclust)
library(ggh4x)
library(ggcorrplot)
library(factoextra)
library(paran)
library(hopkins)
library(cluster)
library(corrr) 
library(FactoMineR)
library(ggpmisc)
library(heatmaply)
library(seqinr)
library(scales)
library(ggdendroplot)
library(ggimage)
library(ggtree)
library(ape)
library(ggnewscale)
library(gdata)
library(phylotools)
library(phylobase)
library(cowplot)

#use this script to calculate cutoffs for the luminex


#first record your wd for this repo (home) 
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/"

#now load in the merged luminex data from the "data" folder

#make sure this is latest catching file
dat <- read.csv(paste0(homewd, "/serology/data/lum_catch_merge.csv"), header = T, stringsAsFactors = F)
head(dat)


length(unique(dat$sampleid))
dat$sampleid[duplicated(dat$sampleid)]

test <- ddply(dat, .(bat_species, bat_sex), summarise, Ntot=length(unique(sampleid)))
test

#now for cutoffs, just slim to antigen data only and convert to long format
lumdat <- dplyr::select(dat, sampleid, bat_species, names(dat)[2:7])
head(lumdat)
lumdat.long <- melt(lumdat)
head(lumdat.long)
names(lumdat.long)[names(lumdat.long)=="variable"] <- "antigen"
names(lumdat.long)[names(lumdat.long)=="value"] <- "MFI"
unique(lumdat.long$antigen)
sort(unique(lumdat.long$MFI))
class(lumdat.long$MFI)
lumdat.long <- subset(lumdat.long, !is.na(MFI))

# # #convert any negative values to 1
# lumdat.long$MFI[lumdat.long$MFI<0] <- 1
# lumdat.long$log10MFI <- log10(lumdat.long$MFI)

lum.split <- dlply(lumdat.long,.(bat_species, antigen))

correct.negs <- function(df){
  abs.min <- abs(min(df$MFI))
  df$MFI <- df$MFI + abs.min + 1
  
  #then, also take the log10
  df$log10MFI <- log10(df$MFI)
  return(df)
}

lumdat.long <- data.table::rbindlist(lapply(lum.split, correct.negs))


#add family
lumdat.long$virus_family <- "Paramyxovirus"
unique(lumdat.long$virus_family)

#and genus information
lumdat.long$virus_genus <- "Henipavirus"
lumdat.long$virus_genus[lumdat.long$antigen=="MojV"] <- "Parahenipavirus"
unique(lumdat.long$virus_genus)

#factor antigens by relatedness
lumdat.long$antigen <- factor(lumdat.long$antigen, levels=c( "NiV", "HeV","CedV", "GhV", "AngV", "MojV"))

lumdat.long$bat_species <- factor(lumdat.long$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))
 
lumdat.long.plot <- lumdat.long

#and heatmap MFI by species/antigen combination 
library(ggh4x)
p1 <- ggplot(lumdat.long.plot) + 
  geom_tile(aes(x=sampleid, y=antigen, fill=MFI)) + theme_bw()+
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(size=12), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        legend.key.width =  unit(c(1.5), "cm"),
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        legend.position = "bottom") +
  scale_fill_viridis_c(trans="log10") + facet_nested(virus_family+virus_genus~bat_species, scales = "free", space = "free");p1

# and also add a mountain plot for additional support
# here, with antigens arranged according to phylogenetic relatedness... 
# makes it a bit easier to see the relatedness between bats
# first pass is all antigens in order of true relatedness

# load aa alignment
henipa.align <- read.alignment(file = paste0(homewd,"serology/data/Henipavirus-Alignment.fasta"), format = "fasta")
# replace names
new.names <- strsplit(henipa.align$nam,"_")
length.names <- lapply(new.names, length)
new.names <- mapply(FUN='[', new.names, length.names)
henipa.align$nam <- new.names

# make distance matrix
dist.henipa <- dist.alignment(henipa.align)

colorz = c('AngV' = "#F8766D", 'HeV' ="#A58AFF", 
           'NiV' = "#DB72FB", 'GhV' = "#00C094", 
           'CedV'= "#53B400", "MojV" = "#619CFF")


#and put them together
image.df2 <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"), antigen=rep("NiV", 3), MFI=rep(21000,3), image=c(paste0(homewd,"/serology/images/eidolon_dupreanum_headshot.png"),paste0(homewd,"/serology/images/rousettus_madagascariensis_headshot.png"),paste0(homewd,"/serology/images/pteropus_rufus_headshot.png")))
image.df2$bat_species <- factor(image.df2$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))
image.df2$virus_family="Paramyxovirus"
image.df <- image.df2

lumdat.long$bat_species <- factor(lumdat.long$bat_species,  levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))

#some cutoffs were established by mixture modeling
cutoff.df <- cbind.data.frame(bat_species=c(rep("Eidolon dupreanum", 4), rep("Rousettus madagascariensis",3)), antigen=c("NiV", "HeV", "GhV", "AngV", "GhV", "BdbV", "RavV"))
cutoff.df$virus_family = "Paramyxovirus"
cutoff.df$MFI = c(936,2851,1141,1763,2642,5647,6617)
mean.cut <- ddply(cutoff.df,.(bat_species, virus_family), summarise, meanMFI = mean(MFI))

plot.cut <- cbind.data.frame(bat_species=rep(c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"),each=2),
                             virus_family=rep(c("Filovirus", "Paramyxovirus"), 3),
                             MFI= c(rep(1672,2), 2642, 6132, rep(1672,2))) #might come back and mixture model all of these


plot.cut$bat_species <- factor(plot.cut$bat_species, levels=c("Eidolon dupreanum", "Rousettus madagascariensis", "Pteropus rufus"))

hc.henipa <- hclust(dist.henipa)

#edit these a bit for ease of viewing
hc.henipa$height[2]<-hc.henipa$height[2]-.3
hc.henipa$height[3]<-hc.henipa$height[3]-.2
hc.henipa$height[4]<-hc.henipa$height[4]-.1

p.h.1 <- ggplot() + geom_dendro(hc.henipa, size=1) + theme_bw() +
      scale_y_reverse() + scale_x_reverse()+ 
      theme(panel.grid = element_blank(), 
            panel.border = element_blank(),
            axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), plot.margin = unit(c(0,1,0,3.3), "cm")) 
p.h.1
p.h.2 <- ggplot() + geom_dendro(hc.henipa, size=1) + theme_bw() +
  scale_y_reverse() + scale_x_reverse()+
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), plot.margin = unit(c(0,1,0,1), "cm")) 
p.h.3 <- ggplot() + geom_dendro(hc.henipa, size=1) + theme_bw() +
  scale_y_reverse() + scale_x_reverse()+
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), plot.margin = unit(c(0,1,0,1), "cm")) 

p.h.row <- cowplot::plot_grid(p.h.1,p.h.2,p.h.3, ncol = 3, nrow = 1, rel_widths = c(1.22,1,1)) + theme(plot.margin = unit(c(0,0,0,0), "cm")) 

p2bsave <- ggplot(subset(lumdat.long)) + coord_cartesian(ylim = c(0,27000)) +
  facet_nested(.~bat_species, scales = "free", space = "free") +
  geom_line(aes(x=antigen, y=MFI, group=sampleid), show.legend=F) + theme_bw()+
  geom_point(aes(x=antigen, y=MFI, color=antigen), show.legend=F) + 
  geom_hline(data = subset(plot.cut, virus_family=="Paramyxovirus"), aes(yintercept = MFI), color="red", linetype=2, size=.8) +
  scale_color_manual(values=colorz) + ylab("mean fluorescence intensity (MFI)\n") +
  #geom_dendro(rowclus) +
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size=16), axis.title.y = element_text(size=16), 
        axis.text.y = element_text(size=16), 
        strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  geom_image(data=subset(image.df, virus_family=="Paramyxovirus"), aes(x=antigen, y=c(23000, 24000,23000), image = image), size=.4, nudge_x = .7); p2bsave
mountain <- cowplot::plot_grid(p2bsave, p.h.row, ncol = 1,nrow = 2, rel_heights = c(1,.2), rel_widths = c(1,1));mountain
mountain<-as.ggplot(mountain)

mountain <- mountain + theme(plot.margin = unit(c(0.1,0,0,2), "cm"))
mountain

##Now plot the phylogenetic tree to go with it

##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/phylo_trees/"
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
tree2 <-  read.tree("ref_right_align.fasta.treefile") 
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
  geom_cladelabel(node=357, label='HeV/NiV sister clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=354, label='AndV clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=113, label='AndV-2 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=350, label='AndV-3 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=250, label='AngV clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=247, label='AngV-2 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=31, label='AchPV2-like virus', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=268, label='ThkPV1-like viruses', offset=8.55, fontsize=4, color="black", align=TRUE)+
  
  geom_cladelabel(node = 239, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 344, label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 293, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 288, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 301, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 255, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 256, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 315, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 312, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 280, label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 307, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 325, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 302, label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  xlim(c(0,14))

p1

p1.2<-ggtree::collapse(p1, 239)+geom_point2(aes(subset=(node==239)), size=4, shape=22, fill="firebrick1")
#p1.3<-ggtree::collapse(p1.2, 344)+geom_point2(aes(subset=(node==344)), size=4, shape=22, fill="firebrick4")
p1.4<-ggtree::collapse(p1.2, 293)+geom_point2(aes(subset=(node==293)), size=4, shape=22, fill="darkgoldenrod1")
p1.5<-ggtree::collapse(p1.4, 288)+geom_point2(aes(subset=(node==288)), size=4, shape=22, fill="darkorange")
p1.6<-ggtree::collapse(p1.5, 301)+geom_point2(aes(subset=(node==301)), size=4, shape=22, fill="darkorange3")
p1.7<-ggtree::collapse(p1.6, 255)+geom_point2(aes(subset=(node==255)), size=4, shape=22, fill="lightpink2")
p1.8<-ggtree::collapse(p1.7, 256)+geom_point2(aes(subset=(node==256)), size=4, shape=22, fill="deeppink1")
p1.9<-ggtree::collapse(p1.8, 315)+geom_point2(aes(subset=(node==315)), size=4, shape=22, fill="coral")
p1.10<-ggtree::collapse(p1.9, 312)+geom_point2(aes(subset=(node==312)), size=4, shape=22, fill="lightsalmon3")
p1.11<-ggtree::collapse(p1.10, 280)+geom_point2(aes(subset=(node==280)), size=4, shape=22, fill="indianred3")
p1.17<-ggtree::collapse(p1.11, 307)+geom_point2(aes(subset=(node==307)), size=4, shape=22, fill="coral3")
p1.18<-ggtree::collapse(p1.17, 325)+geom_point2(aes(subset=(node==325)), size=4, shape=22, fill="deeppink4")
p1.1<-ggtree::collapse(p1.18, 302)+geom_point2(aes(subset=(node==302)), size=4, shape=22, fill="lightpink3")
p1.1

p1.1<- p1.1 +  
  geom_hilight(node=357, fill="deeppink3", alpha=0.3, extendto=12.8)+ #HeV/NiV sister
  geom_hilight(node=354, fill="darkcyan", alpha=0.3, extendto=12.8)+ #AndV
  geom_hilight(node=113, fill="darkturquoise", alpha=0.3, extendto=12.8)+ #AndV2
  geom_hilight(node=350, fill="dodgerblue", alpha=0.3, extendto=12.8)+ #AndV3
  geom_hilight(node=250, fill="springgreen3", alpha=0.3, extendto=12.8)+ #AngV
  geom_hilight(node=247, fill="chartreuse", alpha=0.3, extendto=12.8)+ #AngV2
  geom_hilight(node=31, fill="grey40", alpha=0.3, extendto=12.8)+ #AngV
  geom_hilight(node=268, fill="grey40", alpha=0.3, extendto=12.8) #AngV2
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
        plot.title = element_text(hjust = -0.1),
        plot.margin = unit(c(0,0,0,2), "cm"))
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

###Put together final figure
phylo<-plot_grid(tree1, tree2, labels=c("B","C"),
                 rel_widths=c(1,1), rel_heights = c(1,1),
                 ncol=2,align="l", axis="l", label_size = 23)
phylo


phylo<-as.ggplot(phylo)

Fig1<-plot_grid(mountain, NULL, tree1, labels=c("A","","B"),
                rel_widths=c(1,1,1), rel_heights = c(0.7,0.03,1),
                ncol=1,align="h", axis="l", label_size = 23) +
                theme(plot.background = element_rect(fill = "white"))
Fig1

Fig1<-plot_grid(tree1, NULL, mountain, labels=c("A","","B"),
                rel_widths=c(1,1,1), rel_heights = c(1,0.03,0.7),
                ncol=1,align="h", axis="l", label_size = 23) +
  theme(plot.background = element_rect(fill = "white"))
Fig1


#save trees
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/final_figures_tables/main/"
ggsave(file = paste0(homewd, "Fig1_mountainplot_polymerase_phylogeny.png"),
       plot= Fig1,
       units="mm",  
       width=200, 
       height=160, 
       scale=2, 
       dpi=500)

homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/final_figures_tables/supplemental/"
ggsave(file = paste0(homewd, "FigX_africanbat_polymerase_phylogeny.png"),
       plot= tree2,
       units="mm",  
       width=200, 
       height=120, 
       scale=2, 
       dpi=500)

