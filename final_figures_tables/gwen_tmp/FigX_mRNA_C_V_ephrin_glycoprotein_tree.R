#Script to make genome characterization figure

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
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/genome_characterization/residue_alignments"
setwd(paste0(homewd))

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
        plot.margin = margin(0, 2, 0, 2, unit = "cm"))+
  #legend.position = "bottom", legend.direction = "horizontal")+
  #ggtitle("mRNA editing site")+
  #xlab("amino acid position\n(relative to AndV)")+ 
  scale_x_continuous(breaks = c(2,8), labels = c(1214, 1218),
                     sec.axis = dup_axis(breaks=c(8), labels=c("mRNA editing\n site")))
mRNAmsa.plot

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
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
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
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
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
                ncol=1,align="hv", axis="l", label_size = 23) +theme(plot.margin = margin(0, 0, 0, 10, "pt"));pair

pair<-compare_V/compare_W; pair
pair<-as.ggplot(pair)

#start with the residue alignments
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/genome_characterization/"

#load the msa of ephrin binding conserved sites
ephrinmsa = tidy_msa(msa="ephrin_region_alignment.fasta")
head(ephrinmsa)

#assign new names to each
ephrinmsa$name <- as.character(ephrinmsa$name)
unique(ephrinmsa$name)
ephrinmsa$name[ephrinmsa$name=="KEL86spleen"] <- "PX127395 GhV-like"
ephrinmsa$name[ephrinmsa$name=="CAT017_Plate_42_A10_NODE_34_czid2023"] <- "PX127396 AngV-2"
ephrinmsa$name[ephrinmsa$name=="KEL234_NODE_79_czid2019"] <- "PX127397 AngV-2"
ephrinmsa$name[ephrinmsa$name=="KEL281B_NODE_113_czid2019"] <- "PX127399 AndV-2"
ephrinmsa$name[ephrinmsa$name=="KEL380_Plate_42_B8_NODE_44_czid2023"] <- "PX127400 AngV"
ephrinmsa$name[ephrinmsa$name=="ROU125_NODE_1_czid2023"] <- "PX127342 AndV"
ephrinmsa$name[ephrinmsa$name=="ROU148_manual_assembly"] <- "PX127343 AndV"
ephrinmsa$name[ephrinmsa$name=="ROU143_NODE_287_czid2023"] <- "PX127401 AndV-3"

unique(ephrinmsa$name)
names(ephrinmsa)[names(ephrinmsa)=="character"] <- "aa"

# colz = scales::hue_pal()(length(unique(ephrinmsa$aa)))
# names(colz) <- unique(ephrinmsa$aa)
# names(colz)[names(colz)=="-"] <- "white"

colz = c("A" = "paleturquoise1", "C"="palevioletred1","D"="plum1","E"="purple1","F"="orange",
         "G"="goldenrod1","H"="orchid1","I"="slateblue1","K"="sienna3","L"="royalblue2",
         "N"="yellow1","Q"="tomato1","R"="turquoise1","S"="thistle1","T"="chocolate", "M"="plum4",
         "V"="deeppink4","W"="firebrick3","Y"="lightsteelblue1","-"="white")

colz2 = c("A" = "goldenrod1", "C"="white","D"="white","E"="goldenrod1","F"="white",
          "G"="white","H"="white","I"="white","K"="white","L"="white",
          "N"="white","Q"="goldenrod1","R"="white","S"="white","T"="goldenrod1",
          "V"="white","W"="goldenrod1","Y"="white","-"="white")

#and control the order

ephrinmsa$name <- factor(ephrinmsa$name, levels = rev(c("HeV", "NiV", "YuV1","CedV", "GhV","PX127395 GhV-like",
                                                        "YuV2","PX127342 AndV","PX127343 AndV",
                                                        "PX127399 AndV-2", "PX127401 AndV-3", "SGV","AngV", "PX127400 AngV", 
                                                        "PX127396 AngV-2","PX127397 AngV-2",                         
                                                        "MojV","GAKV","DARV","LayV")))

ephrinmsa.plot <- ggplot(data=ephrinmsa) + 
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
        axis.text.y = element_text(colour=c("black","black","black","black","darkorange1","darkorange1",
                                            "#F8766D","#F8766D","#F8766D","dodgerblue4","dodgerblue1", "cyan4","cyan4","cyan4",
                                            "#00C094","#00C094","black","black", "black","black"), face="bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.text = element_text(size=10),
        plot.margin = margin(0, 0.1, 0, 0.4, unit = "cm"))+
  #legend.position = "bottom", legend.direction = "horizontal")+
  #ggtitle("Ephrin binding residues")+
  scale_x_continuous(breaks = c(1,10,17,23), labels = c(502, 529,556,580), 
                     sec.axis = dup_axis(breaks=c(3,4,11,12,13,14,18,24), labels=c("*","*","*","*","*","*","*","*")))

ephrinmsa.plot


#now add glycoprotein tree

##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/phylo_trees/"
setwd(paste0(homewd,"/glycoprotein"))
tree1 <-  read.tree("pmxGright.treefile")
#drop root for aesthetics
rooted.tree1 <- root(tree1, which(tree1$tip.label == "YP_009094050.1"))
rooted.tree1<-drop.tip(rooted.tree1, "YP_009094050.1")
#load tree data prepared from elsewhere
dat1 <- read.csv(("G_refseq_metadata.csv"), header = T, stringsAsFactors = F)

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

colz = c("Orthoavulavirus"="black","Paraavulavirus"="black","Metaavulavirus"="black","Jeilongvirus"="black","Parajeilongvirus"="black",
         "Morbillivirus"="black","Respirovirus"="black", "Parahenipavirus"="black",
         "Henipavirus"="deepskyblue2","Cynoglossusvirus"="black","Pararubulavirus"="goldenrod","Salemvirus"="black",
         "Ferlavirus"="black","Orthorubulavirus"="black",
         "Hoplichthysvirus"="black","Narmovirus"="black","Aquaparamyxovirus"="black","Scoliodonvirus"="black",
         "Synodonvirus"="black","Unclassified"="black", "Unclassified henipavirus"="black", "Tupaivirus"="black", "Narmovirus"="black",
         "Hippocavirus"="black", "Bovinavirus"="black","Sunshinevirus"="black")

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

#Refseq alignment
p1 <- ggtree(rooted.tree1) %<+% tree.dat1 + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0.5,show.legend = T) +
  #scale_fill_manual(values=colz) +
  guides(color="none")+
  scale_color_manual(values=colz, guide =
                       guide_legend(label.theme = element_text(angle = 0, face = "italic")))+
  scale_shape_manual(values=shapez, guide =
                       guide_legend(label.theme = element_text(angle = 0, face = "italic"), nrow=3)) +
  new_scale_color() +
  geom_tiplab(aes(color = novel, show.legend=F), geom = "text", Genus="Helvetica", size=4, nudge_x=0.05) +
  guides(color="none")+#
  #ggtitle(title)+
  scale_color_manual(values=c("0"="black","1"="red")) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  geom_cladelabel(node=142, label='AndV clade', offset=4.6, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=139, label='AndV-2 clade', offset=4.6, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=138, label='AndV-3 clade', offset=4.6, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=150, label='AngV clade', offset=4.6, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=149, label='AngV-2 clade', offset=4.6, fontsize=4, color="black", align=TRUE)+
  # 
  #geom_cladelabel(node = , label = "Hendra/Nipah virus clade (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  #geom_cladelabel(node = , label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 191, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 181, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 178, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 270, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 259, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 220, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 210, label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 201, label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black") +
  geom_cladelabel(node = 271, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  geom_cladelabel(node = 154, label = "Narmovirus/Salemvirus/Morbillivirus/Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4,fontface=3, color="black")+
  xlim(c(0,12))

p1

p1<-flip(p1, 151, 173)
#p1<-flip(p1, 178, 180)

#p1.2<-ggtree::collapse(p1, 236)+geom_point2(aes(subset=(node==236)), size=4, shape=25, fill="deepskyblue2")
#p1.3<-ggtree::collapse(p1.2, 344)+geom_point2(aes(subset=(node==344)), size=4, shape=25, fill="firebrick4")
p1.4<-ggtree::collapse(p1, 191)+geom_point2(aes(subset=(node==191)), size=4, shape=25, fill="black")
p1.5<-ggtree::collapse(p1.4, 181)+geom_point2(aes(subset=(node==181)), size=4, shape=25, fill="black")
p1.6<-ggtree::collapse(p1.5, 178)+geom_point2(aes(subset=(node==178)), size=4, shape=25, fill="black")
p1.7<-ggtree::collapse(p1.6, 270)+geom_point2(aes(subset=(node==270)), size=4, shape=25, fill="black")
p1.8<-ggtree::collapse(p1.7, 259)+geom_point2(aes(subset=(node==259)), size=4, shape=25, fill="black")
p1.9<-ggtree::collapse(p1.8, 220)+geom_point2(aes(subset=(node==220)), size=4, shape=25, fill="black")
p1.10<-ggtree::collapse(p1.9, 210)+geom_point2(aes(subset=(node==210)), size=4, shape=25, fill="black")
p1.11<-ggtree::collapse(p1.10, 201)+geom_point2(aes(subset=(node==201)), size=4, shape=25, fill="black")
p1.17<-ggtree::collapse(p1.11, 271)+geom_point2(aes(subset=(node==271)), size=4, shape=25, fill="black")
p1.1<-ggtree::collapse(p1.17, 154)+geom_point2(aes(subset=(node==154)), size=4, shape=25, fill="black")
p1.1

p1.1<- p1.1 +  
  geom_hilight(node=142, fill="darkslategray1", alpha=0.3, extendto=9)+ #AndV
  geom_hilight(node=138, fill="dodgerblue1", alpha=0.3, extendto=9)+ #AndV3
  geom_hilight(node=139, fill="dodgerblue4", alpha=0.3, extendto=9)+ #AndV2
  geom_hilight(node=150, fill="#F8766D", alpha=0.3, extendto=9)+ #AngV
  geom_hilight(node=149, fill="darkorange1", alpha=0.3, extendto=9) #AngV2
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
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1),
        plot.margin = unit(c(0.5,0,0,0.5), "cm"))
polytree
























#put everything together
library(patchwork)




part1<-cowplot::plot_grid(mRNAmsa.plot, pair, 
                         labels=c("A","B"),
                         rel_widths=c(1,4), rel_heights = c(1,2),
                         ncol=1,nrow=2,align="hv", axis="lr", label_size = 23)+
  theme(plot.background = element_rect(fill = "white", colour = NA)); part1
part1<-as.ggplot(part1)

fig<-(part1|ephrinmsa.plot)+ plot_layout(widths = c(2, 1))+plot_annotation(tag_levels = list(c('','C'))) &
  theme(plot.tag = element_text(size = 23, face="bold"), 
        plot.background = element_rect(fill = "white", colour = NA));fig





part1<-cowplot::plot_grid(mRNAmsa.plot, ephrinmsa.plot, 
                          labels=c("A","C"),
                          rel_widths=c(1,1), rel_heights = c(1,1),
                          ncol=2,nrow=1,align="hv", axis="lr", label_size = 23)+
  theme(plot.background = element_rect(fill = "white", colour = NA)); part1
part1<-as.ggplot(part1)

fig<-cowplot::plot_grid(part1, pair, 
                        labels=c("","B"),
                        rel_widths=c(1,1), rel_heights = c(1,1),
                        ncol=1,nrow=2,align="hv", axis="lr", label_size = 23)+
  theme(plot.background = element_rect(fill = "white", colour = NA)); fig



#save trees
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/final_figures_tables/gwen_tmp/"
ggsave(file = paste0(homewd, "FigX_mRNA_C_V_ephrin.png"),
       plot= fig,
       units="mm",
       width=200,
       height=100,
       scale=2,
       dpi=500)







