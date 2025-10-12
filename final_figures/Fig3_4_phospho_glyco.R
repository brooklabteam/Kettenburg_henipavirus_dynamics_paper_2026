#Script to make figure 3

rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(ggtree)
library(seqinr)
library(ggmsa)
library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
#install.packages('gdata')
library(gdata)
library(phylotools)
library(phylobase)
library(cowplot)

#start with the residue alignments
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_circulating_clades/genome_characterization/"

#load the msa of mRNA editing sites
setwd(paste0(homewd, "residue_alignments")) 
mRNAmsa = tidy_msa(msa="mRNA_alignment.fasta")
head(mRNAmsa)

#assign new names to each
mRNAmsa$name <- as.character(mRNAmsa$name)
unique(mRNAmsa$name)
mRNAmsa$name[mRNAmsa$name=="ROU125_NODE_1_czid2023"] <- "PX127342 AndV"
mRNAmsa$name[mRNAmsa$name=="ROU148_manual_assembly"] <- "PX127343 AndV"
mRNAmsa$name[mRNAmsa$name=="ROU173_manual_assembly"] <- "PX127344 AndV"
mRNAmsa$name[mRNAmsa$name=="KEL153_NODE_1162_czid2019"] <- "PX127365"

unique(mRNAmsa$name)
names(mRNAmsa)[names(mRNAmsa)=="character"] <- "aa"

colz = c("A" = "tomato", "G"="goldenrod1", "T"="white")

#and control the order

#subset without partial sequence
#mRNAmsa = subset(mRNAmsa, name !="KEL153-virus")

mRNAmsa$name <- factor(mRNAmsa$name, levels = rev(c("PX127342 AndV","PX127343 AndV","PX127344 AndV","PX127365","AngV","SGV","NiV","HeV","MojV",        
                                                    "GhV","GAKV","DARV","LayV","YuV1","YuV2","CedV")))
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
        axis.text = element_text(size=14),
        title = element_text(size=16),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.text = element_text(size=10),
        #plot.background = element_rect(fill="transparent"),
        plot.margin = margin(0, 0, 0, 0, unit = "cm"))+
  #legend.position = "bottom", legend.direction = "horizontal")+
  #ggtitle("mRNA editing site")+
  #xlab("amino acid position\n(relative to AndV)")+ 
  scale_x_continuous(breaks = c(2,8), labels = c(1214, 1218))

mRNAmsa.plot

# ggsave(file = paste0("mRNAplot.png"),
#        plot = mRNAmsa.plot,
#        units="mm",
#        width=60,
#        height=70,
#        scale=2,
#        dpi=300)


#load the msa of c-terminal region of the v protein
Vmsa = tidy_msa(msa="v_alignment.fasta")
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

colz2 = c("A" = "white", "C"="yellow1","D"="white","E"="purple1","F"="white",
          "G"="green1","H"="orange1","I"="white","K"="dodgerblue","L"="white","M"="white",
          "N"="turquoise","P"="yellow","Q"="white","R"="turquoise1","S"="white","T"="white",
          "V"="white","W"="palevioletred1","Y"="white","-"="white")


#and control the order

Vmsa$name <- factor(Vmsa$name, levels = rev(c("PX127342 AndV_add1", "PX127343 AndV_add1", "PX127344 AndV_add1", "AngV_add1", "SGV_add1",  "NiV_add1",    "HeV_add1",    "MojV_add1",   "GhV_add1",   
                                              "GAKV_add1",   "DARV_add1",   "LayV_add1",   "YuV1_add1",   "YuV2_add1",   "CedV_add1" )))
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
        axis.text = element_text(size=14),
        title = element_text(size=16),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.text = element_text(size=10),
        plot.margin = margin(0, 0, 0, 2, unit = "cm"))+
  #legend.position = "bottom", legend.direction = "horizontal")+
  #ggtitle("V protein alignment")+
  scale_x_continuous(breaks = c(1,20,40,60), labels = c(409,429,449,469), 
                     sec.axis = dup_axis(breaks=c(4,23,27,39,41,44,48,53), labels=c("*","*","*","*","*","*","*","*")))

Vmsa.plot

# ggsave(file = paste0("Vplot.png"),
#        plot = Vmsa.plot,
#         units="mm",  
#         width=140, 
#         height=60, 
#         scale=2, 
#         dpi=300)


#load the msa of c-terminal region of the W protein
Wmsa = tidy_msa(msa="w_alignment.fasta")
head(Wmsa)

#assign new names to each
Wmsa$name <- as.character(Wmsa$name)
unique(Wmsa$name)
Wmsa$name[Wmsa$name=="ROU125_NODE_1_czid2023_add2"] <- "PX127342 AndV add2"
Wmsa$name[Wmsa$name=="ROU148_manual_assembly_add2"] <- "PX127343 AndV add2"
Wmsa$name[Wmsa$name=="ROU173_manual_assembly_add2"] <- "PX127344 AndV add2"

unique(Wmsa$name)
names(Wmsa)[names(Wmsa)=="character"] <- "aa"

# colz = scales::hue_pal()(length(unique(Wmsa$aa)))
# names(colz) <- unique(Wmsa$aa)
# names(colz)[names(colz)=="-"] <- "white"

colz = c("A" = "paleturquoise1", "C"="palevioletred1","D"="plum1","E"="purple1","F"="orange",
         "G"="goldenrod1","H"="orchid1","I"="slateblue1","K"="sienna3","L"="royalblue2","M"="grey",
         "N"="yellow1","P"="khaki1","Q"="tomato1","R"="turquoise1","S"="thistle1","T"="chocolate",
         "V"="deeppink4","W"="firebrick3","Y"="lightsteelblue1","-"="white")

colz2 = c("A" = "white", "C"="yellow1","D"="white","E"="purple1","F"="white",
          "G"="green1","H"="orange1","I"="white","K"="dodgerblue","L"="white","M"="white",
          "N"="turquoise","P"="yellow","Q"="white","R"="turquoise1","S"="white","T"="white",
          "V"="white","W"="palevioletred1","Y"="white","-"="white")

#and control the order

Wmsa$name <- factor(Wmsa$name,levels = rev(c("PX127342 AndV add2", "PX127343 AndV add2", "PX127344 AndV add2", "AngV_add2" , "SGV_add2",  "NiV_add2",    "HeV_add2" ,   "MojV_add2",   "GAKV_add2" , 
                                             "DARV_add2" ,  "LayV_add2",   "YuV1_add2" , "YuV2_add2" ,  "GhV_add2" )))
Wmsa.plot <- ggplot(data=Wmsa) + 
  geom_tile(aes(x=position, y=name, fill=aa, color=aa),width=0.9, height=0.9, show.legend = F) +
  geom_text(aes(x=position, y=name, label=aa), show.legend = F)+
  theme_bw() + 
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz) +
  theme(axis.title.y = element_blank(), 
        #panel.background = element_rect(fill="black"),
        #panel.border = element_line(color="black", size=1),
        panel.grid = element_blank(),
        axis.text = element_text(size=14),
        title = element_text(size=16),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.text = element_text(size=10),
        plot.margin = margin(0, 0, 0,0, unit = "cm"))+
#legend.position = "bottom", legend.direction = "horizontal")+
  #ggtitle("W protein alignment")+
  scale_x_continuous(breaks = c(1,18,40,60), labels = c(332,350,372,392))

Wmsa.plot

# ggsave(file = paste0("Wplot.png"),
#        plot = Wmsa.plot,
#        units="mm",  
#        width=170, 
#        height=70, 
#        scale=2, 
#        dpi=300)


#load the msa of ephrin binding conserved sites
ephrinmsa = tidy_msa(msa="ephrin_region_alignment.fasta")
head(ephrinmsa)

#assign new names to each
ephrinmsa$name <- as.character(ephrinmsa$name)
unique(ephrinmsa$name)
ephrinmsa$name[ephrinmsa$name=="CAT017_Plate_42_A10_NODE_34_czid2023"] <- "PX127396"
ephrinmsa$name[ephrinmsa$name=="KEL86spleen_Plate_22_A11_NODE_439_czid2023"] <- "PX127395 GhV-like"
ephrinmsa$name[ephrinmsa$name=="KEL234_NODE_79_czid2019"] <- "PX127397"
ephrinmsa$name[ephrinmsa$name=="KEL281B_NODE_113_czid2019"] <- "PX127399"
ephrinmsa$name[ephrinmsa$name=="KEL380_Plate_42_B8_NODE_44_czid2023"] <- "PX127400"
ephrinmsa$name[ephrinmsa$name=="ROU125_NODE_1_czid2023"] <- "PX127342 AndV"
ephrinmsa$name[ephrinmsa$name=="ROU148_manual_assembly"] <- "PX127343 AndV"
ephrinmsa$name[ephrinmsa$name=="ROU143_NODE_287_czid2023"] <- "PX127401"

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

ephrinmsa$name <- factor(ephrinmsa$name, levels = rev(c("PX127342 AndV",       "PX127343 AndV",   "PX127395 GhV-like",    "PX127396",   "PX127397",      "PX127399",     "PX127400",     
                                                        "PX127401",   "AngV",    "SGV",          "NiV",               "HeV",              
                                                        "MojV",              "GhV",               "GAKV",              "DARV",              "LayV",             
                                                        "YuV1",              "YuV2", "CedV")))
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
        axis.text = element_text(size=14),
        title = element_text(size=16),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.text = element_text(size=10),
        plot.margin = margin(0, 0, 0, 0, unit = "cm"))+
#legend.position = "bottom", legend.direction = "horizontal")+
  #ggtitle("Ephrin binding residues")+
  scale_x_continuous(breaks = c(1,10,17,23), labels = c(502, 529,556,580), 
                     sec.axis = dup_axis(breaks=c(3,4,11,12,13,14,18,24), labels=c("*","*","*","*","*","*","*","*")))

ephrinmsa.plot

# ggsave(file = paste0("ephrinplot.png"),
#        plot = ephrinmsa.plot,
#        units="mm",
#        width=170,
#        height=70,
#        scale=2,
#        dpi=300)

#Now do the phylogenies

##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_circulating_clades/phylo_trees/"

#fusion files first
setwd(paste0(homewd,"/pmx_gene_trees/glycoprotein"))

#load the trees and root them
tree1 <-  read.tree("glycoprotein_align.fasta.treefile")
rooted.tree1 <- root(tree1, which(tree1$tip.label == "YP_009094050.1"))

#Remove root from displaying, still calculates changes correctly without it
rooted.tree1<-drop.tip(rooted.tree1, "YP_009094050.1")

#load tree data prepared from elsewhere
dat1 <- read.csv(("G_metadata_all.csv"), header = T, stringsAsFactors = F)

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
  geom_cladelabel(node = 241, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 346, label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 331, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 338, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 344, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 302, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 303, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 309, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 360, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 321, label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 100, label = "Ferlavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 161, label = "Synodonvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 159, label = "Cynoglossusvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 158, label = "Hoplichthysvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 160, label = "Scoliodonvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  #geom_cladelabel(node = 162, label = "Salemvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 363, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 280, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  xlim(c(0,14))

p1

p1.2<-collapse(p1, 241)+geom_point2(aes(subset=(node==241)), size=4, shape=22, fill="firebrick1")
p1.3<-collapse(p1.2, 346)+geom_point2(aes(subset=(node==346)), size=4, shape=22, fill="firebrick4")
p1.4<-collapse(p1.3, 331)+geom_point2(aes(subset=(node==331)), size=4, shape=22, fill="darkgoldenrod1")
p1.5<-collapse(p1.4, 338)+geom_point2(aes(subset=(node==338)), size=4, shape=22, fill="darkorange")
p1.6<-collapse(p1.5, 344)+geom_point2(aes(subset=(node==344)), size=4, shape=22, fill="darkorange3")
p1.7<-collapse(p1.6, 302)+geom_point2(aes(subset=(node==302)), size=4, shape=22, fill="lightpink2")
p1.8<-collapse(p1.7, 303)+geom_point2(aes(subset=(node==303)), size=4, shape=22, fill="deeppink1")
p1.9<-collapse(p1.8, 309)+geom_point2(aes(subset=(node==309)), size=4, shape=22, fill="coral")
p1.10<-collapse(p1.9, 360)+geom_point2(aes(subset=(node==360)), size=4, shape=22, fill="lightsalmon3")
p1.11<-collapse(p1.10, 321)+geom_point2(aes(subset=(node==321)), size=4, shape=22, fill="indianred3")
#p1.12<-collapse(p1.10, 161)+geom_point2(aes(subset=(node==161)), size=4, shape=22, fill="pink1")
#p1.13<-collapse(p1.10, 159)+geom_point2(aes(subset=(node==159)), size=4, shape=22, fill="grey50")
#p1.14<-collapse(p1.10, 158)+geom_point2(aes(subset=(node==158)), size=4, shape=22, fill="lightsalmon1")
#p1.15<-collapse(p1.10, 160)+geom_point2(aes(subset=(node==160)), size=4, shape=22, fill="lightpink3")
#p1.16<-collapse(p1.10, 162)+geom_point2(aes(subset=(node==162)), size=4, shape=22, fill="indianred1")
p1.17<-collapse(p1.11, 363)+geom_point2(aes(subset=(node==363)), size=4, shape=22, fill="coral3")
p1.1<-collapse(p1.17, 280)+geom_point2(aes(subset=(node==280)), size=4, shape=22, fill="deeppink4")
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
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = -0.1))
tree1


#Put figures together

library(patchwork)

Fig3.1<-mRNAmsa.plot+Vmsa.plot + plot_layout(widths = c(1, 6))
Fig3.1

Fig3<-Fig3.1/Wmsa.plot+ plot_annotation(tag_levels = c('A', 'B','C')) &
  theme(plot.tag = element_text(size = 23, face="bold"))
Fig3


#put figure together: 
Fig4<-plot_grid(tree1, ephrinmsa.plot, labels=c("A","B"),
                rel_widths=c(1.5,1), rel_heights = c(1,1),
                ncol=2,align="hv", axis="l", label_size = 23)
Fig4


#save fig
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_circulating_clades/final_figures/"
ggsave(file = paste0(homewd, "Fig3_phosphoprotein.pdf"),
       plot= Fig3,
       units="mm",
       width=300,
       height=160,
       scale=2,
       dpi=500)

ggsave(file = paste0(homewd, "Fig4_glycoprotein.pdf"),
       plot= Fig4,
       units="mm",
       width=270,
       height=100,
       scale=2,
       dpi=500)


