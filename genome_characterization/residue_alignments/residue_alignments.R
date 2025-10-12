rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(ggtree)
library(seqinr)
library(ggmsa)

#homewd = "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_circulating_clades/genome_characterization/"
homewd = "/Users/gwenddolenkettenburg/Desktop/brook_lab/henipavirus_shedding_cross_correlation_project/genome_characterization/"

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
        plot.margin = unit(c(.1, 2,.1, .1), "cm"))+
        #legend.position = "bottom", legend.direction = "horizontal")+
        ggtitle("mRNA editing site")+
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
Vmsa$name[Vmsa$name=="ROU125_NODE_1_czid2023_add1"] <- "PX127342 AndV add1"
Vmsa$name[Vmsa$name=="ROU148_manual_assembly_add1"] <- "PX127343 AndV add1"
Vmsa$name[Vmsa$name=="ROU173_manual_assembly_add1"] <- "PX127344 AndV add1"

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

Vmsa$name <- factor(Vmsa$name, levels = rev(c("PX127342 AndV add1", "PX127343 AndV add1", "PX127344 AndV add1", "AngV_add1", "SGV_add1",  "NiV_add1",    "HeV_add1",    "MojV_add1",   "GhV_add1",   
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
        plot.margin = unit(c(.1, 2,.1, .1), "cm"))+
  #legend.position = "bottom", legend.direction = "horizontal")+
  ggtitle("V protein alignment")+
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
        plot.margin = unit(c(.1, 2,.1, .1), "cm"))+
  #legend.position = "bottom", legend.direction = "horizontal")+
  ggtitle("W protein alignment")+
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
        plot.margin = unit(c(.1, 2,.1, .1), "cm"))+
  #legend.position = "bottom", legend.direction = "horizontal")+
  ggtitle("Ephrin binding residues")+
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



#Put the whole figure together 

library(patchwork)

Fig2<-(mRNAmsa.plot | ephrinmsa.plot)/Vmsa.plot
Fig2

Fig2<- Fig2+ plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 23))


ggsave(file = paste0(homewd, "residue_alignments/mRNA_V_ephrin.pdf"),
       plot= Fig2,
       units="mm",  
       width=200, 
       height=180, 
       scale=2, 
       dpi=500)


