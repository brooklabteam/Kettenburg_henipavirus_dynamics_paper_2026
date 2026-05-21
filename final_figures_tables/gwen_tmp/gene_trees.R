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

#full genome tree first
setwd(paste0(homewd,"/full_genome_tree"))
tree1 <-  read.tree("fullgenome_align.fasta.treefile")
#drop root for aesthetics
rooted.tree1 <- root(tree1, which(tree1$tip.label == "NC_025345.1"))
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
  # geom_cladelabel(node=357, label='HeV/NiV sister clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=425, label='AndV clade', offset=8, fontsize=4, color="black", align=TRUE)+
  geom_cladelabel(node=213, label='GhV-like', offset=8, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=113, label='AndV-2 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=350, label='AndV-3 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=250, label='AngV clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=247, label='AngV-2 clade', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=31, label='AchPV2-like virus', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # geom_cladelabel(node=268, label='ThkPV1-like viruses', offset=8.55, fontsize=4, color="black", align=TRUE)+
  # 
  geom_cladelabel(node = 275, label = "Hendra/Nipah henipavirus clade",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 362, label = "Pararubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 351, label = "Metaavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 341, label = "Orthoavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 338, label = "Paraavulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 333, label = "Aquaparamyxovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 324, label = "Respirovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 394, label = "Jeilongvirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 391, label = "Narmovirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 372, label = "Orthorubulavirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 385, label = "Morbillivirus (collapsed)",offset=0.1, fontsize = 4, color="black") +
  geom_cladelabel(node = 404, label = "Parahenipavirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  geom_cladelabel(node = 380, label = "Scoliodonvirus/Synodonvirus (collapsed)",offset=0.1, fontsize = 4, color="black")+
  xlim(c(0,12))

p1

p1.2<-ggtree::collapse(p1, 275)+geom_point2(aes(subset=(node==275)), size=4, shape=22, fill="firebrick1")
p1.3<-ggtree::collapse(p1.2, 362)+geom_point2(aes(subset=(node==362)), size=4, shape=22, fill="firebrick4")
p1.4<-ggtree::collapse(p1.3, 351)+geom_point2(aes(subset=(node==351)), size=4, shape=22, fill="darkgoldenrod1")
p1.5<-ggtree::collapse(p1.4, 341)+geom_point2(aes(subset=(node==341)), size=4, shape=22, fill="darkorange")
p1.6<-ggtree::collapse(p1.5, 338)+geom_point2(aes(subset=(node==338)), size=4, shape=22, fill="darkorange3")
p1.7<-ggtree::collapse(p1.6, 333)+geom_point2(aes(subset=(node==333)), size=4, shape=22, fill="lightpink2")
p1.8<-ggtree::collapse(p1.7, 324)+geom_point2(aes(subset=(node==324)), size=4, shape=22, fill="deeppink1")
p1.9<-ggtree::collapse(p1.8, 394)+geom_point2(aes(subset=(node==394)), size=4, shape=22, fill="coral")
p1.10<-ggtree::collapse(p1.9, 391)+geom_point2(aes(subset=(node==391)), size=4, shape=22, fill="lightsalmon3")
p1.11<-ggtree::collapse(p1.10, 372)+geom_point2(aes(subset=(node==372)), size=4, shape=22, fill="indianred3")
p1.17<-ggtree::collapse(p1.11, 385)+geom_point2(aes(subset=(node==385)), size=4, shape=22, fill="coral3")
p1.18<-ggtree::collapse(p1.17, 404)+geom_point2(aes(subset=(node==404)), size=4, shape=22, fill="deeppink4")
p1.1<-ggtree::collapse(p1.18, 380)+geom_point2(aes(subset=(node==380)), size=4, shape=22, fill="lightpink3")
p1.1

p1.1<- p1.1 +  
  #geom_hilight(node=357, fill="deeppink3", alpha=0.3, extendto=12.8)+ #HeV/NiV sister
  geom_hilight(node=425, fill="darkcyan", alpha=0.3, extendto=12.8)+ #AndV
  # geom_hilight(node=113, fill="darkturquoise", alpha=0.3, extendto=12.8)+ #AndV2
  # geom_hilight(node=350, fill="dodgerblue", alpha=0.3, extendto=12.8)+ #AndV3
  # geom_hilight(node=250, fill="springgreen3", alpha=0.3, extendto=12.8)+ #AngV
  # geom_hilight(node=247, fill="chartreuse", alpha=0.3, extendto=12.8)+ #AngV2
  geom_hilight(node=213, fill="orange", alpha=0.3, extendto=12.8) #GhVlike
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




##Now add the seasonal shedding data and GAMs

homewd = "/Users/gwenddolenkettenburg/Desktop/Developer/Kettenburg_henipavirus_dynamics_paper_2026"
homewd="/Users/carabrook/Developer/Kettenburg_henipavirus_dynamics_paper_2026/"
batwd = "/Users/gwenddolenkettenburg/Desktop/Developer/Madagascar-Bat-Data"
batwd= "/Users/carabrook/Developer/Madagascar-Bat-Data/"

# First, load the shedding data 
shed.dat <- read.csv(file = paste0(homewd, "/virus_shedding/data/henipavirus_shedding_data.csv"))
head(shed.dat)
names(shed.dat)

shed.dat <- dplyr::select(shed.dat, sampleid, GhV_shedding_ur, AngV_shedding_ur)

names(shed.dat)[names(shed.dat) == "GhV_shedding_ur"] <- "GhV_pos"
names(shed.dat)[names(shed.dat) == "AngV_shedding_ur"] <- "AngV_pos"

sum(shed.dat$AngV_pos)#45
sum(shed.dat$GhV_pos)#127

# link this to the catching data to plot by sampling session
catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)#Historic-Catching-Data-Through2020/catching_data_historic.csv"))
#catch.dat <- read.csv(file = paste0(batwd, "/Catching-Data/Historic-Catching-Data-Through2020/catching_data_historic.csv"))

head(catch.dat)
names(catch.dat)

# drop most columns and name those you want
catch.dat <-dplyr::select(catch.dat, sampling_session, roost_site, bat_species, processing_date, sampleid, barcode_id, bat_sex, bat_age_class, young_of_year, pupmom_pair, sampleid_pupmom_pair, bat_weight_g, bat_forearm_mm, mass_forearm_residual, recap_new)

#and merge all three datasets
setdiff(shed.dat$sampleid, catch.dat$sampleid)

merge.dat <- merge(shed.dat, catch.dat, by="sampleid")
head(merge.dat)

merge.dat$processing_date <- as.Date(merge.dat$processing_date)

#change roost site to just Moramanga region or Ankarana
unique(merge.dat$roost_site)
merge.dat$roost_site[merge.dat$roost_site=="Ambakoana_Pruf" | merge.dat$roost_site=="Marotsipohy_Pruf" | merge.dat$roost_site=="Analambotaka_Pruf" | merge.dat$roost_site=="Morarano_Pruf" | merge.dat$roost_site=="Marovitsika_Pruf" |  merge.dat$roost_site=="Mahialambo_Pruf"  | merge.dat$roost_site=="Lakato_Edup" | merge.dat$roost_site=="Angavobe_Edup" |merge.dat$roost_site=="Angavokely_Edup"| merge.dat$roost_site=="Maromizaha_Rmad"] <- "Moramanga"
merge.dat$roost_site[merge.dat$roost_site=="Mahabo_Edup" | merge.dat$roost_site=="Mahabo_Pruf"] <- "Mahabo"
merge.dat$roost_site[merge.dat$roost_site=="Makira_Rmad" | merge.dat$roost_site=="Makira_Pruf"] <- "Makira"
merge.dat$roost_site[merge.dat$roost_site!="Moramanga" & merge.dat$roost_site!="Makira" & merge.dat$roost_site!="Mahabo"] <- "Ankarana"

unique(merge.dat$roost_site) #mahabo and makira detract from seasonal analysis, so take them out
single.season.dat = subset(merge.dat, roost_site=="Mahabo" | roost_site == "Makira")

merge.dat = subset(merge.dat, roost_site!="Mahabo" & roost_site != "Makira" )


merge.dat$processing_date <- as.Date(merge.dat$processing_date)
merge.dat$doy <- yday(merge.dat$processing_date)

head(merge.dat)

#and sampling session
#summarise more broadly
sampling.date <- ddply(merge.dat, .(sampling_session), summarise, min_date = min(processing_date), max_date = max(processing_date))
sampling.date$mid_date <- (sampling.date$max_date-sampling.date$min_date)/2 + sampling.date$min_date

# #and slim to merge
sampling.merge <- dplyr::select(sampling.date, sampling_session, mid_date)
# 
merge.dat <- merge(merge.dat, sampling.merge, by = c( "sampling_session"))
merge.dat$year <- year(merge.dat$processing_date)
merge.dat$month <- as.character(month(merge.dat$processing_date))
merge.dat$month[merge.dat$month!="10"& merge.dat$month!="11" &merge.dat$month!="12"]<-paste0("0", merge.dat$month[merge.dat$month!="10"& merge.dat$month!="11" &merge.dat$month!="12"])
merge.dat$mid_date <- as.Date(paste0(merge.dat$year, "-", merge.dat$month, "-15"))
head(merge.dat)

sampling.date <- ddply(merge.dat, .(bat_species, roost_site, mid_date), summarise)
sampling.date# examine
#

tot.merge = subset(merge.dat, year<2022)

head(tot.merge)
tot.merge$processing_date <- as.Date(tot.merge$processing_date)
tot.merge$month <- month(tot.merge$processing_date)
tot.merge$year <- year(tot.merge$processing_date)
tot.merge = subset(tot.merge, bat_species=="Eidolon dupreanum" | bat_species =="Rousettus madagascariensis" )
#tot.merge = subset(tot.merge, bat_species=="Eidolon dupreanum" | bat_species =="Rousettus madagascariensis" | bat_species =="Pteropus rufus")

# look at staggered prevalence
prev.dat <- ddply(subset(tot.merge, !is.na(GhV_pos)),.(bat_species, roost_site, mid_date), summarise, GhVpos = sum(GhV_pos, na.rm = T), AngVpos=sum(AngV_pos, na.rm=T), Ntot = length(sampleid))
prev.dat.sex <- ddply(subset(tot.merge, !is.na(GhV_pos)),.(bat_species, bat_sex, roost_site, mid_date), summarise, GhVpos = sum(GhV_pos, na.rm = T), AngVpos=sum(AngV_pos, na.rm=T), Ntot = length(sampleid))
prev.dat$prevGhV <- prev.dat$GhVpos/prev.dat$Ntot
prev.dat$prevAngV <- prev.dat$AngVpos/prev.dat$Ntot

prev.dat.sex$prevGhV <- prev.dat.sex$GhVpos/prev.dat.sex$Ntot
prev.dat.sex$prevAngV <- prev.dat.sex$AngVpos/prev.dat.sex$Ntot

# prev.dat$prevGhV_lci <- binom.exact(x=prev.dat$GhVpos, n=prev.dat$Ntot)$lower
# prev.dat$prevGhV_uci <- binom.exact(prev.dat$GhVpos, prev.dat$Ntot)$upper
# prev.dat$prevAngV_lci <- binom.exact(prev.dat$AngVpos, prev.dat$Ntot)$lower
# prev.dat$prevAngV_uci <- binom.exact(prev.dat$AngVpos, prev.dat$Ntot)$upper

#and plot prev with CT
p1 <- ggplot(prev.dat) + 
  ylab("prevalence GhV") + theme_bw() + facet_grid(~bat_species) +
  geom_point(aes(x=mid_date, y=prevGhV, size=Ntot, color=roost_site), pch=21);p1

#include era block to correct for the gaps in the time series
prev.dat$era <- "II"
prev.dat$era[prev.dat$mid_date<"2017-01-01"] <- "I"
prev.dat$era <- as.factor(prev.dat$era)

#and by month
prev.dat$month <- month(prev.dat$mid_date)
prev.dat$year <- year(prev.dat$mid_date)
prev.dat <- prev.dat[complete.cases(prev.dat),]
prev.dat$doy <- yday(prev.dat$mid_date)
prev.dat$year <- as.factor(prev.dat$year)

#and by month
prev.dat.sex$month <- month(prev.dat.sex$mid_date)
prev.dat.sex$year <- year(prev.dat.sex$mid_date)
prev.dat.sex <- prev.dat.sex[complete.cases(prev.dat.sex),]
prev.dat.sex$doy <- yday(prev.dat.sex$mid_date)
prev.dat.sex$year <- as.factor(prev.dat.sex$year)

prev.dat.sex$era <- "II"
prev.dat.sex$era[prev.dat.sex$mid_date<"2017-01-01"] <- "I"
prev.dat.sex$era <- as.factor(prev.dat.sex$era)


p2 <- ggplot(prev.dat) + facet_grid(~bat_species) +
  ylab("prevalence GhV") + theme_bw() + #scale_shape_manual(values = c(21,24)) +
  geom_point(aes(x=month, y=prevGhV, size=Ntot, shape=roost_site, color=as.character(year)));p2

# and gam to identify seasonal signatures
prev.dat$GhVneg <- prev.dat$Ntot-prev.dat$GhVpos
prev.dat.sex$GhVneg <- prev.dat.sex$Ntot-prev.dat.sex$GhVpos
prev.dat$roost_site <- as.factor(prev.dat$roost_site)

gamGhVnoSex <- gam(cbind(GhVpos, GhVneg) ~ roost_site +
                     s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                     s(as.numeric(mid_date), bs="tp", k=6, by=interaction(roost_site, era)),
                   family = "binomial",
                   method = "REML", 
                   gamma = 1,
                   select = TRUE,
                   data= subset(prev.dat, bat_species =="Eidolon dupreanum"))
summary(gamGhVnoSex) #moramanga interannual in both eras, mahabo interannual in era I
#plot(gamGhVnoSex)

gamGhVSex <- gam(cbind(GhVpos, GhVneg) ~ roost_site + bat_sex +
                   s(as.numeric(doy), bs="cc", k=10, by = interaction(roost_site, bat_sex)) +
                   s(as.numeric(mid_date), bs="tp", k=6, by = interaction(roost_site, bat_sex, era)),
                 family = "binomial",
                 method = "REML", 
                 gamma = 1,
                 select = TRUE,
                 data= subset(prev.dat.sex, bat_species =="Eidolon dupreanum"))
summary(gamGhVSex) #moramanga female interannual both eras, moramanga males almost sig in era I

AIC(gamGhVnoSex, gamGhVSex) #better with no sex

#and AngV
prev.dat$AngVneg <- prev.dat$Ntot-prev.dat$AngVpos
prev.dat.sex$AngVneg <- prev.dat.sex$Ntot-prev.dat.sex$AngVpos
prev.dat$roost_site <- as.factor(prev.dat$roost_site)

gamAngVnoSex <- gam(cbind(AngVpos, AngVneg) ~ roost_site +
                      s(as.numeric(doy), bs="cc", k=10, by=interaction(roost_site, era)) +
                      s(as.numeric(mid_date), bs="tp", k=6, by=interaction(roost_site, era)),
                    family = "binomial",
                    method = "REML", 
                    gamma = 1,
                    select = TRUE,
                    data= subset(prev.dat, bat_species =="Eidolon dupreanum"))
summary(gamAngVnoSex) #moramanga interannual era II, Mahabo I intraannual almost sig
#plot(gamAngVnoSex)

gamAngVSex <- gam(cbind(AngVpos, AngVneg) ~ roost_site + bat_sex +
                    s(as.numeric(doy), bs="cc", k=10, by = interaction(roost_site, bat_sex, era)) +
                    s(as.numeric(mid_date), bs="tp", k=6, by = interaction(roost_site, bat_sex, era)),
                  family = "binomial",
                  method = "REML", 
                  gamma = 1,
                  select = TRUE,
                  data= subset(prev.dat.sex, bat_species =="Eidolon dupreanum"))
summary(gamAngVSex) #Moramanga female II interannual...I almost sig, Moramanga female intraannual I
AIC(gamAngVnoSex, gamAngVSex) #better with no sex

pred.df <- cbind.data.frame(mid_date = seq(from= as.Date("2015-01-01"), to = as.Date("2015-12-31"), by=1), doy=1:365)
pred.df$bat_species <- "Eidolon dupreanum"
pred.df$antigen <- "GhV"
pred.df2 <- pred.df
pred.df2$antigen <- "AngV"
pred.df <- rbind(pred.df, pred.df2)
pred.df$roost_site <- "Moramanga"

pred.df$era <- "I"
pred.df$era[pred.df$mid_date < "2017-07-01"] <- "II"
pred.df$era <- as.factor(pred.df$era)

pred.df$prevalence <- NA
pred.df$prev_lci <- NA
pred.df$prev_uci <- NA

pred.df$prevalence[pred.df$antigen=="GhV"] <- predict(gamGhVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(mid_date)"))
pred.df$prev_lci[pred.df$antigen=="GhV"] <- predict(gamGhVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(roost_site)", "s(mid_date)")) - 1.96*predict(gamGhVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(roost_site)", "s(mid_date)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="GhV"] <- predict(gamGhVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(roost_site)", "s(mid_date)")) + 1.96*predict(gamGhVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="GhV",], exclude = c("s(roost_site)", "s(mid_date)"), se.fit = TRUE)$se

pred.df$prevalence[pred.df$antigen=="AngV"] <- predict(gamAngVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(roost_site)", "s(mid_date)"))
pred.df$prev_lci[pred.df$antigen=="AngV"] <- predict(gamAngVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(roost_site)", "s(mid_date)")) - 1.96*predict(gamAngVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(roost_site)", "s(mid_date)"), se.fit = TRUE)$se
pred.df$prev_uci[pred.df$antigen=="AngV"] <- predict(gamAngVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(roost_site)", "s(mid_date)")) + 1.96*predict(gamAngVnoSex, type = "response", newdata = pred.df[pred.df$antigen=="AngV",], exclude = c("s(roost_site)", "s(mid_date)"), se.fit = TRUE)$se

pred.df$prev_lci[pred.df$prev_lci<0] <- 0
pred.df$prev_uci[pred.df$prev_uci>1] <- 1

gest.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum",  "Eidolon dupreanum",  "Eidolon dupreanum",  "Eidolon dupreanum", "Rousettus madagascariensis", "Rousettus madagascariensis",  "Rousettus madagascariensis"), start = c(as.Date("2013-08-03"), as.Date("2014-08-03"), as.Date("2015-08-03"), as.Date("2018-08-03"), as.Date("2019-08-03"),  as.Date("2014-09-11"),as.Date("2018-09-11"), as.Date("2019-09-11")), stop=c(as.Date("2013-11-15"), as.Date("2014-11-15"), as.Date("2015-11-15"), as.Date("2018-11-15"), as.Date("2019-11-15"),  as.Date("2014-12-11"), as.Date("2018-12-11"), as.Date("2019-12-11")))
gest.dat$season = "gestation"
lact.dat <- cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis","Rousettus madagascariensis","Rousettus madagascariensis", "Rousettus madagascariensis"), start = c(as.Date("2013-11-16"), as.Date("2014-11-16"), as.Date("2015-11-16"), as.Date("2018-01-01"), as.Date("2018-11-16"), as.Date("2019-11-16"), as.Date("2014-12-12"),  as.Date("2018-01-01"), as.Date("2018-12-12"), as.Date("2019-12-12")), stop=c(as.Date("2014-02-02"), as.Date("2015-02-02"), as.Date("2016-02-02"), as.Date("2018-02-02"), as.Date("2019-02-02"), as.Date("2020-02-02"), as.Date("2015-02-19"), as.Date("2018-02-19"), as.Date("2019-02-19"), as.Date("2020-02-19")))
lact.dat$season = "lactation"

#and winter
cal.dat <-  cbind.data.frame(bat_species=c("Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Eidolon dupreanum", "Rousettus madagascariensis", "Rousettus madagascariensis", "Rousettus madagascariensis"), start = c(as.Date("2013-07-01"), as.Date("2014-05-15"),as.Date("2015-05-15"), as.Date("2018-05-15"), as.Date("2019-05-15"),  as.Date("2014-05-15"), as.Date("2018-05-15"), as.Date("2019-05-15")), stop=c(as.Date("2013-10-10"),as.Date("2014-10-10"),as.Date("2015-10-10"), as.Date("2018-10-10"), as.Date("2019-10-10"),as.Date("2014-10-10"), as.Date("2018-10-10"), as.Date("2019-10-10")))
cal.dat$season <- "winter"
cal.dat <- rbind(cal.dat, gest.dat, lact.dat)
cal.dat$year <- year(cal.dat$stop)
cal.dat$era <- "II"
cal.dat$era[cal.dat$year<2017] <- "I"
cal.dat.save=cal.dat

cal.sub = subset(cal.dat, bat_species =="Eidolon dupreanum" & year==2019)
cal.sub$doy_start <- yday(cal.sub$start)
cal.sub$doy_stop <- yday(cal.sub$stop)
cal.sub$year <- year(cal.sub$stop)
cal.sub$doy_stop[cal.sub$season=="lactation"] <- 365
cal.sub <- rbind(cal.sub, cal.sub[nrow(cal.sub),])
cal.sub$doy_start[4] <- 1
cal.sub$doy_stop[4] <- yday("2019-02-02")

prev.dat$prevGhV_lci <- binom.exact(prev.dat$GhVpos, prev.dat$Ntot)$lower
prev.dat$prevGhV_uci <- binom.exact(prev.dat$GhVpos, prev.dat$Ntot)$upper
prev.dat$prevAngV_lci <- binom.exact(prev.dat$AngVpos, prev.dat$Ntot)$lower
prev.dat$prevAngV_uci <- binom.exact(prev.dat$AngVpos, prev.dat$Ntot)$upper
prev.dat$year <- as.factor(prev.dat$year)

colorz = c('AngV' = "#F8766D", 'HeV' ="#A58AFF", 
           'NiV' = "#DB72FB", 'GhV' = "#00C094", 
           'CedV'= "#53B400", "MojV" = "#619CFF")


filz <- c('gestation' = "yellow", 'lactation' = "greenyellow", 'winter' = "skyblue")

yearz <- c('2013'='3', '2014' = '7', '2015'='8', '2016'='18', '2018' = '16', '2019' = '17', '2020' = '15')
yearz <- c('2013'=3, '2014' = 7, '2015'=8, '2016'=18, '2018' = 16, '2019' = 17, '2020' = 15)


###rousettus

gamGhVnoSexRou <- gam(cbind(GhVpos, GhVneg) ~ roost_site +
                        s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                        s(as.numeric(mid_date), bs="tp", k=6, by=roost_site),
                      family = "binomial",
                      method = "REML", 
                      gamma = 1,
                      select = TRUE,
                      data= subset(prev.dat, bat_species =="Rousettus madagascariensis"))
summary(gamGhVnoSexRou) #moramanga interannual
#plot(gamGhVnoSex)

gamGhVSexRou <- gam(cbind(GhVpos, GhVneg) ~ roost_site + bat_sex +
                      s(as.numeric(doy), bs="cc", k=10, by = interaction(roost_site, bat_sex)) +
                      s(as.numeric(mid_date), bs="tp", k=6, by = interaction(roost_site, bat_sex)),
                    family = "binomial",
                    method = "REML", 
                    gamma = 1,
                    select = TRUE,
                    data= subset(prev.dat.sex, bat_species =="Rousettus madagascariensis"))
summary(gamGhVSexRou) #moramanga male interannual and female borderline

AIC(gamGhVnoSexRou, gamGhVSexRou) #better with no sex

#and AngV

gamAngVnoSexRou <- gam(cbind(AngVpos, AngVneg) ~ roost_site +
                         s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                         s(as.numeric(mid_date), bs="tp", k=6, by=roost_site),
                       family = "binomial",
                       method = "REML", 
                       gamma = 1,
                       select = TRUE,
                       data= subset(prev.dat, bat_species =="Rousettus madagascariensis"))
summary(gamAngVnoSexRou) #nothing
#plot(gamAngVnoSex)

gamAngVSexRou <- gam(cbind(AngVpos, AngVneg) ~ roost_site + bat_sex +
                       s(as.numeric(doy), bs="cc", k=10, by = interaction(roost_site, bat_sex)) +
                       s(as.numeric(mid_date), bs="tp", k=6, by = interaction(roost_site, bat_sex)),
                     family = "binomial",
                     method = "REML", 
                     gamma = 1,
                     select = TRUE,
                     data= subset(prev.dat.sex, bat_species =="Rousettus madagascariensis"))
summary(gamAngVSexRou) #nothing

AIC(gamAngVnoSexRou, gamAngVSexRou) #better with no sex

pred.df.Rou <- cbind.data.frame(mid_date = seq(from= as.Date("2015-01-01"), to = as.Date("2015-12-31"), by=1), doy=1:365)
pred.df.Rou$bat_species <- "Rousettus madagascariensis"
pred.df.Rou$antigen <- "GhV"
pred.df.Rou2 <- pred.df.Rou
pred.df.Rou2$antigen <- "AngV"
pred.df.Rou <- rbind(pred.df.Rou, pred.df.Rou2)
pred.df.Rou$roost_site <- "Moramanga"

pred.df.Rou$era <- "I"
pred.df.Rou$era[pred.df.Rou$mid_date < "2017-07-01"] <- "II"
pred.df.Rou$era <- as.factor(pred.df.Rou$era)


pred.df.Rou$prevalence <- NA
pred.df.Rou$prev_lci <- NA
pred.df.Rou$prev_uci <- NA

pred.df.Rou$prevalence[pred.df.Rou$antigen=="GhV"] <- predict(gamGhVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="GhV",], exclude = c("s(roost_site)", "s(mid_date)"))
pred.df.Rou$prev_lci[pred.df.Rou$antigen=="GhV"] <- predict(gamGhVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="GhV",], exclude = c("s(roost_site)", "s(mid_date)")) - 1.96*predict(gamGhVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="GhV",], exclude = c("s(roost_site)", "s(mid_date)"), se.fit = TRUE)$se
pred.df.Rou$prev_uci[pred.df.Rou$antigen=="GhV"] <- predict(gamGhVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="GhV",], exclude = c("s(roost_site)", "s(mid_date)")) + 1.96*predict(gamGhVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="GhV",], exclude = c("s(roost_site)", "s(mid_date)"), se.fit = TRUE)$se

pred.df.Rou$prevalence[pred.df.Rou$antigen=="AngV"] <- predict(gamAngVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="AngV",], exclude = c("s(roost_site)", "s(mid_date)"))
pred.df.Rou$prev_lci[pred.df.Rou$antigen=="AngV"] <- predict(gamAngVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="AngV",], exclude = c("s(roost_site)", "s(mid_date)")) - 1.96*predict(gamAngVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="AngV",], exclude = c("s(roost_site)", "s(mid_date)"), se.fit = TRUE)$se
pred.df.Rou$prev_uci[pred.df.Rou$antigen=="AngV"] <- predict(gamAngVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="AngV",],exclude = c("s(roost_site)", "s(mid_date)")) + 1.96*predict(gamAngVnoSexRou, type = "response", newdata = pred.df.Rou[pred.df.Rou$antigen=="AngV",], exclude = c("s(roost_site)", "s(mid_date)"), se.fit = TRUE)$se


cal.sub.Rou = subset(cal.dat, bat_species =="Rousettus madagascariensis" & year==2019)
cal.sub.Rou$doy_start <- yday(cal.sub.Rou$start)
cal.sub.Rou$doy_stop <- yday(cal.sub.Rou$stop)
cal.sub.Rou$year <- year(cal.sub.Rou$stop)
cal.sub.Rou$doy_stop[cal.sub.Rou$season=="lactation"] <- 365
cal.sub.Rou <- rbind(cal.sub.Rou, cal.sub.Rou[nrow(cal.sub.Rou),])
cal.sub.Rou$doy_start[4] <- 1
cal.sub.Rou$doy_stop[4] <- yday("2019-02-19")

cal.sub <- rbind(cal.sub, cal.sub.Rou)


prev.dat$prevGhV_lci <- binom.exact(prev.dat$GhVpos, prev.dat$Ntot)$lower
prev.dat$prevGhV_uci <- binom.exact(prev.dat$GhVpos, prev.dat$Ntot)$upper
prev.dat$prevAngV_lci <- binom.exact(prev.dat$AngVpos, prev.dat$Ntot)$lower
prev.dat$prevAngV_uci <- binom.exact(prev.dat$AngVpos, prev.dat$Ntot)$upper
prev.dat$year <- as.factor(prev.dat$year)

colorz = c('AngV' = "#F8766D", 'HeV' ="#A58AFF", 
           'NiV' = "#DB72FB", 'GhV' = "#00C094", 
           'CedV'= "#53B400", "MojV" = "#619CFF")


filz <- c('gestation' = "yellow", 'lactation' = "greenyellow", 'winter' = "skyblue")

yearz <- c('2013'='3', '2014' = '7', '2015'='8', '2016'='18', '2018' = '16', '2019' = '17', '2020' = '15')
yearz <- c('2013'=3, '2014' = 7, '2015'=8, '2016'=18, '2018' = 16, '2019' = 17, '2020' = 15)

pred.df <- rbind(pred.df, pred.df.Rou)

#and plot the interannual for both

#just need new predictions
summary(gamGhVnoSex) # interannual moramanga, mahabo I
summary(gamAngVnoSex)# interannual moramanga II
summary(gamGhVnoSexRou) # interannual moramanga
summary(gamAngVnoSexRou) #nothing


#rousettus 
inter.pred <- cbind.data.frame(mid_date = seq(from =as.Date("2013-08-01"), to=as.Date("2020-03-31"), by=1))
inter.pred$doy <- yday(inter.pred$mid_date)
inter.pred$roost_site <- "Moramanga"
inter.pred$bat_species <- "Eidolon dupreanum"
inter.pred$antigen <- "GhV"
inter.pred2 <- inter.pred
inter.pred2$antigen <- "AngV"
inter.pred <- rbind(inter.pred, inter.pred2)
inter.pred2 <- inter.pred
inter.pred2$bat_species <- "Rousettus madagascariensis"
inter.pred <- rbind(inter.pred, inter.pred2)

inter.pred <- subset(inter.pred, mid_date <"2016-01-31" | mid_date>"2018-01-01")
inter.pred$keep = 1
inter.pred$keep[inter.pred$bat_species=="Rousettus madagascariensis" & inter.pred$mid_date<"2014-07-01"] <- 0
inter.pred$keep[inter.pred$bat_species=="Rousettus madagascariensis" & inter.pred$mid_date>"2015-05-15" &inter.pred$mid_date<"2018-01-01" ] <- 0
inter.pred = subset(inter.pred, keep ==1)


inter.pred$era <- "II"
inter.pred$era[inter.pred$mid_date<"2017-01-01"] <- "I"
inter.pred$era <- as.factor(inter.pred$era)

inter.pred$prevalence <- NA
inter.pred$prev_lci <- NA
inter.pred$prev_uci <- NA


inter.pred$prevalence[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Eidolon dupreanum"] <- predict(gamGhVnoSex, type="response", newdata = inter.pred[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Eidolon dupreanum",], exclude = c("s(roost_site)"))
inter.pred$prev_lci[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Eidolon dupreanum"] <- predict(gamGhVnoSex, type="response", newdata = inter.pred[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Eidolon dupreanum",], exclude = c("s(roost_site)")) - 1.96*predict(gamGhVnoSex, type="response", newdata = inter.pred[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Eidolon dupreanum",], exclude = c("s(roost_site)"), se.fit = T)$se
inter.pred$prev_uci[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Eidolon dupreanum"] <- predict(gamGhVnoSex, type="response", newdata = inter.pred[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Eidolon dupreanum",], exclude = c("s(roost_site)")) + 1.96*predict(gamGhVnoSex, type="response", newdata = inter.pred[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Eidolon dupreanum",], exclude = c("s(roost_site)"), se.fit = T)$se

inter.pred$prevalence[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Eidolon dupreanum"] <- predict(gamAngVnoSex, type="response", newdata = inter.pred[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Eidolon dupreanum",], exclude = c("s(roost_site)"))
inter.pred$prev_lci[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Eidolon dupreanum"] <- predict(gamAngVnoSex, type="response", newdata = inter.pred[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Eidolon dupreanum",], exclude = c("s(roost_site)")) - 1.96*predict(gamAngVnoSex, type="response", newdata = inter.pred[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Eidolon dupreanum",], exclude = c("s(roost_site)"), se.fit = T)$se
inter.pred$prev_uci[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Eidolon dupreanum"] <- predict(gamAngVnoSex, type="response", newdata = inter.pred[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Eidolon dupreanum",], exclude = c("s(roost_site)")) + 1.96*predict(gamAngVnoSex, type="response", newdata = inter.pred[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Eidolon dupreanum",], exclude = c("s(roost_site)"), se.fit = T)$se


inter.pred$prevalence[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Rousettus madagascariensis"] <- predict(gamGhVnoSexRou, type="response", newdata = inter.pred[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Rousettus madagascariensis",], exclude = c("s(roost_site)"))
inter.pred$prev_lci[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Rousettus madagascariensis"] <- predict(gamGhVnoSexRou, type="response", newdata = inter.pred[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Rousettus madagascariensis",], exclude = c("s(roost_site)")) - 1.96*predict(gamGhVnoSexRou, type="response", newdata = inter.pred[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Rousettus madagascariensis",], exclude = c("s(roost_site)"), se.fit = T)$se
inter.pred$prev_uci[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Rousettus madagascariensis"] <- predict(gamGhVnoSexRou, type="response", newdata = inter.pred[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Rousettus madagascariensis",], exclude = c("s(roost_site)")) + 1.96*predict(gamGhVnoSexRou, type="response", newdata = inter.pred[inter.pred$antigen=="GhV" & inter.pred$bat_species=="Rousettus madagascariensis",], exclude = c("s(roost_site)"), se.fit = T)$se

inter.pred$prevalence[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Rousettus madagascariensis"] <- predict(gamAngVnoSexRou, type="response", newdata = inter.pred[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Rousettus madagascariensis",], exclude = c("s(roost_site)"))
inter.pred$prev_lci[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Rousettus madagascariensis"] <- predict(gamAngVnoSexRou, type="response", newdata = inter.pred[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Rousettus madagascariensis",], exclude = c("s(roost_site)")) - 1.96*predict(gamAngVnoSexRou, type="response", newdata = inter.pred[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Rousettus madagascariensis",], exclude = c("s(roost_site)"), se.fit = T)$se
inter.pred$prev_uci[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Rousettus madagascariensis"] <- predict(gamAngVnoSexRou, type="response", newdata = inter.pred[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Rousettus madagascariensis",], exclude = c("s(roost_site)")) + 1.96*predict(gamAngVnoSexRou, type="response", newdata = inter.pred[inter.pred$antigen=="AngV" & inter.pred$bat_species=="Rousettus madagascariensis",], exclude = c("s(roost_site)"), se.fit = T)$se

inter.pred$prev_lci[inter.pred$prev_lci<0] <- 0
inter.pred$prev_uci[inter.pred$prev_uci>1] <- 1

inter.pred$year <- year(inter.pred$mid_date)
inter.pred$group <- paste0(inter.pred$bat_species, "-", inter.pred$antigen, "-", inter.pred$era, "-", inter.pred$year)

# pInterannualShed.All <- ggplot(data=prev.dat) +
#   geom_rect(data=cal.dat.save, aes(xmin=start, xmax=stop, ymin=0, ymax=1, fill=season), alpha=.3) +
#   geom_point(aes(x=as.Date(mid_date), y=prevGhV, size=Ntot, shape=year), color="#00C094") + 
#   geom_point(aes(x=as.Date(mid_date), y=prevAngV, size=Ntot, shape=year), color="#F8766D") +
#   guides(color = guide_legend(title = "virus"))+  
#   scale_fill_manual(values=filz) +
#   ggnewscale::new_scale_fill() +
#   geom_linerange(aes(x=as.Date(mid_date), ymin=prevGhV_lci, ymax=prevGhV_uci),color="#00C094", linetype=2, alpha=.3, show.legend = F) +
#   geom_linerange(aes(x=as.Date(mid_date), ymin=prevAngV_lci, ymax=prevAngV_uci),color="#F8766D", linetype=2, alpha=.3, show.legend = F) +
#   geom_ribbon(data=inter.pred, aes(x=as.Date(mid_date), ymin= prev_lci, ymax=prev_uci, fill=antigen), alpha=.3, show.legend = F) +
#   geom_line(data=inter.pred, aes(x=as.Date(mid_date), y= prevalence, color=antigen, group=group)) +
#   ggh4x::facet_nested(bat_species~era, scales="free", space="free" ) + theme_bw() +
#   theme(axis.text = element_text(size=14), axis.title.y = element_text(size=16),
#         axis.title.x = element_blank(),
#         strip.background = element_rect(fill="white"), 
#         strip.text.y = element_text(size=16, face="italic"), strip.text.x = element_text(size=16, face="italic"),
#         panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
#         plot.margin = unit(c(.5,.5,.5,1), "lines")) +
#   scale_shape_manual(values=yearz) +
#   coord_cartesian(expand = F) + ylab("Virus prevalence") + xlab("Sampling date")+
#   scale_fill_manual(values=colorz) + scale_color_manual(values=colorz) + 
#   guides(size=guide_legend(override.aes = list(color="black")), shape=guide_legend(override.aes = list(color="black")))
# pInterannualShed.All

pInterannualShed.era1 <- ggplot(data=subset(prev.dat, era=="I")) +
  geom_rect(data=subset(cal.dat.save, era=="I"), aes(xmin=start, xmax=stop, ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=as.Date(mid_date), y=prevGhV, size=Ntot, shape=year), color="#00C094") + 
  geom_point(aes(x=as.Date(mid_date), y=prevAngV, size=Ntot, shape=year), color="#F8766D") + 
  scale_fill_manual(values=filz) +
  guides(color = guide_legend(title = "virus"))+  
  guides(fill="none")+
  ggnewscale::new_scale_fill() +
  #geom_linerange(aes(x=as.Date(mid_date), ymin=prevGhV_lci, ymax=prevGhV_uci),color="#00C094", linetype=2, alpha=.3, show.legend = F) +
  #geom_linerange(aes(x=as.Date(mid_date), ymin=prevAngV_lci, ymax=prevAngV_uci),color="#F8766D", linetype=2, alpha=.3, show.legend = F) +
  geom_ribbon(data=subset(inter.pred, era=="I"), aes(x=as.Date(mid_date), ymin= prev_lci, ymax=prev_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(inter.pred, era=="I"), aes(x=as.Date(mid_date), y= prevalence, color=antigen, group=group)) +
  ggh4x::facet_nested(.~bat_species, scales="free", space="free" ) + theme_bw() +
  theme(axis.text = element_text(size=14), axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16, face="italic"), strip.text.x = element_text(size=16, face="italic"),
        plot.title = element_text(size=16),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  scale_shape_manual(values=yearz) +
  coord_cartesian(expand = F) + ylab("Virus prevalence") + xlab("Sampling date")+
  scale_fill_manual(values=colorz) + scale_color_manual(values=colorz) + 
  guides(color="none")+
  ggtitle("Interannual virus prevalence (2013-2016 dataset)")+
  guides(size=guide_legend(override.aes = list(color="black"), title="N bats"), 
         shape=guide_legend(override.aes = list(color="black")))
pInterannualShed.era1

pInterannualShed.era2 <- ggplot(data=subset(prev.dat, era=="II")) +
  geom_rect(data=subset(cal.dat.save, era=="II"), aes(xmin=start, xmax=stop, ymin=0, ymax=1, fill=season), alpha=.3) +
  geom_point(aes(x=as.Date(mid_date), y=prevGhV, size=Ntot, shape=year), color="#00C094") + 
  geom_point(aes(x=as.Date(mid_date), y=prevAngV, size=Ntot, shape=year), color="#F8766D") + 
  scale_fill_manual(values=filz) +
  guides(color = guide_legend(title = "virus"))+  
  guides(fill="none")+
  ggnewscale::new_scale_fill() +
  #geom_linerange(aes(x=as.Date(mid_date), ymin=prevGhV_lci, ymax=prevGhV_uci),color="#00C094", linetype=2, alpha=.3, show.legend = F) +
  #geom_linerange(aes(x=as.Date(mid_date), ymin=prevAngV_lci, ymax=prevAngV_uci),color="#F8766D", linetype=2, alpha=.3, show.legend = F) +
  geom_ribbon(data=subset(inter.pred, era=="II"), aes(x=as.Date(mid_date), ymin= prev_lci, ymax=prev_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(inter.pred, era=="II"), aes(x=as.Date(mid_date), y= prevalence, color=antigen, group=group)) +
  ggh4x::facet_nested(.~bat_species, scales="free", space="free" ) + theme_bw() +
  theme(axis.text = element_text(size=14), axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16, face="italic"), strip.text.x = element_text(size=16, face="italic"),
        plot.title = element_text(size=16),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  scale_shape_manual(values=yearz) +
  coord_cartesian(expand = F) + ylab("Virus prevalence") + xlab("Sampling date")+
  scale_fill_manual(values=colorz) + scale_color_manual(values=colorz) + 
  guides(color="none")+
  ggtitle("Interannual virus prevalence (2018-2020 dataset)")+
  guides(size=guide_legend(override.aes = list(color="black"), title="N bats"), 
         shape=guide_legend(override.aes = list(color="black")))
pInterannualShed.era2


###Put together final figure

library(patchwork)

#v1
inter<-pInterannualShed.era1/pInterannualShed.era2+ plot_annotation(tag_levels = list(c('B','C'))) &
  theme(plot.tag = element_text(size = 23, face="bold"))
inter
inter<-as.ggplot(inter)

Fig1<-plot_grid(tree1, inter, labels=c("A",""),
                rel_widths=c(1.2,1), rel_heights = c(1,1),
                ncol=2,align="hV", axis="l", label_size = 23)+
  theme(plot.background = element_rect(fill="white"))

Fig1


#save figures
homewd="/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_henipavirus_dynamics_paper_2026/final_figures_tables/main/"
ggsave(file = paste0(homewd, "Fig1_polymerase_phylogeny_seasonal_shedding.png"),
       plot= Fig1,
       units="mm",  
       width=360, 
       height=130, 
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

