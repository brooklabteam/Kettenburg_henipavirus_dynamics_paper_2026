rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)

homewd = "/Users/carabrook/Developer/luminex-serology-madabat"

setwd(homewd)

#load raw data
lum.dat <- read.csv(file = paste0(homewd,"/data/lum_raw.csv"), header = T, stringsAsFactors = F)
head(lum.dat)


#clean

sort(unique(lum.dat$Sample))

lum.dat = subset(lum.dat, Sample!="Neg 2.5 µg/ml" & Sample!="Pos 2.5 µg/ml" & Sample!="α-NiV 0.625 µg/ml")

#average over replicates where there are some
names(lum.dat)[names(lum.dat)=="Sample"] <- "sampleid"

#convert to numeric
lum.dat$HeV = as.numeric(as.character(lum.dat$HeV))
lum.dat$NiV =  as.numeric(as.character(lum.dat$NiV))
lum.dat$CedV = as.numeric(as.character(lum.dat$CedV))
lum.dat$MojV = as.numeric(as.character(lum.dat$MojV))
lum.dat$GhV =  as.numeric(as.character(lum.dat$GhV))
lum.dat$AngV = as.numeric(as.character(lum.dat$AngV))
lum.dat$EboV = as.numeric(as.character(lum.dat$EboV))
lum.dat$BdbV = as.numeric(as.character(lum.dat$BdbV))
lum.dat$BomV = as.numeric(as.character(lum.dat$BomV))
lum.dat$SudV = as.numeric(as.character(lum.dat$SudV))
lum.dat$RestV = as.numeric(as.character(lum.dat$RestV))
lum.dat$RavV = as.numeric(as.character(lum.dat$RavV))
lum.dat$MlaV =  as.numeric(as.character(lum.dat$MlaV))
lum.dat$LloV =  as.numeric(as.character(lum.dat$LloV))
lum.dat$Mock =  as.numeric(as.character(lum.dat$Mock))



lum.dat$sampleid[duplicated(lum.dat$sampleid)]

#a few of these are problem samples in the dataset that are indistinguishable. remove here.
 lum.dat$keep = 1 #
 lum.dat$keep[lum.dat$sampleid=="MIZ305"|
                lum.dat$sampleid=="MIZ306"|
                lum.dat$sampleid=="MIZ307"|
                lum.dat$sampleid=="MIZ545"|
                lum.dat$sampleid=="MIZ546"|
                lum.dat$sampleid=="MIZ547"|
                lum.dat$sampleid=="MIZ548" |
                lum.dat$sampleid=="MIZ549"] <- 0

lum.dat = subset(lum.dat, keep==1)
lum.dat <- dplyr::select(lum.dat, -(keep))

# now split by plate and control for binding negative control beads
lum.split <- dlply(lum.dat, .(Plate))

correct.binding <- function(df){
  control =subset(df, sampleid=="Blank") 
  out =subset(df, sampleid!="Blank") 
  control.avg <- ddply(control,.(sampleid), summarise, HeV = mean(HeV, na.rm=T), NiV = mean(NiV, na.rm=T), CedV = mean(CedV, na.rm=T), MojV = mean(MojV, na.rm=T), GhV = mean(GhV, na.rm=T), AngV = mean(AngV, na.rm=T), EboV = mean(EboV, na.rm=T), BdbV = mean(BdbV, na.rm=T), BomV = mean(BomV, na.rm=T), SudV = mean(SudV, na.rm=T), RestV = mean(RestV, na.rm=T), RavV = mean(RavV, na.rm=T), MlaV = mean(MlaV, na.rm=T), LloV = mean(LloV, na.rm=T), Mock = mean(Mock, na.rm=T) )
  
  #and subtract from all values 
 
  
  out$HeV =  out$HeV - control.avg$HeV 
  out$NiV = out$NiV - control.avg$NiV 
  out$CedV = out$CedV - control.avg$CedV 
  out$MojV = out$MojV - control.avg$MojV 
  out$GhV = out$GhV  - control.avg$GhV  
  out$AngV = out$AngV - control.avg$AngV 
  out$EboV = out$EboV - control.avg$EboV 
  out$BdbV = out$BdbV - control.avg$BdbV 
  out$BomV = out$BomV - control.avg$BomV 
  out$SudV = out$SudV - control.avg$SudV 
  out$RestV = out$RestV - control.avg$RestV 
  out$RavV = out$RavV - control.avg$RavV 
  out$MlaV = out$MlaV - control.avg$MlaV 
  out$LloV =  out$LloV - control.avg$LloV 
  out$Mock = out$Mock -control.avg$Mock 
  
    return(out)
    
    
    
    
    
  
}

lum.dat.new <- data.table::rbindlist(lapply(lum.split, correct.binding))





#now, correct by mock
lum.dat.new$HeV = lum.dat.new$HeV - lum.dat.new$Mock
lum.dat.new$NiV =  lum.dat.new$NiV -lum.dat.new$Mock
lum.dat.new$CedV = lum.dat.new$CedV - lum.dat.new$Mock
lum.dat.new$MojV = lum.dat.new$MojV -lum.dat.new$Mock
lum.dat.new$GhV =  lum.dat.new$GhV - lum.dat.new$Mock
lum.dat.new$AngV = lum.dat.new$AngV -lum.dat.new$Mock
lum.dat.new$EboV = lum.dat.new$EboV -lum.dat.new$Mock
lum.dat.new$BdbV = lum.dat.new$BdbV -lum.dat.new$Mock
lum.dat.new$BomV = lum.dat.new$BomV -lum.dat.new$Mock
lum.dat.new$SudV = lum.dat.new$SudV -lum.dat.new$Mock
lum.dat.new$RestV = lum.dat.new$RestV -lum.dat.new$Mock
lum.dat.new$RavV = lum.dat.new$RavV -lum.dat.new$Mock
lum.dat.new$MlaV =  lum.dat.new$MlaV -lum.dat.new$Mock
lum.dat.new$LloV =  lum.dat.new$LloV -lum.dat.new$Mock


#there are some samples for which there were two aliquots
unique(lum.dat.new$sampleid[duplicated(lum.dat.new$sampleid)])

#average those
lum.avg <- ddply(lum.dat.new, .(sampleid), summarise, HeV = mean(HeV), NiV = mean(NiV), CedV = mean(CedV), MojV = mean(MojV), GhV = mean(GhV), AngV = mean(AngV), EboV = mean(EboV), BdbV = mean(BdbV), BomV = mean(BomV), SudV = mean(SudV), RestV = mean(RestV), RavV = mean(RavV), MlaV = mean(MlaV), LloV = mean(LloV), Mock = mean(Mock) )

head(lum.avg) #1178


#merge with catching data
#link to wherever this private repo is cloned on your home computer
datawd = "/Users/carabrook/Developer/Madagascar-Bat-Data/"

catch.dat <- read.csv(file = paste0(datawd,"/Catching-Data/catching_data_clean.csv"), header = T, stringsAsFactors = F)
head(catch.dat)
names(catch.dat)

#collect only the relevant columns for merging
catch.dat <- dplyr::select(catch.dat, 1:33, mass_forearm_residual, recap_new)

head(catch.dat)

#there are still some names that don't show up in the catching data
setdiff(unique(lum.avg$sampleid),unique(catch.dat$sampleid)) #these are fine to let go


lum.merge <- merge(lum.avg, catch.dat, by="sampleid") #1172

head(lum.merge) 

unique(lum.merge$bat_species) #"Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis"
 
unique(lum.merge$duplicated_sampleid)
lum.merge <- dplyr::select(lum.merge, -(duplicated_sampleid))
#hist(lum.merge$GhV[lum.merge$bat_species=="Rousettus madagascariensis"])

#save to data folder. later you will need to melt it to determine seropositives
write.csv(lum.merge, file = paste0(homewd,"/data/lum_catch_merge.csv"), row.names = F)

