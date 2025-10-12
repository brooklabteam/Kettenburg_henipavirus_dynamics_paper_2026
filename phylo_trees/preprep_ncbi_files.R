rm(list=ls())

library(plyr)
library(dplyr)


#set wd
homewd = "/Users/gwenddolenkettenburg/Desktop/brook_lab/henipavirus_shedding_cross_correlation_project/phylo_trees"

#Just change out each folder throughout the document
setwd(paste0(homewd, "/ncbi_csv/"))

#load the dataset and query
dat <- read.csv(file = "ncbi_full_genome_refseq.csv", header = T, stringsAsFactors = F)
head(dat)

accession_num <- paste(c(dat$Accession), collapse = ",")

#now put this into your webbrowser to download
text.for.NCBI <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&retmode=text&id=",accession_num)
text.for.NCBI
#once downloaded, send to MAFFT for alignment