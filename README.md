# Repo: Kettenburg_dynamics_paper_2026 for paper linked [here]()

This is a GitHub repository for the project to describe novel circulating henipavirus clades, and antibody/virus shedding dynamics
in Malagasy fruit bats. The github is organized by folder, with each analysis corresponding to the names of each folder. The R 
script to generate each figure is within its respective analysis folder, and final figures and tables are within the 'final_figures_tables' folder. 
A brief description of each anlysis is outlined under each subheading below. All R scripts are annotated, and non-truncated versions of analyses are
hosted within each analysis folder linked below: 

---
## Final figures and tables

The exploratory scripts for the analyses presented in the figures are housed within scripts in the respective analysis folders listed below. 

-All scripts to make final figures for main text, extended data, and supplemental figures is here [here](https://github.com/brooklabteam/Kettenburg_henipavirus_dynamics_paper_2026/tree/main/final_figures_tables)


---
## Genome characterization and annotation

Raw reads from sequenced samples were run on the CZID (now SeqtoID) platform, excluding "bat" host reads. Assembled contigs mapping to paramyxoviruses 
were manually BLAST'd, classified into genes, then further processed phylogenetically (see next section)

All genome sequences were submitted to NCBI and assigned accession numbers PX095965-PX095971 (Sanger sequences), PX127342-PX127344 (AndV full and nearly-full genomes), 
PX127345-PX127510 (all other partial sequences). 

Specific analyses (alignments with MAFFT within Geneious, figures made with ggplot in R): 

-Read mapping - folder for this analysis is [here](https://github.com/brooklabteam/Kettenburg_henipavirus_dynamics_paper_2026/tree/main/genome_characterization)

-mRNA editing site - folder for this analysis is [here](https://github.com/brooklabteam/Kettenburg_henipavirus_dynamics_paper_2026/tree/main/genome_characterization/residue_alignments)

-Ephrin binding sites - folder for this analysis is [here](https://github.com/brooklabteam/Kettenburg_henipavirus_dynamics_paper_2026/tree/main/genome_characterization/residue_alignments)


---
## Phylogenetic tree analysis

Using either full genome or partial/full polymerase sequences, we phylogenetically classified novel henipaviruses and pararubulaviruses. 

Specific analyses (alignment with MAFFT within Geneious, substitution model selection and tree building done in IQ-TREE, figures made with ggtree/ggplot in R): 

-Paramyxovirus trees - folder for this analysis is [here](https://github.com/brooklabteam/Kettenburg_henipavirus_dynamics_paper_2026/tree/main/phylo_trees)


---
## Virus shedding

Using full genome sequences for GhV, AngV, and AndV, we designed multiplex qPCR panels for genes L, M, and N. The two most sensitive assays were used to confirm
positivity. 

Specific analyses (statistical analysis of seasonality, body condition correlation with viral shedding, coinfection, and sampling session prevalence correlation 
with virus shedding using GAMs/GLMs/LMs run in R, figures made with ggplot in R):

-qPCR data and scripts - folder for this analysis is [here](https://github.com/brooklabteam/Kettenburg_henipavirus_dynamics_paper_2026/tree/main/virus_shedding)


---
## Serology

Using a previously designed multi-plex panel of henipavirus glycoproteins, we added AngV and ran the assay on sera samples. 

Specific analyses (statistical analysis of seasonality, body condition correlation with MFI, using GAMs/GLMs/LMs run in R, age-seroprevalence curves, estimation 
of duration of maternally-inhereted immunity and adult immunity, figures made with ggplot in R)

-Seropositive data and scripts - folder for this analysis is [here](https://github.com/brooklabteam/Kettenburg_henipavirus_dynamics_paper_2026/tree/main/serology)

---
## Mechanistic modeling

Fit MSIR models using joint inference from three datasets: age-seroprevalence data, shedding data by biweek, and seropositive data by biweek

Specific anlyses (figures made with ggplot in R)

-MSIR model fitting - folder for this analysis is [here](https://github.com/brooklabteam/Kettenburg_henipavirus_dynamics_paper_2026/tree/main/model_fitting)






