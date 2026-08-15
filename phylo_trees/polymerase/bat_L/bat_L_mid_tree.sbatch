#!/bin/bash
#SBATCH --account=ac_brooklab2
#SBATCH --job-name=b_L_mid_tree
#SBATCH --partition=savio3
#SBATCH --output=bat_L_mid_tree
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gkettenburg@berkeley.edu

module load iqtree

iqtree3 -s bat_align_L_mid.fasta -m MFP -bb 1000 --prefix batLmid -T 8