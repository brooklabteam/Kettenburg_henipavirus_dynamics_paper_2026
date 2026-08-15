#!/bin/bash
#SBATCH --account=ac_brooklab2
#SBATCH --job-name=b_G_right_tree
#SBATCH --partition=savio3
#SBATCH --output=bat_G_right_tree
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gkettenburg@berkeley.edu

module load iqtree

iqtree3 -s bat_align_G_right.fasta -m MFP -bb 1000 --prefix batGright -T 8