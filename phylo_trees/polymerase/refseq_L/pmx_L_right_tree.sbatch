#!/bin/bash
#SBATCH --account=ac_brooklab2
#SBATCH --job-name=p_L_right_tree
#SBATCH --partition=savio3
#SBATCH --output=pmx_L_right_tree
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gkettenburg@berkeley.edu

module load iqtree

iqtree3 -s refseq_align_L_right.fasta -m MFP -bb 1000 --prefix pmxLright -T 8