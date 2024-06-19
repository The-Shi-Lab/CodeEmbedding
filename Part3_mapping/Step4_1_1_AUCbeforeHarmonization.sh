#!/bin/bash

#SBATCH --job-name=burden
#SBATCH --error=/net/snowwhite/home/jiacong/EHRmapping/mapping2/slurm/%A_%a.err.txt
#SBATCH --output=/net/snowwhite/home/jiacong/EHRmapping/mapping2/slurm/%A_%a.out
#SBATCH --mem=120G
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1

Rscript /net/snowwhite/home/jiacong/EHRmapping/mapping2/scripts/Step4_1_1_AUCbeforeHarmonization.R