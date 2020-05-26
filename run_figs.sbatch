#!/bin/bash
#!SBATCH -n 4 
#!SBATCH --mem=32G
#! --time=04:00:00

module load r/3.6.0

Rscript --vanilla code/01_DE_RNA.R
Rscript --vanilla code/02_Heatmap_DE_rna.R
Rscript --vanilla code/03_signatures.R


