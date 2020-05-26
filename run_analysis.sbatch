#!/bin/bash
#SBATCH -n 1 
#SBATCH --time=08:00:00

nextflow run code/rna_pipe.nf -work-dir /pine/scr/j/r/jraab/agracz/work/  -resume



