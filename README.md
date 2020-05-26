[![DOI](https://zenodo.org/badge/267152732.svg)](https://zenodo.org/badge/latestdoi/267152732)


# Code associated with Tulasi et. al 2020

The code was run using nextflow using the Slurm workload manager.
Assuming you have the fastq files from the manuscript, you could re-reun this code using 
`sbatch run_analysis.sh`

You may  need to change `nextflow.config` to match your system

In `code/rna_pipe.nf` there are a few paths that should be altered to match your desired output locations

In `data/rna/biliar_sox9_agracz_samples.csv` you can change the location of the fastq files in the sample sheet to match where you have the fastq files. 


Please contact jraab@med.unc.edu with any questions

# Software/Versions used in this repo
* fastqc/0.11.8
* multiqc/1.7
* Salmon/1.1.0
* R/3.6.0

* R packages
    * DESeq2/1.26.0
    * tximport/1.14.2
    * tidyverse/1.3.0
    * biomaRt/2.42.1
    * UpSetR/1.4.0
    * ComplexHeatmap/2.2.0
    * singscore/1.6.0
    * broom/0.5.6



