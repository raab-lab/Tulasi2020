#!/usr/bin/env nextflow 
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/
params.outdir = '/proj/jraablab/users/jraab/collabs/agracz/sox9_biliary/results/'
params.scratch = 'pine/scr/j/r/jraab/agracz/'
params.transcripts = '/proj/jraablab/users/share/resources/transcripts/gencode.vM24.transcripts.fa.gz' 


ch_genome_fa = file('/proj/jraablab/users/share/resources/genomes/GRCm38.primary_assembly.genome.fa')
ch_star_gtf  = file('/proj/jraablab/users/share/resources/transcripts/gencode.vM24.annotation.gtf')


Channel.fromPath(params.transcripts)
   .into { index_transcripts; alignment_transcripts} 

designFile = file("data/rna/biliary_sox9_agracz_samples.csv")
Channel.from(designFile)
       .splitCsv(header: true)
       .map {row -> tuple(row.Sample, file(row.File), row.Group ) } 
       .into{ ch_fastq_quant; ch_fastq_align; ch_fastq_qc } 
//ch_fastq_qc.subscribe { println it } 

process unzip_transcripts { 
   input: 
   file(fa) from alignment_transcripts

   output: 
   file('unzipped.fa') into ch_unzipped

   script:
   """
   gunzip -c ${fa}  > unzipped.fa
   """
} 

// QC
process fastc { 
   label 'small'
   module 'fastqc/0.11.8'
   
   input: 
   set val(name), file(fastq), val(group) from ch_fastq_qc

   output: 
   file("*.{zip,html}") into ch_qc_files

   script: 
   """
   fastqc $fastq  -t ${task.cpus} && touch tmp.txt
   """
 }

// Summarise QC
process multiqc { 
   label 'small'
   module 'multiqc/1.7'
 
   publishDir "${params.outdir}/qc", mode: 'copy'
   input: 
   file(html) from ch_qc_files.collect()

   output: 
   file('multiqc_report.html') into ch_multiqc

   script: 
   """
   multiqc .
   """
}

// Salmon Index
process  salmon_index {  
   label 'medium'
   module 'salmon/1.1.0' 
   
   input: 
   file(transcripts)  from index_transcripts
   
   output: 
   file("${transcripts.baseName}.idx") into salmon_index

   script: 
   """
   salmon index -t ${transcripts} -i ${transcripts.baseName}.idx
   """

} 


// Salmon Count
process salmon_count{ 

   tag {name}
   module 'salmon/1.1.0'
   label 'medium'

   publishDir "${params.outdir}/quant/", mode: 'copy'

   input: 
   set val(name), file(fn), val(group) from ch_fastq_quant
   file(idx) from salmon_index.collect()

   output: 
   file("${name}/*") into ch_quant

   script: 
   """
   salmon quant -i ${idx} -l A -r ${fn} -p ${task.cpus} -o ${name} --writeUnmappedNames --seqBias --gcBias --validateMappings
   """
} 

  
