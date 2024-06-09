#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastqDir = "$launchDir/data/bam/"
params.gene_bed="$launchDir/bed/carlos.format.bed"

process Sort_Samtools {
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
    publishDir './results/minimap2'
    cpus 3

    input:
    tuple val (sample_id), path (bam)
    
    output:
    val sample_id
    path "${sample_id}.bam"
    path "${sample_id}.bam.bai"

    script:
    """
    samtools index $bam
    """
}

process CoverageGeneDepth_Mosdepth {
    container "quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0"
    publishDir './results/mosdepth', mode: 'copy'
    cpus 4

    input:
    val sample_id 
    path bam
    path bai
    path bed
    
    output:
    path "mosdepth.${sample_id}.regions.bed"

    script:
    """
    mosdepth -F 3844 \\
    -t $task.cpus \\
    -b $bed \\
    "mosdepth.${sample_id}" \\
    $bam

    gunzip mosdepth.${sample_id}.regions.bed.gz
    """
}


workflow {
    channel_fastq = channel.fromFilePairs("$params.fastqDir/Barcode*.bam",size:1)
    //channel_fastq.view()
    Sort_Samtools(channel_fastq)
    CoverageGeneDepth_Mosdepth(Sort_Samtools.out,params.gene_bed)

}
