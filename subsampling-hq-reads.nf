#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastqDir = "$launchDir/data/20240319-amplicones-KV-CB"
params.ref="$launchDir/reference/reference.fasta"
params.scriptsDir="$launchDir/scripts"

process ConcatFastq {
    publishDir 'results/fastq_cat', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz")
    
    script:
    """
    cat ${sample_id}/*.fastq.gz > ${sample_id}.fastq.gz

    """
}

process QualityControl_Nanoplot {
    container "quay.io/biocontainers/nanoplot:1.42.0--pyhdfd78af_0"
    publishDir './results/nanoplot', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)

    output:
    path "${sample_id}"

    script:
    """
    NanoPlot --fastq $query_file -o $sample_id

    """
}

process Trimming_Nanofilt {
    container "quay.io/biocontainers/nanofilt:2.8.0--py_0"

    conda '/home/crowfoot2/anaconda3/envs/nanofilt'
    publishDir 'results/nanofilt', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)

    output:
    tuple val(sample_id), path ("${sample_id}_hq/trimmed.${sample_id}.fastq.gz")


    script:
    """
    mkdir ${sample_id}_hq

    gunzip -c $query_file \\
    | NanoFilt -q 15 -l 1000 \\
    | gzip > ${sample_id}_hq/trimmed.${sample_id}.fastq.gz
    """
}

process Trimmed_Nanoplot {
    container "quay.io/biocontainers/nanoplot:1.42.0--pyhdfd78af_0"
    publishDir './results/nanoplot_trimmed', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)

    output:
    path "${sample_id}"

    script:
    """
    NanoPlot --fastq $query_file -o $sample_id

    """
}

workflow {
    channel_fastq = channel.fromFilePairs("$params.fastqDir/*",type:"dir",size:1)
    ConcatFastq(channel_fastq)
    QualityControl_Nanoplot(ConcatFastq.out)
    Trimming_Nanofilt(ConcatFastq.out)
    Trimmed_Nanoplot(Trimming_Nanofilt.out)
    
}
