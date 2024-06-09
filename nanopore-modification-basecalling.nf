#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastqDir = "$launchDir/results/pod5"
params.ref="$launchDir/data/M_tuberculosis_H37Rv_genome_NC_000962.3.fasta"
params.scriptsDir="$launchDir/scripts"

process ConcatPod5 {
    publishDir 'results/pod5_cat', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)

    output:
    tuple val(sample_id), path("${sample_id}.pod5")
    
    script:
    """
    pod5 merge ${sample_id}/*.pod5 -o ${sample_id}.pod5

    """
}

process Basecall_Methylation {
    publishDir 'results/basecalled_meth', mode: 'copy'
    maxForks 1

    input:
    tuple val(sample_id), path(query_file)

    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    shell:
    '''
    export PATH=$PATH:/data/programs/dorado-0.6.0-linux-x64/bin

    dorado basecaller sup@v3.3,5mCG_5hmCG !{sample_id}.pod5 > !{sample_id}.bam
    '''
}

workflow {
    channel_fastq = channel.fromFilePairs("$params.fastqDir/barcode*",type:"dir",size:1)
    ConcatPod5(channel_fastq)
    Basecall_Methylation(ConcatPod5.out)
}
