#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastqDir = "$launchDir/data" 

process ConcatFastq {
    conda "/home/crowfoot2/anaconda3/envs/tbp"
    publishDir './results/concat', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz")
    
    script:
    """
    cat ${sample_id}/*.fastq.gz > ${sample_id}.fastq.gz

    """
}

process Tbprofiler {
    conda "/home/crowfoot2/anaconda3/envs/tbp"
    publishDir './results/tbprofiler', mode: 'copy'
    cpus 30
    maxForks 1

    input:
    tuple val(sample_id), path(query_file)
    
    output:
    val sample_id
    path "${sample_id}.tbprofiler"


    script:
    """
    mkdir -p ${sample_id}.tbprofiler

    tb-profiler profile \\
    --read1 $query_file/${sample_id}_d_s.fastq \\
    --platform nanopore \\
    -t $task.cpus \\
    -p ${sample_id}.tbprofiler \\
    --txt

    mv results ${sample_id}.tbprofiler
    mv bam ${sample_id}.tbprofiler
    mv vcf ${sample_id}.tbprofiler
    """
}

process Tbprofiler_jsons {
    conda "/home/crowfoot2/anaconda3/envs/tbp"
    publishDir './results/tbprofiler-jsons/results', mode: 'copy'
    cpus 30
    maxForks 1

    input:
    val sample_id
    path "${sample_id}.tbprofiler"
    
    output:
    path "${sample_id}.tbprofiler.results.json"


    script:
    """
    cp "${sample_id}.tbprofiler/results/${sample_id}.tbprofiler.results.json" .
    """
}

process Tbprofiler_collate {
    conda "/home/crowfoot2/anaconda3/envs/tbp"
    publishDir './results/tbprofiler-collate', mode: 'copy'
    cpus 1
    maxForks 1

    input:
    path jsons
    
    output:
    path "tbprofiler.*"


    script:
    """
    mkdir -p results
    cp *.json results
    tb-profiler collate
    """
}

workflow {
    channel_fastq = channel.fromFilePairs("$params.fastqDir/barcode*",size:1,type:'dir')
    //channel_fastq = channel.fromFilePairs("$params.fastqDir/barcode*.fastq.gz",size:1)
    //channel_fastq.view()
    ConcatFastq(channel_fastq)
    Tbprofiler(channel_fastq)
    Tbprofiler_jsons(Tbprofiler.out)   
    Tbprofiler_collate(Tbprofiler_jsons.out.collect())

}