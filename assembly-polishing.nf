#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.fastqDir = "$launchDir/combined/results/extract_reads/filtered"

process Polishing_medaka {
    conda '/home/tao/anaconda3/envs/medaka'
    cpus=6
    publishDir 'results/medaka', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)
    val sample_id2
    path assembly
    each taxid
    output:
    val sample_id
    path "medaka.${taxid}"
    val taxid 
    script:
    """
    export CUDA_VISIBLE_DEVICES=""
    medaka_consensus -i $query_file \\
    -d $assembly/assembly.fasta \\
    -o medaka.${taxid} \\
    -m r941_min_sup_g507 \\
    -t ${task.cpus} -b 150
    """
}
process Checkm2_polished {
    conda '/home/tao/anaconda3/envs/checkm2'
    cpus=10
    publishDir 'results/checkm2_polished', mode: 'copy'
    input:
    val sample_id
    path query_file
    val taxid
    output:
    val sample_id
    path "polished_${sample_id}_${taxid}"
    val taxid
    script:
    """
    checkm2 predict --threads ${task.cpus} -x .fasta\\
     --input ${query_file}/consensus.fasta\\
     --output-directory polished_${sample_id}_${taxid}
    """
}
workflow {
    channel_fastq = channel.fromFilePairs("${params.fastqDir}/*.fastq.gz", size:1)
    channel_fastq.view()
    taxids=['562','573','470','1352','1496','28901']
    Mapping_Minimap2(channel_fastq, taxids, params.ref)
    FromBAM_toFastq(Mapping_Minimap2.out)
    Assembly_flye(FromBAM_toFastq.out)
    Contam_checkm2(Assembly_flye.out)
    Polishing_medaka(channel_fastq, Assembly_flye.out)
    Checkm2_polished(Polishing_medaka.out)

}