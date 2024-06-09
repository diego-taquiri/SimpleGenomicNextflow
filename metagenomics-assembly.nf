#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastqDir = "$launchDir/data"
params.ref="$launchDir/ref_genomes"

process Assembly_flye {
    conda '/home/tao/anaconda3/envs/flye'
    cpus=3
    publishDir 'results/flye', mode: 'copy'

    input:
    val sample_id
    path query_file
    val taxid
    output:
    val sample_id
    path "assembly.${taxid}"
    val taxid
    script:
    """
    flye --nano-hq $query_file --iterations 1 --out-dir assembly.${taxid} \\
    -t ${task.cpus}
    """

}

process Contam_checkm2 {
    conda '/home/tao/anaconda3/envs/checkm2'
    cpus=6
    publishDir 'results/checkm2', mode: 'copy'

    input:
    val sample_id
    path query_file
    val taxid
    output:
    path "checkm2_${sample_id}_${taxid}"
    script:
    """
    checkm2 predict --threads ${task.cpus} -x .fasta\\
     --input ${query_file}/assembly.fasta\\
     --output-directory checkm2_${sample_id}_${taxid}
    """
} 

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
    Assembly_flye(FromBAM_toFastq.out)
    Contam_checkm2(Assembly_flye.out)
    Polishing_medaka(channel_fastq, Assembly_flye.out)
    Checkm2_polished(Polishing_medaka.out)

}