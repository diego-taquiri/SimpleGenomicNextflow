#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.dataDir = "$launchDir/data/20240314_Mmar_KV"
params.ref="$launchDir/reference/reference.fasta"
params.scriptsDir="$launchDir/scripts"
params.checkm2DB="/data/data_analysis/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd"

process ConcatFastq {
    publishDir 'results/fastq_cat', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)

    output:
    tuple val(sample_id), path("${sample_id}_combined.fastq.gz")
    
    shell:
    '''    
    find ./!{sample_id}/ -type f -name "*.fastq.gz" -print0 | xargs -0 cat > "!{sample_id}_combined.fastq.gz"

    '''
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

process Assembly_flye {
    container "quay.io/biocontainers/flye:2.9.3--py310h2b6aa90_0"
    cpus=5
    publishDir 'results/flye', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)

    output:
    tuple val(sample_id), path ("assembly.${sample_id}")

    script:
    """
    flye --nano-hq $query_file \\
    --iterations 1 \\
    --out-dir assembly.${sample_id} \\
    -t ${task.cpus}
    """
}

process Contam_checkm2 {
    container "quay.io/biocontainers/checkm2:1.0.1--pyh7cba7a3_0"
    cpus=6
    publishDir 'results/checkm2', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)
    path checkm2db

    output:
    path "checkm2.${sample_id}"
    
    script:
    """
    checkm2 predict --threads ${task.cpus} -x .fasta\\
     --input ${query_file}/assembly.fasta\\
     --output-directory checkm2.${sample_id} \\
     --database_path $checkm2db
    """

}

process Mapping_Minimap2 {
    container "quay.io/biocontainers/minimap2:2.28--he4a0461_0"
    publishDir './results/minimap2', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)
    path ref

    output:
    path "${sample_id}.bam"
    path "${sample_id}.bam.bai"
    path "${sample_id}_flagstat.txt"
    val sample_id


    script:
    """
    minimap2 -ax map-ont $ref $query_file \\
        | samtools sort \\
        | samtools view -F 2048 -b > ${sample_id}.bam
	samtools index ${sample_id}.bam
    samtools flagstat ${sample_id}.bam > "${sample_id}_flagstat.txt"

    """
}

workflow {
    channel_data = channel.fromFilePairs("$params.dataDir/muestra*",type:"dir",size:1)
    ConcatFastq(channel_data)
    QualityControl_Nanoplot(ConcatFastq.out)
    Assembly_flye(ConcatFastq.out)
    Contam_checkm2(Assembly_flye.out,params.checkm2DB)

}
