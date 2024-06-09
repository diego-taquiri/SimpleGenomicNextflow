#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastqDir = "$launchDir/data/TB_Run5_23/fastq_pass"
params.ref="$launchDir/data/M_tuberculosis_H37Rv_genome_NC_000962.3.fasta"
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

/*
process Ref_index_Samtools {
    conda '/Users/diegotaquiridiaz/miniconda3/envs/sam_tools'
    publishDir "results/ref_genomes", mode: 'copy'

    input:
    path ref
    each taxid

    output:
    path "ref_genome_taxid${taxid}.fasta"
    path "ref_genome_taxid${taxid}.fasta.fai"

    script:
    """
    cp ${ref}/ref_genome_taxid${taxid}.fasta ref_genome_taxid${taxid}.fasta
    samtools faidx ref_genome_taxid${taxid}.fasta > ref_genome_taxid${taxid}.fasta.fai

    """
}
*/
process Mapping_Minimap2 {
    container "quay.io/biocontainers/minimap2:2.27--he4a0461_1"
    publishDir './results/minimap2', mode: 'copy'
    cpus 14

    input:
    tuple val(sample_id), path(query_file)
    path ref

    output:
    path "${sample_id}.sam"
    val sample_id

    script:
    """
    minimap2 -ax map-ont \\
    -t $task.cpus \\
    $ref \\
    $query_file \\
    -o ${sample_id}.sam

    """
}

process Sort_Samtools {
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
    publishDir './results/minimap2', mode: 'copy'
    cpus 14

    input:
    path sam
    val sample_id
    
    output:
    path "${sample_id}.bam"
    path "${sample_id}.bam.bai"
    path "${sample_id}_flagstat.txt"
    val sample_id

    script:
    """
    samtools sort $sam \\
        --threads $task.cpus \\
        | samtools view -F 2048 -b > ${sample_id}.bam
    samtools index ${sample_id}.bam
    samtools flagstat ${sample_id}.bam > "${sample_id}_flagstat.txt"
    """
}


process CoverageDepth_Samtools {
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
    publishDir './results/minimap2', mode: 'copy'

    input:
    path bam
    path bam_bai
    path flagstat
    val sample_id
    
    output:
    path "${sample_id}_coverage.txt"

    script:
    """
    samtools coverage $bam > "${sample_id}_coverage.txt"

    """
}

process SummaryStats {
    publishDir '.', mode: 'copy'

    input:
    path nanoplot
    path coverage

    output:
    path './results/summary/final_table.txt'
    
    shell:
    '''
    bash !{params.scriptsDir}/summary_table_0.2.sh 
    
    '''
}

workflow {
    channel_fastq = channel.fromFilePairs("$params.fastqDir/barcode_unclassified",type:"dir",size:1)
    ConcatFastq(channel_fastq)
    QualityControl_Nanoplot(ConcatFastq.out)
    Mapping_Minimap2(ConcatFastq.out, params.ref)
    Sort_Samtools(Mapping_Minimap2.out)
    CoverageDepth_Samtools(Sort_Samtools.out)
    SummaryStats(QualityControl_Nanoplot.out.collect(),CoverageDepth_Samtools.out.collect())
    
}
