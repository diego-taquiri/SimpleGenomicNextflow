#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastqDir = "$launchDir/data"
params.ref="$launchDir/reference/reference.fasta"
params.scriptsDir="$launchDir/scripts"

process ConcatFastq {
    conda '/home/crowfoot2/anaconda3/envs/kraken'
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
    conda '/home/crowfoot2/anaconda3/envs/filt'
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
    conda '/home/crowfoot2/anaconda3/envs/nanofilt'
    publishDir 'results/nanofilt', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)

    output:
    tuple val(sample_id), path ("trimmed.${sample_id}.fastq.gz")


    script:
    """
    gunzip -c $query_file \\
    | NanoFilt -q 18 -l 4000 \\
    | gzip > trimmed.${sample_id}.fastq.gz
    """
}

process Trimmed_Nanoplot {
    conda '/home/crowfoot2/anaconda3/envs/filt'
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

process Mapping_Minimap2 {
    conda '/home/crowfoot2/anaconda3/envs/medaka'
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
        | samtools view -F 2304 -b > ${sample_id}.bam
	samtools index ${sample_id}.bam
    samtools flagstat ${sample_id}.bam > "${sample_id}_flagstat.txt"

    """
}

process CoverageDepth_Samtools {
    conda '/home/crowfoot2/anaconda3/envs/medaka'
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

process Consensus_Samtools {
    conda '/home/crowfoot2/anaconda3/envs/medaka'
    publishDir './results/consensus', mode: 'copy'

    input:
    path bam
    path bam_bai
    path flagstat
    val sample_id
    
    output:
    path "consensus_${sample_id}.fasta"

    script:
    """
    samtools consensus --show-del yes \\
     $bam \\
     -o  consensus_${sample_id}.fasta

    """
}

/*
process SummaryStats {
    conda '/home/crowfoot2/anaconda3/envs/medaka'
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
*/
workflow {
    channel_fastq = channel.fromFilePairs("$params.fastqDir/*",type:"dir",size:1)
    ConcatFastq(channel_fastq)
    QualityControl_Nanoplot(ConcatFastq.out)
    Trimming_Nanofilt(ConcatFastq.out)
    Trimmed_Nanoplot(Trimming_Nanofilt.out)
    Mapping_Minimap2(Trimming_Nanofilt.out, params.ref)
    CoverageDepth_Samtools(Mapping_Minimap2.out)
    Consensus_Samtools(Mapping_Minimap2.out)
    /*SummaryStats(QualityControl_Nanoplot.out.collect(),CoverageDepth_Samtools.out.collect())*/
}
