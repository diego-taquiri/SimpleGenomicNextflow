#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastqDir='/Users/diegotaquiridiaz/Desktop/Proyectos/Proyectos_2023-1/Mirko_Lab/Wastewater/fastq'
params.ref="$launchDir/data/"
params.taxdb='/Users/diegotaquiridiaz/Desktop/Proyectos/Proyectos_2023-1/Mirko_Lab/Wastewater/database'
params.krakenU='/Users/diegotaquiridiaz/Desktop/Proyectos/Proyectos_2023-1/Mirko_Lab/Wastewater/krakenUniq/reports/krakenU'

process ExtractReads_KrakenUniq {
    conda '/Users/diegotaquiridiaz/miniconda3/envs/kraken'
    publishDir 'results/extract_reads', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)
    each taxid
    path taxdb
    path krakenU

    output:
    val sample_id
    path "${sample_id}.taxid${taxid}.fastq"
    val taxid

    script:
    """
    krakenuniq-extract-reads -t $taxdb/TAXDB $taxid \\
    $krakenU/${sample_id}.krakenU.tsv \\
    $query_file > ${sample_id}.taxid${taxid}.fastq

    """
}

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

process Mapping_Minimap2 {
    conda '/Users/diegotaquiridiaz/miniconda3/envs/sam_tools'
    publishDir "results/minimap2", mode: 'copy'

    input:
    val sample_id
    path query_file
    val taxid
    path ref

    output:
    path "${sample_id}.taxid${taxid}.bam"
    path "${sample_id}.taxid${taxid}.bam.bai"
    path "${sample_id}.taxid${taxid}.flagstat.txt"
    val sample_id
    val taxid


    script:
    """
    minimap2 -ax map-ont ${ref}/ref_genome_taxid${taxid}.fasta $query_file \\
        | samtools sort \\
        | samtools view -F 2048 -b > ${sample_id}.taxid${taxid}.bam
	samtools index ${sample_id}.taxid${taxid}.bam
    samtools flagstat ${sample_id}.taxid${taxid}.bam > "${sample_id}.taxid${taxid}.flagstat.txt"

    """
}

process CoverageDepth_Samtools {
    conda '/Users/diegotaquiridiaz/miniconda3/envs/sam_tools'
    publishDir "results/minimap2", mode: 'copy'

    input:
    path bam
    path bam_bai
    path flagstat
    val sample_id
    val taxid
    
    output:
    path "${sample_id}.taxid${taxid}.coverage.txt"

    script:
    """
    samtools coverage $bam > "${sample_id}.taxid${taxid}.coverage.txt"

    """
}

process ExtractMappedReads {
    conda '/Users/diegotaquiridiaz/miniconda3/envs/sam_tools'
    publishDir 'results/mapped_reads', mode: 'copy'

    input:
    path bam
    path bam_bai
    path flagstat
    val sample_id
    val taxid

    output:
    path "${sample_id}.taxid${taxid}.mappedReads.fasta"

    script:
    """
    samtools fasta -F 0x104 $bam > "${sample_id}.taxid${taxid}.mappedReads.fasta"

    """
}

process QCExtractedReads_Nanoplot {
    conda '/Users/diegotaquiridiaz/miniconda3/envs/sam_tools'
    publishDir "results/nanoplot"

    input:
    path bam
    path bam_bai
    path flagstat
    val sample_id
    val taxid

    output:
    path "${sample_id}.taxid${taxid}"

    script:
    """
    NanoPlot --bam $bam -o ${sample_id}.taxid${taxid}

    """
}


workflow {
    ch_fastq = channel.fromFilePairs("${params.fastqDir}/*.fastq.gz", size: 1)
    taxids = ['561','573','470','1352','1496','28901','287','1280','446','487','550','1313','546','622','615','1639','584']
    ExtractReads_KrakenUniq(ch_fastq, taxids, params.taxdb, params.krakenU)
    Ref_index_Samtools(params.ref, taxids)
    Mapping_Minimap2(ExtractReads_KrakenUniq.out, params.ref)
    CoverageDepth_Samtools(Mapping_Minimap2.out)
    ExtractMappedReads(Mapping_Minimap2.out)
    

}
