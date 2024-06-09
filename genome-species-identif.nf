#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.metabatDir="$launchDir/data/metabat2_result/"
params.ref="$launchDir/genomes/ZymoBIOMICS.STD.refseq.v3/Genomes/zymo_all_ref_genomes.fasta"


process GTDBtk {
    conda '/home/crowfoot2/anaconda3/envs/gtdbtk-2.1.1'
    cpus=32
    publishDir './results/gtdb', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)
    
    output:
    path "gtdb_${sample_id}/gtdb_identify_${sample_id}"
    path "gtdb_${sample_id}/gtdb_align_${sample_id}"
    path "gtdb_${sample_id}/gtdb_classify_${sample_id}"



    script:
    """
    mkdir gtdb_${sample_id}

    gtdbtk identify --genome_dir ${query_file}/consensus.fasta.metabat-bins*/ \\
    --out_dir gtdb_${sample_id}/gtdb_identify_${sample_id} \\
    -x fa

    gtdbtk align --identify_dir gtdb_${sample_id}/gtdb_identify_${sample_id} \\
    --out_dir gtdb_${sample_id}/gtdb_align_${sample_id}

    gtdbtk classify --genome_dir ${query_file}/consensus.fasta.metabat-bins*/ \\
    --align_dir gtdb_${sample_id}/gtdb_align_${sample_id} \\
    --out_dir gtdb_${sample_id}/gtdb_classify_${sample_id} \\
    --cpus ${task.cpus} \\
    -x fa
    """
}
workflow {
    ch_bins=channel.fromFilePairs("${params.metabatDir}/barcode0*",size:1,type:'dir')
    GTDBtk(ch_bins)
}