params.fastqDir = "$launchDir/data"
params.kraken2DB = "$launchDir/DBTITICACA" 

process ConcatFastq {
    conda '/home/tao/anaconda3/envs/kraken2'
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

process Kraken2 {
    conda '/home/tao/anaconda3/envs/kraken2'
    publishDir 'results/kraken2', mode: 'copy' 
    cpus 16
    input:
    tuple val(sample_id), path(query_file)

    output:
    path "${sample_id}_kraken2_report.tsv"
    path "${sample_id}_kraken2_classification.tsv"
    
    script:
    """
    kraken2 --db $params.kraken2DB \\
    --threads ${task.cpus} \\
    --report ${sample_id}_kraken2_report.tsv \\
    --output ${sample_id}_kraken2_classification.tsv \\
    $query_file
    """
}

workflow {
    channel_fastq = channel.fromFilePairs("$params.fastqDir/*", type:"dir", size:1)
    channel_fastq.view()
    //ConcatFastq(channel_fastq)
    //Kraken2(ConcatFastq.out)
}

