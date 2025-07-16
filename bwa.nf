nextflow.enable.dsl=2

process BWA_MEM {
    conda "${params.bwamem2_env}"
    scratch true
    label 'BWA_MEM'
    if (params.publish.tokenize(',').contains('raw_bam')) {
        publishDir("${params.bam_dir}", mode: 'copy')
    }
    maxRetries 3

    input:
    val(meta)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    """
    bwa-mem2 mem -T 0 -t 16 -R "@RG\\tID:${meta.patient}\\tSM:${meta.patient}_${meta.sample}_${meta.status}\\tPL:ILLUMINA" ${params.ref} ${meta.fastq_1} ${meta.fastq_2} | samtools view -bh --input-fmt-option nthreads=16 - | samtools sort -@ 16 -T ${meta.patient}_${meta.sample}_${meta.status}  -o ${meta.patient}_${meta.sample}_${meta.status}_raw.bam -
    """
}
