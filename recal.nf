nextflow.enable.dsl=2

process RECALIBRATE_BaseRecal {
    conda "${params.java_env}"
    scratch true
    label 'RECALIBRATE'
    errorStrategy = 'retry'
    maxRetries 3

    input:
    tuple val(patient), val(meta), path(bam), path(bai)
    each chunk

    output:
    tuple val(patient), val(meta), path(bam), path("*table"), path(bai), emit: MergeReport_input

    script:
    """
    ${params.database_path}/EVC_nextflow/gatk-4.6.0.0/gatk BaseRecalibrator \
    -I ${bam} \
    -R ${params.ref} \
    --known-sites ${params.database_dir}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites ${params.database_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -L ${params.database_dir}/interval_list_20/${chunk}-scattered.interval_list \
    -O ${meta.patient}_${meta.sample}_${meta.status}_${chunk}_recal.table
    """

}
