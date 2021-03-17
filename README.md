# RNA-seq_variant-calling
This is the workflow for RNA-seq germline variant calling based on [GATK RNAseq short variant discovery workflows](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) 
## Pipeline workflow
![image](https://github.com/Tina04021997/RNA-seq_variant-calling/blob/main/RNA-seq%20variant%20calling%20workflow.jpg)
## Environments
- STAR v2.7.8a
- Picard v2.23.4
- GATK v4.1.9.0
- VEP v101
- VAtools v4.1.0
## Input data
- Paired-end fastq files
- Reference: mm10 resource bundle see [Create GATK mm10 resource bundle](https://github.com/igordot/genomics/blob/master/workflows/gatk-mouse-mm10.md)

## Output data
- Tsv files (SNP/INDEL) transformed from annotated vcf files
