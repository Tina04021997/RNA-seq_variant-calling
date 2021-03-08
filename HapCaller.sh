#!/bin/bash
# Author:Tina Yang
# Date: Mar 5,2021
# Single sample variant calling by GATK HaplotypeCaller
# GATK v4.1.9.0
  
#Dataset
REF_PATH="/LVM_data/tina/0220_variant_calling/mm10"
DATA_PATH="/LVM_data/tina/0220_variant_calling"
SNP_PATH="/LVM_data/tina/0220_variant_calling/mm10_SNP"
 
 
# Variant calling for germline SNPs and indels with bamout to show realigned reads
#If need to specify specific intervals to call variant, use --intervals
 
gatk --java-options "-Xmx4g" HaplotypeCaller \
     --dbsnp $SNP_PATH/dbSNP.vcf.gz \
     -R $REF_PATH/GRCm38.primary_assembly.genome.fa \
     -I $DATA_PATH/RECAL.bam \
     -O $DATA_PATH/RawVariant.vcf.gz \
     -L $REF_PATH/interval.list \
     -bamout $DATA_PATH/HapBamOut.bam \
     --native-pair-hmm-threads 20
