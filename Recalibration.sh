#!/bin/bash
# Author:Tina Yang
# Date: Mar 5,2021
# Base quality recalibration
# GATK v4.1.9.0

#Dataset
REF_PATH="/LVM_data/tina/0220_variant_calling/mm10"
DATA_PATH="/LVM_data/tina/0220_variant_calling"
INDEL_PATH="/LVM_data/tina/0220_variant_calling/mm10_INDEL"
SNP_PATH="/LVM_data/tina/0220_variant_calling/mm10_SNP"
 
 
# Steps:
  #1. Generates recalibration table for Base Quality Score Recalibration (BQSR): BaseRecalibrator
  #2. Apply base quality score recalibration: ApplyBQSR
  #3. Evaluate and compare base quality score recalibration (BQSR) tables: AnalyzeCovariates
 
# Create index
tabix $SNP_PATH/dbSNP.vcf.gz
gatk IndexFeatureFile -I $INDEL_PATH/indels.vcf
samtools index /LVM_data/tina/0220_variant_calling/BAMforRECAL.bam
 
#Add Read Group information to bam file so that GATK can work
picard AddOrReplaceReadGroups \
  I=$DATA_PATH/STAR_PASS2_resultAligned.sortedByCoord.out.mark.SPLIT.bam \
  O=$DATA_PATH/BAMforRECAL.bam \
  RGID=4 \
  RGLB=lib1 \
  RGPL=ILLUMINA \
  RGPU=unit1 \
  RGSM=MTCQ1


# Step1 Base Recalibrate
gatk BaseRecalibrator \
     -I $DATA_PATH/BAMforRECAL.bam  \
     -R $REF_PATH/GRCm38.primary_assembly.genome.fa \
     --known-sites $SNP_PATH/dbSNP.vcf.gz \
     --known-sites $INDEL_PATH/indels.vcf \
     -L $REF_PATH/interval.list \
     -O $DATA_PATH/recal.table

# Step2 Recalibrates the base qualities of the input reads by base quality score recalibration
gatk ApplyBQSR \
     -R $REF_PATH/GRCm38.primary_assembly.genome.fa \
     -I $DATA_PATH/BAMforRECAL.bam \
     --bqsr-recal-file $DATA_PATH/recal.table \
     -L $REF_PATH/interval.list  \
     -O $DATA_PATH/RECAL.bam

# Step3 Generates plots to assess the quality of a recalibration run
# Get analysis-ready RNAseq Reads
gatk AnalyzeCovariates \
     -bqsr $DATA_PATH/recal.table \
     -plots $DATA_PATH/AnalyzeCovariates.pdf
