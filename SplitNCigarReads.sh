#!/bin/bash
# Author:Tina Yang
# Date: Mar 2,2021
# GATK v4.1.9.0

# SplitNCigarReads, reformat RNA alignments that span introns for HaplotypeCaller
# Splits reads with N in the cigar into multiple supplementary alignments and hard clips mismatching overhangs
  
#Dataset
REF_PATH="/LVM_data/tina/0220_variant_calling/mm10"
DATA_PATH="/LVM_data/tina/0220_variant_calling"

#remember to create fa.fai & fa.dict files first via $samtools faidx ref.fa & gatk CreateSequenceDictionary -R ref.fa
gatk SplitNCigarReads \
-R $REF_PATH/GRCm38.primary_assembly.genome.fa \
-I $DATA_PATH/STAR_PASS2_resultAligned.sortedByCoord.out.mark.bam \
-O $DATA_PATH/STAR_PASS2_resultAligned.sortedByCoord.out.mark.SPLIT.bam
