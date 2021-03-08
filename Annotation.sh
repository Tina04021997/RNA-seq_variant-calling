#!/bin/bash
# Author:Tina Yang
# Date: Mar 6,2021
# By using Ensembl VEP to determine the effect of the variants
# According to Ensembl, the VEP Cache version should be the same as your Ensembl VEP installation
# VEP v101
  
#Dataset
REF_PATH="/LVM_data/tina/0220_variant_calling/mm10"
DATA_PATH="/LVM_data/tina/0220_variant_calling"

#Running VEP
vep -i $DATA_PATH/FilteredVariant.vcf.gz \
    --cache \
    --cache_version 101 \
    --dir_cache $DATA_PATH/.vep \
    --species mus_musculus \
    --everything \
    --fork 8 \
    --offline \
    --fasta  $REF_PATH/GRCm38.primary_assembly.genome.fa \
    -o $DATA_PATH/VEPannotated.vcf
