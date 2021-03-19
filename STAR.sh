 #!/bin/bash
 # Author:Tina Yang
 # Date: Mar 2,2021
 # STAR 2-pass mapping script for PE
 # STAR v2.7.8a
 
 # Reference: GENCODE mm10
 # Steps (based on Rajesh Detroja's lines on rna-star@googlegroups.com):
     #1. Indexing genome with annotations
     #2. 1-pass mapping with indexed genome
     #3. Indexing genome with annotations and SJ.out.tab files
     #4. 2-pass mapping with new indexed genome with annotations and SJ.out.tab files
 
 #Dataset
 REF_PATH="/LVM_data/tina/0220_variant_calling/mm10"
 ANNOTATION_PATH="/LVM_data/tina/0220_variant_calling/mm10"
 DATA_PATH="/LVM_data/tina/0220_variant_calling"
 
 
 # Step1 Indexing
 STAR \
 --runMode genomeGenerate \
 --genomeDir $REF_PATH \
 --genomeFastaFiles $REF_PATH/GRCm38.primary_assembly.genome.fa \
 --sjdbGTFfile $ANNOTATION_PATH/gencode.vM22.annotation.gtf \
 --runThreadN 20
 
 # Step2 Mapping-1-pass
 STAR \
 --genomeDir $REF_PATH \
 --readFilesIn $DATA_PATH/MTCQ_HLGV2DRXX_L1_R1.clean.fastq $DATA_PATH/MTCQ_HLGV2DRXX_L1_R2.clean.fastq \
 --outFileNamePrefix STAR_PASS1_result \
 --runThreadN 20
 
 # Step3 Indexing genome with annotations and SJ.out.tab file
 STAR \
 --runMode genomeGenerate \
 --genomeDir $REF_PATH/SJ_index \
 --genomeFastaFiles $REF_PATH/GRCm38.primary_assembly.genome.fa \
 --sjdbGTFfile $ANNOTATION_PATH/gencode.vM22.annotation.gtf \
 --runThreadN 20 \
 --sjdbFileChrStartEnd $DATA_PATH/STAR_PASS1_resultSJ.out.tab
 
 
 # Step 4 2-pass mapping output alignments in BAM format and sort BAM by coordinate
 STAR \
 --genomeDir $REF_PATH/SJ_index \
 --readFilesIn $DATA_PATH/MTCQ_HLGV2DRXX_L1_R1.clean.fastq $DATA_PATH/MTCQ_HLGV2DRXX_L1_R2.clean.fastq \
 --outFileNamePrefix STAR_PASS2_result \
 --runThreadN 20 \
 --outSAMtype BAM SortedByCoordinate
