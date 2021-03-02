#!/bin/bash                                                                                     
# Author:Tina Yang                                                                                
# Date: Mar 2,2021                                                                                
# Picrad mark duplicates                                                                          
# Picard v2.23.4                                                                                 

picard MarkDuplicates \                                                                            
    I=/LVM_data/tina/0220_variant_calling/STAR_PASS2_resultAligned.sortedByCoord.out.bam \        
    O=/LVM_data/tina/0220_variant_calling/STAR_PASS2_resultAligned.sortedByCoord.out.mark.bam \        
    M=/LVM_data/tina/0220_variant_calling/STAR_PASS2_resultAligned.sortedByCoord.out.mark.txt 
