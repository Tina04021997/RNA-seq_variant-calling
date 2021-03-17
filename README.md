# RNA-seq_variant-calling
This is the workflow for RNA-seq germline variant calling based on [GATK RNAseq short variant discovery workflows](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) and [VAtools](https://vatools.readthedocs.io/en/latest/) for data cleaning.
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
## Output data
-  tsv file transformed from annotated vcf files

## Reference data 
**mm10 resource bundle** see [Create GATK mm10 resource bundle](https://github.com/igordot/genomics/blob/master/workflows/gatk-mouse-mm10.md)

- If you encounter problems while concatenating dbSNP VCF files, try this:
```
bgzip -c vcf_chr_number.vcf > vcf_chr_number.vcf.gz
tabix  vcf_chr_number.vcf.gz
bcftools concat vcf_chr_number.vcf.gz vcf_chr_number.vcf.gz -Oz -o dbSNP.vcf.gz 
```
- If you encounter problems while sorting MGP indels VCF file, try this:
```
grep "^#" mgp.v5.indels.pass.chr.vcf > indels.vcf && grep -v "^#" mgp.v5.indels.pass.chr.vcf | \sort -V -k1,1 -k2,2n >> indels.vcf
```
**GENCODE mm10 fasta (PRI) & GTF files**
```
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/GRCm38.primary_assembly.genome.fa.gz
gunzip GRCm38.primary_assembly.genome.fa.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.annotation.gtf.gz
gunzip gencode.vM22.annotation.gtf.gz
```
