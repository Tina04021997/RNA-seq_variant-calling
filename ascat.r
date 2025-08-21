## PrepareHTS and Running ASCAT for all samples ##
R version 4.2.2
ASCAT 3.2.0

Add rho_manual = xxx, psi_manual = xxx when running ascat.runAscat (e.g., ascat.runAscat(ascat.bc, img.prefix = "first",gamma=GAMMA, rho_manual = 0.95, psi_manual = 1.9))
Where rho is purity and psi is ploidy

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
gender <- args[2]

library(ASCAT)

# # PrepareHTS: Extracting logR and BAF from HTS data (bam files)
ascat.prepareHTS(
  tumourseqfile = paste0("PATH//",sample,"_tumor.bam"),
  normalseqfile = paste0("/PATH/",sample,"_normal.bam"),
  tumourname = paste0(sample,"_tumor"),
  normalname = paste0(sample,"_normal"),
  allelecounter_exe = "/PATH/bin/alleleCounter",
  alleles.prefix = "/PATH/hg38/Alleles/G1000_alleles_hg38_chr",
  loci.prefix = "/PATH/hg38/Loci/G1000_loci_hg38_chr",
  gender = paste0(gender),
  genomeVersion = "hg38")

# Running ASCAT
# For HTS data (WGS, WES and targeted sequencing), gamma must be set to 1 in ascat.runASCAT
ascat.bc = ascat.loadData(Tumor_LogR_file = paste0(sample,"_tumor_tumourLogR.txt"), 
                          Tumor_BAF_file = paste0(sample,"_tumor_tumourBAF.txt"), 
                          Germline_LogR_file = paste0(sample,"_tumor_normalLogR.txt"), 
                          Germline_BAF_file = paste0(sample,"_tumor_normalBAF.txt"), 
                          gender = paste0(gender), 
                          genomeVersion = "hg38")

ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "/PATH/hg38/GC_Correction/GC_G1000_hg38.txt", replictimingfile = "/PATH/hg38/RT_Correction/RT_G1000_hg38.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc, out.dir = NA) # penalty=25 for targeted sequencing data
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma = 1, write_segments = T) # gamma=1 for HTS data
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc,ascat.output, QC, file= 'ASCAT_objects.Rdata')
