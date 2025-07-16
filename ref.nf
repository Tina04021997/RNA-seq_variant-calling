nextflow.enable.dsl=2

params.database_path = "/tscc/projects/ps-lalexandrov/shared"

params.sample = "sample.csv"
params.ref="${params.database_path}/EVC_nextflow/GRCh38_ref/GRCh38.d1.vd1.fa"
params.bam_dir="$projectDir/RESULTS/BAM"
params.report_dir="$projectDir/RESULTS/REPORT"
params.bed="${params.database_path}/EVC_nextflow/GRCh38_ref/hg38_chr.bed.gz"
params.manta_bed="${params.database_path}/EVC_nextflow/GRCh38_ref/manta.bed.gz"
params.database_dir="${params.database_path}/EVC_nextflow/Databases/GRCh38"
params.mkdup_temp_dir="$projectDir/mkdup_tmp"
params.FASTQC="${params.database_path}/EVC_nextflow/FastQC"
params.FASTQC_dir="$projectDir/RESULTS/FASTQC"
params.mkdup_dir="$projectDir/RESULTS/MKDUP"
params.recal_dir="$projectDir/RESULTS/RECALIBRATE"
params.SAGE_java="${params.database_path}/EVC_nextflow/SAGE/sage_v3.3.jar"
params.SAGE_ref_dir="${params.database_path}/EVC_nextflow/SAGE/v5_34/ref/38"
params.SAGE_dir="$projectDir/RESULTS/SAGE"
params.strelka_dir="$projectDir/RESULTS/STRELKA"
params.MuSE2="${params.database_path}/EVC_nextflow/MuSE/MuSE"
params.muse2_dir="$projectDir/RESULTS/MuSE2"
params.MUTECT2_dir="$projectDir/RESULTS/Mutect2"
params.mosdepth_dir="$projectDir/RESULTS/mosdepth"
params.conpair_dir="$projectDir/RESULTS/Conpair"
params.cnvkit_dir="$projectDir/RESULTS/CNVkit"
params.Delly_dir="$projectDir/RESULTS/Delly"
params.ascat_dir="$projectDir/RESULTS/ASCAT"
params.conpair="${params.database_path}/EVC_nextflow/Conpair-0.2"
params.jre="${params.database_path}/EVC_nextflow/jre1.8.0_401"
params.manta_dir="$projectDir/RESULTS/MANTA"
params.postevc_dir="$projectDir/RESULTS/POST"

params.bwamem2_env = "${params.database_path}/EVC_nextflow/yml/bwamem2.yml"
params.mkdup_env = "${params.database_path}/EVC_nextflow/yml/mkdup.yml"
params.conpair_env = "${params.database_path}/EVC_nextflow/yml/conpair.yml"
params.samtools_env = "${params.database_path}/EVC_nextflow/yml/samtools.yml"
params.strelka_env = "${params.database_path}/EVC_nextflow/yml/strelka_env.yml"
params.mosdepth_env = "${params.database_path}/EVC_nextflow/yml/mosdepth_env.yml"
params.summary_env = "${params.database_path}/EVC_nextflow/yml/py_summary.yml"
params.java_env = "${params.database_path}/EVC_nextflow/yml/java.yml"
params.muse2_env = "${params.database_path}/EVC_nextflow/yml/muse2.yml"
params.fastqc_env = "${params.database_path}/EVC_nextflow/yml/fastqc_env.yml"
params.cnvkit_env = "${params.database_path}/EVC_nextflow/yml/cnvkit.yml"
params.delly_env = "${params.database_path}/EVC_nextflow/yml/delly.yml"
params.ascat_env = "${params.database_path}/EVC_nextflow/yml/ascat.yml"
params.alleleCounter_env = "${params.database_path}/EVC_nextflow/yml/alleleCounter.yml"
params.tmp_dir = '${workflow.workDir}/ascat_tmp'
params.manta_env = "${params.database_path}/EVC_nextflow/yml/manta.yml"

include { FASTQC } from './Modules/FASTQC'
include { BWA_MEM } from './Modules/BWA_MEM'
include { CHECK_BAM_BWA } from './Modules/CHECK_BAM_BWA'
include { COMBINE_REPORTS_BWA } from './Modules/COMBINE_REPORTS_BWA'
include { RECALIBRATE_BaseRecal } from './Modules/RECALIBRATE_BaseRecal.nf'
include { RECALIBRATE_BaseRecal_exome } from './Modules/RECALIBRATE_BaseRecal_exome.nf'
include { RECALIBRATE_MergeReport } from './Modules/RECALIBRATE_MergeReport.nf'
include { RECALIBRATE_BQSR } from './Modules/RECALIBRATE_BQSR.nf'
include { RECALIBRATE_BQSR_exome } from './Modules/RECALIBRATE_BQSR_exome.nf'
include { RECALIBRATE_MergeBam } from './Modules/RECALIBRATE_MergeBam.nf'
include { RECALIBRATE_SortBam } from './Modules/RECALIBRATE_SortBam.nf'
include { CHECK_BAM_RECAL } from './Modules/CHECK_BAM_RECAL'
include { COMBINE_REPORTS_RECAL } from './Modules/COMBINE_REPORTS_RECAL'
include { MKDUP } from './Modules/MKDUP'
include { CHECK_BAM_MKDUP } from './Modules/CHECK_BAM_MKDUP'
include { COMBINE_REPORTS_MKDUP } from './Modules/COMBINE_REPORTS_MKDUP'
include { MOSDEPTH } from './Modules/MOSDEPTH'
include { CONPAIR } from './Modules/CONPAIR'
include { SAGE } from './Modules/SAGE'
include { STRELKA } from './Modules/STRELKA'
include { MuSE2 } from './Modules/MuSE2'

include { CNVkit_buildcnn } from './Modules/CNVkit_buildcnn'
include { CNVkit } from './Modules/CNVkit'

include { Create_sample_tsv } from './Modules/Delly/Create_sample_tsv'
include { Delly_SVcalling } from './Modules/Delly/Delly_SVcalling'
include { Delly_Prefiltering } from './Modules/Delly/Delly_Prefiltering'
include { Delly_Filtering } from './Modules/Delly/Delly_Filtering'
include { Delly_Filtering_final } from './Modules/Delly/Delly_Filtering_final'

include { MANTA } from './Modules/MANTA'

include { ASCAT } from './Modules/ASCAT'
include { ASCAT_allelecount } from './Modules/ASCAT_exome/ASCAT_allelecount'
include { ASCAT_logrbaf } from './Modules/ASCAT_exome/ASCAT_logrbaf'
include { ASCAT_exome } from './Modules/ASCAT_exome/ASCAT_exome'

include { MUTECT2_CALLING } from './Modules/Mutect2/MUTECT2_CALLING'
include { MUTECT2_CALLING_exome } from './Modules/Mutect2/MUTECT2_CALLING_exome'
include { GETpileUP } from './Modules/Mutect2/GETpileUP'
include { GETpileUP_exome } from './Modules/Mutect2/GETpileUP_exome'
include { GETpileUP_Merge } from './Modules/Mutect2/GETpileUP_Merge'
include { LearnReadOrientationModel } from './Modules/Mutect2/LearnReadOrientationModel'
include { LearnReadOrientationModel_exome } from './Modules/Mutect2/LearnReadOrientationModel_exome'
include { MergeMutectStats } from './Modules/Mutect2/MergeMutectStats'
include { MergeVcfs } from './Modules/Mutect2/MergeVcfs'
include { CalculateContamination } from './Modules/Mutect2/CalculateContamination'
include { FilterMutectCalls } from './Modules/Mutect2/FilterMutectCalls'

include { SUMMARY } from './Modules/SUMMARY.nf'

include { POSTEVC } from './Modules/POSTEVC.nf'

include { SAVE_CSV_FASTQC } from './Modules/SAVE_CSV/SAVE_CSV_FASTQC'
include { SAVE_CSV_BWA_MEM } from './Modules/SAVE_CSV/SAVE_CSV_BWA_MEM'
include { SAVE_CSV_MKDUP } from './Modules/SAVE_CSV/SAVE_CSV_MKDUP'
include { SAVE_CSV_RECAL } from './Modules/SAVE_CSV/SAVE_CSV_RECAL'
include { SAVE_CSV_SAGE } from './Modules/SAVE_CSV/SAVE_CSV_SAGE'
include { SAVE_CSV_STRELKA } from './Modules/SAVE_CSV/SAVE_CSV_STRELKA'
include { SAVE_CSV_MuSE2 } from './Modules/SAVE_CSV/SAVE_CSV_MuSE2'
include { SAVE_CSV_Mutect2 } from './Modules/SAVE_CSV/SAVE_CSV_Mutect2'
include { SAVE_CSV_CONPAIR } from './Modules/SAVE_CSV/SAVE_CSV_CONPAIR'
include { SAVE_CSV_MOSDEPTH } from './Modules/SAVE_CSV/SAVE_CSV_MOSDEPTH'

def isToolSelected(String tool) {
    if (params.tool) {
        def tools = params.tool.split(',').collect { it.trim() }
        return tools.contains(tool)
    }
    return false
}

workflow {
    chunk = Channel.of(1..20)
    
    // Sanity check for ascat input format
    if (params.tool.toString().contains('ascat')) {
        def sampleFile = file(params.sample)
        def hasSexColumn = false
        
        if (sampleFile.exists()) {
            def headerLine = sampleFile.readLines()[0]
            def headers = headerLine.split(',')
            hasSexColumn = headers.contains("sex")
            
            if (!hasSexColumn) {
                log.error "ERROR: ASCAT tool requires 'sex' column in sample.csv file"
                exit 1
            }
        } else {
            log.error "ERROR: Cannot find input file ${params.sample}"
            exit 1
        }
    }

    // starts from scratch - mapping
    if (params.step in ['mapping']) {
        // Sanity check for sample sheet
        def requiredColumns = ['patient', 'sample', 'status', 'fastq_1', 'fastq_2']
        if (params.tool && isToolSelected('ascat')) {
            requiredColumns.add('sex')
        }
        
        // Check if sample sheet exists
        def sampleFile = file(params.sample)
        if (!sampleFile.exists()) {
            log.error "ERROR: Sample sheet file '${params.sample}' does not exist"
            exit 1
        }
        
        // Check for required columns in the sample sheet
        def headerLine = sampleFile.readLines()[0]
        def headers = headerLine.split(',').collect { it.trim() }
        def missingColumns = requiredColumns.findAll { !headers.contains(it) }
        
        if (missingColumns) {
            def errorMsg = "ERROR: Sample sheet is missing required column(s): ${missingColumns.join(', ')}"
            log.error errorMsg
            exit 1
        }

        if (params.tool && params.tool.contains('ascat')) {
            Channel.fromPath(params.sample)
            | splitCsv(header:true)
            | map { row ->
                meta = row.subMap('patient', 'sample', 'status', 'fastq_1', 'fastq_2', 'sex')
                return meta
            } | set { sample_sheet }
        } else {
            Channel.fromPath(params.sample)
            | splitCsv(header:true)
            | map { row ->
                meta = row.subMap('patient', 'sample', 'status', 'fastq_1', 'fastq_2')
                return meta
            } | set { sample_sheet }
        }
        FASTQC(sample_sheet)
        BWA_MEM(sample_sheet).set { BWA_MEM_out }
        // CHECK_BAM_BWA(BWA_MEM_out.bam).set{ CHECK_BAM_BWA_out }
        // COMBINE_REPORTS_BWA(CHECK_BAM_BWA_out.individual_reports_bwa.collect())
    }

    // starts from markduplicate
    if (params.step in ['mapping', 'markdup']) {
        if (params.step == "markdup"){
            // Sanity check for sample sheet
            def requiredColumns = ['patient', 'sample', 'status', 'bam', 'bai']
            if (params.tool && isToolSelected('ascat')) {
                requiredColumns.add('sex')
            }
            
            // Check if sample sheet exists
            def sampleFile = file(params.sample)
            if (!sampleFile.exists()) {
                log.error "ERROR: Sample sheet file '${params.sample}' does not exist"
                exit 1
            }
            
            // Check for required columns in the sample sheet
            def headerLine = sampleFile.readLines()[0]
            def headers = headerLine.split(',').collect { it.trim() }
            def missingColumns = requiredColumns.findAll { !headers.contains(it) }
            
            if (missingColumns) {
                def errorMsg = "ERROR: Sample sheet is missing required column(s): ${missingColumns.join(', ')}"
                log.error errorMsg
                exit 1
            }

            if (params.tool && params.tool.contains('ascat')) {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map { row ->
                    meta = row.subMap('patient', 'sample', 'status', 'sex')
                    return [meta, file(row.bam)]
                } | set { sample_sheet }
            } else {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map { row ->
                    meta = row.subMap('patient', 'sample', 'status')
                    return [meta, file(row.bam)]
                } | set { sample_sheet }
            }
            MKDUP(sample_sheet).set{ MKDUP_out }      
        } else {
            MKDUP(BWA_MEM_out).set{ MKDUP_out }
            // CHECK_BAM_MKDUP(MKDUP_out.mkdup_bam).set{ CHECK_BAM_MKDUP_out }
            // COMBINE_REPORTS_MKDUP(CHECK_BAM_MKDUP_out.individual_reports_mkdup.collect())
        }
    }

    // starts from recalibration
    if (params.step in ['mapping', 'markdup', 'recalibration']) {
        if (params.step == "recalibration"){
            // Sanity check for sample sheet
            def requiredColumns = ['patient', 'sample', 'status', 'bam', 'bai']
            if (params.tool && isToolSelected('ascat')) {
                requiredColumns.add('sex')
            }
            
            // Check if sample sheet exists
            def sampleFile = file(params.sample)
            if (!sampleFile.exists()) {
                log.error "ERROR: Sample sheet file '${params.sample}' does not exist"
                exit 1
            }
            
            // Check for required columns in the sample sheet
            def headerLine = sampleFile.readLines()[0]
            def headers = headerLine.split(',').collect { it.trim() }
            def missingColumns = requiredColumns.findAll { !headers.contains(it) }
            
            if (missingColumns) {
                def errorMsg = "ERROR: Sample sheet is missing required column(s): ${missingColumns.join(', ')}"
                log.error errorMsg
                exit 1
            }

            if (params.tool && params.tool.contains('ascat')) {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map {row ->
                    meta = row.subMap('patient', 'sample', 'status', 'sex')
                    return [meta.patient, meta, file(row.bam), file(row.bai)]
                } | set {sample_sheet}
            } else {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map { row ->
                    meta = row.subMap('patient', 'sample', 'status')
                    return [meta.patient, meta, file(row.bam), file(row.bai)]
                } | set { sample_sheet }
            }

            if (params.type == "exome") {
                RECALIBRATE_BaseRecal_exome(sample_sheet).set{ RECALIBRATE_BaseRecal_out }
                RECALIBRATE_BQSR_exome(RECALIBRATE_BaseRecal_out.BQSR_input).set { RECALIBRATE_BQSR_out }
                RECALIBRATE_SortBam(RECALIBRATE_BQSR_out.SortBam_input).set{ RECALIBRATE_out }
            } else if (params.type == "genome") {
                RECALIBRATE_BaseRecal(sample_sheet, chunk).set{ RECALIBRATE_BaseRecal_out }

                RECALIBRATE_BaseRecal_out.MergeReport_input.groupTuple(by:[0,1]).set { BaseRecal_out_pair }
                BaseRecal_out_pair.map{
                    [patient:it[0], status:it[1].status, meta:it[1], bam:it[2][0], table:it[3], bai:it[4][0]]
                }.set{ BaseRecal_out_MAP }

                RECALIBRATE_MergeReport(BaseRecal_out_MAP).set{ RECALIBRATE_MergeReport_out }

                RECALIBRATE_BQSR(RECALIBRATE_MergeReport_out.BQSR_input, chunk).set { RECALIBRATE_BQSR_out }

                RECALIBRATE_BQSR_out.MergeBam_input.groupTuple(by:[0,1]).set { BaseRecal_BQSR_out_pair }
                BaseRecal_BQSR_out_pair.map{
                    [patient:it[0], meta:it[1], bam:it[2]]
                }.set{BaseRecal_BQSR_out_MAP}

                RECALIBRATE_MergeBam(BaseRecal_BQSR_out_MAP).set{ RECALIBRATE_MergeBam_out }
                RECALIBRATE_SortBam(RECALIBRATE_MergeBam_out.MergeBam_input).set{ RECALIBRATE_out }
            }
        } else {
            if (params.type == "exome") {
                RECALIBRATE_BaseRecal_exome(MKDUP_out.pair_mutect).set{ RECALIBRATE_BaseRecal_out }
                RECALIBRATE_BQSR_exome(RECALIBRATE_BaseRecal_out.BQSR_input).set { RECALIBRATE_BQSR_out }
                RECALIBRATE_SortBam(RECALIBRATE_BQSR_out.SortBam_input).set{ RECALIBRATE_out }

            } else if (params.type == "genome") {
                RECALIBRATE_BaseRecal(MKDUP_out.pair_mutect, chunk).set{ RECALIBRATE_BaseRecal_out }

                RECALIBRATE_BaseRecal_out.MergeReport_input.groupTuple(by:[0,1]).set { BaseRecal_out_pair }
                BaseRecal_out_pair.map{
                    [patient:it[0], status:it[1].status, meta:it[1], bam:it[2][0], table:it[3], bai:it[4][0]]
                }.set{ BaseRecal_out_MAP }

                RECALIBRATE_MergeReport(BaseRecal_out_MAP).set{ RECALIBRATE_MergeReport_out }

                RECALIBRATE_BQSR(RECALIBRATE_MergeReport_out.BQSR_input, chunk).set { RECALIBRATE_BQSR_out }

                RECALIBRATE_BQSR_out.MergeBam_input.groupTuple(by:[0,1]).set { BaseRecal_BQSR_out_pair }
                BaseRecal_BQSR_out_pair.map{
                    [patient:it[0], meta:it[1], bam:it[2]]
                }.set{BaseRecal_BQSR_out_MAP}

                RECALIBRATE_MergeBam(BaseRecal_BQSR_out_MAP).set{ RECALIBRATE_MergeBam_out }
                RECALIBRATE_SortBam(RECALIBRATE_MergeBam_out.MergeBam_input).set{ RECALIBRATE_out }
            }
        }
        // CHECK_BAM_RECAL(RECALIBRATE_out.pair_recal).set{ CHECK_BAM_RECAL_out }
        // COMBINE_REPORTS_RECAL(CHECK_BAM_RECAL_out.individual_reports_recal.collect())
    }

    // starts from variant calling
    if (params.step in ['mapping', 'markdup', 'recalibration', 'variant_calling']) {
        if (params.step == "variant_calling"){
            // Sanity check for sample sheet
            def requiredColumns = ['patient', 'sample', 'status', 'bam', 'bai']
            if (params.tool && isToolSelected('ascat')) {
                requiredColumns.add('sex')
            }
            
            // Check if sample sheet exists
            def sampleFile = file(params.sample)
            if (!sampleFile.exists()) {
                log.error "ERROR: Sample sheet file '${params.sample}' does not exist"
                exit 1
            }
            
            // Check for required columns in the sample sheet
            def headerLine = sampleFile.readLines()[0]
            def headers = headerLine.split(',').collect { it.trim() }
            def missingColumns = requiredColumns.findAll { !headers.contains(it) }
            
            if (missingColumns) {
                def errorMsg = "ERROR: Sample sheet is missing required column(s): ${missingColumns.join(', ')}"
                log.error errorMsg
                exit 1
            }

            if (params.tool && params.tool.contains('ascat')) {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)  
                | map {row ->
                    meta = row.subMap('patient', 'sample', 'status', 'sex')
                    return [meta.patient, meta, file(row.bam), file(row.bai)]
                } | set {sample_sheet}
            } else {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map { row ->
                    meta = row.subMap('patient', 'sample', 'status')
                    return [meta.patient, meta, file(row.bam), file(row.bai)]
                } | set { sample_sheet }
            }

            if (params.type == "exome") {
                GETpileUP_exome(sample_sheet).set { GETpileUP_out }
            } else if (params.type == "genome") {
                GETpileUP(sample_sheet, chunk).set { GETpileUP_out }
            }

            sample_sheet.filter{it[1].status == 'normal'}.set{normal}
            sample_sheet.filter{it[1].status == 'tumor'}.set{tumor}
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2], tumor:tumor[2], tumor_meta:tumor[1], normal_meta:normal[1]]
            }.set{ RECALIBRATE_out_MAP }

            // for CN analysis
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], gender:normal[1].sex, normal:normal[2],tumor:tumor[2], sample:tumor[1].sample]
            }.set{ RECALIBRATE_out_MAP_CN }

            // Collect all normal bams for CNVkit and Delly
            normal.map { patient, meta, bam, bai -> [bam, bai] }.flatten().collect().set { normal_bams }

            SAGE(RECALIBRATE_out_MAP)
            STRELKA(RECALIBRATE_out_MAP)
            MuSE2(RECALIBRATE_out_MAP)
            
            CONPAIR(RECALIBRATE_out_MAP)
            MOSDEPTH(RECALIBRATE_out_MAP)

            if (params.type == "exome") {
                MUTECT2_CALLING_exome(RECALIBRATE_out_MAP).set { MUTECT2_CALLING_out }
                LearnReadOrientationModel_exome(MUTECT2_CALLING_out.LearnReadOrientationModel_input).set { LearnReadOrientationModel_out }

                GETpileUP_out.CalculateContamination_input.filter{it[1].status== 'normal'}.set{ GETpileUP_normal }
                GETpileUP_out.CalculateContamination_input.filter{it[1].status== 'tumor'}.set{ GETpileUP_tumor }
                GETpileUP_normal.cross(GETpileUP_tumor){it[0]}.map{
                    normal, tumor ->
                    [patient:normal[0], normal:normal[2],tumor:tumor[2], tumor_meta:tumor[1]]
                }.set{ GETpileUP_out_MAP }

                CalculateContamination(GETpileUP_out_MAP).set { CalculateContamination_out }
                

                def vcfChannel = MUTECT2_CALLING_out.MUTECT2_vcf
                    .map { map, vcf -> 
                        ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
                    }
                def contaminationChannel = CalculateContamination_out.MUTECT2_contamination_table
                    .map { map, cont -> 
                        ["${map.patient}_${map.tumor_meta.sample}", [map, cont]]
                    }
                def segmentsChannel = CalculateContamination_out.MUTECT2_segments_table
                    .map { map, seg -> 
                        ["${map.patient}_${map.tumor_meta.sample}", [map, seg]]
                    }
                def orientationChannel = LearnReadOrientationModel_out.MUTECT2_read_orientation
                    .map { map, ori -> 
                        ["${map.patient}_${map.tumor_meta.sample}", [map, ori]]
                    }
                def statsChannel = MUTECT2_CALLING_out.MUTECT2_stats
                    .map { map, stats -> 
                        ["${map.patient}_${map.tumor_meta.sample}", [map, stats]]
                    }

                vcfChannel
                    .join(contaminationChannel)
                    .join(segmentsChannel)
                    .join(orientationChannel)
                    .join(statsChannel)
                    .map { patient_sample, vcfData, contData, segData, oriData, statsData ->
                        [vcfData[0], vcfData[1], contData[1], segData[1], oriData[1], statsData[1]]
                    }
                    .set { filterInput }
                FilterMutectCalls(filterInput).set { FILTER_OUT }
		
            //SAVE_CSV_Mutect2(FilterMutectCalls.out.Mutect2_out,params.report_dir,params.MUTECT2_dir)

            } else if (params.type == "genome") {
                MUTECT2_CALLING(RECALIBRATE_out_MAP, chunk).set { MUTECT2_CALLING_out }
                MUTECT2_CALLING_out.LearnReadOrientationModel_input.groupTuple(by:0).set { LearnReadOrientationModel_input_pair }
                MUTECT2_CALLING_out.MergeMutectStats_input.groupTuple(by:0).set { MergeMutectStats_input_pair }
                MUTECT2_CALLING_out.MergeVcfs_input.groupTuple(by:0).set { MergeVcfs_input_pair }
                LearnReadOrientationModel(LearnReadOrientationModel_input_pair).set { LearnReadOrientationModel_out }
                MergeMutectStats(MergeMutectStats_input_pair).set { MergeMutectStats_out }
                MergeVcfs(MergeVcfs_input_pair).set { MergeVcfs_out }

                GETpileUP_out.GETpileUP_Merge_input.map{patient, meta, table, chunk -> [patient, meta, [table, chunk]]}.groupTuple(by:[0,1],sort: {it[1]}).set { GETpileUP_out_pair }
                GETpileUP_out_pair.map{
                [patient:it[0], meta:it[1], mix:it[2]]
                }.set{ GETpileUP_out_pair_MAP }

                GETpileUP_Merge(GETpileUP_out_pair_MAP).set { GETpileUP_Merge_out }

                GETpileUP_Merge_out.CalculateContamination_input.filter{it[1].status== 'normal'}.set{ GETpileUP_normal }
                GETpileUP_Merge_out.CalculateContamination_input.filter{it[1].status== 'tumor'}.set{ GETpileUP_tumor }
                GETpileUP_normal.cross(GETpileUP_tumor){it[0]}.map{
                    normal, tumor ->
                    [patient:normal[0], normal:normal[2],tumor:tumor[2], tumor_meta:tumor[1]]
                }.set{ GETpileUP_out_MAP }

                CalculateContamination(GETpileUP_out_MAP).set { CalculateContamination_out }
                def vcfChannel = MergeVcfs_out.MUTECT2_vcf
                    .map{ map, vcf ->
                        ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
                    }
                def contaminationChannel = CalculateContamination_out.MUTECT2_contamination_table
                    .map{ map, table ->
                        ["${map.patient}_${map.tumor_meta.sample}", [map, table]]
                    }
                def segmentsChannel = CalculateContamination_out.MUTECT2_segments_table
                    .map{ map, segments ->
                        ["${map.patient}_${map.tumor_meta.sample}", [map, segments]]
                    }
                def orientationChannel = LearnReadOrientationModel_out.MUTECT2_read_orientation
                    .map{ map, orientation ->
                        ["${map.patient}_${map.tumor_meta.sample}", [map, orientation]]
                    }
                def statsChannel = MergeMutectStats_out.MUTECT2_stats
                    .map{ map, stats ->
                        ["${map.patient}_${map.tumor_meta.sample}", [map, stats]]
                    }

                // Join by the composite key (patient_sample)
                vcfChannel
                    .join(contaminationChannel)
                    .join(segmentsChannel)
                    .join(orientationChannel)
                    .join(statsChannel)
                    .map{ patient_sample, vcfData, contData, segData, oriData, statsData ->
                        [vcfData[0], vcfData[1], contData[1], segData[1], oriData[1], statsData[1]]
                    }
                    .set{ filterInput }

                FilterMutectCalls(filterInput).set{ FILTER_OUT }
            // SAVE_CSV_Mutect2(FilterMutectCalls.out.Mutect2_out,params.report_dir,params.MUTECT2_dir)
            }

            // TOOL-SPECIFIC BLOCKS - MOVED INSIDE variant_calling STEP
            // ASCAT
            if (params.tool && params.tool.contains('ascat')) {
                chromosomes = Channel.of( *(1..22).collect { it.toString() } + ['X'] )
                if (params.type == "exome") {
                    ASCAT_allelecount(RECALIBRATE_out_MAP_CN, chromosomes)

                    ASCAT_allelecount.out.allelecount
                        .map { map, chr, file ->
                            def patient = map.patient
                            def sample = map.sample
                            def gender = map.gender
                            [patient, sample, gender, file]
                        }
                        .groupTuple(by: [0, 1, 2])  // Group by patient, sample, and gender
                        // Restructure from pairs to separate normal and tumor lists
                        .map { patient, sample, gender, files ->
                            // Extract normal and tumor files from the nested structure
                            def normalFiles = []
                            def tumorFiles = []

                            // Flatten any nested structure if present
                            def flatFiles = files.flatten()

                            // Separate normal and tumor files
                            flatFiles.each { file ->
                                if (file.toString().contains("_normal_")) {
                                    normalFiles << file
                                } else if (file.toString().contains("_tumor_")) {
                                    tumorFiles << file
                                }
                            }

                            [patient, sample, gender, normalFiles, tumorFiles]
                        }
                        .set { ascat_logrbaf_input }

                    ASCAT_logrbaf(ascat_logrbaf_input)
                    ASCAT_exome(ASCAT_logrbaf.out.ascat_input)

                } else if (params.type == "genome") {
                    ASCAT(RECALIBRATE_out_MAP_CN)
                }
            }

            // CNVkit
            if (params.tool && params.tool.contains('cnvkit')) {
                CNVkit_buildcnn(normal_bams).set { CNVkit_buildcnn_out }
                CNVkit(RECALIBRATE_out_MAP_CN, CNVkit_buildcnn_out.CNVkit_ref_cnn).set { CNVkit_out }
            }

            // Delly
            if (params.tool && params.tool.contains('delly')) {
                if (params.type == 'exome') {
                    log.error "ERROR: Delly is not recommended for exome data. Please use a different SV caller (manta) for exome data"
                    exit 1  
                } else {
                    Create_sample_tsv(Channel.fromPath(params.sample, checkIfExists: true))
                    sample_tsv_collected = Create_sample_tsv.out.samples_tsv.collect()

                    Delly_SVcalling(RECALIBRATE_out_MAP_CN).set { SVcalling_out }
                    SVcalling_out.SVcalling_bcf
                    .combine(sample_tsv_collected)
                    .set { prefiltering_input }

                    Delly_Prefiltering(prefiltering_input).set { Prefiltering_out }
                    Delly_Filtering(Prefiltering_out.Delly_pre_bcf, normal_bams).set { Filtering_out }

                    Filtering_out.Delly_geno_bcf
                    .combine(sample_tsv_collected)
                    .set { filtering_final_input }

                    Delly_Filtering_final(filtering_final_input)
                }
            }

            // Manta
            if (params.tool && params.tool.contains('manta')) {
                MANTA(RECALIBRATE_out_MAP_CN)
            }

        } else {

            RECALIBRATE_out.pair_recal.filter{it[1].status == 'normal'}.set{normal}
            RECALIBRATE_out.pair_recal.filter{it[1].status == 'tumor'}.set{tumor}
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2], tumor:tumor[2], tumor_meta:tumor[1], normal_meta:normal[1]]
            }.set{ RECALIBRATE_out_MAP }

            // for CN analysis
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], gender:normal[1].sex, normal:normal[2],tumor:tumor[2], sample:tumor[1].sample]
            }.set{ RECALIBRATE_out_MAP_CN }

            // Collect all normal bams for CNVkit and Delly
            normal.map { patient, meta, bam, bai -> [bam, bai] }.flatten().collect().set { normal_bams }

            SAGE(RECALIBRATE_out_MAP)
            STRELKA(RECALIBRATE_out_MAP)
            MuSE2(RECALIBRATE_out_MAP)
            
            CONPAIR(RECALIBRATE_out_MAP)
            MOSDEPTH(RECALIBRATE_out_MAP)

            if (params.type == "exome") {
                GETpileUP_exome(RECALIBRATE_out.pair_recal).set { GETpileUP_out }
                MUTECT2_CALLING_exome(RECALIBRATE_out_MAP).set { MUTECT2_CALLING_out }
                LearnReadOrientationModel_exome(MUTECT2_CALLING_out.LearnReadOrientationModel_input).set { LearnReadOrientationModel_out }

                GETpileUP_out.CalculateContamination_input.filter{it[1].status== 'normal'}.set{ GETpileUP_normal }
                GETpileUP_out.CalculateContamination_input.filter{it[1].status== 'tumor'}.set{ GETpileUP_tumor }
                GETpileUP_normal.cross(GETpileUP_tumor){it[0]}.map{
                    normal, tumor ->
                    [patient:normal[0], normal:normal[2],tumor:tumor[2], tumor_meta:tumor[1]]
                }.set{ GETpileUP_out_MAP }

                CalculateContamination(GETpileUP_out_MAP).set { CalculateContamination_out }
                

                def vcfChannel = MUTECT2_CALLING_out.MUTECT2_vcf
                    .map { map, vcf -> 
                        ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
                    }
                def contaminationChannel = CalculateContamination_out.MUTECT2_contamination_table
                    .map { map, cont -> 
                        ["${map.patient}_${map.tumor_meta.sample}", [map, cont]]
                    }
                def segmentsChannel = CalculateContamination_out.MUTECT2_segments_table
                    .map { map, seg -> 
                        ["${map.patient}_${map.tumor_meta.sample}", [map, seg]]
                    }
                def orientationChannel = LearnReadOrientationModel_out.MUTECT2_read_orientation
                    .map { map, ori -> 
                        ["${map.patient}_${map.tumor_meta.sample}", [map, ori]]
                    }
                def statsChannel = MUTECT2_CALLING_out.MUTECT2_stats
                    .map { map, stats -> 
                        ["${map.patient}_${map.tumor_meta.sample}", [map, stats]]
                    }

                vcfChannel
                    .join(contaminationChannel)
                    .join(segmentsChannel)
                    .join(orientationChannel)
                    .join(statsChannel)
                    .map { patient_sample, vcfData, contData, segData, oriData, statsData ->
                        [vcfData[0], vcfData[1], contData[1], segData[1], oriData[1], statsData[1]]
                    }
                    .set { filterInput }
                FilterMutectCalls(filterInput).set { FILTER_OUT }
		
            //SAVE_CSV_Mutect2(FilterMutectCalls.out.Mutect2_out,params.report_dir,params.MUTECT2_dir)

            } else if (params.type == "genome") {
                GETpileUP(RECALIBRATE_out.pair_recal, chunk).set { GETpileUP_out }
                MUTECT2_CALLING(RECALIBRATE_out_MAP, chunk).set { MUTECT2_CALLING_out }
                MUTECT2_CALLING_out.LearnReadOrientationModel_input.groupTuple(by:0).set { LearnReadOrientationModel_input_pair }
                MUTECT2_CALLING_out.MergeMutectStats_input.groupTuple(by:0).set { MergeMutectStats_input_pair }
                MUTECT2_CALLING_out.MergeVcfs_input.groupTuple(by:0).set { MergeVcfs_input_pair }
                LearnReadOrientationModel(LearnReadOrientationModel_input_pair).set { LearnReadOrientationModel_out }
                MergeMutectStats(MergeMutectStats_input_pair).set { MergeMutectStats_out }
                MergeVcfs(MergeVcfs_input_pair).set { MergeVcfs_out }

                GETpileUP_out.GETpileUP_Merge_input.map{patient, meta, table, chunk -> [patient, meta, [table, chunk]]}.groupTuple(by:[0,1],sort: {it[1]}).set { GETpileUP_out_pair }
                GETpileUP_out_pair.map{
                [patient:it[0], meta:it[1], mix:it[2]]
                }.set{ GETpileUP_out_pair_MAP }

                GETpileUP_Merge(GETpileUP_out_pair_MAP).set { GETpileUP_Merge_out }

                GETpileUP_Merge_out.CalculateContamination_input.filter{it[1].status== 'normal'}.set{ GETpileUP_normal }
                GETpileUP_Merge_out.CalculateContamination_input.filter{it[1].status== 'tumor'}.set{ GETpileUP_tumor }
                GETpileUP_normal.cross(GETpileUP_tumor){it[0]}.map{
                    normal, tumor ->
                    [patient:normal[0], normal:normal[2],tumor:tumor[2], tumor_meta:tumor[1]]
                }.set{ GETpileUP_out_MAP }

                CalculateContamination(GETpileUP_out_MAP).set { CalculateContamination_out }
                def vcfChannel = MergeVcfs_out.MUTECT2_vcf
                    .map{ map, vcf ->
                        ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
                    }
                def contaminationChannel = CalculateContamination_out.MUTECT2_contamination_table
                    .map{ map, table ->
                        ["${map.patient}_${map.tumor_meta.sample}", [map, table]]
                    }
                def segmentsChannel = CalculateContamination_out.MUTECT2_segments_table
                    .map{ map, segments ->
                        ["${map.patient}_${map.tumor_meta.sample}", [map, segments]]
                    }
                def orientationChannel = LearnReadOrientationModel_out.MUTECT2_read_orientation
                    .map{ map, orientation ->
                        ["${map.patient}_${map.tumor_meta.sample}", [map, orientation]]
                    }
                def statsChannel = MergeMutectStats_out.MUTECT2_stats
                    .map{ map, stats ->
                        ["${map.patient}_${map.tumor_meta.sample}", [map, stats]]
                    }

                // Join by the composite key (patient_sample)
                vcfChannel
                    .join(contaminationChannel)
                    .join(segmentsChannel)
                    .join(orientationChannel)
                    .join(statsChannel)
                    .map{ patient_sample, vcfData, contData, segData, oriData, statsData ->
                        [vcfData[0], vcfData[1], contData[1], segData[1], oriData[1], statsData[1]]
                    }
                    .set{ filterInput }

                FilterMutectCalls(filterInput).set{ FILTER_OUT }
            // SAVE_CSV_Mutect2(FilterMutectCalls.out.Mutect2_out,params.report_dir,params.MUTECT2_dir)
            }

            // TOOL-SPECIFIC BLOCKS - ALSO FOR NON variant_calling START STEPS
	    // ASCAT
	    if (params.tool in ['ascat']) {
	        chromosomes = Channel.of( *(1..22).collect { it.toString() } + ['X'] )
	        if (params.type == "exome") {
	            ASCAT_allelecount(RECALIBRATE_out_MAP_CN, chromosomes)
	
	            ASCAT_allelecount.out.allelecount
	                .map { map, chr, file ->
	                    def patient = map.patient
	                    def sample = map.sample
	                    def gender = map.gender
	                    [patient, sample, gender, file]
	                }
	                .groupTuple(by: [0, 1, 2])  // Group by patient, sample, and gender
	                // Restructure from pairs to separate normal and tumor lists
	                .map { patient, sample, gender, files ->
	                    // Extract normal and tumor files from the nested structure
	                    def normalFiles = []
	                    def tumorFiles = []
	
	                    // Flatten any nested structure if present
	                    def flatFiles = files.flatten()
	
	                    // Separate normal and tumor files
	                    flatFiles.each { file ->
	                        if (file.toString().contains("_normal_")) {
	                            normalFiles << file
	                        } else if (file.toString().contains("_tumor_")) {
	                            tumorFiles << file
	                        }
	                    }
	
	                    [patient, sample, gender, normalFiles, tumorFiles]
	                }
	                .set { ascat_logrbaf_input }
	
	            ASCAT_logrbaf(ascat_logrbaf_input)
	            ASCAT_exome(ASCAT_logrbaf.out.ascat_input)
	
	        } else if (params.type == "genome") {
	            ASCAT(RECALIBRATE_out_MAP_CN)
	        }
	    }
	
	    // CNVkit
	    if (params.tool in ['cnvkit']) {
	        CNVkit_buildcnn(normal_bams).set { CNVkit_buildcnn_out }
	        CNVkit(RECALIBRATE_out_MAP_CN, CNVkit_buildcnn_out.CNVkit_ref_cnn).set { CNVkit_out }
	    }
	
	    // Delly
	    if (params.tool in ['delly']) {
	        if (params.type == 'exome' && params.tool.toString().contains('delly')) {
	            log.error "ERROR: Delly is not recommended for exome data. Please use a different SV caller (manta) for exome data"
	            exit 1  
	        } else {
	            Create_sample_tsv(Channel.fromPath(params.sample, checkIfExists: true))
	            sample_tsv_collected = Create_sample_tsv.out.samples_tsv.collect()
	
	            Delly_SVcalling(RECALIBRATE_out_MAP_CN).set { SVcalling_out }
	            SVcalling_out.SVcalling_bcf
	            .combine(sample_tsv_collected)
	            .set { prefiltering_input }
	
	            Delly_Prefiltering(prefiltering_input).set { Prefiltering_out }
	            Delly_Filtering(Prefiltering_out.Delly_pre_bcf, normal_bams).set { Filtering_out }
	
	            Filtering_out.Delly_geno_bcf
	            .combine(sample_tsv_collected)
	            .set { filtering_final_input }
	
	            Delly_Filtering_final(filtering_final_input)
	        }
	    }
	
	    // Manta
	    if (params.tool in ['manta']) {
	        MANTA(RECALIBRATE_out_MAP_CN)
	    }
        }
    }

    // POST - filtering and consensus calling
    if (params.step in ['mapping', 'markdup', 'recalibration', 'variant_calling']) {
        
        // Use the outputs that are guaranteed to exist from variant calling steps
        def mutect2Channel = FilterMutectCalls.out.Mutect2_out
            .map { map, vcf ->
                ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
            }
        
        def muse2Channel = MuSE2.out.MuSE2_out
            .map { map, vcf ->
                ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
            }
        
        def strelkaChannel = STRELKA.out.STRELKA_out
            .map { map, snv, indel ->
                ["${map.patient}_${map.tumor_meta.sample}", [map, snv, indel]]
            }
        
        def sageChannel = SAGE.out.SAGE_out
            .map { map, vcf ->
                ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
            }
        
        // Recreate recalChannel from the available data
        // This needs to be derived from the same source that created RECALIBRATE_out_MAP
        def recalChannel
        
        if (params.step == "variant_calling") {
            sample_sheet.filter{it[1].status == 'normal'}.set{normal_post}
            sample_sheet.filter{it[1].status == 'tumor'}.set{tumor_post}
            normal_post.cross(tumor_post){it[0]}.map{
                normal, tumor ->
                ["${tumor[1].patient}_${tumor[1].sample}", [tumor[1], tumor[2]]]
            }.set { recalChannel }
        } else {
            RECALIBRATE_out.pair_recal.filter{it[1].status == 'normal'}.set{normal_post}
            RECALIBRATE_out.pair_recal.filter{it[1].status == 'tumor'}.set{tumor_post}
            normal_post.cross(tumor_post){it[0]}.map{
                normal, tumor ->
                ["${tumor[1].patient}_${tumor[1].sample}", [tumor[1], tumor[2]]]
            }.set { recalChannel }
        }
        
        mutect2Channel
            .join(muse2Channel)
            .join(strelkaChannel)
            .join(sageChannel)
            .join(recalChannel)
            .map { patient_sample, mutect2_vcf, muse2_vcf, strelka_vcf, sage_vcf, recal_bam ->
                [patient_sample, mutect2_vcf[1], muse2_vcf[1], strelka_vcf[1], strelka_vcf[2], sage_vcf[1], recal_bam[1]]
            }
            .set { postevcInput }
        
        POSTEVC(postevcInput)
    }
}
