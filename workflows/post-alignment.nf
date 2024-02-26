/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowConfigureinput.initialise(params, log)

// Initialize sample id, mapped bam and bam index channels
ch_input = Channel.fromSamplesheet("input")
            .map{ meta, mapped, index ->

            if (meta.filetype != mapped.getExtension().toString()) {
                error('The file extension does not fit the specified file_type.\n' + mapped.toString() )
            }

            meta.index  = index ? true : false

            return [meta, mapped, index]

            }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOWS
//

include { PREPARE_INDICES                             } from '../subworkflows/local/prepare_indices'
include { PREPARE_GENOME                              } from '../subworkflows/local/prepare_genome/main'
include { BAM_MARKDUPLICATES                          } from '../subworkflows/local/bam_markduplicates/main'
include { BAM_BASERECALIBRATOR                        } from '../subworkflows/local/bam_baserecalibrator/main'
include { BAM_APPLYBQSR                               } from '../subworkflows/local/bam_applybqsr/main'
include { QC_BAM                                      } from '../subworkflows/local/qc_bam'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow POSTALIGNMENT {

    ch_versions = Channel.empty()

    // Initialize file channels for PREPARE_GENOME subworkflow
    ch_genome_fasta             = Channel.fromPath(params.fasta).map { it -> [[id:it[0].simpleName], it] }.collect()
    ch_genome_fai               = params.fai            ? Channel.fromPath(params.fai).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                        : Channel.empty()
    known_indels                = params.known_indels   ? Channel.fromPath(params.known_indels).collect()            : Channel.value([])
    dbsnp                       = params.dbsnp          ? Channel.fromPath(params.dbsnp).collect()                   : Channel.value([])
    ch_intervals_wgs            = params.intervals_wgs  ? Channel.fromPath(params.intervals_wgs).collect()
                                                        : Channel.empty()
    ch_intervals_y              = params.intervals_y    ? Channel.fromPath(params.intervals_y).collect()
                                                        : Channel.empty()

    // Prepare references and indices.
    PREPARE_GENOME(
        ch_genome_fasta,
        dbsnp,
        known_indels
    )

    // Gather built indices or get them from the params
    known_indels_tbi            = params.known_indels   ? params.known_indels_tbi      ? Channel.fromPath(params.known_indels_tbi).collect()      : PREPARE_GENOME.out.known_indels_tbi      : Channel.value([])
    dbsnp_tbi                   = params.dbsnp          ? params.dbsnp_tbi             ? Channel.fromPath(params.dbsnp_tbi).collect()             : PREPARE_GENOME.out.dbsnp_tbi             : Channel.value([])
    dict                        = params.dict           ? Channel.fromPath(params.dict).map{ it -> [ [id:'dict'], it ] }.collect()
                                                        : PREPARE_GENOME.out.dict
    ch_genome_chrsizes          = PREPARE_GENOME.out.genome_chrom_sizes


    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    known_sites_indels     = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()
    
    // PREPARE INDICES IF NOT PROVIDED
    PREPARE_INDICES(
        ch_input
    )

    ch_versions = ch_versions.mix(PREPARE_INDICES.out.versions)
    genome_bam_bai  = PREPARE_INDICES.out.genome_bam_bai

    //MARKDUPLICATES AND REMOVE DUPLICATES
    BAM_MARKDUPLICATES(
                genome_bam_bai,
                ch_genome_fasta,
                ch_genome_fai,
    )
    .set {ch_mapped}

    ch_versions = ch_versions.mix(BAM_MARKDUPLICATES.out.versions)
    mapped_bam_bai = BAM_MARKDUPLICATES.out.a_genome_bam_bai

    //BASERECALIBRATOR
    BAM_BASERECALIBRATOR(
                mapped_bam_bai,
                dict,
                ch_genome_fasta,
                ch_genome_fai,
                known_sites_indels,
                known_sites_indels_tbi
    )
    
    ch_versions = ch_versions.mix(BAM_BASERECALIBRATOR.out.versions)
    ch_table_bqsr = BAM_BASERECALIBRATOR.out.table_bqsr
    bam_applybqsr = mapped_bam_bai.join(ch_table_bqsr)

    //APPLYBQSR
    BAM_APPLYBQSR(
                bam_applybqsr,
                dict,
                ch_genome_fasta,
                ch_genome_fai,
    )
    .set{ch_recalibrated}

    ch_versions = ch_versions.mix(BAM_APPLYBQSR.out.versions)

    // BAM QUALITY CHECK
    QC_BAM (
        ch_recalibrated.bqsr_bam,
        ch_recalibrated.bqsr_bai,
        ch_recalibrated.a_bqsr_bam_bai,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_chrsizes,
        ch_intervals_wgs,
        ch_intervals_y
    )

    ch_versions = ch_versions.mix(QC_BAM.out.versions)

    //
    // MODULE: Pipeline reporting
    //

    // The template v2.7.1 template update introduced: ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml')
    // This caused the pipeline to stall
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    ///
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowConfigureinput.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowConfigureinput.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.multiple_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.qualimap_results.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.global_dist.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.cov.map{it[1]}.collect().ifEmpty([]))
//    ch_multiqc_files = ch_multiqc_files.mix(PEDDY_CHECK.out.ped.map{it[1]}.collect().ifEmpty([]))
//    ch_multiqc_files = ch_multiqc_files.mix(PEDDY_CHECK.out.csv.map{it[1]}.collect().ifEmpty([]))


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
