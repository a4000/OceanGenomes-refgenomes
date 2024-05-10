/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HIFIADAPTERFILT        } from '../modules/local/hifiadapterfilt/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MERYL_MERYL            } from '../modules/local/meryl/meryl/main'
include { MERYL_HISTOGRAM        } from '../modules/nf-core/meryl/histogram/main'
include { GENOMESCOPE2           } from '../modules/nf-core/genomescope2/main'
include { HIFIASM                } from '../modules/nf-core/hifiasm/main'
include { GFASTATS               } from '../modules/nf-core/gfastats/main'
include { BUSCO_BUSCO            } from '../modules/nf-core/busco/busco/main'
include { BUSCO_GENERATEPLOT     } from '../modules/nf-core/busco/generateplot/main'
include { MERQURY                } from '../modules/nf-core/merqury/main'
include { OMNIC                  } from '../subwokflows/local/omnic/main'
include { YAHS                   } from '../modules/nf-core/yahs/main'
include { FCS_FCSGX              } from '../modules/nf-core/fcs/fcsgx/main'
include { TIARA_TIARA            } from '../modules/nf-core/tiara/tiara/main'
include { BBMAP_FILTERBYNAME     } from '../modules/nf-core/bbmap/filterbyname/main'
include { RENAMECATSCAFF         } from '../modules/local/custom/renamecatscaff/main'
include { GFASTATS as GFASTATS2  } from '../modules/nf-core/gfastats/main'
include { BUSCO_BUSCO as BUSCO2  } from '../modules/nf-core/busco/busco/main'
include { MERQURY as MERQURY2    } from '../modules/nf-core/merqury/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_oceangenomesrefgenomes_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow OCEANGENOMESREFGENOMES {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run HiFiAdapterFilt
    //
    HIFIADAPTERFILT (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(HIFIADAPTERFILT.out.versions.first())

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        HIFIADAPTERFILT.out.reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run Meryl
    //
    MERYL_MERYL (

    )
    ch_versions = ch_versions.mix(MERYL_MERYL.out.versions.first())

    //
    // MODULE: Run Meryl histogram
    //
    MERYL_HISTOGRAM (

    )
    ch_versions = ch_versions.mix(MERYL_HISTOGRAM.out.versions.first())

    //
    // MODULE: Run Genomescope2
    //
    GENOMESCOPE2 (

    )
    ch_versions = ch_versions.mix(GENOMESCOPE2.out.versions.first())

    //
    // MODULE: Run Hifiasm
    //
    HIFIASM (

    )
    ch_versions = ch_versions.mix(HIFIASM.out.versions.first())

    //
    // MODULE: Run Gfastats
    //
    GFASTATS (

    )
    ch_versions = ch_versions.mix(GFASTATS.out.versions.first())

    //
    // MODULE: Run Busco
    //
    BUSCO_BUSCO (

    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    //
    // MODULE: Run Busco generate_plot
    //
    BUSCO_GENERATEPLOT (

    )
    ch_versions = ch_versions.mix(BUSCO_GENERATEPLOT.out.versions.first())

    //
    // MODULE: Run Merqury
    //
    MERQURY (

    )
    ch_versions = ch_versions.mix(MERQURY.out.versions.first())

    //
    // SUBWORKFLOW: Run omnic workflow
    //
    OMNIC (

    )
    ch_versions = ch_versions.mix(OMNIC.out.versions.first())

    //
    // MODULE: Run Yahs
    //
    YAHS (

    )
    ch_versions = ch_versions.mix(YAHS.out.versions.first())

    //
    // MODULE: Run Fcsgx
    //
    FCS_FCSGX (

    )
    ch_versions = ch_versions.mix(FCS_FCSGX.out.versions.first())

    //
    // MODULE: Run Tiara
    //
    TIARA_TIARA (

    )
    ch_versions = ch_versions.mix(TIARA_TIARA.out.versions.first())

    //
    // MODULE: Run BBmap filterbyname
    //
    BBMAP_FILTERBYNAME (

    )
    ch_versions = ch_versions.mix(BBMAP_FILTERBYNAME.out.versions.first())

    //
    // MODULE: Rename, and concatenate scaffolds
    //
    RENAMECATSCAFF (

    )
    ch_versions = ch_versions.mix(RENAMECATSCAFF.out.versions.first())

    //
    // MODULE: Run Gfastats again
    //
    GFASTATS2 (

    )
    ch_versions = ch_versions.mix(GFASTATS2.out.versions.first())

    //
    // MODULE: Run Busco again
    //
    BUSCO2 (

    )
    ch_versions = ch_versions.mix(BUSCO2.out.versions.first())

    //
    // MODULE: Run Merqury again
    //
    MERQURY2 (

    )
    ch_versions = ch_versions.mix(MERQURY2.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
