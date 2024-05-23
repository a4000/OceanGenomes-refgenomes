/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HIFIADAPTERFILT              } from '../modules/local/hifiadapterfilt/main'
include { FASTQC                       } from '../modules/nf-core/fastqc/main'
include { MERYL_COUNT                  } from '../modules/nf-core/meryl/count/main'
include { MERYL_HISTOGRAM              } from '../modules/nf-core/meryl/histogram/main'
include { GENOMESCOPE2                 } from '../modules/nf-core/genomescope2/main'
include { CAT_HIC                      } from '../modules/local/cat_hic/main'
include { HIFIASM                      } from '../modules/nf-core/hifiasm/main'
include { GFASTATS as GFASTATS_TOFASTA } from '../modules/nf-core/gfastats/main'
include { GFASTATS                     } from '../modules/nf-core/gfastats/main'
include { BUSCO_BUSCO                  } from '../modules/nf-core/busco/busco/main'
include { BUSCO_GENERATEPLOT           } from '../modules/nf-core/busco/generateplot/main'
include { MERQURY                      } from '../modules/nf-core/merqury/main'
include { OMNIC                        } from '../subworkflows/local/omnic/main'
include { YAHS                         } from '../modules/nf-core/yahs/main'
include { FCS_FCSGX                    } from '../modules/nf-core/fcs/fcsgx/main'
include { TIARA_TIARA                  } from '../modules/nf-core/tiara/tiara/main'
include { BBMAP_FILTERBYNAME           } from '../modules/local/bbmap/filterbyname/main'
include { CAT                          } from '../modules/local/cat/main'
include { GFASTATS as GFASTATS2        } from '../modules/nf-core/gfastats/main'
include { BUSCO_BUSCO as BUSCO_BUSCO2  } from '../modules/nf-core/busco/busco/main'
include { MERQURY as MERQURY2          } from '../modules/nf-core/merqury/main'
include { RCLONE                       } from '../modules/local/rclone/main'
include { MULTIQC                      } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap             } from 'plugin/nf-validation'
include { paramsSummaryMultiqc         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText       } from '../subworkflows/local/utils_nfcore_oceangenomesrefgenomes_pipeline'

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

    ch_hifi = ch_samplesheet.map {
        meta ->
            return [ meta[0], meta.hifi_dir[0] ]
    }

    ch_hic = ch_samplesheet.map {
        meta ->
            if (meta.hic_dir[0] != null) {
                return [ meta[0], meta.hic_dir[0] ]
            }
    }

    //
    // MODULE: Run HiFiAdapterFilt
    //
    HIFIADAPTERFILT (
        ch_hifi
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
    MERYL_COUNT (
        HIFIADAPTERFILT.out.reads,
        params.kvalue
    )
    ch_versions = ch_versions.mix(MERYL_COUNT.out.versions.first())

    //
    // MODULE: Run Meryl histogram
    //
    MERYL_HISTOGRAM (
        MERYL_COUNT.out.meryl_db,
        params.kvalue
    )
    ch_versions = ch_versions.mix(MERYL_HISTOGRAM.out.versions.first())

    //
    // MODULE: Run Genomescope2
    //
    GENOMESCOPE2 (
        MERYL_HISTOGRAM.out.hist
    )
    ch_versions = ch_versions.mix(GENOMESCOPE2.out.versions.first())

    //
    // MODULE: Concatenate Hi-C files together for cases when there is multiple R1 and multiple R2 files
    //
    CAT_HIC (
        ch_hic
    )

    //
    // MODULE: Run Hifiasm
    //
    ch_hifiasm_in = HIFIADAPTERFILT.out.reads.join(CAT_HIC.out.cat_files)
    .map {
        meta, hifi, hic ->
            return [ meta, hifi, hic[0], hic[1] ]
    }

    HIFIASM (
        ch_hifiasm_in,
        [],
        []
    )
    ch_versions = ch_versions.mix(HIFIASM.out.versions.first())

    //
    // MODULE: Run Gfastats
    //
    ch_gfastats_in = HIFIASM.out.paternal_contigs.join(GENOMESCOPE2.out.summary)

    GFASTATS (
        ch_gfastats_in,
        "fasta",
        "",
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GFASTATS.out.versions.first())

    //
    // MODULE: Run Busco
    //
    BUSCO_BUSCO (
        GFASTATS.out.assembly,
        params.busco_mode,
        params.busco_db,
        []
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    //
    // MODULE: Run Busco generate_plot
    //
    BUSCO_GENERATEPLOT (
        BUSCO_BUSCO.out.short_summaries_txt
    )
    ch_versions = ch_versions.mix(BUSCO_GENERATEPLOT.out.versions.first())

    //
    // MODULE: Run Merqury
    //
    ch_merqury_in = MERYL_COUNT.out.meryl_db.join(GFASTATS.out.assembly)

    MERQURY (
        ch_merqury_in
    )
    ch_versions = ch_versions.mix(MERQURY.out.versions.first())

    //
    // SUBWORKFLOW: Run omnic workflow
    //
    ch_omnic_in = CAT_HIC.out.cat_files.join(GFASTATS.out.assembly)

    OMNIC (
        ch_omnic_in
    )
    ch_versions = ch_versions.mix(OMNIC.out.versions.first())

    //
    // MODULE: Run Yahs
    //
    ch_yahs_in = OMNIC.out.bam.join(GFASTATS.out.assembly).join(OMNIC.out.bai)

    YAHS (
        ch_yahs_in
    )
    ch_versions = ch_versions.mix(YAHS.out.versions.first())

    //
    // MODULE: Run Fcsgx
    //
    FCS_FCSGX (
        YAHS.out.scaffolds_fasta,
        params.ncbi_db
    )
    ch_versions = ch_versions.mix(FCS_FCSGX.out.versions.first())

    //
    // MODULE: Run Tiara
    //
    TIARA_TIARA (
        YAHS.out.scaffolds_fasta
    )
    ch_versions = ch_versions.mix(TIARA_TIARA.out.versions.first())

    //
    // MODULE: Run BBmap filterbyname
    //
    ch_bbmap_filterbyname_in = YAHS.out.scaffolds_fasta.join(TIARA_TIARA.out.classifications)

    BBMAP_FILTERBYNAME (
        ch_bbmap_filterbyname_in
    )
    ch_versions = ch_versions.mix(BBMAP_FILTERBYNAME.out.versions.first())

    //
    // MODULE: Rename, and concatenate scaffolds
    //
    //CAT (
    //    TIARA_TIARA.out.fasta,
    //    "fasta"
    //)
    //ch_versions = ch_versions.mix(CAT.out.versions.first())

    //
    // MODULE: Run Gfastats again
    //
    //GFASTATS2 (
    //    BBMAP_FILTERBYNAME.out.scaffolds,
    //    "fasta",
    //    [], // Get genome size from this file
    //    [],
    //    [],
    //    [],
    //    []
    //)
    //ch_versions = ch_versions.mix(GFASTATS2.out.versions.first())

    //
    // MODULE: Run Busco again
    //
    //BUSCO_BUSCO2 (
    //    BBMAP_FILTERBYNAME.out.scaffolds,
    //    params.busco_mode,
    //    params.busco_db,
    //    []
    //)
    //ch_versions = ch_versions.mix(BUSCO_BUSCO2.out.versions.first())

    //
    // MODULE: Run Busco generate_plot again
    //
    //BUSCO_GENERATEPLOT2 (
    //    BUSCO_BUSCO2.out.short_summaries_txt
    //)
    //ch_versions = ch_versions.mix(BUSCO_GENERATEPLOT2.out.versions.first())

    //
    // MODULE: Run Merqury again
    //
    //ch_merqury2_in = MERYL_COUNT.out.meryl_db.join(GFASTATS2.out.assembly)
    //MERQURY2 (
    //    ch_merqury2_in
    //)
    //ch_versions = ch_versions.mix(MERQURY2.out.versions.first())

    //
    // MODULE: Run Rclone
    //
    //RCLONE (
    //    ch_all_out_files,
    //    params.rclone_dest
    //)
    //ch_versions = ch_versions.mix(RCLONE.out.versions.first())

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
