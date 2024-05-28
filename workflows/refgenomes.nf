/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HIFIADAPTERFILT                                   } from '../modules/local/hifiadapterfilt/main'
include { FASTQC as FASTQC_HIFI                             } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_HIC                              } from '../modules/nf-core/fastqc/main'
include { MERYL_COUNT                                       } from '../modules/nf-core/meryl/count/main'
include { MERYL_HISTOGRAM                                   } from '../modules/nf-core/meryl/histogram/main'
include { GENOMESCOPE2                                      } from '../modules/nf-core/genomescope2/main'
include { CAT_HIC                                           } from '../modules/local/cat_hic/main'
include { HIFIASM                                           } from '../modules/nf-core/hifiasm/main'
include { GFASTATS as GFASTATS_PATERNAL                     } from '../modules/nf-core/gfastats/main'
include { GFASTATS as GFASTATS_MATERNAL                     } from '../modules/nf-core/gfastats/main'
include { BUSCO_BUSCO                                       } from '../modules/nf-core/busco/busco/main'
include { BUSCO_GENERATEPLOT                                } from '../modules/nf-core/busco/generateplot/main'
include { MERQURY                                           } from '../modules/nf-core/merqury/main'
include { OMNIC as OMNIC_PATERNAL                           } from '../subworkflows/local/omnic/main'
include { OMNIC as OMNIC_MATERNAL                           } from '../subworkflows/local/omnic/main'
include { YAHS as YAHS_PATERNAL                             } from '../modules/nf-core/yahs/main'
include { YAHS as YAHS_MATERNAL                             } from '../modules/nf-core/yahs/main'
include { FCS_FCSGX as FCS_FCSGX_PATERNAL                   } from '../modules/nf-core/fcs/fcsgx/main'
include { FCS_FCSGX as FCS_FCSGX_MATERNAL                   } from '../modules/nf-core/fcs/fcsgx/main'
include { TIARA_TIARA as TIARA_TIARA_PATERNAL               } from '../modules/nf-core/tiara/tiara/main'
include { TIARA_TIARA as TIARA_TIARA_MATERNAL               } from '../modules/nf-core/tiara/tiara/main'
include { BBMAP_FILTERBYNAME as BBMAP_FILTERBYNAME_PATERNAL } from '../modules/local/bbmap/filterbyname/main'
include { BBMAP_FILTERBYNAME as BBMAP_FILTERBYNAME_MATERNAL } from '../modules/local/bbmap/filterbyname/main'
include { CAT_SCAFFOLDS                                     } from '../modules/local/cat_scaffolds/main'
include { GFASTATS as GFASTATS_FINAL                        } from '../modules/nf-core/gfastats/main'
include { BUSCO_BUSCO as BUSCO_BUSCO_FINAL                  } from '../modules/nf-core/busco/busco/main'
include { BUSCO_GENERATEPLOT as BUSCO_GENERATEPLOT_FINAL    } from '../modules/nf-core/busco/generateplot/main'
include { MERQURY as MERQURY_FINAL                          } from '../modules/nf-core/merqury/main'
include { TAR                                               } from '../modules/local/tar/main'
include { RCLONE                                            } from '../modules/local/rclone/main'
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                  } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                            } from '../subworkflows/local/utils_oceangenomesrefgenomes_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow REFGENOMES {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_rclone_in = Channel.empty()

    ch_hifi = ch_samplesheet
        .map {
            meta ->
                return [ meta[0], meta.hifi_dir[0] ]
        }

    ch_hic = ch_samplesheet
        .map {
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
    // MODULE: Run FastQC on HiFi fastqc files
    //
    FASTQC_HIFI (
        HIFIADAPTERFILT.out.reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_HIFI.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_HIFI.out.versions.first())

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
    // MODULE: Run FastQC on Hi-C fastqc files
    //
    FASTQC_HIC (
        CAT_HIC.out.cat_files
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_HIC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_HIC.out.versions.first())

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
    ch_gfastats_pat_in = HIFIASM.out.paternal_contigs.join(GENOMESCOPE2.out.summary)
    ch_gfastats_mat_in = HIFIASM.out.maternal_contigs.join(GENOMESCOPE2.out.summary)

    GFASTATS_PATERNAL (
        ch_gfastats_pat_in,
        "fasta",
        "",
        "hap1",
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GFASTATS_PATERNAL.out.versions.first())

    GFASTATS_MATERNAL (
        ch_gfastats_mat_in,
        "fasta",
        "",
        "hap2",
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GFASTATS_MATERNAL.out.versions.first())

    ch_gfastats_assemblies = GFASTATS_PATERNAL.out.assembly.join(GFASTATS_MATERNAL.out.assembly)
        .map {
            meta, paternal_contigs, maternal_contigs ->
                return [ meta, [ paternal_contigs, maternal_contigs ] ]
        }

    //
    // MODULE: Run Busco
    //
    BUSCO_BUSCO (
        ch_gfastats_assemblies,
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
    ch_merqury_in = MERYL_COUNT.out.meryl_db.join(ch_gfastats_assemblies)

    MERQURY (
        ch_merqury_in
    )
    ch_versions = ch_versions.mix(MERQURY.out.versions.first())

    //
    // SUBWORKFLOW: Run omnic workflow
    //
    ch_omnic_pat_in = CAT_HIC.out.cat_files.join(ch_gfastats_assemblies)
        .map {
            meta, reads, assemblies ->
                return [ meta, reads, assemblies[0] ]
        }

    ch_omnic_mat_in = CAT_HIC.out.cat_files.join(ch_gfastats_assemblies)
        .map {
            meta, reads, assemblies ->
                return [ meta, reads, assemblies[1] ]
        }

    OMNIC_PATERNAL (
        ch_omnic_pat_in,
        "hap1"
    )
    ch_versions = ch_versions.mix(OMNIC_PATERNAL.out.versions.first())

    OMNIC_MATERNAL (
        ch_omnic_mat_in,
        "hap2"
    )
    ch_versions = ch_versions.mix(OMNIC_MATERNAL.out.versions.first())

    //
    // MODULE: Run Yahs
    //
    ch_yahs_pat_in = OMNIC_PATERNAL.out.bam.join(GFASTATS_PATERNAL.out.assembly).join(OMNIC_PATERNAL.out.fai)
    ch_yahs_mat_in = OMNIC_MATERNAL.out.bam.join(GFASTATS_MATERNAL.out.assembly).join(OMNIC_MATERNAL.out.fai)

    YAHS_PATERNAL (
        ch_yahs_pat_in,
        "hap1"
    )
    ch_versions = ch_versions.mix(YAHS_PATERNAL.out.versions.first())

    YAHS_MATERNAL (
        ch_yahs_pat_in,
        "hap2"
    )
    ch_versions = ch_versions.mix(YAHS_MATERNAL.out.versions.first())

    //
    // MODULE: Run Fcsgx
    //
    FCS_FCSGX_PATERNAL (
        YAHS_PATERNAL.out.scaffolds_fasta,
        params.gx_db,
        "hap1"
    )
    ch_versions = ch_versions.mix(FCS_FCSGX_PATERNAL.out.versions.first())

    FCS_FCSGX_MATERNAL (
        YAHS_MATERNAL.out.scaffolds_fasta,
        params.gx_db,
        "hap2"
    )
    ch_versions = ch_versions.mix(FCS_FCSGX_MATERNAL.out.versions.first())

    //
    // MODULE: Run Tiara
    //
    TIARA_TIARA_PATERNAL (
        YAHS_PATERNAL.out.scaffolds_fasta,
        "hap1"
    )
    ch_versions = ch_versions.mix(TIARA_TIARA_PATERNAL.out.versions.first())

    TIARA_TIARA_MATERNAL (
        YAHS_MATERNAL.out.scaffolds_fasta,
        "hap2"
    )
    ch_versions = ch_versions.mix(TIARA_TIARA_MATERNAL.out.versions.first())

    //
    // MODULE: Run BBmap filterbyname
    //
    ch_bbmap_filterbyname_pat_in = YAHS_PATERNAL.out.scaffolds_fasta.join(TIARA_TIARA_PATERNAL.out.classifications)
    ch_bbmap_filterbyname_mat_in = YAHS_MATERNAL.out.scaffolds_fasta.join(TIARA_TIARA_MATERNAL.out.classifications)

    BBMAP_FILTERBYNAME_PATERNAL (
        ch_bbmap_filterbyname_pat_in,
        "hap1_filtered_scaffolds.fa"
    )
    ch_versions = ch_versions.mix(BBMAP_FILTERBYNAME_PATERNAL.out.versions.first())

    BBMAP_FILTERBYNAME_MATERNAL (
        ch_bbmap_filterbyname_mat_in,
        "hap2_filtered_scaffolds.fa"
    )
    ch_versions = ch_versions.mix(BBMAP_FILTERBYNAME_MATERNAL.out.versions.first())

    //
    // MODULE: Rename, and concatenate scaffolds
    //
    ch_filtered_scaffolds = BBMAP_FILTERBYNAME_PATERNAL.out.scaffolds.join(BBMAP_FILTERBYNAME_MATERNAL.out.scaffolds)
        .map {
            meta, paternal_scaffolds, maternal_scaffolds ->
                return [ meta, [ paternal_scaffolds, maternal_scaffolds ] ]
        }

    CAT_SCAFFOLDS (
        ch_filtered_scaffolds
    )
    ch_versions = ch_versions.mix(CAT_SCAFFOLDS.out.versions.first())

    //
    // MODULE: Run Gfastats again
    //
    ch_gfastats_fin_in = CAT_SCAFFOLDS.out.cat_file.join(GENOMESCOPE2.out.summary)

    GFASTATS_FINAL (
        ch_gfastats_fin_in,
        "fasta",
        "",
        "final",
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GFASTATS_FINAL.out.versions.first())

    //
    // MODULE: Run Busco again
    //
    BUSCO_BUSCO_FINAL (
        CAT_SCAFFOLDS.out.cat_file,
        params.busco_mode,
        params.busco_db,
        []
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO_FINAL.out.versions.first())

    //
    // MODULE: Run Busco generate_plot again
    //
    BUSCO_GENERATEPLOT_FINAL (
        BUSCO_BUSCO_FINAL.out.short_summaries_txt
    )
    ch_versions = ch_versions.mix(BUSCO_GENERATEPLOT_FINAL.out.versions.first())

    //
    // MODULE: Run Merqury again
    //
    ch_merqury_fin_in = MERYL_COUNT.out.meryl_db.join(CAT_SCAFFOLDS.out.paternal_scaffold).join(CAT_SCAFFOLDS.out.maternal_scaffold)
        .map {
            meta, meryl_db, paternal_scaffold, maternal_scaffold ->
                return [ meta, meryl_db, [ paternal_scaffold, maternal_scaffold ] ]
        }

    MERQURY_FINAL (
        ch_merqury_fin_in
    )
    ch_versions = ch_versions.mix(MERQURY_FINAL.out.versions.first())

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
    ch_multiqc_files                      = ch_multiqc_files.mix(BUSCO_GENERATEPLOT.out.png.map { meta, png -> return [ png ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(GENOMESCOPE2.out.linear_plot_png.map { meta, linear_plot_png -> return [ linear_plot_png ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(GENOMESCOPE2.out.transformed_linear_plot_png.map { meta, transformed_linear_plot_png -> return [ transformed_linear_plot_png ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(GENOMESCOPE2.out.log_plot_png.map { meta, log_plot_png -> return [ log_plot_png ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(GENOMESCOPE2.out.transformed_log_plot_png.map { meta, transformed_log_plot_png -> return [ transformed_log_plot_png ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_cn_fl_png.map { meta, spectra_cn_fl_png -> return [ spectra_cn_fl_png[0] ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_cn_ln_png.map { meta, spectra_cn_ln_png -> return [ spectra_cn_ln_png[0] ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_cn_st_png.map { meta, spectra_cn_st_png -> return [ spectra_cn_st_png[0] ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_asm_fl_png.map { meta, spectra_asm_fl_png -> return [ spectra_asm_fl_png[0] ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_asm_ln_png.map { meta, spectra_asm_ln_png -> return [ spectra_asm_ln_png[0] ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_asm_st_png.map { meta, spectra_asm_st_png -> return [ spectra_asm_st_png[0] ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY_FINAL.out.spectra_cn_fl_png.map { meta, spectra_cn_fl_png -> return [ spectra_cn_fl_png[0] ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY_FINAL.out.spectra_cn_ln_png.map { meta, spectra_cn_ln_png -> return [ spectra_cn_ln_png[0] ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY_FINAL.out.spectra_cn_st_png.map { meta, spectra_cn_st_png -> return [ spectra_cn_st_png[0] ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY_FINAL.out.spectra_asm_fl_png.map { meta, spectra_asm_fl_png -> return [ spectra_asm_fl_png[0] ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY_FINAL.out.spectra_asm_ln_png.map { meta, spectra_asm_ln_png -> return [ spectra_asm_ln_png[0] ] })
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY_FINAL.out.spectra_asm_st_png.map { meta, spectra_asm_st_png -> return [ spectra_asm_st_png[0] ] })

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    //
    // Collect HIFIADAPTERFILT files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        HIFIADAPTERFILT.out.reads
            .map {
                meta, reads ->
                    return [ meta, reads, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/hifiadapterfilt" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIADAPTERFILT.out.stats
            .map {
                meta, stats ->
                    return [ meta, stats, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/hifiadapterfilt/qc_stats" ]
            }
    )

    //
    // Collect MERYL files for rclone
    //
    TAR (
        MERYL_COUNT.out.meryl_db,
        "meryldb.tar.gz"
    )
    ch_versions = ch_versions.mix(TAR.out.versions.first())

    ch_rclone_in = ch_rclone_in.mix(
        TAR.out.tar_file
            .map {
                meta, meryl_db ->
                    return [ meta, meryl_db, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/meryl" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERYL_HISTOGRAM.out.hist
            .map {
                meta, hist ->
                    return [ meta, hist, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/meryl" ]
            }
    )

    //
    // Collect GENOMESCOPE2 files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        GENOMESCOPE2.out.linear_plot_png
            .map {
                meta, linear_plot_png ->
                    return [ meta, linear_plot_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/genomescope2" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GENOMESCOPE2.out.transformed_linear_plot_png
            .map {
                meta, transformed_linear_plot_png ->
                    return [ meta, transformed_linear_plot_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/genomescope2" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GENOMESCOPE2.out.log_plot_png
            .map {
                meta, log_plot_png ->
                    return [ meta, log_plot_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/genomescope2" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GENOMESCOPE2.out.transformed_log_plot_png
            .map {
                meta, transformed_log_plot_png ->
                    return [ meta, transformed_log_plot_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/genomescope2" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GENOMESCOPE2.out.model
            .map {
                meta, model ->
                    return [ meta, model, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/genomescope2" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GENOMESCOPE2.out.summary
            .map {
                meta, summary ->
                    return [ meta, summary, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/genomescope2" ]
            }
    )

    //
    // Collect HIFIASM files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.raw_unitigs
            .map {
                meta, raw_unitigs ->
                    return [ meta, raw_unitigs, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.corrected_reads
            .map {
                meta, corrected_reads ->
                    return [ meta, corrected_reads, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.source_overlaps
            .map {
                meta, source_overlaps ->
                    return [ meta, source_overlaps, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.reverse_overlaps
            .map {
                meta, reverse_overlaps ->
                    return [ meta, reverse_overlaps, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.processed_unitigs
            .map {
                meta, processed_unitigs ->
                    return [ meta, processed_unitigs, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.paternal_contigs
            .map {
                meta, paternal_contigs ->
                    return [ meta, paternal_contigs, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.maternal_contigs
            .map {
                meta, maternal_contigs ->
                    return [ meta, maternal_contigs, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.log
            .map {
                meta, log ->
                    return [ meta, log, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/hifiasm" ]
            }
    )

    //
    // Collect GFASTATS files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        GFASTATS_PATERNAL.out.assembly
            .map {
                meta, assembly ->
                    return [ meta, assembly, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/gfastats" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GFASTATS_MATERNAL.out.assembly
            .map {
                meta, assembly ->
                    return [ meta, assembly, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/gfastats" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GFASTATS_PATERNAL.out.assembly_summary
            .map {
                meta, assembly_summary ->
                    return [ meta, assembly_summary, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/gfastats" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GFASTATS_MATERNAL.out.assembly_summary
            .map {
                meta, assembly_summary ->
                    return [ meta, assembly_summary, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/gfastats" ]
            }
    )

    //
    // Collect BUSCO files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO.out.batch_summary
            .map {
                meta, batch_summary ->
                    return [ meta, batch_summary, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/busco" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO.out.short_summaries_txt
            .map {
                meta, short_summaries_txt ->
                    return [ meta, short_summaries_txt, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/busco" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO.out.short_summaries_json
            .map {
                meta, short_summaries_json ->
                    return [ meta, short_summaries_json, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/busco" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO.out.busco_dir
            .map {
                meta, busco_dir ->
                    return [ meta, busco_dir, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/busco" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_GENERATEPLOT.out.png
            .map {
                meta, png ->
                    return [ meta, png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/busco" ]
            }
    )

    //
    // Collect MERQURY files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.assembly_only_kmers_bed
            .map {
                meta, assembly_only_kmers_bed ->
                    return [ meta, assembly_only_kmers_bed, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.assembly_only_kmers_wig
            .map {
                meta, assembly_only_kmers_wig ->
                    return [ meta, assembly_only_kmers_wig, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.stats
            .map {
                meta, stats ->
                    return [ meta, stats, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.dist_hist
            .map {
                meta, dist_hist ->
                    return [ meta, dist_hist, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_cn_fl_png
            .map {
                meta, spectra_cn_fl_png ->
                    return [ meta, spectra_cn_fl_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_cn_hist
            .map {
                meta, spectra_cn_hist ->
                    return [ meta, spectra_cn_hist, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_cn_ln_png
            .map {
                meta, spectra_cn_ln_png ->
                    return [ meta, spectra_cn_ln_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_cn_st_png
            .map {
                meta, spectra_cn_st_png ->
                    return [ meta, spectra_cn_st_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_asm_fl_png
            .map {
                meta, spectra_asm_fl_png ->
                    return [ meta, spectra_asm_fl_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_asm_hist
            .map {
                meta, spectra_asm_hist ->
                    return [ meta, spectra_asm_hist, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_asm_ln_png
            .map {
                meta, spectra_asm_ln_png ->
                    return [ meta, spectra_asm_ln_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_asm_st_png
            .map {
                meta, spectra_asm_st_png ->
                    return [ meta, spectra_asm_st_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.assembly_qv
            .map {
                meta, assembly_qv ->
                    return [ meta, assembly_qv, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.scaffold_qv
            .map {
                meta, scaffold_qv ->
                    return [ meta, scaffold_qv, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.read_ploidy
            .map {
                meta, read_ploidy ->
                    return [ meta, read_ploidy, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury/${meta.tolid}_png" ]
            }
    )

    //
    // Collect OMNIC files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        OMNIC_PATERNAL.out.fai
            .map {
                meta, fai ->
                    return [ meta, fai, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/omnic" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        OMNIC_MATERNAL.out.fai
            .map {
                meta, fai ->
                    return [ meta, fai, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/omnic" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        OMNIC_PATERNAL.out.bam
            .map {
                meta, bam ->
                    return [ meta, bam, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/omnic" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        OMNIC_MATERNAL.out.bam
            .map {
                meta, bam ->
                    return [ meta, bam, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/omnic" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        OMNIC_PATERNAL.out.bai
            .map {
                meta, bai ->
                    return [ meta, bai, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/omnic" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        OMNIC_MATERNAL.out.bai
            .map {
                meta, bai ->
                    return [ meta, bai, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/omnic" ]
            }
    )

    //
    // Collect YAHS files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        YAHS_PATERNAL.out.scaffolds_fasta
            .map {
                meta, scaffolds_fasta ->
                    return [ meta, scaffolds_fasta, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/yahs" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        YAHS_MATERNAL.out.scaffolds_fasta
            .map {
                meta, scaffolds_fasta ->
                    return [ meta, scaffolds_fasta, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/yahs" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        YAHS_PATERNAL.out.scaffolds_agp
            .map {
                meta, scaffolds_agp ->
                    return [ meta, scaffolds_agp, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/yahs" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        YAHS_MATERNAL.out.scaffolds_agp
            .map {
                meta, scaffolds_agp ->
                    return [ meta, scaffolds_agp, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/yahs" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        YAHS_PATERNAL.out.binary
            .map {
                meta, binary ->
                    return [ meta, binary, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/yahs" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        YAHS_MATERNAL.out.binary
            .map {
                meta, binary ->
                    return [ meta, binary, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/yahs" ]
            }
    )

    //
    // Collect FCSGX files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        FCS_FCSGX_PATERNAL.out.fcs_gx_report
            .map {
                meta, fcs_gx_report ->
                    return [ meta, fcs_gx_report, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/fcsgx" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        FCS_FCSGX_MATERNAL.out.fcs_gx_report
            .map {
                meta, fcs_gx_report ->
                    return [ meta, fcs_gx_report, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/fcsgx" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        FCS_FCSGX_PATERNAL.out.taxonomy_report
            .map {
                meta, taxonomy_report ->
                    return [ meta, taxonomy_report, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/fcsgx" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        FCS_FCSGX_MATERNAL.out.taxonomy_report
            .map {
                meta, taxonomy_report ->
                    return [ meta, taxonomy_report, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/fcsgx" ]
            }
    )

    //
    // Collect TIARA files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        TIARA_TIARA_PATERNAL.out.classifications
            .map {
                meta, classifications ->
                    return [ meta, classifications, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/tiara" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        TIARA_TIARA_MATERNAL.out.classifications
            .map {
                meta, classifications ->
                    return [ meta, classifications, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/tiara" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        TIARA_TIARA_PATERNAL.out.log
            .map {
                meta, log ->
                    return [ meta, log, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/tiara" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        TIARA_TIARA_MATERNAL.out.log
            .map {
                meta, log ->
                    return [ meta, log, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/tiara" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        TIARA_TIARA_PATERNAL.out.fasta
            .map {
                meta, fasta ->
                    return [ meta, fasta, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/tiara" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        TIARA_TIARA_MATERNAL.out.fasta
            .map {
                meta, fasta ->
                    return [ meta, fasta, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/tiara" ]
            }
    )

    //
    // Collect concatenated BBMAP_FILTERBYNAME files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        CAT_SCAFFOLDS.out.cat_file
            .map {
                meta, cat_file ->
                    return [ meta, cat_file, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/bbmap_filterbyname_cat" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        CAT_SCAFFOLDS.out.paternal_scaffold
            .map {
                meta, paternal_scaffold ->
                    return [ meta, paternal_scaffold, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/bbmap_filterbyname_cat" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        CAT_SCAFFOLDS.out.maternal_scaffold
            .map {
                meta, maternal_scaffold ->
                    return [ meta, maternal_scaffold, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/bbmap_filterbyname_cat" ]
            }
    )

    //
    // Collect final GFASTATS files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        GFASTATS_FINAL.out.assembly
            .map {
                meta, assembly ->
                    return [ meta, assembly, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/gfastats_final" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GFASTATS_FINAL.out.assembly_summary
            .map {
                meta, assembly_summary ->
                    return [ meta, assembly_summary, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/gfastats_final" ]
            }
    )

    //
    // Collect final BUSCO files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO_FINAL.out.batch_summary
            .map {
                meta, batch_summary ->
                    return [ meta, batch_summary, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/busco_final" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO_FINAL.out.short_summaries_txt
            .map {
                meta, short_summaries_txt ->
                    return [ meta, short_summaries_txt, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/busco_final" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO_FINAL.out.short_summaries_json
            .map {
                meta, short_summaries_json ->
                    return [ meta, short_summaries_json, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/busco_final" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO_FINAL.out.busco_dir
            .map {
                meta, busco_dir ->
                    return [ meta, busco_dir, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/busco_final" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_GENERATEPLOT_FINAL.out.png
            .map {
                meta, png ->
                    return [ meta, png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/busco_final" ]
            }
    )

    //
    // Collect final MERQURY files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.assembly_only_kmers_bed
            .map {
                meta, assembly_only_kmers_bed ->
                    return [ meta, assembly_only_kmers_bed, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.assembly_only_kmers_wig
            .map {
                meta, assembly_only_kmers_wig ->
                    return [ meta, assembly_only_kmers_wig, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.stats
            .map {
                meta, stats ->
                    return [ meta, stats, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.dist_hist
            .map {
                meta, dist_hist ->
                    return [ meta, dist_hist, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.spectra_cn_fl_png
            .map {
                meta, spectra_cn_fl_png ->
                    return [ meta, spectra_cn_fl_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.spectra_cn_hist
            .map {
                meta, spectra_cn_hist ->
                    return [ meta, spectra_cn_hist, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.spectra_cn_ln_png
            .map {
                meta, spectra_cn_ln_png ->
                    return [ meta, spectra_cn_ln_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.spectra_cn_st_png
            .map {
                meta, spectra_cn_st_png ->
                    return [ meta, spectra_cn_st_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.spectra_asm_hist
            .map {
                meta, spectra_asm_hist ->
                    return [ meta, spectra_asm_hist, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.spectra_asm_ln_png
            .map {
                meta, spectra_asm_ln_png ->
                    return [ meta, spectra_asm_ln_png, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.assembly_qv
            .map {
                meta, assembly_qv ->
                    return [ meta, assembly_qv, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.scaffold_qv
            .map {
                meta, scaffold_qv ->
                    return [ meta, scaffold_qv, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.read_ploidy
            .map {
                meta, read_ploidy ->
                    return [ meta, read_ploidy, "${params.rclone_dest}/${meta.id}/${meta.id}_${meta.version}/merqury_final/${meta.tolid}_png" ]
            }
    )

    //
    // MODULE: Run rclone
    //
    RCLONE (
        ch_rclone_in
    )
    ch_versions = ch_versions.mix(RCLONE.out.versions.first())

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
