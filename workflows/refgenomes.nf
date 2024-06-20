/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HIFIADAPTERFILT                                } from '../modules/local/hifiadapterfilt/main'
include { FASTQC as FASTQC_HIFI                          } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_HIC                           } from '../modules/nf-core/fastqc/main'
include { MERYL_COUNT                                    } from '../modules/nf-core/meryl/count/main'
include { MERYL_HISTOGRAM                                } from '../modules/nf-core/meryl/histogram/main'
include { GENOMESCOPE2                                   } from '../modules/nf-core/genomescope2/main'
include { CAT_HIC                                        } from '../modules/local/cat_hic/main'
include { HIFIASM                                        } from '../modules/nf-core/hifiasm/main'
include { GFASTATS as GFASTATS_HAP1                      } from '../modules/nf-core/gfastats/main'
include { GFASTATS as GFASTATS_HAP2                      } from '../modules/nf-core/gfastats/main'
include { BUSCO_BUSCO                                    } from '../modules/nf-core/busco/busco/main'
include { BUSCO_GENERATEPLOT                             } from '../modules/nf-core/busco/generateplot/main'
include { MERQURY                                        } from '../modules/nf-core/merqury/main'
include { OMNIC as OMNIC_HAP1                            } from '../subworkflows/local/omnic/main'
include { OMNIC as OMNIC_HAP2                            } from '../subworkflows/local/omnic/main'
include { YAHS as YAHS_HAP1                              } from '../modules/nf-core/yahs/main'
include { YAHS as YAHS_HAP2                              } from '../modules/nf-core/yahs/main'
include { FCS_FCSGX as FCS_FCSGX_HAP1                    } from '../modules/nf-core/fcs/fcsgx/main'
include { FCS_FCSGX as FCS_FCSGX_HAP2                    } from '../modules/nf-core/fcs/fcsgx/main'
include { TIARA_TIARA as TIARA_TIARA_HAP1                } from '../modules/nf-core/tiara/tiara/main'
include { TIARA_TIARA as TIARA_TIARA_HAP2                } from '../modules/nf-core/tiara/tiara/main'
include { BBMAP_FILTERBYNAME as BBMAP_FILTERBYNAME_HAP1  } from '../modules/local/bbmap/filterbyname/main'
include { BBMAP_FILTERBYNAME as BBMAP_FILTERBYNAME_HAP2  } from '../modules/local/bbmap/filterbyname/main'
include { CAT_SCAFFOLDS                                  } from '../modules/local/cat_scaffolds/main'
include { GFASTATS as GFASTATS_FINAL                     } from '../modules/nf-core/gfastats/main'
include { BUSCO_BUSCO as BUSCO_BUSCO_FINAL               } from '../modules/nf-core/busco/busco/main'
include { BUSCO_GENERATEPLOT as BUSCO_GENERATEPLOT_FINAL } from '../modules/nf-core/busco/generateplot/main'
include { MERQURY as MERQURY_FINAL                       } from '../modules/nf-core/merqury/main'
include { TAR                                            } from '../modules/local/tar/main'
include { MD5SUM as MD5SUM_OMNICS_HAP1                   } from '../modules/local/md5sum/main'
include { MD5SUM as MD5SUM_OMNICS_HAP2                   } from '../modules/local/md5sum/main'
include { MD5SUM as MD5SUM_YAHS_HAP1                     } from '../modules/local/md5sum/main'
include { MD5SUM as MD5SUM_YAHS_HAP2                     } from '../modules/local/md5sum/main'
include { RCLONE                                         } from '../modules/local/rclone/main'
include { MULTIQC                                        } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                               } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                         } from '../subworkflows/local/utils_oceangenomesrefgenomes_pipeline'

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
                meta = meta[0]
                meta.id = meta.sample + "_" + meta.version + "_" + meta.date
                return [ meta, meta.hifi_dir ]
        }

    ch_hic = ch_samplesheet
        .map {
            meta ->
                if (meta.hic_dir[0] != null) {
                    meta = meta[0]
                    meta.id = meta.sample + "_" + meta.version + "_" + meta.date
                    return [ meta, meta.hic_dir ]
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
        HIFIADAPTERFILT.out.reads,
        "hifi"
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_HIFI.out.zip)
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
        CAT_HIC.out.cat_files,
        "hic"
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_HIC.out.zip)
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
    ch_gfastats_hap1_in = HIFIASM.out.hap1_contigs.join(GENOMESCOPE2.out.summary)
    ch_gfastats_hap2_in = HIFIASM.out.hap2_contigs.join(GENOMESCOPE2.out.summary)

    GFASTATS_HAP1 (
        ch_gfastats_hap1_in,
        "fasta",
        "",
        "hap1",
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GFASTATS_HAP1.out.versions.first())

    GFASTATS_HAP2 (
        ch_gfastats_hap2_in,
        "fasta",
        "",
        "hap2",
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GFASTATS_HAP2.out.versions.first())

    ch_contig_assemblies = GFASTATS_HAP1.out.assembly.join(GFASTATS_HAP2.out.assembly)
        .map {
            meta, hap1_contigs, hap2_contigs ->
                return [ meta, [ hap1_contigs, hap2_contigs ] ]
        }

    //
    // MODULE: Run Busco
    //
    BUSCO_BUSCO (
        ch_contig_assemblies,
        params.buscomode,
        params.buscodb,
        []
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    //
    // MODULE: Run Busco generate_plot
    //
    BUSCO_GENERATEPLOT (
        BUSCO_BUSCO.out.short_summaries_txt,
        "contigs"
    )
    ch_versions = ch_versions.mix(BUSCO_GENERATEPLOT.out.versions.first())

    //
    // MODULE: Run Merqury
    //
    ch_merqury_in = MERYL_COUNT.out.meryl_db.join(ch_contig_assemblies)

    MERQURY (
        ch_merqury_in,
        "contigs"
    )
    ch_versions = ch_versions.mix(MERQURY.out.versions.first())

    //
    // SUBWORKFLOW: Run omnic workflow
    //
    ch_omnic_hap1_in = CAT_HIC.out.cat_files.join(ch_contig_assemblies)
        .map {
            meta, reads, assemblies ->
                return [ meta, reads, assemblies[0] ]
        }

    ch_omnic_hap2_in = CAT_HIC.out.cat_files.join(ch_contig_assemblies)
        .map {
            meta, reads, assemblies ->
                return [ meta, reads, assemblies[1] ]
        }

    OMNIC_HAP1 (
        ch_omnic_hap1_in,
        "hap1"
    )
    ch_versions = ch_versions.mix(OMNIC_HAP1.out.versions.first())

    OMNIC_HAP2 (
        ch_omnic_hap2_in,
        "hap2"
    )
    ch_versions = ch_versions.mix(OMNIC_HAP2.out.versions.first())

    //
    // MODULE: Run Yahs
    //
    ch_yahs_hap1_in = OMNIC_HAP1.out.bam.join(GFASTATS_HAP1.out.assembly).join(OMNIC_HAP1.out.fai)
    ch_yahs_hap2_in = OMNIC_HAP2.out.bam.join(GFASTATS_HAP2.out.assembly).join(OMNIC_HAP2.out.fai)

    YAHS_HAP1 (
        ch_yahs_hap1_in,
        "hap1"
    )
    ch_versions = ch_versions.mix(YAHS_HAP1.out.versions.first())

    YAHS_HAP2 (
        ch_yahs_hap1_in,
        "hap2"
    )
    ch_versions = ch_versions.mix(YAHS_HAP2.out.versions.first())

    //
    // MODULE: Run Fcsgx
    //
    FCS_FCSGX_HAP1 (
        YAHS_HAP1.out.scaffolds_fasta,
        params.gxdb,
        "hap1"
    )
    ch_versions = ch_versions.mix(FCS_FCSGX_HAP1.out.versions.first())

    FCS_FCSGX_HAP2 (
        YAHS_HAP2.out.scaffolds_fasta,
        params.gxdb,
        "hap2"
    )
    ch_versions = ch_versions.mix(FCS_FCSGX_HAP2.out.versions.first())

    //
    // MODULE: Run Tiara
    //
    TIARA_TIARA_HAP1 (
        YAHS_HAP1.out.scaffolds_fasta,
        "hap1"
    )
    ch_versions = ch_versions.mix(TIARA_TIARA_HAP1.out.versions.first())

    TIARA_TIARA_HAP2 (
        YAHS_HAP2.out.scaffolds_fasta,
        "hap2"
    )
    ch_versions = ch_versions.mix(TIARA_TIARA_HAP2.out.versions.first())

    //
    // MODULE: Run BBmap filterbyname
    //
    ch_bbmap_filterbyname_hap1_in = YAHS_HAP1.out.scaffolds_fasta.join(TIARA_TIARA_HAP1.out.classifications)
    ch_bbmap_filterbyname_hap2_in = YAHS_HAP2.out.scaffolds_fasta.join(TIARA_TIARA_HAP2.out.classifications)

    BBMAP_FILTERBYNAME_HAP1 (
        ch_bbmap_filterbyname_hap1_in,
        "hap1_filtered_scaffolds.fa"
    )
    ch_versions = ch_versions.mix(BBMAP_FILTERBYNAME_HAP1.out.versions.first())

    BBMAP_FILTERBYNAME_HAP2 (
        ch_bbmap_filterbyname_hap2_in,
        "hap2_filtered_scaffolds.fa"
    )
    ch_versions = ch_versions.mix(BBMAP_FILTERBYNAME_HAP2.out.versions.first())

    //
    // MODULE: Rename, and concatenate scaffolds
    //
    ch_filtered_scaffolds = BBMAP_FILTERBYNAME_HAP1.out.scaffolds.join(BBMAP_FILTERBYNAME_HAP2.out.scaffolds)
        .map {
            meta, hap1_scaffolds, hap2_scaffolds ->
                return [ meta, [ hap1_scaffolds, hap2_scaffolds ] ]
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
        params.buscomode,
        params.buscodb,
        []
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO_FINAL.out.versions.first())

    //
    // MODULE: Run Busco generate_plot again
    //
    BUSCO_GENERATEPLOT_FINAL (
        BUSCO_BUSCO_FINAL.out.short_summaries_txt,
        "scaffolds"
    )
    ch_versions = ch_versions.mix(BUSCO_GENERATEPLOT_FINAL.out.versions.first())

    //
    // MODULE: Run Merqury again
    //
    ch_merqury_fin_in = MERYL_COUNT.out.meryl_db.join(CAT_SCAFFOLDS.out.hap1_scaffold).join(CAT_SCAFFOLDS.out.hap2_scaffold)
        .map {
            meta, meryl_db, hap1_scaffold, hap2_scaffold ->
                return [ meta, meryl_db, [ hap1_scaffold, hap2_scaffold ] ]
        }

    MERQURY_FINAL (
        ch_merqury_fin_in,
        "scaffolds"
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

    ch_multiqc_files                      = ch_multiqc_files.mix(BUSCO_GENERATEPLOT.out.png)
    ch_multiqc_files                      = ch_multiqc_files.mix(BUSCO_GENERATEPLOT_FINAL.out.png)
    ch_multiqc_files                      = ch_multiqc_files.mix(GENOMESCOPE2.out.linear_plot_png)
    ch_multiqc_files                      = ch_multiqc_files.mix(GENOMESCOPE2.out.transformed_linear_plot_png)
    ch_multiqc_files                      = ch_multiqc_files.mix(GENOMESCOPE2.out.log_plot_png)
    ch_multiqc_files                      = ch_multiqc_files.mix(GENOMESCOPE2.out.transformed_log_plot_png)
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_cn_fl_png)
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_cn_ln_png)
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_cn_st_png)
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_asm_fl_png)
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_asm_ln_png)
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_asm_st_png)
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY_FINAL.out.spectra_cn_fl_png)
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY_FINAL.out.spectra_cn_ln_png)
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY_FINAL.out.spectra_cn_st_png)
    ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY_FINAL.out.spectra_asm_ln_png)
    ch_multiqc_wf_summary                 = ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    ch_multiqc_versions                   = ch_collated_versions
    ch_multiqc_method_desc                = ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false)

    ch_multiqc_files = ch_multiqc_files
        .groupTuple()
        .map {
            meta, files ->
                return [ meta, files[0], files[1], files[2], files[3], files[4], files[5], files[6], files[7], files[8], files[9], files[10], files[11], files[12], files[13], files[14], files[15], files[16], files[17] ]
        }

    MULTIQC (
        ch_multiqc_files,
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        ch_multiqc_wf_summary.first(),
        ch_multiqc_versions.first(),
        ch_multiqc_method_desc.first()
    )

    //
    // Collect HIFIADAPTERFILT files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        HIFIADAPTERFILT.out.reads
            .map {
                meta, reads ->
                    return [ meta, "hifiadapterfilt_reads", reads, "${params.rclonedest}/${meta.sample}/${meta.id}/hifiadapterfilt" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIADAPTERFILT.out.stats
            .map {
                meta, stats ->
                    return [ meta, "hifiadapterfilt_stats", stats, "${params.rclonedest}/${meta.sample}/${meta.id}/hifiadapterfilt/qc_stats" ]
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
                    return [ meta, "meryl_tar_file", meryl_db, "${params.rclonedest}/${meta.sample}/${meta.id}/meryl" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERYL_HISTOGRAM.out.hist
            .map {
                meta, hist ->
                    return [ meta, "meryl_hist", hist, "${params.rclonedest}/${meta.sample}/${meta.id}/meryl" ]
            }
    )

    //
    // Collect GENOMESCOPE2 files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        GENOMESCOPE2.out.linear_plot_png
            .map {
                meta, linear_plot_png ->
                    return [ meta, "genomescope2_linear_plot_png", linear_plot_png, "${params.rclonedest}/${meta.sample}/${meta.id}/genomescope2" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GENOMESCOPE2.out.transformed_linear_plot_png
            .map {
                meta, transformed_linear_plot_png ->
                    return [ meta, "genomescope2_transformed_linear_plot_png", transformed_linear_plot_png, "${params.rclonedest}/${meta.sample}/${meta.id}/genomescope2" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GENOMESCOPE2.out.log_plot_png
            .map {
                meta, log_plot_png ->
                    return [ meta, "genomescope2_log_plot_png", log_plot_png, "${params.rclonedest}/${meta.sample}/${meta.id}/genomescope2" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GENOMESCOPE2.out.transformed_log_plot_png
            .map {
                meta, transformed_log_plot_png ->
                    return [ meta, "genomescope2_transformed_log_plot_png", transformed_log_plot_png, "${params.rclonedest}/${meta.sample}/${meta.id}/genomescope2" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GENOMESCOPE2.out.model
            .map {
                meta, model ->
                    return [ meta, "genomescope2_model", model, "${params.rclonedest}/${meta.sample}/${meta.id}/genomescope2" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GENOMESCOPE2.out.summary
            .map {
                meta, summary ->
                    return [ meta, "genomescope2_summary", summary, "${params.rclonedest}/${meta.sample}/${meta.id}/genomescope2" ]
            }
    )

    //
    // Collect HIFIASM files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.raw_unitigs
            .map {
                meta, raw_unitigs ->
                    return [ meta, "hifiasm_raw_unitigs", raw_unitigs, "${params.rclonedest}/${meta.sample}/${meta.id}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.corrected_reads
            .map {
                meta, corrected_reads ->
                    return [ meta, "hifiasm_corrected_reads", corrected_reads, "${params.rclonedest}/${meta.sample}/${meta.id}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.source_overlaps
            .map {
                meta, source_overlaps ->
                    return [ meta, "hifiasm_source_overlaps", source_overlaps, "${params.rclonedest}/${meta.sample}/${meta.id}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.reverse_overlaps
            .map {
                meta, reverse_overlaps ->
                    return [ meta, "hifiasm_reverse_overlaps", reverse_overlaps, "${params.rclonedest}/${meta.sample}/${meta.id}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.processed_unitigs
            .map {
                meta, processed_unitigs ->
                    return [ meta, "hifiasm_processed_untigs", processed_unitigs, "${params.rclonedest}/${meta.sample}/${meta.id}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.hap1_contigs
            .map {
                meta, hap1_contigs ->
                    return [ meta, "hifiasm_hap1_logs", hap1_contigs, "${params.rclonedest}/${meta.sample}/${meta.id}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.hap2_contigs
            .map {
                meta, hap2_contigs ->
                    return [ meta, "hifiasm_hap2_logs", hap2_contigs, "${params.rclonedest}/${meta.sample}/${meta.id}/hifiasm" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIASM.out.log
            .map {
                meta, log ->
                    return [ meta, "hifiasm_log", log, "${params.rclonedest}/${meta.sample}/${meta.id}/hifiasm" ]
            }
    )

    //
    // Collect GFASTATS files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        GFASTATS_HAP1.out.assembly
            .map {
                meta, assembly ->
                    return [ meta, "gfastats_hap1_assembly", assembly, "${params.rclonedest}/${meta.sample}/${meta.id}/gfastats" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GFASTATS_HAP2.out.assembly
            .map {
                meta, assembly ->
                    return [ meta, "gfastats_hap2_assembly", assembly, "${params.rclonedest}/${meta.sample}/${meta.id}/gfastats" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GFASTATS_HAP1.out.assembly_summary
            .map {
                meta, assembly_summary ->
                    return [ meta, "gfastats_hap1_assembly_summary", assembly_summary, "${params.rclonedest}/${meta.sample}/${meta.id}/gfastats" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GFASTATS_HAP2.out.assembly_summary
            .map {
                meta, assembly_summary ->
                    return [ meta, "gfastats_hap2_assembly_summary", assembly_summary, "${params.rclonedest}/${meta.sample}/${meta.id}/gfastats" ]
            }
    )

    //
    // Collect BUSCO files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO.out.batch_summary
            .map {
                meta, batch_summary ->
                    return [ meta, "busco_batch_summary", batch_summary, "${params.rclonedest}/${meta.sample}/${meta.id}/busco" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO.out.short_summaries_txt
            .map {
                meta, short_summaries_txt ->
                    return [ meta, "busco_short_summaries_txt", short_summaries_txt, "${params.rclonedest}/${meta.sample}/${meta.id}/busco" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO.out.short_summaries_json
            .map {
                meta, short_summaries_json ->
                    return [ meta, "busco_short_summaries_json", short_summaries_json, "${params.rclonedest}/${meta.sample}/${meta.id}/busco" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO.out.busco_dir
            .map {
                meta, busco_dir ->
                    return [ meta, "busco_busco_dir", busco_dir, "${params.rclonedest}/${meta.sample}/${meta.id}/busco" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_GENERATEPLOT.out.png
            .map {
                meta, png ->
                    return [ meta, "busco_png", png, "${params.rclonedest}/${meta.sample}/${meta.id}/busco" ]
            }
    )

    //
    // Collect MERQURY files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.assembly_only_kmers_bed
            .map {
                meta, assembly_only_kmers_bed ->
                    return [ meta, "merqury_assembly_only_kmers_bed", assembly_only_kmers_bed, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.assembly_only_kmers_wig
            .map {
                meta, assembly_only_kmers_wig ->
                    return [ meta, "merqury_assembly_only_kmers_wig", assembly_only_kmers_wig, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.stats
            .map {
                meta, stats ->
                    return [ meta, "merqury_stats", stats, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.dist_hist
            .map {
                meta, dist_hist ->
                    return [ meta, "merqury_dist_hist", dist_hist, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_cn_fl_png
            .map {
                meta, spectra_cn_fl_png ->
                    return [ meta, "merqury_spectra_cn_fl_png", spectra_cn_fl_png, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_cn_hist
            .map {
                meta, spectra_cn_hist ->
                    return [ meta, "merqury_spectra_cn_hist", spectra_cn_hist, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_cn_ln_png
            .map {
                meta, spectra_cn_ln_png ->
                    return [ meta, "merqury_spectra_cn_ln_png", spectra_cn_ln_png, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_cn_st_png
            .map {
                meta, spectra_cn_st_png ->
                    return [ meta, "merqury_spectra_cn_st_png", spectra_cn_st_png, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_asm_fl_png
            .map {
                meta, spectra_asm_fl_png ->
                    return [ meta, "merqury_spectra_asm_fl_png", spectra_asm_fl_png, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_asm_hist
            .map {
                meta, spectra_asm_hist ->
                    return [ meta, "merqury_spectra_asm_hist", spectra_asm_hist, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_asm_ln_png
            .map {
                meta, spectra_asm_ln_png ->
                    return [ meta, "merqury_spectra_asm_ln_png", spectra_asm_ln_png, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.spectra_asm_st_png
            .map {
                meta, spectra_asm_st_png ->
                    return [ meta, "merqury_spectra_asm_st_png", spectra_asm_st_png, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.assembly_qv
            .map {
                meta, assembly_qv ->
                    return [ meta, "merqury_assembly_qv", assembly_qv, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.scaffold_qv
            .map {
                meta, scaffold_qv ->
                    return [ meta, "merqury_scaffold_qv", scaffold_qv, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY.out.read_ploidy
            .map {
                meta, read_ploidy ->
                    return [ meta, "merqury_read_ploidy", read_ploidy, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury/${meta.tolid}_png" ]
            }
    )

    //
    // Collect OMNIC files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        OMNIC_HAP1.out.fai
            .map {
                meta, fai ->
                    return [ meta, "omnic_hap1_fai", fai, "${params.rclonedest}/${meta.sample}/${meta.id}/omnic" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        OMNIC_HAP2.out.fai
            .map {
                meta, fai ->
                    return [ meta, "omnic_hap2_fai", fai, "${params.rclonedest}/${meta.sample}/${meta.id}/omnic" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        OMNIC_HAP1.out.bam
            .map {
                meta, bam ->
                    return [ meta, "omnic_hap1_bam", bam, "${params.rclonedest}/${meta.sample}/${meta.id}/omnic" ]
            }
    )
    MD5SUM_OMNICS_HAP1(
        OMNIC_HAP1.out.bam
            .map {
                meta, bam ->
                    return [ meta, "omnic_hap1_bam", bam ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MD5SUM_OMNICS_HAP1.out.txt
            .map {
                meta, txt ->
                    return [ meta, "omnic_hap1_bam_md5sum", txt, "${params.rclonedest}/${meta.sample}/${meta.id}/omnic" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        OMNIC_HAP2.out.bam
            .map {
                meta, bam ->
                    return [ meta, "omnic_hap2_bam", bam, "${params.rclonedest}/${meta.sample}/${meta.id}/omnic" ]
            }
    )
    MD5SUM_OMNICS_HAP2(
        OMNIC_HAP2.out.bam
            .map {
                meta, bam ->
                    return [ meta, "omnic_hap2_bam", bam ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MD5SUM_OMNICS_HAP2.out.txt
            .map {
                meta, txt ->
                    return [ meta, "omnic_hap2_bam_md5sum", txt, "${params.rclonedest}/${meta.sample}/${meta.id}/omnic" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        OMNIC_HAP1.out.bai
            .map {
                meta, bai ->
                    return [ meta, "omnic_hap1_bai", bai, "${params.rclonedest}/${meta.sample}/${meta.id}/omnic" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        OMNIC_HAP2.out.bai
            .map {
                meta, bai ->
                    return [ meta, "omnic_hap2_bai", bai, "${params.rclonedest}/${meta.sample}/${meta.id}/omnic" ]
            }
    )

    //
    // Collect YAHS files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        YAHS_HAP1.out.scaffolds_fasta
            .map {
                meta, scaffolds_fasta ->
                    return [ meta, "yahs_hap1_scaffolds_fasta", scaffolds_fasta, "${params.rclonedest}/${meta.sample}/${meta.id}/yahs" ]
            }
    )
    MD5SUM_YAHS_HAP1(
        YAHS_HAP1.out.scaffolds_fasta
            .map {
                meta, scaffolds_fasta ->
                    return [ meta, "yahs_hap1_scaffolds_fasta", scaffolds_fasta ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MD5SUM_YAHS_HAP1.out.txt
            .map {
                meta, txt ->
                    return [ meta, "yahs_hap1_scaffolds_fasta_md5sum", txt, "${params.rclonedest}/${meta.sample}/${meta.id}/yahs" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        YAHS_HAP2.out.scaffolds_fasta
            .map {
                meta, scaffolds_fasta ->
                    return [ meta, "yahs_hap2_scaffolds_fasta", scaffolds_fasta, "${params.rclonedest}/${meta.sample}/${meta.id}/yahs" ]
            }
    )
    MD5SUM_YAHS_HAP2(
        YAHS_HAP2.out.scaffolds_fasta
            .map {
                meta, scaffolds_fasta ->
                    return [ meta, "yahs_hap2_scaffolds_fasta", scaffolds_fasta ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MD5SUM_YAHS_HAP2.out.txt
            .map {
                meta, txt ->
                    return [ meta, "yahs_hap2_scaffolds_fasta_md5sum", txt, "${params.rclonedest}/${meta.sample}/${meta.id}/yahs" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        YAHS_HAP1.out.scaffolds_agp
            .map {
                meta, scaffolds_agp ->
                    return [ meta, "yahs_hap1_scaffolds_agp", scaffolds_agp, "${params.rclonedest}/${meta.sample}/${meta.id}/yahs" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        YAHS_HAP2.out.scaffolds_agp
            .map {
                meta, scaffolds_agp ->
                    return [ meta, "yahs_hap2_scaffolds_agp", scaffolds_agp, "${params.rclonedest}/${meta.sample}/${meta.id}/yahs" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        YAHS_HAP1.out.binary
            .map {
                meta, binary ->
                    return [ meta, "yahs_hap1_binary", binary, "${params.rclonedest}/${meta.sample}/${meta.id}/yahs" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        YAHS_HAP2.out.binary
            .map {
                meta, binary ->
                    return [ meta, "yahs_hap2_binary", binary, "${params.rclonedest}/${meta.sample}/${meta.id}/yahs" ]
            }
    )

    //
    // Collect FCSGX files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        FCS_FCSGX_HAP1.out.fcs_gx_report
            .map {
                meta, fcs_gx_report ->
                    return [ meta, "fcsgx_hap1_gx_report", fcs_gx_report, "${params.rclonedest}/${meta.sample}/${meta.id}/fcsgx" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        FCS_FCSGX_HAP2.out.fcs_gx_report
            .map {
                meta, fcs_gx_report ->
                    return [ meta, "fcsgx_hap2_gx_report", fcs_gx_report, "${params.rclonedest}/${meta.sample}/${meta.id}/fcsgx" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        FCS_FCSGX_HAP1.out.taxonomy_report
            .map {
                meta, taxonomy_report ->
                    return [ meta, "fcsgx_hap1_taxonomy_report", taxonomy_report, "${params.rclonedest}/${meta.sample}/${meta.id}/fcsgx" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        FCS_FCSGX_HAP2.out.taxonomy_report
            .map {
                meta, taxonomy_report ->
                    return [ meta, "fcsgx_hap2_taxonomy_report", taxonomy_report, "${params.rclonedest}/${meta.sample}/${meta.id}/fcsgx" ]
            }
    )

    //
    // Collect TIARA files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        TIARA_TIARA_HAP1.out.classifications
            .map {
                meta, classifications ->
                    return [ meta, "tiara_hap1_classifications", classifications, "${params.rclonedest}/${meta.sample}/${meta.id}/tiara" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        TIARA_TIARA_HAP2.out.classifications
            .map {
                meta, classifications ->
                    return [ meta, "tiara_hap2_classifications", classifications, "${params.rclonedest}/${meta.sample}/${meta.id}/tiara" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        TIARA_TIARA_HAP1.out.log
            .map {
                meta, log ->
                    return [ meta, "tiara_hap1_log", log, "${params.rclonedest}/${meta.sample}/${meta.id}/tiara" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        TIARA_TIARA_HAP2.out.log
            .map {
                meta, log ->
                    return [ meta, "tiara_hap2_log", log, "${params.rclonedest}/${meta.sample}/${meta.id}/tiara" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        TIARA_TIARA_HAP1.out.fasta
            .map {
                meta, fasta ->
                    return [ meta, "tiara_hap1_fasta", fasta, "${params.rclonedest}/${meta.sample}/${meta.id}/tiara" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        TIARA_TIARA_HAP2.out.fasta
            .map {
                meta, fasta ->
                    return [ meta, "tiara_hap2_fasta", fasta, "${params.rclonedest}/${meta.sample}/${meta.id}/tiara" ]
            }
    )

    //
    // Collect concatenated BBMAP_FILTERBYNAME files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        CAT_SCAFFOLDS.out.cat_file
            .map {
                meta, cat_file ->
                    return [ meta, "cat_scaffolds_hap1_hap2_scaffold", cat_file, "${params.rclonedest}/${meta.sample}/${meta.id}/bbmap" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        CAT_SCAFFOLDS.out.hap1_scaffold
            .map {
                meta, hap1_scaffold ->
                    return [ meta, "cat_scaffolds_hap1_scaffold", hap1_scaffold, "${params.rclonedest}/${meta.sample}/${meta.id}/bbmap" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        CAT_SCAFFOLDS.out.hap2_scaffold
            .map {
                meta, hap2_scaffold ->
                    return [ meta, "cat_scaffolds_hap2_scaffold", hap2_scaffold, "${params.rclonedest}/${meta.sample}/${meta.id}/bbmap" ]
            }
    )

    //
    // Collect final GFASTATS files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        GFASTATS_FINAL.out.assembly
            .map {
                meta, assembly ->
                    return [ meta, "gfastats_final_assembly", assembly, "${params.rclonedest}/${meta.sample}/${meta.id}/gfastats_final" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        GFASTATS_FINAL.out.assembly_summary
            .map {
                meta, assembly_summary ->
                    return [ meta, "gfastats_final_assembly_summary", assembly_summary, "${params.rclonedest}/${meta.sample}/${meta.id}/gfastats_final" ]
            }
    )

    //
    // Collect final BUSCO files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO_FINAL.out.batch_summary
            .map {
                meta, batch_summary ->
                    return [ meta, "busco_final_batch_summary", batch_summary, "${params.rclonedest}/${meta.sample}/${meta.id}/busco_final" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO_FINAL.out.short_summaries_txt
            .map {
                meta, short_summaries_txt ->
                    return [ meta, "busco_final_short_summaries_txt", short_summaries_txt, "${params.rclonedest}/${meta.sample}/${meta.id}/busco_final" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO_FINAL.out.short_summaries_json
            .map {
                meta, short_summaries_json ->
                    return [ meta, "busco_final_short_summaries_json", short_summaries_json, "${params.rclonedest}/${meta.sample}/${meta.id}/busco_final" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_BUSCO_FINAL.out.busco_dir
            .map {
                meta, busco_dir ->
                    return [ meta, "busco_final_busco_dir", busco_dir, "${params.rclonedest}/${meta.sample}/${meta.id}/busco_final" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        BUSCO_GENERATEPLOT_FINAL.out.png
            .map {
                meta, png ->
                    return [ meta, "busco_final_png", png, "${params.rclonedest}/${meta.sample}/${meta.id}/busco_final" ]
            }
    )

    //
    // Collect final MERQURY files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.assembly_only_kmers_bed
            .map {
                meta, assembly_only_kmers_bed ->
                    return [ meta, "merqury_final_assembly_only_kmers_bed", assembly_only_kmers_bed, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.assembly_only_kmers_wig
            .map {
                meta, assembly_only_kmers_wig ->
                    return [ meta, "merqury_final_assembly_only_kmers_wig", assembly_only_kmers_wig, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.stats
            .map {
                meta, stats ->
                    return [ meta, "merqury_final_stats", stats, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.dist_hist
            .map {
                meta, dist_hist ->
                    return [ meta, "merqury_final_dist_hist", dist_hist, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.spectra_cn_fl_png
            .map {
                meta, spectra_cn_fl_png ->
                    return [ meta, "merqury_final_spectra_cn_fl_png", spectra_cn_fl_png, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.spectra_cn_hist
            .map {
                meta, spectra_cn_hist ->
                    return [ meta, "merqury_spectra_cn_hist", spectra_cn_hist, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.spectra_cn_ln_png
            .map {
                meta, spectra_cn_ln_png ->
                    return [ meta, "merqury_final_spectra_cn_ln_png", spectra_cn_ln_png, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.spectra_cn_st_png
            .map {
                meta, spectra_cn_st_png ->
                    return [ meta, "merqury_final_spectra_cn_st_png", spectra_cn_st_png, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.spectra_asm_hist
            .map {
                meta, spectra_asm_hist ->
                    return [ meta, "merqury_final_spectra_hist", spectra_asm_hist, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.spectra_asm_ln_png
            .map {
                meta, spectra_asm_ln_png ->
                    return [ meta, "merqury_final_spectra_asm_ln_png", spectra_asm_ln_png, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.assembly_qv
            .map {
                meta, assembly_qv ->
                    return [ meta, "merqury_final_assembly_qv", assembly_qv, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.scaffold_qv
            .map {
                meta, scaffold_qv ->
                    return [ meta, "merqury_final_scaffold_qv", scaffold_qv, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MERQURY_FINAL.out.read_ploidy
            .map {
                meta, read_ploidy ->
                    return [ meta, "merqury_final_read_ploidy", read_ploidy, "${params.rclonedest}/${meta.sample}/${meta.id}/merqury_final/${meta.tolid}_png" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        MULTIQC.out.report
            .map {
                meta, report ->
                    return [ meta, "multiqc", report, "${params.rclonedest}/${meta.sample}/${meta.id}/multiqc" ]
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
