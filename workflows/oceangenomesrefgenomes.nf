/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HIFIADAPTERFILT                                   } from '../modules/local/hifiadapterfilt/main'
include { FASTQC                                            } from '../modules/nf-core/fastqc/main'
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
include { RCLONE                                            } from '../modules/local/rclone/main'
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                  } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                            } from '../subworkflows/local/utils_nfcore_oceangenomesrefgenomes_pipeline'

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
    ch_gfastats_pat_in = HIFIASM.out.paternal_contigs.join(GENOMESCOPE2.out.summary)
    ch_gfastats_mat_in = HIFIASM.out.maternal_contigs.join(GENOMESCOPE2.out.summary)

    GFASTATS_PATERNAL (
        ch_gfastats_pat_in,
        "hap1.fasta",
        "",
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GFASTATS_PATERNAL.out.versions.first())

    GFASTATS_MATERNAL (
        ch_gfastats_mat_in,
        "hap2.fasta",
        "",
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
        ch_omnic_pat_in
    )
    ch_versions = ch_versions.mix(OMNIC_PATERNAL.out.versions.first())

    OMNIC_MATERNAL (
        ch_omnic_mat_in
    )
    ch_versions = ch_versions.mix(OMNIC_MATERNAL.out.versions.first())

    //
    // MODULE: Run Yahs
    //
    ch_yahs_pat_in = OMNIC_PATERNAL.out.bam.join(GFASTATS_PATERNAL.out.assembly).join(OMNIC_PATERNAL.out.fai)
    ch_yahs_mat_in = OMNIC_MATERNAL.out.bam.join(GFASTATS_MATERNAL.out.assembly).join(OMNIC_MATERNAL.out.fai)

    YAHS_PATERNAL (
        ch_yahs_pat_in
    )
    ch_versions = ch_versions.mix(YAHS_PATERNAL.out.versions.first())

    YAHS_MATERNAL (
        ch_yahs_pat_in
    )
    ch_versions = ch_versions.mix(YAHS_MATERNAL.out.versions.first())

    //
    // MODULE: Run Fcsgx
    //
    FCS_FCSGX_PATERNAL (
        YAHS_PATERNAL.out.scaffolds_fasta,
        params.gx_db
    )
    ch_versions = ch_versions.mix(FCS_FCSGX_PATERNAL.out.versions.first())

    FCS_FCSGX_MATERNAL (
        YAHS_MATERNAL.out.scaffolds_fasta,
        params.gx_db
    )
    ch_versions = ch_versions.mix(FCS_FCSGX_MATERNAL.out.versions.first())

    //
    // MODULE: Run Tiara
    //
    TIARA_TIARA_PATERNAL (
        YAHS_PATERNAL.out.scaffolds_fasta
    )
    ch_versions = ch_versions.mix(TIARA_TIARA_PATERNAL.out.versions.first())

    TIARA_TIARA_MATERNAL (
        YAHS_MATERNAL.out.scaffolds_fasta
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
    ch_merqury_fin_in = MERYL_COUNT.out.meryl_db.join(GFASTATS_FINAL.out.assembly)

    MERQURY_FINAL (
        ch_merqury_fin_in
    )
    ch_versions = ch_versions.mix(MERQURY_FINAL.out.versions.first())

    //
    // Collect files for rclone
    //
    ch_rclone_in = ch_rclone_in.mix(
        HIFIADAPTERFILT.out.reads
            .map {
                meta, reads ->
                    return [ meta, reads, "${params.rclone_dest}/${meta.id}/hifi" ]
            }
    )
    ch_rclone_in = ch_rclone_in.mix(
        HIFIADAPTERFILT.out.stats
            .map {
                meta, stats ->
                    return [ meta, stats, "${params.rclone_dest}/${meta.id}/hifi/qc_stats" ]
            }
    )

    //meryl
    //tar -czvf "${processed_dir}/${tolid}_${seq_date}_meryldb.tar.gz" "${processed_dir}/${sample}.meryl" && rclone copy "${processed_dir}/${tolid}_${seq_date}_meryldb.tar.gz" pawsey0812:oceanomics-assemblies/${sample}/meryl
    //rclone copy "${processed_dir}/${sample}.meryl.hist" pawsey0812:oceanomics-assemblies/${sample}/meryl

    //genomescope2
    //rclone copy "genomescope/" "pawsey0812:oceanomics-assemblies/${sample}/genomescope" --checksum

    //gfastats
    //rclone move "${sample}_${ver}.0.hifiasm.p_ctg.assembly.summary.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/gfastats
//rclone move "${sample}_${ver}.0.hifiasm.a_ctg.assembly.summary.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/gfastats
//rclone move "${sample}_${ver}.0.hifiasm.hap1.p_ctg.assembly.summary.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/gfastats
//rclone move "${sample}_${ver}.0.hifiasm.hap2.p_ctg.assembly.summary.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/gfastats
//rclone move "${sample}_${ver}.0.hifiasm.p_ctg.fasta" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/assembly
//rclone move "${sample}_${ver}.0.hifiasm.a_ctg.fasta" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/assembly
//rclone move "${sample}_${ver}.0.hifiasm.hap1.p_ctg.fasta" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/assembly
//rclone move "${sample}_${ver}.0.hifiasm.hap2.p_ctg.fasta" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/assembly
//rclone move "${sample}_${ver}.0.hifiasm.p_ctg.gfa" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/assembly
//rclone move "${sample}_${ver}.0.hifiasm.a_ctg.gfa" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/assembly
//rclone move "${sample}_${ver}.0.hifiasm.hap1.p_ctg.gfa" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/assembly
//rclone move "${sample}_${ver}.0.hifiasm.hap2.p_ctg.gfa" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/assembly

//busco
//rclone move "${sample}_${ver}.0.hifiasm.busco_figure.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.batch_summary.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.hap1_short_summary.specific.${busco_db}.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.hap1_short_summary.specific.${busco_db}.json" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.hap1_full_table.tsv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.hap1_missing_busco_list.tsv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.hap2_short_summary.specific.${busco_db}.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.hap2_short_summary.specific.${busco_db}.json" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.hap2_full_table.tsv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.hap2_missing_busco_list.tsv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.pri_short_summary.specific.${busco_db}.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.pri_short_summary.specific.${busco_db}.json" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.pri_full_table.tsv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.pri_missing_busco_list.tsv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.alt_short_summary.specific.${busco_db}.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.alt_short_summary.specific.${busco_db}.json" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.alt_full_table.tsv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.0.hifiasm.alt_missing_busco_list.tsv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco

//merqury
//rclone move "${sample}_${ver}.0.hifiasm.hap1.p_ctg.spectra-cn.fl.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.hap1.p_ctg.spectra-cn.ln.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.hap1.p_ctg.spectra-cn.st.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.hap2.p_ctg.spectra-cn.fl.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.hap2.p_ctg.spectra-cn.ln.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.hap2.p_ctg.spectra-cn.st.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.spectra-asm.fl.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.spectra-asm.ln.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.spectra-asm.st.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.spectra-cn.fl.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.spectra-cn.ln.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.spectra-cn.st.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.completeness.stats" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_stats
//rclone move "${sample}_${ver}.0.hifiasm.summary.qv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_qv
//rclone move "${sample}_${ver}.0.hifiasm.hap1.p_ctg.qv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_qv
//rclone move "${sample}_${ver}.0.hifiasm.hap2.p_ctg.qv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_qv

//yahs
//rclone move "${sample}_${ver}.1.yahs.hap1.tiara_filter_summary.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/decontam
//rclone move "${sample}_${ver}.1.yahs.hap1.tiara.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/decontam
//rclone move "${sample}_${ver}.1.yahs.hap2.tiara_filter_summary.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/decontam
//rclone move "${sample}_${ver}.1.yahs.hap2.tiara.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/decontam

//cat
//rclone copy "${sample}_${ver}.2.tiara.hap1_scaffolds.fa" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/assembly
//rclone copy "${sample}_${ver}.2.tiara.hap2_scaffolds.fa" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/assembly
//rclone move "${sample}_${ver}.2.tiara.hap1_scaffolds_renamed.fa" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/assembly
//rclone move "${sample}_${ver}.2.tiara.hap2_scaffolds_renamed.fa" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/assembly
//rclone move "${sample}_${ver}.2.tiara.hap1.hap2_combined_scaffolds.fa" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/assembly

//gfastats_final
//rclone copy "${sample}_${ver}.2.tiara.hap1.scaffolds.summary.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/gfastats
//rclone copy "${sample}_${ver}.2.tiara.hap2.scaffolds.summary.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/gfastats

//busco_final
//rclone move "${sample}_${ver}.2.tiara.busco_figure.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.2.tiara.batch_summary.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.2.tiara.hap1_short_summary.specific.${busco_db}.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.2.tiara.hap1_short_summary.specific.${busco_db}.json" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.2.tiara.hap1_full_table.tsv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.2.tiara.hap1_missing_busco_list.tsv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.2.tiara.hap2_short_summary.specific.${busco_db}.txt" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.2.tiara.hap2_short_summary.specific.${busco_db}.json" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.2.tiara.hap2_full_table.tsv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco
//rclone move "${sample}_${ver}.2.tiara.hap2_missing_busco_list.tsv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/busco

//merqury_final
//rclone move "${sample}_${ver}.0.hifiasm.hap1.p_ctg.spectra-cn.fl.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.hap1.p_ctg.spectra-cn.ln.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.hap1.p_ctg.spectra-cn.st.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.hap2.p_ctg.spectra-cn.fl.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.hap2.p_ctg.spectra-cn.ln.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.hap2.p_ctg.spectra-cn.st.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.spectra-asm.fl.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.spectra-asm.ln.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.spectra-asm.st.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.spectra-cn.fl.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.spectra-cn.ln.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.spectra-cn.st.png" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_png
//rclone move "${sample}_${ver}.0.hifiasm.completeness.stats" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_stats
//rclone move "${sample}_${ver}.0.hifiasm.summary.qv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_qv
//rclone move "${sample}_${ver}.0.hifiasm.hap1.p_ctg.qv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_qv
//rclone move "${sample}_${ver}.0.hifiasm.hap2.p_ctg.qv" pawsey0812:oceanomics-assemblies/${sample}/${sample}_${ver}/merqury/${tolid}_qv


    //
    // MODULE: Run rclone
    //
    RCLONE (
        ch_rclone_in
    )
    ch_versions = ch_versions.mix(RCLONE.out.versions.first())

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
