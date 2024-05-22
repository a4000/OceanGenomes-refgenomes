include { SAMTOOLS_FAIDX         } from '../../../modules/nf-core/samtools/faidx/main'
include { CUT                    } from '../../../modules/local/cut/main'
include { BWA_INDEX              } from '../../../modules/nf-core/bwa/index/main'
include { BWA_MEM                } from '../../../modules/nf-core/bwa/mem/main'
include { PAIRTOOLS_PARSE        } from '../../../modules/nf-core/pairtools/parse/main'
include { PAIRTOOLS_SORT         } from '../../../modules/nf-core/pairtools/sort/main'
include { PAIRTOOLS_DEDUP        } from '../../../modules/nf-core/pairtools/dedup/main'
include { PAIRTOOLS_SPLIT        } from '../../../modules/local/pairtools/split/main'
include { SAMTOOLS_SORT          } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX         } from '../../../modules/nf-core/samtools/index/main'

workflow OMNIC {

    take:
    ch_hic_read     // channel: [ val(meta), [ reads ] ]
    ch_fasta        // channel: [ val(meta), fasta ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Create an index for the reference
    //
    SAMTOOLS_FAIDX (
        ch_fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    //
    // MODULE: Use the index to generate a genome file
    //
    CUT (
        SAMTOOLS_FAIDX.out.fai
    )
    ch_versions = ch_versions.mix(CUT.out.versions.first())

    //
    // MODULE: Generate bwa index files
    //
    BWA_INDEX (
        ch_fasta
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

    //
    // MODULE: Run Alignment
    //
    BWA_MEM (
        ch_fasta,
        ch_hic_read
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    //
    // MODULE: Record valid ligation events
    //
    PAIRTOOLS_PARSE (
        CUT.out.genome,
        BWA_MEM.out.sam
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_PARSE.out.versions.first())

    //
    // MODULE: Sort the pairsam file
    //
    PAIRTOOLS_SORT (
        PAIRTOOLS_PARSE.out.pairsam
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_SORT.out.versions.first())

    //
    // MODULE: Remove PCR duplicates
    //
    PAIRTOOLS_DEDUP (
        PAIRTOOLS_SORT.out.pairsam
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_DEDUP.out.versions.first())

    //
    // MODULE: Generate .pairs and bam files
    //
    PAIRTOOLS_SPLIT (
        PAIRTOOLS_DEDUP.out.pairsam
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_SPLIT.out.versions.first())

    //
    // MODULE: Generate the final bam file
    //
    SAMTOOLS_SORT (
        PAIRTOOLS_SPLIT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    //
    // MODULE: Index the bam file
    //
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    bam      = SAMTOOLS_SORT.out.bam    // channel: [ val(meta), bam ]
    bai      = SAMTOOLS_INDEX.out.bai   // channel: [ val(meta), bai ]
    versions = ch_versions              // channel: [ versions.yml ]
}
