//
// Control oligo QC subworkflow: Map reads to viral database, generate BAM statistics, and visualize control oligo performance
//

include { MINIMAP2_ALIGN      } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX      } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_IDXSTATS   } from '../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_STATS      } from '../../modules/nf-core/samtools/stats/main'

include { PLOT_OLIGO_QC       } from '../../modules/local/plot_oligo_qc'
include { SAMTOOLS_VIEW_COUNT } from '../../modules/local/samtools_view_count'

workflow CONTROL_OLIGO_QC {
    take:
    ch_reads    // channel: [meta, reads]
    ch_viral_db // channel: [meta2, viral_database.fa]

    main:
    ch_versions = Channel.empty()

    // Map reads to viral database containing control oligo sequences using minimap2
    MINIMAP2_ALIGN(
        ch_reads,
        ch_viral_db,
        true,
        "bai",
        false,
        false,
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    // Index BAM files
    SAMTOOLS_INDEX(MINIMAP2_ALIGN.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // Combine BAM and BAI files for downstream processes
    ch_bam_bai = MINIMAP2_ALIGN.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0])
        .map { meta, bam, bai -> [meta, bam, bai] }

    // Get BAM statistics using idxstats
    SAMTOOLS_IDXSTATS(ch_bam_bai)
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

    // Get comprehensive BAM statistics
    SAMTOOLS_STATS(
        ch_bam_bai,
        [[], []],
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    // Extract control oligo read counts using samtools view
    SAMTOOLS_VIEW_COUNT(
        ch_bam_bai,
        "CTRL_OLIGO",
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_COUNT.out.versions)

    // Concatenate all idxstats files with sample names for visualization
    ch_concatenated_counts = SAMTOOLS_IDXSTATS.out.idxstats
        .map { meta, idxstats ->
            def sample_name = meta.id
            [sample_name, idxstats]
        }
        .collectFile(name: "control_qc_counts.tsv", newLine: false) { sample_name, idxstats ->
            def content = idxstats.text
                .split('\n')
                .findAll { it.trim() }
                .collect { line -> "${sample_name}.bam\t${line}" }
                .join('\n')
            return content + '\n'
        }


    // Generate visualization
    PLOT_OLIGO_QC(ch_concatenated_counts)
    ch_versions = ch_versions.mix(PLOT_OLIGO_QC.out.versions)

    emit:
    bam      = MINIMAP2_ALIGN.out.bam // channel: [meta, bam]
    bai      = SAMTOOLS_INDEX.out.bai // channel: [meta, bai]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [meta, idxstats]
    stats    = SAMTOOLS_STATS.out.stats // channel: [meta, stats]
    versions = ch_versions // channel: versions
}
