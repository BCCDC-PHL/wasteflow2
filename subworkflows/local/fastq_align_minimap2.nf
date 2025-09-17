//
// Alignment with Bowtie2
//

include { MINIMAP2_ALIGN          } from '../../modules/nf-core/minimap2/align'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_MINIMAP2 {
    take:
    ch_reads // channel: [ val(meta), [ reads ] ]
    ch_fasta // channel: [ val(meta), [ ref ] ]

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with Minimap2
    //
    bam_format = 'true'
    bam_format_ext = 'bai'
    cigar_paf_format = 'false'
    cigar_bam = 'false'
    MINIMAP2_ALIGN(ch_reads, ch_fasta, bam_format, bam_format_ext, cigar_paf_format, cigar_bam)
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS(MINIMAP2_ALIGN.out.bam, ch_fasta)
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    paf_orig       = MINIMAP2_ALIGN.out.paf // channel: [ val(meta), paf     ]
    bam_orig       = MINIMAP2_ALIGN.out.bam // channel: [ val(meta), aligned ]
    minimap2_index = MINIMAP2_ALIGN.out.index // channel: [ val(meta), [ index ] ]
    bam            = BAM_SORT_STATS_SAMTOOLS.out.bam // channel: [ val(meta), [ bam ] ]
    bai            = BAM_SORT_STATS_SAMTOOLS.out.bai // channel: [ val(meta), [ bai ] ]
    csi            = BAM_SORT_STATS_SAMTOOLS.out.csi // channel: [ val(meta), [ csi ] ]
    stats          = BAM_SORT_STATS_SAMTOOLS.out.stats // channel: [ val(meta), [ stats ] ]
    flagstat       = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats       = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    versions       = ch_versions // channel: [ versions.yml ]
}
