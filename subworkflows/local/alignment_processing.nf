#!/usr/bin/env nextflow

// Import nf-core subworkflows
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_RSV_A     } from '../nf-core/bam_sort_stats_samtools'
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_RSV_B     } from '../nf-core/bam_sort_stats_samtools'
// Import local subworkflows
include { FASTQ_ALIGN_MINIMAP2 as FASTQ_ALIGN_MINIMAP2_SARS_COV2       } from './fastq_align_minimap2'
include { FASTQ_ALIGN_MINIMAP2 as FASTQ_ALIGN_MINIMAP2_RSV_COMBINED    } from './fastq_align_minimap2'
include { FASTQ_ALIGN_MINIMAP2 as FASTQ_ALIGN_MINIMAP2_H1N1            } from './fastq_align_minimap2'
include { FASTQ_ALIGN_MINIMAP2 as FASTQ_ALIGN_MINIMAP2_H3N2            } from './fastq_align_minimap2'
include { FASTQ_ALIGN_MINIMAP2 as FASTQ_ALIGN_MINIMAP2_H5N1            } from './fastq_align_minimap2'

// Import modules
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_RSV_A                         } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_RSV_B                         } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_DICT as SAMTOOLS_DICT_RSV_A                         } from '../../modules/nf-core/samtools/dict/main'
include { SAMTOOLS_DICT as SAMTOOLS_DICT_RSV_B                         } from '../../modules/nf-core/samtools/dict/main'
include { SAMTOOLS_REHEADER_VIA_SAM as SAMTOOLS_REHEADER_VIA_SAM_RSV_A } from '../../modules/local/samtools_reheader'
include { SAMTOOLS_REHEADER_VIA_SAM as SAMTOOLS_REHEADER_VIA_SAM_RSV_B } from '../../modules/local/samtools_reheader'

workflow ALIGNMENT {
    take:
    ch_all_extracted_reads // channel: [meta, reads] - all extracted reads mixed together
    sars_cov2_fasta        // channel: SARS-CoV-2 fasta
    rsv_combined_fasta     // channel: RSV combined fasta
    rsv_a_fasta            // channel: RSV-A fasta
    rsv_b_fasta            // channel: RSV-B fasta
    h1n1_fasta             // channel: H1N1 fasta
    h3n2_fasta             // channel: H3N2 fasta
    h5n1_fasta             // channel: H5N1 fasta

    main:
    ch_versions = Channel.empty()
    ch_bam = Channel.empty()
    ch_bai = Channel.empty()
    ch_minimap2_multiqc = Channel.empty()
    ch_minimap2_flagstat_multiqc = Channel.empty()

    // Create dictionary files from individual references
    SAMTOOLS_DICT_RSV_A(
        rsv_a_fasta.map { fasta -> [[id: 'PP109421.1'], fasta] }
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DICT_RSV_A.out.versions)

    SAMTOOLS_DICT_RSV_B(
        rsv_b_fasta.map { fasta -> [[id: 'OP975389.1'], fasta] }
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DICT_RSV_B.out.versions)

    // Filter and prepare SARS-CoV-2 samples (taxid: 2697049)
    ch_sars_cov2_samples = ch_all_extracted_reads
        .filter { meta, reads ->
            def is_sars_cov2 = meta.taxid == '2697049'
            return is_sars_cov2
        }
        .map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.genome = 'MN908947.3'
            [new_meta, reads]
        }

    // Filter and prepare RSV samples (taxid: 3049954) - using COMBINED RSV reference
    ch_rsv_samples = ch_all_extracted_reads
        .filter { meta, reads ->
            def is_rsv = meta.taxid == '3049954'
            return is_rsv
        }
        .map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.genome = 'RSV_COMBINED'
            [new_meta, reads]
        }

    // Filter and prepare RSV samples (serotype: H1N1) - using H1N1 reference
    ch_h1n1_samples = ch_all_extracted_reads
        .filter { meta, reads ->
            def is_h1n1 = meta.serotype == 'H1N1'
            return is_h1n1
        }
        .map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.genome = 'H1N1'
            [new_meta, reads]
        }
    
    // Filter and prepare RSV samples (serotype: H3N2) - using COMBINED RSV reference
    ch_h3n2_samples = ch_all_extracted_reads
        .filter { meta, reads ->
            def is_h3n2 = meta.serotype == 'H3N2'
            return is_h3n2
        }
        .map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.genome = 'H3N2'
            [new_meta, reads]
        }
    // Filter and prepare RSV samples (serotype: H5N1) - using COMBINED RSV reference
    ch_h5n1_samples = ch_all_extracted_reads
        .filter { meta, reads ->
            def is_h5n1 = meta.serotype == 'H5N1'
            return is_h5n1
        }
        .map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.genome = 'H5N1'
            [new_meta, reads]
        }
    // Run SARS-CoV-2 alignment
    FASTQ_ALIGN_MINIMAP2_SARS_COV2(
        ch_sars_cov2_samples,
        sars_cov2_fasta.map { fasta -> [[id: 'MN908947.3'], fasta] },
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_MINIMAP2_SARS_COV2.out.versions)

    FASTQ_ALIGN_MINIMAP2_RSV_COMBINED(
        ch_rsv_samples,
        rsv_combined_fasta.map { fasta -> [[id: 'RSV_COMBINED'], fasta] },
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_MINIMAP2_RSV_COMBINED.out.versions)

    FASTQ_ALIGN_MINIMAP2_H1N1(
        ch_h1n1_samples,
        h1n1_fasta.map { fasta -> [[id: 'H1N1'], fasta] },
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_MINIMAP2_H1N1.out.versions)

    FASTQ_ALIGN_MINIMAP2_H3N2(
        ch_h3n2_samples,
        h3n2_fasta.map { fasta -> [[id: 'H3N2'], fasta] },
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_MINIMAP2_H3N2.out.versions)
    
    FASTQ_ALIGN_MINIMAP2_H5N1(
        ch_h5n1_samples,
        h5n1_fasta.map { fasta -> [[id: 'H5N1'], fasta] },
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_MINIMAP2_H5N1.out.versions)

    // Mix SARS-CoV-2 and Influenza outputs
    ch_bam = ch_bam.mix(FASTQ_ALIGN_MINIMAP2_SARS_COV2.out.bam)
    ch_bam = ch_bam.mix(FASTQ_ALIGN_MINIMAP2_H1N1.out.bam)
    ch_bam = ch_bam.mix(FASTQ_ALIGN_MINIMAP2_H3N2.out.bam)
    ch_bam = ch_bam.mix(FASTQ_ALIGN_MINIMAP2_H5N1.out.bam)
    ch_bai = ch_bai.mix(FASTQ_ALIGN_MINIMAP2_SARS_COV2.out.bai)
    ch_bai = ch_bai.mix(FASTQ_ALIGN_MINIMAP2_H1N1.out.bai)
    ch_bai = ch_bai.mix(FASTQ_ALIGN_MINIMAP2_H3N2.out.bai)
    ch_bai = ch_bai.mix(FASTQ_ALIGN_MINIMAP2_H5N1.out.bai)
    ch_minimap2_flagstat_multiqc = ch_minimap2_flagstat_multiqc.mix(FASTQ_ALIGN_MINIMAP2_SARS_COV2.out.flagstat)
    ch_minimap2_flagstat_multiqc = ch_minimap2_flagstat_multiqc.mix(FASTQ_ALIGN_MINIMAP2_H1N1.out.flagstat)
    ch_minimap2_flagstat_multiqc = ch_minimap2_flagstat_multiqc.mix(FASTQ_ALIGN_MINIMAP2_H3N2.out.flagstat)
    ch_minimap2_flagstat_multiqc = ch_minimap2_flagstat_multiqc.mix(FASTQ_ALIGN_MINIMAP2_H5N1.out.flagstat)

    // Prepare RSV separation inputs
    ch_rsv_bam_bai = FASTQ_ALIGN_MINIMAP2_RSV_COMBINED.out.bam.join(FASTQ_ALIGN_MINIMAP2_RSV_COMBINED.out.bai, by: [0])

    // Prepare RSV-A separation inputs
    ch_rsv_a_input = ch_rsv_bam_bai.map { meta, bam, bai ->
        def meta_rsv_a = meta.clone()
        meta_rsv_a.genome = 'PP109421.1'
        [meta_rsv_a, bam, bai]
    }

    // Prepare RSV-B separation inputs
    ch_rsv_b_input = ch_rsv_bam_bai.map { meta, bam, bai ->
        def meta_rsv_b = meta.clone()
        meta_rsv_b.genome = 'OP975389.1'
        [meta_rsv_b, bam, bai]
    }

    // Separate RSV-A reads
    SAMTOOLS_VIEW_RSV_A(
        ch_rsv_a_input,
        rsv_combined_fasta.map { fasta -> [[id: 'RSV_COMBINED'], fasta] },
        [],
        'bai',
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_RSV_A.out.versions)

    // Separate RSV-B reads
    SAMTOOLS_VIEW_RSV_B(
        ch_rsv_b_input,
        rsv_combined_fasta.map { fasta -> [[id: 'RSV_COMBINED'], fasta] },
        [],
        'bai',
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_RSV_B.out.versions)

    // Prepare channels for re-headering
    ch_rsv_a_for_reheader = SAMTOOLS_VIEW_RSV_A.out.bam
        .combine(SAMTOOLS_DICT_RSV_A.out.dict.map { meta, dict -> dict })
        .map { meta, bam, dict ->
            def new_meta = meta.clone()
            [new_meta, bam, dict]
        }

    ch_rsv_b_for_reheader = SAMTOOLS_VIEW_RSV_B.out.bam
        .combine(SAMTOOLS_DICT_RSV_B.out.dict.map { meta, dict -> dict })
        .map { meta, bam, dict ->
            def new_meta = meta.clone()
            [new_meta, bam, dict]
        }

    // Re-header RSV BAMs
    SAMTOOLS_REHEADER_VIA_SAM_RSV_A(ch_rsv_a_for_reheader)
    SAMTOOLS_REHEADER_VIA_SAM_RSV_B(ch_rsv_b_for_reheader)
    ch_versions = ch_versions.mix(SAMTOOLS_REHEADER_VIA_SAM_RSV_A.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_REHEADER_VIA_SAM_RSV_B.out.versions)

    // Run BAM_SORT_STATS_SAMTOOLS on RE-HEADERED RSV BAMs
    BAM_SORT_STATS_SAMTOOLS_RSV_A(
        SAMTOOLS_REHEADER_VIA_SAM_RSV_A.out.bam,
        rsv_a_fasta.map { fasta -> [[id: 'PP109421.1'], fasta] },
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_RSV_A.out.versions)

    BAM_SORT_STATS_SAMTOOLS_RSV_B(
        SAMTOOLS_REHEADER_VIA_SAM_RSV_B.out.bam,
        rsv_b_fasta.map { fasta -> [[id: 'OP975389.1'], fasta] },
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_RSV_B.out.versions)

    // Add RSV outputs to channels
    ch_bam = ch_bam.mix(BAM_SORT_STATS_SAMTOOLS_RSV_A.out.bam)
    ch_bam = ch_bam.mix(BAM_SORT_STATS_SAMTOOLS_RSV_B.out.bam)

    ch_bai = ch_bai.mix(BAM_SORT_STATS_SAMTOOLS_RSV_A.out.bai)
    ch_bai = ch_bai.mix(BAM_SORT_STATS_SAMTOOLS_RSV_B.out.bai)

    ch_minimap2_flagstat_multiqc = ch_minimap2_flagstat_multiqc.mix(FASTQ_ALIGN_MINIMAP2_RSV_COMBINED.out.flagstat)
    ch_minimap2_flagstat_multiqc = ch_minimap2_flagstat_multiqc.mix(BAM_SORT_STATS_SAMTOOLS_RSV_A.out.flagstat)
    ch_minimap2_flagstat_multiqc = ch_minimap2_flagstat_multiqc.mix(BAM_SORT_STATS_SAMTOOLS_RSV_B.out.flagstat)

    emit:
    bam                       = ch_bam
    bai                       = ch_bai
    minimap2_multiqc          = ch_minimap2_multiqc
    minimap2_flagstat_multiqc = ch_minimap2_flagstat_multiqc
    versions                  = ch_versions
}
