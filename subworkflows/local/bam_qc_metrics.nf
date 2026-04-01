#!/usr/bin/env nextflow

// Import subworkflows
include { BAM_MARKDUPLICATES_PICARD as BAM_MARKDUPLICATES_PICARD_SARS_COV2 } from '../nf-core/bam_markduplicates_picard'
include { BAM_MARKDUPLICATES_PICARD as BAM_MARKDUPLICATES_PICARD_RSV_A     } from '../nf-core/bam_markduplicates_picard'
include { BAM_MARKDUPLICATES_PICARD as BAM_MARKDUPLICATES_PICARD_RSV_B     } from '../nf-core/bam_markduplicates_picard'
include { BAM_MARKDUPLICATES_PICARD as BAM_MARKDUPLICATES_PICARD_SEGMENTS  } from '../nf-core/bam_markduplicates_picard'

// Import modules
include { PICARD_COLLECTMULTIPLEMETRICS as PICARD_COLLECTMULTIPLEMETRICS_SARS_COV2 } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS as PICARD_COLLECTMULTIPLEMETRICS_RSV_A     } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS as PICARD_COLLECTMULTIPLEMETRICS_RSV_B     } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS as PICARD_COLLECTMULTIPLEMETRICS_SEGMENTS  } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { MOSDEPTH as MOSDEPTH_GENOME_SARS_COV2    } from '../../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_GENOME_RSV_A        } from '../../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_GENOME_RSV_B        } from '../../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_SEGMENTS            } from '../../modules/nf-core/mosdepth/main'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME_SARS_COV2 } from '../../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME_RSV_A     } from '../../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME_RSV_B     } from '../../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_SEGMENTS_H1N1         } from '../../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_SEGMENTS_H3N2         } from '../../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_SEGMENTS_H5N1         } from '../../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_SEGMENTS_FLU_B_VIC         } from '../../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS_AGG as PLOT_MOSDEPTH_REGIONS_AGG_SARS_COV2 } from '../../modules/local/plot_mosdepth_regions_aggregate'
include { PLOT_MOSDEPTH_REGIONS_AGG as PLOT_MOSDEPTH_REGIONS_AGG_RSV_A     } from '../../modules/local/plot_mosdepth_regions_aggregate'
include { PLOT_MOSDEPTH_REGIONS_AGG as PLOT_MOSDEPTH_REGIONS_AGG_RSV_B     } from '../../modules/local/plot_mosdepth_regions_aggregate'
include { PLOT_MULTIPANEL_COVERAGE_HEATMAP } from '../../modules/local/plot_multilpanel_heatmap'

workflow BAM_QC_METRICS {
    take:
    ch_bam                 // channel: [meta, bam] - includes both full-genome and segment BAMs
    ch_bai                 // channel: [meta, bai]
    sars_cov2_fasta        // channel: fasta
    sars_cov2_fai          // channel: fai
    rsv_a_fasta            // channel: fasta
    rsv_a_fai              // channel: fai
    rsv_b_fasta            // channel: fasta
    rsv_b_fai              // channel: fai
    ch_segment_fasta       // channel: [meta[virus, segment, id], fasta] - all segments
    ch_segment_fai         // channel: [meta[virus, segment, id], fai]

    main:

    ch_versions = Channel.empty()
    ch_markduplicates_flagstat_multiqc = Channel.empty()
    ch_picard_multiqc = Channel.empty()
    ch_mosdepth_multiqc = Channel.empty()
    ch_amplicon_heatmap_multiqc = Channel.empty()

    // Initialize output BAM/BAI channels
    ch_output_bam = ch_bam
    ch_output_bai = ch_bai

    //
    // Branch BAMs by type: full-genome vs segments
    //
    ch_bam
        .branch { meta, bam ->
            sars_cov2: meta.genome == 'MN908947.3'
                return [meta, bam]
            rsv_a: meta.genome == 'PP109421.1'
                return [meta, bam]
            rsv_b: meta.genome == 'OP975389.1'
                return [meta, bam]
            flu_segments: meta.genome in ['H1N1', 'H3N2', 'H5N1', 'FLU-B-VIC'] && meta.containsKey('segment')
                return [meta, bam]
            
        }
        .set { bam_by_type }

    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    if (!params.skip_markduplicates) {
        
        // SARS-CoV-2 mark duplicates
        BAM_MARKDUPLICATES_PICARD_SARS_COV2(
            bam_by_type.sars_cov2,
            sars_cov2_fasta.map { fasta -> [[id: 'MN908947.3'], fasta] },
            sars_cov2_fai
        )
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD_SARS_COV2.out.versions)

        // RSV-A mark duplicates
        BAM_MARKDUPLICATES_PICARD_RSV_A(
            bam_by_type.rsv_a,
            rsv_a_fasta.map { fasta -> [[id: 'PP109421.1'], fasta] },
            rsv_a_fai
        )
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD_RSV_A.out.versions)

        // RSV-B mark duplicates
        BAM_MARKDUPLICATES_PICARD_RSV_B(
            bam_by_type.rsv_b,
            rsv_b_fasta.map { fasta -> [[id: 'OP975389.1'], fasta] },
            rsv_b_fai
        )
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD_RSV_B.out.versions)

        // Influenza segments mark duplicates
        // Join segment BAMs with their corresponding references
        bam_by_type.flu_segments
            .map { meta, bam -> 
                // Create join key: virus + segment
                def key = "${meta.genome}_${meta.segment}"
                [key, meta, bam]
            }
            .combine(
                ch_segment_fasta
                    .map { ref_meta, fasta ->
                        def key = "${ref_meta.virus}_${ref_meta.segment}"
                        [key, ref_meta, fasta]
                    },
                by: 0
            )
            .combine(
                ch_segment_fai
                    .map { ref_meta, fai ->
                        def key = "${ref_meta.virus}_${ref_meta.segment}"
                        [key, fai]
                    },
                by: 0
            )
            .map { key, meta, bam, ref_meta, fasta, fai ->
                [meta, bam, [ref_meta, fasta], fai]
            }
            .set { ch_segments_with_refs }
        
        BAM_MARKDUPLICATES_PICARD_SEGMENTS(
            ch_segments_with_refs.map { meta, bam, fasta_tuple, fai -> [meta, bam] },
            ch_segments_with_refs.map { meta, bam, fasta_tuple, fai -> fasta_tuple },
            ch_segments_with_refs.map { meta, bam, fasta_tuple, fai -> fai }
        )
        
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD_SEGMENTS.out.versions)
        
        // Combine all mark duplicates outputs
        ch_output_bam = BAM_MARKDUPLICATES_PICARD_SARS_COV2.out.bam
            .mix(BAM_MARKDUPLICATES_PICARD_RSV_A.out.bam)
            .mix(BAM_MARKDUPLICATES_PICARD_RSV_B.out.bam)
            .mix(BAM_MARKDUPLICATES_PICARD_SEGMENTS.out.bam)

        ch_output_bai = BAM_MARKDUPLICATES_PICARD_SARS_COV2.out.bai
            .mix(BAM_MARKDUPLICATES_PICARD_RSV_A.out.bai)
            .mix(BAM_MARKDUPLICATES_PICARD_RSV_B.out.bai)
            .mix(BAM_MARKDUPLICATES_PICARD_SEGMENTS.out.bai)

        ch_markduplicates_flagstat_multiqc = BAM_MARKDUPLICATES_PICARD_SARS_COV2.out.flagstat
            .mix(BAM_MARKDUPLICATES_PICARD_RSV_A.out.flagstat)
            .mix(BAM_MARKDUPLICATES_PICARD_RSV_B.out.flagstat)
            .mix(BAM_MARKDUPLICATES_PICARD_SEGMENTS.out.flagstat)
    }

    //
    // MODULE: Picard metrics
    //
    if (!params.skip_picard_metrics) {
        
        // Combine BAM and BAI files
        ch_bam_bai = ch_output_bam.join(ch_output_bai, by: [0])

        // Branch by type
        ch_bam_bai
            .branch { meta, bam, bai ->
                sars_cov2: meta.genome == 'MN908947.3'
                    return [meta, bam, bai]
                rsv_a: meta.genome == 'PP109421.1'
                    return [meta, bam, bai]
                rsv_b: meta.genome == 'OP975389.1'
                    return [meta, bam, bai]
                flu_segments: meta.genome in ['H1N1', 'H3N2', 'H5N1', 'FLU-B-VIC'] && meta.containsKey('segment')
                    return [meta, bam, bai]
            }
            .set { bam_bai_by_type }

        // SARS-CoV-2 Picard metrics
        PICARD_COLLECTMULTIPLEMETRICS_SARS_COV2(
            bam_bai_by_type.sars_cov2,
            sars_cov2_fasta.map { fasta -> [[id: 'MN908947.3'], fasta] },
            [[:], []]
        )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS_SARS_COV2.out.versions.first())

        // RSV-A Picard metrics
        PICARD_COLLECTMULTIPLEMETRICS_RSV_A(
            bam_bai_by_type.rsv_a,
            rsv_a_fasta.map { fasta -> [[id: 'PP109421.1'], fasta] },
            [[:], []]
        )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS_RSV_A.out.versions.first())

        // RSV-B Picard metrics
        PICARD_COLLECTMULTIPLEMETRICS_RSV_B(
            bam_bai_by_type.rsv_b,
            rsv_b_fasta.map { fasta -> [[id: 'OP975389.1'], fasta] },
            [[:], []]
        )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS_RSV_B.out.versions.first())

        // Segments Picard metrics - join with references
        bam_bai_by_type.flu_segments
            .map { meta, bam, bai -> 
                def key = "${meta.genome}_${meta.segment}"
                [key, meta, bam, bai]
            }
            .combine(
                ch_segment_fasta
                    .map { ref_meta, fasta ->
                        def key = "${ref_meta.virus}_${ref_meta.segment}"
                        [key, ref_meta, fasta]
                    },
                by: 0
            )
            .map { key, meta, bam, bai, ref_meta, fasta ->
                [[meta, bam, bai], [ref_meta, fasta]]
            }
            .set { ch_segments_bam_bai_with_refs }

        PICARD_COLLECTMULTIPLEMETRICS_SEGMENTS(
            ch_segments_bam_bai_with_refs.map { bam_tuple, fasta_tuple -> bam_tuple },
            ch_segments_bam_bai_with_refs.map { bam_tuple, fasta_tuple -> fasta_tuple },
            [[:], []]
        )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS_SEGMENTS.out.versions.first())
        
        // Combine Picard outputs
        ch_picard_multiqc = PICARD_COLLECTMULTIPLEMETRICS_SARS_COV2.out.metrics
            .mix(PICARD_COLLECTMULTIPLEMETRICS_RSV_A.out.metrics)
            .mix(PICARD_COLLECTMULTIPLEMETRICS_RSV_B.out.metrics)
            .mix(PICARD_COLLECTMULTIPLEMETRICS_SEGMENTS.out.metrics)
    }

    //
    // MODULE: Coverage QC with MOSDEPTH
    //
    if (!params.skip_mosdepth) {
        
        // Combine BAM and BAI, branch by type
        ch_bam_bai = ch_output_bam.join(ch_output_bai, by: [0])
        ch_bam_bai
            .branch { meta, bam, bai ->
                sars_cov2: meta.genome == 'MN908947.3'
                    return [meta, bam, bai, []]
                rsv_a: meta.genome == 'PP109421.1'
                    return [meta, bam, bai, []]
                rsv_b: meta.genome == 'OP975389.1'
                    return [meta, bam, bai, []]
                flu_segments: meta.genome in ['H1N1', 'H3N2', 'H5N1', 'FLU-B-VIC'] && meta.containsKey('segment')
                    return [meta, bam, bai, []]
            }
            .set { bam_bai_by_type }

        // SARS-CoV-2 MOSDEPTH
        MOSDEPTH_GENOME_SARS_COV2(
            bam_bai_by_type.sars_cov2,
            sars_cov2_fasta.map { fasta -> [[id: 'MN908947.3'], fasta] }
        )
        ch_versions = ch_versions.mix(MOSDEPTH_GENOME_SARS_COV2.out.versions)

        // RSV-A MOSDEPTH
        MOSDEPTH_GENOME_RSV_A(
            bam_bai_by_type.rsv_a,
            rsv_a_fasta.map { fasta -> [[id: 'PP109421.1'], fasta] }
        )
        ch_versions = ch_versions.mix(MOSDEPTH_GENOME_RSV_A.out.versions)

        // RSV-B MOSDEPTH
        MOSDEPTH_GENOME_RSV_B(
            bam_bai_by_type.rsv_b,
            rsv_b_fasta.map { fasta -> [[id: 'OP975389.1'], fasta] }
        )
        ch_versions = ch_versions.mix(MOSDEPTH_GENOME_RSV_B.out.versions)

        // Segments MOSDEPTH - join with references
        bam_bai_by_type.flu_segments
            .map { meta, bam, bai, bed -> 
                def key = "${meta.genome}_${meta.segment}"
                [key, meta, bam, bai, bed]
            }
            .combine(
                ch_segment_fasta
                    .map { ref_meta, fasta ->
                        def key = "${ref_meta.virus}_${ref_meta.segment}"
                        [key, ref_meta, fasta]
                    },
                by: 0
            )
            .map { key, meta, bam, bai, bed, ref_meta, fasta ->
                [[meta, bam, bai, bed], [ref_meta, fasta]]
            }
            .set { ch_segments_bam_bai_bed_with_refs }

        MOSDEPTH_SEGMENTS(
            ch_segments_bam_bai_bed_with_refs.map { bam_tuple, fasta_tuple -> bam_tuple },
            ch_segments_bam_bai_bed_with_refs.map { bam_tuple, fasta_tuple -> fasta_tuple }
        )
        ch_versions = ch_versions.mix(MOSDEPTH_SEGMENTS.out.versions)

        // Combine mosdepth outputs
        ch_mosdepth_multiqc = MOSDEPTH_GENOME_SARS_COV2.out.global_txt
            .mix(MOSDEPTH_GENOME_RSV_A.out.global_txt)
            .mix(MOSDEPTH_GENOME_RSV_B.out.global_txt)
            .mix(MOSDEPTH_SEGMENTS.out.global_txt)

        // Coverage plots
        sars_cov_2_coverage = Channel.empty()
        rsva_coverage = Channel.empty()
        rsvb_coverage = Channel.empty()
        h1n1_coverage = Channel.empty()
        h3n2_coverage = Channel.empty()
        h5n1_coverage = Channel.empty()

        PLOT_MOSDEPTH_REGIONS_GENOME_SARS_COV2(
            MOSDEPTH_GENOME_SARS_COV2.out.regions_bed.collect { it[1] }
        )
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_GENOME_SARS_COV2.out.versions)
        sars_cov_2_coverage = PLOT_MOSDEPTH_REGIONS_GENOME_SARS_COV2.out.all_coverage_tsv.map { tsv -> [[id: 'MN908947.3'], tsv] }

        PLOT_MOSDEPTH_REGIONS_GENOME_RSV_A(
            MOSDEPTH_GENOME_RSV_A.out.regions_bed.collect { it[1] }
        )
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_GENOME_RSV_A.out.versions)
        rsva_coverage = PLOT_MOSDEPTH_REGIONS_GENOME_RSV_A.out.all_coverage_tsv.map { tsv -> [[id: 'PP109421.1'], tsv] }

        PLOT_MOSDEPTH_REGIONS_GENOME_RSV_B(
            MOSDEPTH_GENOME_RSV_B.out.regions_bed.collect { it[1] }
        )
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_GENOME_RSV_B.out.versions)
        rsvb_coverage = PLOT_MOSDEPTH_REGIONS_GENOME_RSV_B.out.all_coverage_tsv.map { tsv -> [[id: 'OP975389.1'], tsv] }

        //
        // Filter MOSDEPTH outputs by genome type
        //
        MOSDEPTH_SEGMENTS.out.regions_bed
            .branch { meta, bed ->
                h1n1: meta.genome =~ /(?i)H1N1/
                h3n2: meta.genome =~ /(?i)H3N2/
                h5n1: meta.genome =~ /(?i)H5N1/
                flu_b_vic: meta.genome =~ /(?i)FLU-B-VIC/
            }
            .set { ch_mosdepth_branched }


        // Segments coverage plot
        PLOT_MOSDEPTH_REGIONS_SEGMENTS_H1N1(
            ch_mosdepth_branched.h1n1.collect { it[1] }
        )
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_SEGMENTS_H1N1.out.versions)

        PLOT_MOSDEPTH_REGIONS_SEGMENTS_H3N2(
            ch_mosdepth_branched.h3n2.collect { it[1] }
        )
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_SEGMENTS_H3N2.out.versions)

        PLOT_MOSDEPTH_REGIONS_SEGMENTS_H5N1(
            ch_mosdepth_branched.h5n1.collect { it[1] }
        )
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_SEGMENTS_H5N1.out.versions)

        PLOT_MOSDEPTH_REGIONS_SEGMENTS_FLU_B_VIC(
            ch_mosdepth_branched.flu_b_vic.collect { it[1] }
        )
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_SEGMENTS_FLU_B_VIC.out.versions)

        h1n1_coverage = PLOT_MOSDEPTH_REGIONS_SEGMENTS_H1N1.out.all_coverage_tsv.map { tsv -> [[id: 'H1N1'], tsv] }
        h3n2_coverage = PLOT_MOSDEPTH_REGIONS_SEGMENTS_H3N2.out.all_coverage_tsv.map { tsv -> [[id: 'H3N2'], tsv] }
        h5n1_coverage = PLOT_MOSDEPTH_REGIONS_SEGMENTS_H5N1.out.all_coverage_tsv.map { tsv -> [[id: 'H5N1'], tsv] }
        flu_b_vic_coverage = PLOT_MOSDEPTH_REGIONS_SEGMENTS_FLU_B_VIC.out.all_coverage_tsv.map { tsv -> [[id: 'FLU-B-VIC'], tsv] }

        // Aggregate plots
        sars_cov_2_coverage_safe = sars_cov_2_coverage.ifEmpty([[id: 'MN908947.3_empty'], []])
        rsva_coverage_safe = rsva_coverage.ifEmpty([[id: 'PP109421.1_empty'], []])
        rsvb_coverage_safe = rsvb_coverage.ifEmpty([[id: 'OP975389.1_empty'], []])
        h1n1_coverage_safe = h1n1_coverage.ifEmpty([[id: 'H1N1_empty'], []])
        h3n2_coverage_safe = h3n2_coverage.ifEmpty([[id: 'H3N2_empty'], []])
        h5n1_coverage_safe = h5n1_coverage.ifEmpty([[id: 'H5N1_empty'], []])
        flu_b_vic_coverage_safe = flu_b_vic_coverage.ifEmpty([[id: 'FLU-B-VIC_empty'], []])

        PLOT_MOSDEPTH_REGIONS_AGG_SARS_COV2(sars_cov_2_coverage, 'sars-cov-2')
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_AGG_SARS_COV2.out.versions)

        PLOT_MOSDEPTH_REGIONS_AGG_RSV_A(rsva_coverage, 'rsva')
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_AGG_RSV_A.out.versions)

        PLOT_MOSDEPTH_REGIONS_AGG_RSV_B(rsvb_coverage, 'rsvb')
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_AGG_RSV_B.out.versions)

        // Multipanel heatmap
        ch_sars_bed = params.sars_cov2_bed ? Channel.fromPath(params.sars_cov2_bed).map { bed -> [[id: 'MN908947.3'], bed] } : Channel.empty()
        ch_rsv_a_bed = params.rsv_a_bed ? Channel.fromPath(params.rsv_a_bed).map { bed -> [[id: 'PP109421.1'], bed] } : Channel.empty()
        ch_rsv_b_bed = params.rsv_b_bed ? Channel.fromPath(params.rsv_b_bed).map { bed -> [[id: 'OP975389.1'], bed] } : Channel.empty()
        ch_h1n1_bed = params.genomes['H1N1'].bed ? Channel.fromPath(params.genomes['H1N1'].bed ).map { bed -> [[id: 'H1N1'], bed] } : Channel.empty() 
        ch_h3n2_bed = params.genomes['H3N2'].bed ? Channel.fromPath(params.genomes['H3N2'].bed ).map { bed -> [[id: 'H3N2'], bed] } : Channel.empty()
        ch_h5n1_bed = params.genomes['H5N1'].bed ? Channel.fromPath(params.genomes['H5N1'].bed ).map { bed -> [[id: 'H5N1'], bed] } : Channel.empty()
        ch_flu_b_vic_bed = params.genomes['FLU-B-VIC'].bed ? Channel.fromPath(params.genomes['FLU-B-VIC'].bed ).map { bed -> [[id: 'FLU-B-VIC'], bed] } : Channel.empty()
        ch_metadata = params.metadata ? Channel.fromPath(params.metadata) : Channel.value([])
        
        PLOT_MULTIPANEL_COVERAGE_HEATMAP(
            sars_cov_2_coverage_safe,
            rsva_coverage_safe,
            rsvb_coverage_safe,
            h1n1_coverage_safe,
            h3n2_coverage_safe,
            h5n1_coverage_safe,
            flu_b_vic_coverage_safe,
            ch_sars_bed,
            ch_rsv_a_bed,
            ch_rsv_b_bed,
            ch_h1n1_bed,
            ch_h3n2_bed,
            ch_h5n1_bed,
            ch_flu_b_vic_bed,
            ch_metadata
        )
        ch_versions = ch_versions.mix(PLOT_MULTIPANEL_COVERAGE_HEATMAP.out.versions)
        
    }

    emit:
    bam                             = ch_output_bam
    bai                             = ch_output_bai
    markduplicates_flagstat_multiqc = ch_markduplicates_flagstat_multiqc
    picard_multiqc                  = ch_picard_multiqc
    mosdepth_multiqc                = ch_mosdepth_multiqc
    amplicon_heatmap_multiqc        = ch_amplicon_heatmap_multiqc
    sars_cov_2_boxplots             = PLOT_MOSDEPTH_REGIONS_AGG_SARS_COV2.out.boxplots
    rsva_boxplots                   = PLOT_MOSDEPTH_REGIONS_AGG_RSV_A.out.boxplots
    rsvb_boxplots                   = PLOT_MOSDEPTH_REGIONS_AGG_RSV_B.out.boxplots
    //plotsars_rsv_multipanel_coverage_heatmap = PLOT_MULTIPANEL_COVERAGE_HEATMAP.out.heatmap_sars
    //influenza_multipanel_coverage_heatmap = PLOT_MULTIPANEL_COVERAGE_HEATMAP.out.heatmap_influenza
    versions                        = ch_versions
}