#!/usr/bin/env nextflow
// Import subworkflows
include { BAM_MARKDUPLICATES_PICARD as BAM_MARKDUPLICATES_PICARD_SARS_COV2         } from '../nf-core/bam_markduplicates_picard'
include { BAM_MARKDUPLICATES_PICARD as BAM_MARKDUPLICATES_PICARD_RSV_A             } from '../nf-core/bam_markduplicates_picard'
include { BAM_MARKDUPLICATES_PICARD as BAM_MARKDUPLICATES_PICARD_RSV_B             } from '../nf-core/bam_markduplicates_picard'

// Import modules
include { PICARD_COLLECTMULTIPLEMETRICS as PICARD_COLLECTMULTIPLEMETRICS_SARS_COV2 } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS as PICARD_COLLECTMULTIPLEMETRICS_RSV_A     } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS as PICARD_COLLECTMULTIPLEMETRICS_RSV_B     } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { MOSDEPTH as MOSDEPTH_GENOME_SARS_COV2                                    } from '../../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_GENOME_RSV_A                                        } from '../../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_GENOME_RSV_B                                        } from '../../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_AMPLICON_SARS_COV2                                  } from '../../modules/nf-core/mosdepth/main'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME_SARS_COV2          } from '../../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME_RSV_A              } from '../../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME_RSV_B              } from '../../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_AMPLICON_SARS_COV2        } from '../../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS_AGG as PLOT_MOSDEPTH_REGIONS_AGG_SARS_COV2         } from '../../modules/local/plot_mosdepth_regions_aggregate'
include { PLOT_MOSDEPTH_REGIONS_AGG as PLOT_MOSDEPTH_REGIONS_AGG_RSV_A             } from '../../modules/local/plot_mosdepth_regions_aggregate'
include { PLOT_MOSDEPTH_REGIONS_AGG as PLOT_MOSDEPTH_REGIONS_AGG_RSV_B             } from '../../modules/local/plot_mosdepth_regions_aggregate'
include { PLOT_MULTIPANEL_COVERAGE_HEATMAP                                         } from '../../modules/local/plot_multilpanel_heatmap'

workflow BAM_QC_METRICS {
    take:
    ch_bam          // channel: [meta, bam]
    ch_bai          // channel: [meta, bai]
    sars_cov2_fasta // channel: fasta
    sars_cov2_fai   // channel: fai
    rsv_a_fasta     // channel: fasta
    rsv_a_fai       // channel: fai
    rsv_b_fasta     // channel: fasta
    rsv_b_fai       // channel: fai

    main:

    ch_versions = Channel.empty()
    ch_markduplicates_flagstat_multiqc = Channel.empty()
    ch_picard_multiqc = Channel.empty()
    ch_mosdepth_multiqc = Channel.empty()
    ch_amplicon_heatmap_multiqc = Channel.empty()

    // Initialize output BAM/BAI channels with input (in case mark duplicates is skipped)
    ch_output_bam = ch_bam
    ch_output_bai = ch_bai

    //
    // SUBWORKFLOW: Mark duplicate reads
    //

    if (!params.skip_markduplicates) {
        // Separate BAM files by genome type for appropriate reference matching
        ch_bam
            .branch { meta, bam ->
                sars_cov2: meta.genome == 'MN908947.3'
                return [meta, bam]
                rsv_a: meta.genome == 'PP109421.1'
                return [meta, bam]
                rsv_b: meta.genome == 'OP975389.1'
                return [meta, bam]
            }
            .set { bam_by_genome }

        sars_cov2_fasta.map { fasta -> [[id: 'MN908947.3'], fasta] }
        rsv_a_fasta.map { fasta -> [[id: 'PP109421.1'], fasta] }
        rsv_b_fasta.map { fasta -> [[id: 'OP975389.1'], fasta] }

        // SARS-CoV-2 mark duplicates
        BAM_MARKDUPLICATES_PICARD_SARS_COV2(
            bam_by_genome.sars_cov2,
            sars_cov2_fasta.map { fasta -> [[id: 'MN908947.3'], fasta] },
            sars_cov2_fai,
        )
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD_SARS_COV2.out.versions)

        // RSV-A mark duplicates
        BAM_MARKDUPLICATES_PICARD_RSV_A(
            bam_by_genome.rsv_a,
            rsv_a_fasta.map { fasta -> [[id: 'PP109421.1'], fasta] },
            rsv_a_fai,
        )
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD_RSV_A.out.versions)

        // RSV-B mark duplicates
        BAM_MARKDUPLICATES_PICARD_RSV_B(
            bam_by_genome.rsv_b,
            rsv_b_fasta.map { fasta -> [[id: 'OP975389.1'], fasta] },
            rsv_b_fai,
        )
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD_RSV_B.out.versions)

        // Combine all mark duplicates outputs
        ch_output_bam = BAM_MARKDUPLICATES_PICARD_SARS_COV2.out.bam
            .mix(BAM_MARKDUPLICATES_PICARD_RSV_A.out.bam)
            .mix(BAM_MARKDUPLICATES_PICARD_RSV_B.out.bam)

        ch_output_bai = BAM_MARKDUPLICATES_PICARD_SARS_COV2.out.bai
            .mix(BAM_MARKDUPLICATES_PICARD_RSV_A.out.bai)
            .mix(BAM_MARKDUPLICATES_PICARD_RSV_B.out.bai)

        ch_markduplicates_flagstat_multiqc = BAM_MARKDUPLICATES_PICARD_SARS_COV2.out.flagstat
            .mix(BAM_MARKDUPLICATES_PICARD_RSV_A.out.flagstat)
            .mix(BAM_MARKDUPLICATES_PICARD_RSV_B.out.flagstat)
    }

    //
    // MODULE: Picard metrics
    //
    if (!params.skip_picard_metrics) {

        // Combine BAM and BAI files
        ch_bam_bai = ch_output_bam.join(ch_output_bai, by: [0])

        // Separate BAM files by genome type for appropriate reference matching
        ch_bam_bai
            .branch { meta, bam, bai ->
                sars_cov2: meta.genome == 'MN908947.3'
                return [meta, bam, bai]
                rsv_a: meta.genome == 'PP109421.1'
                return [meta, bam, bai]
                rsv_b: meta.genome == 'OP975389.1'
                return [meta, bam, bai]
            }
            .set { bam_bai_by_genome }

        // SARS-CoV-2 collect multiple metrics
        PICARD_COLLECTMULTIPLEMETRICS_SARS_COV2(
            bam_bai_by_genome.sars_cov2,
            sars_cov2_fasta.map { fasta -> [[id: 'MN908947.3'], fasta] },
            [[:], []],
        )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS_SARS_COV2.out.versions.first())

        // RSV-A collect multiple metrics
        PICARD_COLLECTMULTIPLEMETRICS_RSV_A(
            bam_bai_by_genome.rsv_a,
            rsv_a_fasta.map { fasta -> [[id: 'PP109421.1'], fasta] },
            [[:], []],
        )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS_RSV_A.out.versions.first())

        // RSV-B collect multiple metrics
        PICARD_COLLECTMULTIPLEMETRICS_RSV_B(
            bam_bai_by_genome.rsv_b,
            rsv_b_fasta.map { fasta -> [[id: 'OP975389.1'], fasta] },
            [[:], []],
        )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS_RSV_B.out.versions.first())

        // Combine Picard outputs for MultiQC
        ch_picard_multiqc = PICARD_COLLECTMULTIPLEMETRICS_SARS_COV2.out.metrics
            .mix(PICARD_COLLECTMULTIPLEMETRICS_RSV_A.out.metrics)
            .mix(PICARD_COLLECTMULTIPLEMETRICS_RSV_B.out.metrics)
    }
    //
    // MODULE: Genome-wide and amplicon-specific coverage QC plots
    //

    if (!params.skip_mosdepth) {
        // Combine BAM and BAI files and separate by genome type
        ch_bam_bai = ch_output_bam.join(ch_output_bai, by: [0])


        ch_bam_bai
            .branch { meta, bam, bai ->
                sars_cov2: meta.genome == 'MN908947.3'
                return [meta, bam, bai, []]
                rsv_a: meta.genome == 'PP109421.1'
                return [meta, bam, bai, []]
                rsv_b: meta.genome == 'OP975389.1'
                return [meta, bam, bai, []]
            }
            .set { bam_bai_by_genome }

        MOSDEPTH_GENOME_SARS_COV2(
            bam_bai_by_genome.sars_cov2,
            sars_cov2_fasta.map { fasta -> [[id: 'MN908947.3'], fasta] },
        )
        ch_versions = ch_versions.mix(MOSDEPTH_GENOME_SARS_COV2.out.versions)

        // RSV-A MOSDEPTH
        MOSDEPTH_GENOME_RSV_A(
            bam_bai_by_genome.rsv_a,
            rsv_a_fasta.map { fasta -> [[id: 'PP109421.1'], fasta] },
        )
        ch_versions = ch_versions.mix(MOSDEPTH_GENOME_RSV_A.out.versions)

        // RSV-B MOSDEPTH
        MOSDEPTH_GENOME_RSV_B(
            bam_bai_by_genome.rsv_b,
            rsv_b_fasta.map { fasta -> [[id: 'OP975389.1'], fasta] },
        )
        ch_versions = ch_versions.mix(MOSDEPTH_GENOME_RSV_B.out.versions)

        // Combine mosdepth outputs for MultiQC
        ch_mosdepth_multiqc = MOSDEPTH_GENOME_SARS_COV2.out.global_txt
            .mix(MOSDEPTH_GENOME_RSV_A.out.global_txt)
            .mix(MOSDEPTH_GENOME_RSV_B.out.global_txt)


        sars_cov_2_coverage = Channel.empty()
        rsva_coverage = Channel.empty()
        rsvb_coverage = Channel.empty()

        PLOT_MOSDEPTH_REGIONS_GENOME_SARS_COV2(
            MOSDEPTH_GENOME_SARS_COV2.out.regions_bed.collect { it[1] }
        )
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_GENOME_SARS_COV2.out.versions)
        sars_cov_2_coverage = PLOT_MOSDEPTH_REGIONS_GENOME_SARS_COV2.out.all_coverage_tsv.map { tsv -> [[id: 'MN908947.3'], tsv] }

        // RSV-A combined plot
        PLOT_MOSDEPTH_REGIONS_GENOME_RSV_A(
            MOSDEPTH_GENOME_RSV_A.out.regions_bed.collect { it[1] }
        )
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_GENOME_RSV_A.out.versions)
        rsva_coverage = PLOT_MOSDEPTH_REGIONS_GENOME_RSV_A.out.all_coverage_tsv.map { tsv -> [[id: 'PP109421.1'], tsv] }

        // RSV-B combined plot
        PLOT_MOSDEPTH_REGIONS_GENOME_RSV_B(
            MOSDEPTH_GENOME_RSV_B.out.regions_bed.collect { it[1] }
        )
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_GENOME_RSV_B.out.versions)
        rsvb_coverage = PLOT_MOSDEPTH_REGIONS_GENOME_RSV_B.out.all_coverage_tsv.map { tsv -> [[id: 'OP975389.1'], tsv] }


        sars_cov_2_coverage_safe = sars_cov_2_coverage.ifEmpty([[id: 'MN908947.3_empty'], []])
        rsva_coverage_safe = rsva_coverage.ifEmpty([[id: 'PP109421.1_empty'], []])
        rsvb_coverage_safe = rsvb_coverage.ifEmpty([[id: 'OP975389.1_empty'], []])

        PLOT_MOSDEPTH_REGIONS_AGG_SARS_COV2(sars_cov_2_coverage, 'sars-cov-2')
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_AGG_SARS_COV2.out.versions)

        PLOT_MOSDEPTH_REGIONS_AGG_RSV_A(rsva_coverage, 'rsva')
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_AGG_RSV_A.out.versions)

        PLOT_MOSDEPTH_REGIONS_AGG_RSV_B(rsvb_coverage, 'rsvb')
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_AGG_RSV_B.out.versions)


        // Multipanel gene level heatmap
        ch_sars_bed = params.sars_cov2_bed ? Channel.fromPath(params.sars_cov2_bed).map { bed -> [[id: 'MN908947.3'], bed] } : Channel.empty()
        ch_rsv_a_bed = params.rsv_a_bed ? Channel.fromPath(params.rsv_a_bed).map { bed -> [[id: 'PP109421.1'], bed] } : Channel.empty()
        ch_rsv_b_bed = params.rsv_b_bed ? Channel.fromPath(params.rsv_b_bed).map { bed -> [[id: 'OP975389.1'], bed] } : Channel.empty()
        ch_metadata = params.metadata ? Channel.fromPath(params.metadata) : Channel.value([])

        PLOT_MULTIPANEL_COVERAGE_HEATMAP(
            sars_cov_2_coverage_safe,
            rsva_coverage_safe,
            rsvb_coverage_safe,
            ch_sars_bed,
            ch_rsv_a_bed,
            ch_rsv_b_bed,
            ch_metadata,
        )
        ch_versions = ch_versions.mix(PLOT_MULTIPANEL_COVERAGE_HEATMAP.out.versions)
    }

    emit:
    bam                             = ch_output_bam // channel: [meta, bam] - processed BAM files (after mark duplicates if enabled)
    bai                             = ch_output_bai // channel: [meta, bai] - processed BAI files (after mark duplicates if enabled)
    markduplicates_flagstat_multiqc = ch_markduplicates_flagstat_multiqc // channel: [meta, flagstat] - flagstat files from mark duplicates for MultiQC
    picard_multiqc                  = ch_picard_multiqc // channel: [meta, metrics] - Picard metrics files for MultiQC
    mosdepth_multiqc                = ch_mosdepth_multiqc // channel: [meta, global_txt] - MOSDEPTH global coverage files for MultiQC
    amplicon_heatmap_multiqc        = ch_amplicon_heatmap_multiqc // channel: heatmap_tsv - amplicon coverage heatmap for MultiQC (SARS-CoV-2 only)
    sars_cov_2_boxplots             = PLOT_MOSDEPTH_REGIONS_AGG_SARS_COV2.out.boxplots // channel: [meta, png] - SARS-CoV-2 genome coverage boxplots
    rsva_boxplots                   = PLOT_MOSDEPTH_REGIONS_AGG_RSV_A.out.boxplots // channel: [meta, png] - RSV-A genome coverage boxplots
    rsvb_boxplots                   = PLOT_MOSDEPTH_REGIONS_AGG_RSV_B.out.boxplots // channel: [meta, png] - RSV-B genome coverage box
    multipanel_coverage_heatmap     = PLOT_MULTIPANEL_COVERAGE_HEATMAP.out.heatmap // channel: [meta, png] - multi-panel coverage heatmap
    versions                        = ch_versions // channel: versions.yml - software versions
}
