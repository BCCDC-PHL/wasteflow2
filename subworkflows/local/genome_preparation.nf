#!/usr/bin/env nextflow

include { PREPARE_GENOME as PREPARE_GENOME_SARS_COV2     } from './prepare_genome'
include { PREPARE_GENOME as PREPARE_GENOME_RSV_A         } from './prepare_genome'
include { PREPARE_GENOME as PREPARE_GENOME_RSV_B         } from './prepare_genome'
include { PREPARE_GENOME as PREPARE_GENOME_RSV_COMBINED  } from './prepare_genome'
include { PREPARE_GENOME as PREPARE_GENOME_PANEL_CONTROL } from './prepare_genome'
include { PREPARE_GENOME as PREPARE_GENOME_H1N1         } from './prepare_genome'
include { PREPARE_GENOME as PREPARE_GENOME_H3N2         } from './prepare_genome'
include { PREPARE_GENOME as PREPARE_GENOME_H5N1         } from './prepare_genome'

workflow GENOME_PREPARATION {
    main:
    ch_versions = Channel.empty()


    // SARS-CoV-2 genome preparation (non-segmented, but has BED for genes)
    PREPARE_GENOME_SARS_COV2(
        'MN908947.3',
        params.genomes['MN908947.3'].fasta,
        params.genomes['MN908947.3'].gff,
        null,                                                    // primer_bed
        null,                                                    // bowtie2_index
        null,                                                    // nextclade_dataset
        params.genomes['MN908947.3'].nextclade_dataset_name,
        params.genomes['MN908947.3'].nextclade_dataset_tag,
        null,                                                    // segments_bed (use null for non-segmented)
        null                                                     // segment_gffs (null = non-segmented)
    )

    // RSV-B genome preparation (non-segmented, but has BED for genes)
    PREPARE_GENOME_RSV_B(
        'OP975389.1',
        params.genomes['OP975389.1'].fasta,
        params.genomes['OP975389.1'].gff,
        null,                                                    // primer_bed
        null,                                                    // bowtie2_index
        null,                                                    // nextclade_dataset
        params.genomes['OP975389.1'].nextclade_dataset_name,
        params.genomes['OP975389.1'].nextclade_dataset_tag,
        null,                                                    // segments_bed (use null for non-segmented)
        null                                                     // segment_gffs (null = non-segmented)
    )

    // RSV-A genome preparation (non-segmented, but has BED for genes)
    PREPARE_GENOME_RSV_A(
        'PP109421.1',
        params.genomes['PP109421.1'].fasta,
        params.genomes['PP109421.1'].gff,
        null,                                                    // primer_bed
        null,                                                    // bowtie2_index
        null,                                                    // nextclade_dataset
        params.genomes['PP109421.1'].nextclade_dataset_name,
        params.genomes['PP109421.1'].nextclade_dataset_tag,
        null,                                                    // segments_bed (use null for non-segmented)
        null                                                     // segment_gffs (null = non-segmented)
    )

    // RSV Combined genome preparation
    PREPARE_GENOME_RSV_COMBINED(
        'RSV_combined',
        params.genomes['RSV_combined'].fasta,
        null,
        null,
        null,
        null,
        null,
        null,
        null,
        null
    )

    // Panel Control genome preparation
    PREPARE_GENOME_PANEL_CONTROL(
        'Panel_Control_DB',
        params.genomes['panel_control_db'].fasta,
        null,
        null,
        null,
        null,
        null,
        null,
        null,
        null
    )

    // H1N1 genome preparation (SEGMENTED - has segment_gffs)
    PREPARE_GENOME_H1N1(
        'H1N1',
        params.genomes['H1N1'].fasta,
        null,                                                    // gff (null for segmented)
        null,                                                    // primer_bed
        null,                                                    // bowtie2_index
        null,                                                    // nextclade_dataset
        null,                                                    // nextclade_dataset_name
        null,                                                    // nextclade_dataset_tag
        params.genomes['H1N1'].bed,                             // segments_bed
        params.genomes['H1N1'].segment_gffs                     // segment_gffs map (NOT null = segmented)
    )

    // H3N2 genome preparation (SEGMENTED - has segment_gffs)
    PREPARE_GENOME_H3N2(
        'H3N2',
        params.genomes['H3N2'].fasta,
        null,                                                    // gff (null for segmented)
        null,                                                    // primer_bed
        null,                                                    // bowtie2_index
        null,                                                    // nextclade_dataset
        null,                                                    // nextclade_dataset_name
        null,                                                    // nextclade_dataset_tag
        params.genomes['H3N2'].bed,                             // segments_bed
        params.genomes['H3N2'].segment_gffs                     // segment_gffs map (NOT null = segmented)
    )

    // H5N1 genome preparation (SEGMENTED - has segment_gffs)
    PREPARE_GENOME_H5N1(
        'H5N1',
        params.genomes['H5N1'].fasta,
        null,                                                    // gff (null for segmented)
        null,                                                    // primer_bed
        null,                                                    // bowtie2_index
        null,                                                    // nextclade_dataset
        null,                                                    // nextclade_dataset_name
        null,                                                    // nextclade_dataset_tag
        params.genomes['H5N1'].bed,                             // segments_bed
        params.genomes['H5N1'].segment_gffs                     // segment_gffs map (NOT null = segmented)
    )

    ch_versions = ch_versions.mix(PREPARE_GENOME_SARS_COV2.out.versions)
    ch_versions = ch_versions.mix(PREPARE_GENOME_RSV_A.out.versions)
    ch_versions = ch_versions.mix(PREPARE_GENOME_RSV_B.out.versions)
    ch_versions = ch_versions.mix(PREPARE_GENOME_RSV_COMBINED.out.versions)
    ch_versions = ch_versions.mix(PREPARE_GENOME_PANEL_CONTROL.out.versions)
    ch_versions = ch_versions.mix(PREPARE_GENOME_H1N1.out.versions)
    ch_versions = ch_versions.mix(PREPARE_GENOME_H3N2.out.versions)
    ch_versions = ch_versions.mix(PREPARE_GENOME_H5N1.out.versions)

    emit:
    sars_cov2_fasta                = PREPARE_GENOME_SARS_COV2.out.fasta
    sars_cov2_fai                  = PREPARE_GENOME_SARS_COV2.out.fai
    sars_cov2_gff                  = PREPARE_GENOME_SARS_COV2.out.gff
    sars_cov2_bowtie2_index        = PREPARE_GENOME_SARS_COV2.out.bowtie2_index
    sars_cov2_chrom_sizes          = PREPARE_GENOME_SARS_COV2.out.chrom_sizes
    sars_cov2_snpeff_db            = PREPARE_GENOME_SARS_COV2.out.snpeff_db
    sars_cov2_snpeff_config        = PREPARE_GENOME_SARS_COV2.out.snpeff_config
    rsv_a_fasta                    = PREPARE_GENOME_RSV_A.out.fasta
    rsv_a_fai                      = PREPARE_GENOME_RSV_A.out.fai
    rsv_a_gff                      = PREPARE_GENOME_RSV_A.out.gff
    rsv_a_bowtie2_index            = PREPARE_GENOME_RSV_A.out.bowtie2_index
    rsv_a_chrom_sizes              = PREPARE_GENOME_RSV_A.out.chrom_sizes
    rsv_a_snpeff_db                = PREPARE_GENOME_RSV_A.out.snpeff_db
    rsv_a_snpeff_config            = PREPARE_GENOME_RSV_A.out.snpeff_config
    rsv_b_fasta                    = PREPARE_GENOME_RSV_B.out.fasta
    rsv_b_fai                      = PREPARE_GENOME_RSV_B.out.fai
    rsv_b_gff                      = PREPARE_GENOME_RSV_B.out.gff
    rsv_b_bowtie2_index            = PREPARE_GENOME_RSV_B.out.bowtie2_index
    rsv_b_chrom_sizes              = PREPARE_GENOME_RSV_B.out.chrom_sizes
    rsv_b_snpeff_db                = PREPARE_GENOME_RSV_B.out.snpeff_db
    rsv_b_snpeff_config            = PREPARE_GENOME_RSV_B.out.snpeff_config
    rsv_combined_fasta             = PREPARE_GENOME_RSV_COMBINED.out.fasta
    rsv_combined_bowtie2_index     = PREPARE_GENOME_RSV_COMBINED.out.bowtie2_index
    panel_control_db_fasta         = PREPARE_GENOME_PANEL_CONTROL.out.fasta
    panel_control_db_bowtie2_index = PREPARE_GENOME_PANEL_CONTROL.out.bowtie2_index
    
    // H1N1 outputs
    h1n1_fasta                     = PREPARE_GENOME_H1N1.out.fasta              // Combined for alignment
    h1n1_bowtie2_index             = PREPARE_GENOME_H1N1.out.bowtie2_index
    h1n1_fai                       = PREPARE_GENOME_H1N1.out.fai
    h1n1_chrom_sizes               = PREPARE_GENOME_H1N1.out.chrom_sizes
    h1n1_segment_fasta             = PREPARE_GENOME_H1N1.out.segment_fasta      // Split segments
    h1n1_segment_gff               = PREPARE_GENOME_H1N1.out.segment_gff
    h1n1_segment_fai               = PREPARE_GENOME_H1N1.out.segment_fai
    h1n1_segment_chrom_sizes       = PREPARE_GENOME_H1N1.out.segment_chrom_sizes
    h1n1_segment_snpeff_db         = PREPARE_GENOME_H1N1.out.segment_snpeff_db
    h1n1_segment_snpeff_config     = PREPARE_GENOME_H1N1.out.segment_snpeff_config
    
    // H3N2 outputs
    h3n2_fasta                     = PREPARE_GENOME_H3N2.out.fasta
    h3n2_bowtie2_index             = PREPARE_GENOME_H3N2.out.bowtie2_index
    h3n2_fai                       = PREPARE_GENOME_H3N2.out.fai
    h3n2_chrom_sizes               = PREPARE_GENOME_H3N2.out.chrom_sizes
    h3n2_segment_fasta             = PREPARE_GENOME_H3N2.out.segment_fasta
    h3n2_segment_gff               = PREPARE_GENOME_H3N2.out.segment_gff
    h3n2_segment_fai               = PREPARE_GENOME_H3N2.out.segment_fai
    h3n2_segment_chrom_sizes       = PREPARE_GENOME_H3N2.out.segment_chrom_sizes
    h3n2_segment_snpeff_db         = PREPARE_GENOME_H3N2.out.segment_snpeff_db
    h3n2_segment_snpeff_config     = PREPARE_GENOME_H3N2.out.segment_snpeff_config
    
    // H5N1 outputs
    h5n1_fasta                     = PREPARE_GENOME_H5N1.out.fasta
    h5n1_bowtie2_index             = PREPARE_GENOME_H5N1.out.bowtie2_index
    h5n1_fai                       = PREPARE_GENOME_H5N1.out.fai
    h5n1_chrom_sizes               = PREPARE_GENOME_H5N1.out.chrom_sizes
    h5n1_segment_fasta             = PREPARE_GENOME_H5N1.out.segment_fasta
    h5n1_segment_gff               = PREPARE_GENOME_H5N1.out.segment_gff
    h5n1_segment_fai               = PREPARE_GENOME_H5N1.out.segment_fai
    h5n1_segment_chrom_sizes       = PREPARE_GENOME_H5N1.out.segment_chrom_sizes
    h5n1_segment_snpeff_db         = PREPARE_GENOME_H5N1.out.segment_snpeff_db
    h5n1_segment_snpeff_config     = PREPARE_GENOME_H5N1.out.segment_snpeff_config
    
    versions                       = ch_versions
}