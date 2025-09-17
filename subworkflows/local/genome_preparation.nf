#!/usr/bin/env nextflow

include { PREPARE_GENOME as PREPARE_GENOME_SARS_COV2     } from './prepare_genome'
include { PREPARE_GENOME as PREPARE_GENOME_RSV_A         } from './prepare_genome'
include { PREPARE_GENOME as PREPARE_GENOME_RSV_B         } from './prepare_genome'
include { PREPARE_GENOME as PREPARE_GENOME_RSV_COMBINED  } from './prepare_genome'
include { PREPARE_GENOME as PREPARE_GENOME_PANEL_CONTROL } from './prepare_genome'

workflow GENOME_PREPARATION {
    main:
    ch_versions = Channel.empty()


    // SARS-CoV-2 genome preparation
    PREPARE_GENOME_SARS_COV2(
        params.genomes['MN908947.3'].fasta,
        params.genomes['MN908947.3'].gff,
        null,
        null,
        null,
        params.genomes['MN908947.3'].nextclade_dataset_name,
        params.genomes['MN908947.3'].nextclade_dataset_tag,
    )

    // RSV-B genome preparation
    PREPARE_GENOME_RSV_B(
        params.genomes['OP975389.1'].fasta,
        params.genomes['OP975389.1'].gff,
        null,
        null,
        null,
        params.genomes['OP975389.1'].nextclade_dataset_name,
        params.genomes['OP975389.1'].nextclade_dataset_tag,
    )

    // RSV-A genome preparation
    PREPARE_GENOME_RSV_A(
        params.genomes['PP109421.1'].fasta,
        params.genomes['PP109421.1'].gff,
        null,
        null,
        null,
        params.genomes['PP109421.1'].nextclade_dataset_name,
        params.genomes['PP109421.1'].nextclade_dataset_tag,
    )

    // RSV Combined genome preparation
    PREPARE_GENOME_RSV_COMBINED(
        params.genomes['RSV_combined'].fasta,
        null,
        null,
        null,
        null,
        null,
        null,
    )

    // Panel Control genome preparation
    PREPARE_GENOME_PANEL_CONTROL(
        params.genomes['panel_control_db'].fasta,
        null,
        null,
        null,
        null,
        null,
        null,
    )

    ch_versions = ch_versions.mix(PREPARE_GENOME_SARS_COV2.out.versions)
    ch_versions = ch_versions.mix(PREPARE_GENOME_RSV_A.out.versions)
    ch_versions = ch_versions.mix(PREPARE_GENOME_RSV_B.out.versions)
    ch_versions = ch_versions.mix(PREPARE_GENOME_RSV_COMBINED.out.versions)
    ch_versions = ch_versions.mix(PREPARE_GENOME_PANEL_CONTROL.out.versions)

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
    versions                       = ch_versions
}
