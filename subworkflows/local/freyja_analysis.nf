#!/usr/bin/env nextflow

include { BAM_VARIANT_DEMIX_BOOT_FREYJA as BAM_VARIANT_DEMIX_BOOT_FREYJA_SARS_COV2 } from '../nf-core/bam_variant_demix_boot_freyja/main'
include { BAM_VARIANT_DEMIX_BOOT_FREYJA as BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_A     } from '../nf-core/bam_variant_demix_boot_freyja/main'
include { BAM_VARIANT_DEMIX_BOOT_FREYJA as BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_B     } from '../nf-core/bam_variant_demix_boot_freyja/main'
include { BAM_VARIANT_DEMIX_BOOT_FREYJA as BAM_VARIANT_DEMIX_BOOT_FREYJA_H1N1      } from '../nf-core/bam_variant_demix_boot_freyja/main'
include { BAM_VARIANT_DEMIX_BOOT_FREYJA as BAM_VARIANT_DEMIX_BOOT_FREYJA_H3N2      } from '../nf-core/bam_variant_demix_boot_freyja/main'
include { BAM_VARIANT_DEMIX_BOOT_FREYJA as BAM_VARIANT_DEMIX_BOOT_FREYJA_H5N1      } from '../nf-core/bam_variant_demix_boot_freyja/main'
include { BAM_VARIANT_DEMIX_BOOT_FREYJA as BAM_VARIANT_DEMIX_BOOT_FREYJA_FLU_B_VIC      } from '../nf-core/bam_variant_demix_boot_freyja/main'

workflow FREYJA_ANALYSIS {
    take:
    ch_bam          // channel: [meta, bam]
    sars_cov2_fasta // channel: fasta
    rsv_a_fasta     // channel: fasta
    rsv_b_fasta     // channel: fasta
    h1n1_ha_fasta   // channel: fasta
    h3n2_ha_fasta   // channel: fasta
    h5n1_ha_fasta   // channel: fasta
    flu_b_vic_ha_fasta // channel: fasta

    main:
    ch_versions = Channel.empty()
    ch_freyja_multiqc = Channel.empty()

    // Separate BAM files by genome type for Freyja processing
    ch_bam
        .branch { meta, bam ->
            sars_cov2: meta.genome == 'MN908947.3'
                return [meta, bam]
            rsv_a: meta.genome == 'PP109421.1'
                return [meta, bam]
            rsv_b: meta.genome == 'OP975389.1'
                return [meta, bam]
            h1n1_ha: meta.genome == 'H1N1' && meta.segment == 'HA'
                return [meta, bam]
            h3n2_ha: meta.genome == 'H3N2' && meta.segment == 'HA'
                return [meta, bam]
            h5n1_ha: meta.genome == 'H5N1' && meta.segment == 'HA'
                return [meta, bam]
            flu_b_vic_ha: meta.genome == 'FLU-B-VIC' && meta.segment == 'HA'
                return [meta, bam]
        }
        .set { bam_by_genome }

    // SARS-CoV-2 Freyja Analysis
    BAM_VARIANT_DEMIX_BOOT_FREYJA_SARS_COV2(
        bam_by_genome.sars_cov2,
        sars_cov2_fasta,
        params.skip_freyja_boot,
        params.freyja_repeats,
        params.freyja_db_name_sars_cov2,
        params.freyja_barcodes_sars_cov2,
        params.freyja_lineages_sars_cov2,
    )

    // RSV-A Freyja Analysis
    BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_A(
        bam_by_genome.rsv_a,
        rsv_a_fasta,
        params.skip_freyja_boot,
        params.freyja_repeats,
        params.freyja_db_name_rsv_a,
        params.freyja_barcodes_rsv_a,
        [],
    )

    // RSV-B Freyja Analysis
    BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_B(
        bam_by_genome.rsv_b,
        rsv_b_fasta,
        params.skip_freyja_boot,
        params.freyja_repeats,
        params.freyja_db_name_rsv_b,
        params.freyja_barcodes_rsv_b,
        [],
    )

    // H1N1 HA Freyja Analysis
    BAM_VARIANT_DEMIX_BOOT_FREYJA_H1N1(
        bam_by_genome.h1n1_ha,
        h1n1_ha_fasta,
        params.skip_freyja_boot,
        params.freyja_repeats,
        params.freyja_db_name_h1n1,
        params.freyja_barcodes_h1n1,
        [],
    )

    // H3N2 HA Freyja Analysis
    BAM_VARIANT_DEMIX_BOOT_FREYJA_H3N2(
        bam_by_genome.h3n2_ha,
        h3n2_ha_fasta,
        params.skip_freyja_boot,
        params.freyja_repeats,
        params.freyja_db_name_h3n2,
        params.freyja_barcodes_h3n2,
        [],
    )

    // H5N1 HA Freyja Analysis
    BAM_VARIANT_DEMIX_BOOT_FREYJA_H5N1(
        bam_by_genome.h5n1_ha,
        h5n1_ha_fasta,
        params.skip_freyja_boot,
        params.freyja_repeats,
        params.freyja_db_name_h5n1,
        params.freyja_barcodes_h5n1,
        [],
    )

    BAM_VARIANT_DEMIX_BOOT_FREYJA_FLU_B_VIC(
        bam_by_genome.flu_b_vic_ha,
        flu_b_vic_ha_fasta,
        params.skip_freyja_boot,
        params.freyja_repeats,
        params.freyja_db_name_flu_b_vic,
        params.freyja_barcodes_flu_b_vic,
        [],
    )

    // Combine versions from all Freyja processes
    ch_versions = ch_versions.mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_SARS_COV2.out.versions)
    ch_versions = ch_versions.mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_A.out.versions)
    ch_versions = ch_versions.mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_B.out.versions)
    ch_versions = ch_versions.mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_H1N1.out.versions)
    ch_versions = ch_versions.mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_H3N2.out.versions)
    ch_versions = ch_versions.mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_H5N1.out.versions)

    // Combine MultiQC outputs from all Freyja processes
    ch_freyja_multiqc = BAM_VARIANT_DEMIX_BOOT_FREYJA_SARS_COV2.out.demix
        .mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_A.out.demix)
        .mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_B.out.demix)
        .mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_H1N1.out.demix)
        .mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_H3N2.out.demix)
        .mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_H5N1.out.demix)

    // Organize Freyja files for MultiQC
    ch_freyja_organized = ch_freyja_multiqc
        .map { meta, file ->
            def virus_suffix = ""
            def virus_dir = ""

            if (meta.genome == 'MN908947.3') {
                virus_suffix = "sars_cov2"
                virus_dir = "freyja_sars_cov2"
            }
            else if (meta.genome == 'PP109421.1') {
                virus_suffix = "rsv_a"
                virus_dir = "freyja_rsv_a"
            }
            else if (meta.genome == 'OP975389.1') {
                virus_suffix = "rsv_b"
                virus_dir = "freyja_rsv_b"
            }
            else if (meta.genome == 'H1N1' && meta.segment == 'HA') {
                virus_suffix = "h1n1_ha"
                virus_dir = "freyja_h1n1_ha"
            }
            else if (meta.genome == 'H3N2' && meta.segment == 'HA') {
                virus_suffix = "h3n2_ha"
                virus_dir = "freyja_h3n2_ha"
            }
            else if (meta.genome == 'H5N1' && meta.segment == 'HA') {
                virus_suffix = "h5n1_ha"
                virus_dir = "freyja_h5n1_ha"
            }
            else if (meta.genome == 'FLU-B-VIC' && meta.segment == 'HA') {
                virus_suffix = "flu_b_vic_ha"
                virus_dir = "freyja_flu_b_vic_ha"
            }
            else {
                log.warn("Unknown genome type: ${meta.genome} for sample ${meta.id}")
                virus_suffix = "unknown"
                virus_dir = "freyja_unknown"
            }

            // Create new filename and copy file
            def newFileName = "${meta.id}_${virus_suffix}.tsv"

            // Use copyTo - this automatically creates directories
            return file.copyTo("${virus_dir}/${newFileName}")
        }
        .collect()

    emit:
    freyja_organized   = ch_freyja_organized // channel: [meta, file] - Organized Freyja files for MultiQC
    freyja_multiqc     = ch_freyja_multiqc // channel: [meta, demix] - Freyja demix results for MultiQC
    sars_cov2_demix    = BAM_VARIANT_DEMIX_BOOT_FREYJA_SARS_COV2.out.demix // channel: [meta, demix] - SARS-CoV-2 specific demix results
    sars_cov2_variants = BAM_VARIANT_DEMIX_BOOT_FREYJA_SARS_COV2.out.variants // channel: [meta, variants] - SARS-CoV-2 variants
    rsv_a_demix        = BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_A.out.demix // channel: [meta, demix] - RSV-A specific demix results
    rsv_a_variants     = BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_A.out.variants // channel: [meta, variants] - RSV-A variants
    rsv_b_demix        = BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_B.out.demix // channel: [meta, demix] - RSV-B specific demix results
    rsv_b_variants     = BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_B.out.variants // channel: [meta, variants] - RSV-B variants
    h1n1_demix         = BAM_VARIANT_DEMIX_BOOT_FREYJA_H1N1.out.demix // channel: [meta, demix] - H1N1 HA specific demix results
    h1n1_variants      = BAM_VARIANT_DEMIX_BOOT_FREYJA_H1N1.out.variants // channel: [meta, variants] - H1N1 HA variants
    h3n2_demix         = BAM_VARIANT_DEMIX_BOOT_FREYJA_H3N2.out.demix // channel: [meta, demix] - H3N2 HA specific demix results
    h3n2_variants      = BAM_VARIANT_DEMIX_BOOT_FREYJA_H3N2.out.variants // channel: [meta, variants] - H3N2 HA variants
    h5n1_demix         = BAM_VARIANT_DEMIX_BOOT_FREYJA_H5N1.out.demix // channel: [meta, demix] - H5N1 HA specific demix results
    h5n1_variants      = BAM_VARIANT_DEMIX_BOOT_FREYJA_H5N1.out.variants // channel: [meta, variants] - H5N1 HA variants
    flu_b_vic_demix         = BAM_VARIANT_DEMIX_BOOT_FREYJA_FLU_B_VIC.out.demix // channel: [meta, demix] - Flu B VIC HA specific demix results
    flu_b_vic_variants      = BAM_VARIANT_DEMIX_BOOT_FREYJA_FLU_B_VIC.out.variants // channel: [meta, variants] - Flu B VIC HA variants
    versions           = ch_versions // channel: versions.yml - software versions
}