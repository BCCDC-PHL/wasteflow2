#!/usr/bin/env nextflow

include { BAM_VARIANT_DEMIX_BOOT_FREYJA as BAM_VARIANT_DEMIX_BOOT_FREYJA_SARS_COV2 } from '../nf-core/bam_variant_demix_boot_freyja/main'
include { BAM_VARIANT_DEMIX_BOOT_FREYJA as BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_A     } from '../nf-core/bam_variant_demix_boot_freyja/main'
include { BAM_VARIANT_DEMIX_BOOT_FREYJA as BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_B     } from '../nf-core/bam_variant_demix_boot_freyja/main'

workflow FREYJA_ANALYSIS {
    take:
    ch_bam          // channel: [meta, bam]
    sars_cov2_fasta // channel: fasta
    rsv_a_fasta     // channel: fasta
    rsv_b_fasta     // channel: fasta

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

    // Combine versions from all Freyja processes
    ch_versions = ch_versions.mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_SARS_COV2.out.versions)
    ch_versions = ch_versions.mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_A.out.versions)
    ch_versions = ch_versions.mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_B.out.versions)

    // Combine MultiQC outputs from all Freyja processes
    ch_freyja_multiqc = BAM_VARIANT_DEMIX_BOOT_FREYJA_SARS_COV2.out.demix
        .mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_A.out.demix)
        .mix(BAM_VARIANT_DEMIX_BOOT_FREYJA_RSV_B.out.demix)

    // Organize Freyja files for MultiQC - FIXED VERSION using copyTo
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
    versions           = ch_versions // channel: versions.yml - software versions
}
