#!/usr/bin/env nextflow

include { VARIANTS_IVAR as VARIANTS_IVAR_SARS_COV2 } from './variants_ivar'
include { VARIANTS_IVAR as VARIANTS_IVAR_RSV_A     } from './variants_ivar'
include { VARIANTS_IVAR as VARIANTS_IVAR_RSV_B     } from './variants_ivar'
include { VARIANTS_LONG_TABLE                      } from './variants_long_table'

workflow VARIANT_CALLING {
    take:
    ch_bam                  // channel: [meta, bam]
    sars_cov2_fasta         // channel: fasta
    sars_cov2_fai           // channel: fai
    sars_cov2_chrom_sizes   // channel: chrom_sizes
    sars_cov2_gff           // channel: gff
    sars_cov2_snpeff_db     // channel: snpeff_db
    sars_cov2_snpeff_config // channel: snpeff_config
    rsv_a_fasta             // channel: fasta
    rsv_a_fai               // channel: fai
    rsv_a_chrom_sizes       // channel: chrom_sizes
    rsv_a_gff               // channel: gff
    rsv_a_snpeff_db         // channel: snpeff_db
    rsv_a_snpeff_config     // channel: snpeff_config
    rsv_b_fasta             // channel: fasta
    rsv_b_fai               // channel: fai
    rsv_b_chrom_sizes       // channel: chrom_sizes
    rsv_b_gff               // channel: gff
    rsv_b_snpeff_db         // channel: snpeff_db
    rsv_b_snpeff_config     // channel: snpeff_config

    main:
    ch_versions = Channel.empty()
    sars_cov_2_vcf = Channel.empty()
    sars_cov_2_tbi = Channel.empty()
    rsv_a_vcf = Channel.empty()
    rsv_a_tbi = Channel.empty()
    rsv_b_vcf = Channel.empty()
    rsv_b_tbi = Channel.empty()
    sars_cov_2_ivar_counts_multiqc = Channel.empty()
    rsv_a_ivar_counts_multiqc = Channel.empty()
    rsv_b_ivar_counts_multiqc = Channel.empty()
    sars_cov_2_bcftools_stats_multiqc = Channel.empty()
    rsv_a_bcftools_stats_multiqc = Channel.empty()
    rsv_b_bcftools_stats_multiqc = Channel.empty()
    sars_cov_2_snpeff_multiqc = Channel.empty()
    rsv_a_snpeff_multiqc = Channel.empty()
    rsv_b_snpeff_multiqc = Channel.empty()
    sars_cov_2_snpsift_txt = Channel.empty()
    rsv_a_snpsift_txt = Channel.empty()
    rsv_b_snpsift_txt = Channel.empty()

    ch_ivar_variants_header_mqc = file("${projectDir}/assets/headers/ivar_variants_header_mqc.txt", checkIfExists: true)

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

    // SARS-CoV-2 variant calling
    VARIANTS_IVAR_SARS_COV2(
        bam_by_genome.sars_cov2,
        sars_cov2_fasta,
        sars_cov2_fai,
        sars_cov2_chrom_sizes,
        sars_cov2_gff,
        sars_cov2_snpeff_db,
        sars_cov2_snpeff_config,
        ch_ivar_variants_header_mqc,
    )

    // RSV-A variant calling
    VARIANTS_IVAR_RSV_A(
        bam_by_genome.rsv_a,
        rsv_a_fasta,
        rsv_a_fai,
        rsv_a_chrom_sizes,
        rsv_a_gff,
        rsv_a_snpeff_db,
        rsv_a_snpeff_config,
        ch_ivar_variants_header_mqc,
    )

    // RSV-B variant calling
    VARIANTS_IVAR_RSV_B(
        bam_by_genome.rsv_b,
        rsv_b_fasta,
        rsv_b_fai,
        rsv_b_chrom_sizes,
        rsv_b_gff,
        rsv_b_snpeff_db,
        rsv_b_snpeff_config,
        ch_ivar_variants_header_mqc,
    )

    // Combine versions
    ch_versions = ch_versions.mix(VARIANTS_IVAR_SARS_COV2.out.versions)
    ch_versions = ch_versions.mix(VARIANTS_IVAR_RSV_A.out.versions)
    ch_versions = ch_versions.mix(VARIANTS_IVAR_RSV_B.out.versions)


    sars_cov_2_vcf = VARIANTS_IVAR_SARS_COV2.out.vcf
    // channel: [meta, vcf] - VCF files with variants
    sars_cov_2_tbi = VARIANTS_IVAR_SARS_COV2.out.tbi
    // channel: [meta, tbi] - Tabix index files for VCF
    rsv_a_vcf = VARIANTS_IVAR_RSV_A.out.vcf
    // channel: [meta, vcf] - VCF files with variants
    rsv_a_tbi = VARIANTS_IVAR_RSV_A.out.tbi
    // channel: [meta, tbi] - Tabix index files for VCF
    rsv_b_vcf = VARIANTS_IVAR_RSV_B.out.vcf
    // channel: [meta, vcf] - VCF files with variants
    rsv_b_tbi = VARIANTS_IVAR_RSV_B.out.tbi
    // channel: [meta, tbi] - Tabix
    sars_cov_2_ivar_counts_multiqc = VARIANTS_IVAR_SARS_COV2.out.multiqc_tsv
    // channel: [meta, tsv] - iVar variant counts for MultiQC
    rsv_a_ivar_counts_multiqc = VARIANTS_IVAR_RSV_A.out.multiqc_tsv
    // channel: [meta, tsv] - iVar variant counts for MultiQC
    rsv_b_ivar_counts_multiqc = VARIANTS_IVAR_RSV_B.out.multiqc_tsv
    // channel: [meta, tsv] - iVar variant counts for
    sars_cov_2_bcftools_stats_multiqc = VARIANTS_IVAR_SARS_COV2.out.stats
    // channel: [meta, stats] - BCFtools stats for MultiQC
    rsv_a_bcftools_stats_multiqc = VARIANTS_IVAR_RSV_A.out.stats
    // channel: [meta, stats] - BCFtools stats for MultiQC
    rsv_b_bcftools_stats_multiqc = VARIANTS_IVAR_RSV_B.out.stats
    // channel: [meta, stats] - BCFtools stats for MultiQC
    sars_cov_2_snpeff_multiqc = VARIANTS_IVAR_SARS_COV2.out.snpeff_csv
    // channel: [meta, csv] - SnpEff annotation results for MultiQC
    rsv_a_snpeff_multiqc = VARIANTS_IVAR_RSV_A.out.snpeff_csv
    // channel: [meta, csv] - SnpEff annotation results for MultiQC
    rsv_b_snpeff_multiqc = VARIANTS_IVAR_RSV_B.out.snpeff_csv
    // channel: [meta, csv] - SnpEff annotation results for MultiQC
    sars_cov_2_snpsift_txt = VARIANTS_IVAR_SARS_COV2.out.snpsift_txt
    // channel: [meta, txt] - SnpSift filtered variants
    rsv_a_snpsift_txt = VARIANTS_IVAR_RSV_A.out.snpsift_txt
    // channel: [meta, txt] - SnpSift filtered variants
    rsv_b_snpsift_txt = VARIANTS_IVAR_RSV_B.out.snpsift_txt

    emit:
    sars_cov_2_vcf                    // channel: [meta, vcf] - VCF files with variants
    sars_cov_2_tbi                    // channel: [meta, tbi] - Tabix index files for VCF
    rsv_a_vcf                         // channel: [meta, vcf] - VCF files with variants
    rsv_a_tbi                         // channel: [meta, tbi] - Tabix index files for VCF
    rsv_b_vcf                         // channel: [meta, vcf] - VCF files with variants
    rsv_b_tbi                         // channel: [meta, tbi] - Tabix
    sars_cov_2_ivar_counts_multiqc    // channel: [meta, tsv] - iVar variant counts for MultiQC
    rsv_a_ivar_counts_multiqc         // channel: [meta, tsv] - iVar variant counts for MultiQC
    rsv_b_ivar_counts_multiqc         // channel: [meta, tsv] - iVar variant counts for
    sars_cov_2_bcftools_stats_multiqc // channel: [meta, stats] - BCFtools stats for MultiQC
    rsv_a_bcftools_stats_multiqc      // channel: [meta, stats] - BCFtools stats for MultiQC
    rsv_b_bcftools_stats_multiqc      // channel: [meta, stats] - BCFtools stats for MultiQC
    sars_cov_2_snpeff_multiqc         // channel: [meta, csv] - SnpEff annotation results for MultiQC
    rsv_a_snpeff_multiqc              // channel: [meta, csv] - SnpEff annotation results for MultiQC
    rsv_b_snpeff_multiqc              // channel: [meta, csv] - SnpEff annotation results for MultiQC
    sars_cov_2_snpsift_txt            // channel: [meta, txt] - SnpSift filtered variants
    rsv_a_snpsift_txt                 // channel: [meta, txt] - SnpSift filtered variants
    rsv_b_snpsift_txt                 // channel: [meta, txt] - SnpSift filtered variants
    versions                          = ch_versions // channel: versions.yml - software versions
}
