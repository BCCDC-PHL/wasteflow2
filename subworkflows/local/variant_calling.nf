#!/usr/bin/env nextflow

include { VARIANTS_IVAR as VARIANTS_IVAR_SARS_COV2 } from './variants_ivar'
include { VARIANTS_IVAR as VARIANTS_IVAR_RSV_A     } from './variants_ivar'
include { VARIANTS_IVAR as VARIANTS_IVAR_RSV_B     } from './variants_ivar'
include { VARIANTS_IVAR as VARIANTS_IVAR_INFLUENZA_SEGMENTS } from './variants_ivar'
include { VARIANTS_LONG_TABLE                      } from './variants_long_table'

workflow VARIANT_CALLING {
    take:
    ch_bam                  // channel: [meta, bam]
    sars_cov2_fasta         // channel: fasta
    sars_cov2_fai           // channel: fai
    sars_cov2_chrom_sizes   // channel: chrom_sizes
    sars_cov2_gff           // channel: gff
    rsv_a_fasta             // channel: fasta
    rsv_a_fai               // channel: fai
    rsv_a_chrom_sizes       // channel: chrom_sizes
    rsv_a_gff               // channel: gff
    rsv_b_fasta             // channel: fasta
    rsv_b_fai               // channel: fai
    rsv_b_chrom_sizes       // channel: chrom_sizes
    rsv_b_gff               // channel: gff
    ch_segment_fasta        // channel: fasta for all influenza segments
    ch_segment_fai          // channel: fai for all influenza segments
    ch_segment_chrom_sizes  // channel: chrom_sizes for all influenza segments
    ch_segment_gff          // channel: gff for all influenza segments
    

    main:
    ch_versions = channel.empty()

    ch_ivar_variants_header_mqc = file("${projectDir}/assets/headers/ivar_variants_header_mqc.txt", checkIfExists: true)

    // Split BAMs by virus type for processing
    ch_bam
        .branch { meta, bam ->
            sars_cov2: meta.genome == 'MN908947.3'
                return [meta, bam]
            rsv_a: meta.genome == 'PP109421.1'
                return [meta, bam]
            rsv_b: meta.genome == 'OP975389.1'
                return [meta, bam]
            segments: meta.genome in ['H1N1', 'H3N2', 'H5N1'] && meta.containsKey('segment')
                return [meta, bam]
        }
        .set { bam_by_type }

    // SARS-CoV-2 variant calling
    VARIANTS_IVAR_SARS_COV2(
        bam_by_type.sars_cov2,
        sars_cov2_fasta,
        sars_cov2_fai,
        sars_cov2_chrom_sizes,
        sars_cov2_gff,
        ch_ivar_variants_header_mqc,
    )

    // RSV-A variant calling
    VARIANTS_IVAR_RSV_A(
        bam_by_type.rsv_a,
        rsv_a_fasta,
        rsv_a_fai,
        rsv_a_chrom_sizes,
        rsv_a_gff,
        ch_ivar_variants_header_mqc,
    )

    // RSV-B variant calling
    VARIANTS_IVAR_RSV_B(
        bam_by_type.rsv_b,
        rsv_b_fasta,
        rsv_b_fai,
        rsv_b_chrom_sizes,
        rsv_b_gff,
        ch_ivar_variants_header_mqc,
    )

    // Influenza segments variant calling
    // Join segment BAMs with their corresponding references
    bam_by_type.segments
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
        .combine(
            ch_segment_chrom_sizes
                .map { ref_meta, chrom_sizes ->
                    def key = "${ref_meta.virus}_${ref_meta.segment}"
                    [key, chrom_sizes]
                },
            by: 0
        )
        .combine(
            ch_segment_gff
                .map { ref_meta, gff ->
                    def key = "${ref_meta.virus}_${ref_meta.segment}"
                    [key, gff]
                },
            by: 0
        )
        .map { key, meta, bam, ref_meta, fasta, fai, chrom_sizes, gff ->
            [meta, bam, [ref_meta, fasta], fai, chrom_sizes, gff]
        }
        .set { ch_segments_with_refs }
        
    // Enrich meta with matched fasta and fai
    ch_segments_with_refs
        .map { meta, bam, ref_meta_fasta, fai, chrom_sizes, gff ->
            // Store fasta in meta for later retrieval
            def enriched_meta = meta.clone()
            enriched_meta.matched_fasta = ref_meta_fasta[1]
            enriched_meta.matched_fai = fai
            
            [enriched_meta, bam, ref_meta_fasta[1], fai, chrom_sizes, gff]
        }
        .set { ch_segments_enriched }

    // Split for process inputs
    ch_segments_enriched
        .multiMap { meta, bam, fasta, fai, chrom_sizes, gff ->
            bams: [meta, bam]
            fastas: fasta
            fais: fai
            chrom_sizes: chrom_sizes
            gffs: gff
        }
        .set { ch_segments_inputs }

    VARIANTS_IVAR_INFLUENZA_SEGMENTS(
        ch_segments_inputs.bams,
        ch_segments_inputs.fastas,
        ch_segments_inputs.fais,
        ch_segments_inputs.chrom_sizes,
        ch_segments_inputs.gffs,
        ch_ivar_variants_header_mqc
    )

    // Combine versions
    ch_versions = ch_versions.mix(VARIANTS_IVAR_SARS_COV2.out.versions)
    ch_versions = ch_versions.mix(VARIANTS_IVAR_RSV_A.out.versions)
    ch_versions = ch_versions.mix(VARIANTS_IVAR_RSV_B.out.versions)
    ch_versions = ch_versions.mix(VARIANTS_IVAR_INFLUENZA_SEGMENTS.out.versions)

    // Mix all VCFs from different virus types
    def ch_all_vcf = VARIANTS_IVAR_SARS_COV2.out.vcf
        .mix(VARIANTS_IVAR_RSV_A.out.vcf)
        .mix(VARIANTS_IVAR_RSV_B.out.vcf)
        .mix(VARIANTS_IVAR_INFLUENZA_SEGMENTS.out.vcf)

    // Mix all TBIs from different virus types
    def ch_all_tbi = VARIANTS_IVAR_SARS_COV2.out.tbi
        .mix(VARIANTS_IVAR_RSV_A.out.tbi)
        .mix(VARIANTS_IVAR_RSV_B.out.tbi)
        .mix(VARIANTS_IVAR_INFLUENZA_SEGMENTS.out.tbi)

    // Mix all iVar counts for MultiQC
    def ch_all_ivar_counts_multiqc = VARIANTS_IVAR_SARS_COV2.out.multiqc_tsv
        .mix(VARIANTS_IVAR_RSV_A.out.multiqc_tsv)
        .mix(VARIANTS_IVAR_RSV_B.out.multiqc_tsv)
        .mix(VARIANTS_IVAR_INFLUENZA_SEGMENTS.out.multiqc_tsv)

    // Mix all BCFtools stats for MultiQC
    def ch_all_bcftools_stats_multiqc = VARIANTS_IVAR_SARS_COV2.out.stats
        .mix(VARIANTS_IVAR_RSV_A.out.stats)
        .mix(VARIANTS_IVAR_RSV_B.out.stats)
        .mix(VARIANTS_IVAR_INFLUENZA_SEGMENTS.out.stats)

    emit:
    vcf                 = ch_all_vcf                    // channel: [meta, vcf] - All VCF files with variants
    tbi                 = ch_all_tbi                    // channel: [meta, tbi] - All Tabix index files
    ivar_counts_multiqc = ch_all_ivar_counts_multiqc    // channel: [meta, tsv] - All iVar variant counts for MultiQC
    bcftools_stats      = ch_all_bcftools_stats_multiqc // channel: [meta, stats] - All BCFtools stats for MultiQC
    versions            = ch_versions                   // channel: versions.yml - software versions
}