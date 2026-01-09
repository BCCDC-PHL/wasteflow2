#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog        } from 'plugin/nf-schema'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_wasteflow2_pipeline'

// Import subworkflows
include { GENOME_PREPARATION      } from '../subworkflows/local/genome_preparation'
include { READ_PROCESSING         } from '../subworkflows/local/read_processing'
include { TAXONOMY_CLASSIFICATION } from '../subworkflows/local/taxonomy_classification'
include { ALIGNMENT               } from '../subworkflows/local/alignment_processing'
include { BAM_QC_METRICS          } from '../subworkflows/local/bam_qc_metrics'
include { VARIANT_CALLING         } from '../subworkflows/local/variant_calling'
include { FREYJA_ANALYSIS         } from '../subworkflows/local/freyja_analysis'
include { INFLUENZA_SEROTYPING    } from '../subworkflows/local/influenza_serotyping'
include { CONTROL_OLIGO_QC        } from '../subworkflows/local/control_oligo_qc'

// Import remaining modules
include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
include { BAM_TRIM_PRIMERS_IVAR   } from '../subworkflows/local/bam_trim_primers_ivar'

// Import local modules
include { SPLIT_BAM_BY_SEGMENT     } from '../modules/local/split_bam_by_segment'
include { REHEADER_SEGMENT_BAM } from '../modules/local/reheader_segment_bam'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WASTEFLOW {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    def pass_mapped_reads = [:]
    def fail_mapped_reads = [:]

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    multiqc_report = Channel.empty()

    //
    // SUBWORKFLOW: Prepare reference genomes
    //
    if (!params.skip_genome_preparation) {
        GENOME_PREPARATION()
        
        // SARS-CoV-2
        ch_sars_cov2_fasta = GENOME_PREPARATION.out.sars_cov2_fasta
        ch_sars_cov2_fai = GENOME_PREPARATION.out.sars_cov2_fai
        ch_sars_cov2_chrom_sizes = GENOME_PREPARATION.out.sars_cov2_chrom_sizes
        ch_sars_cov2_gff = GENOME_PREPARATION.out.sars_cov2_gff
        ch_sars_cov2_snpeff_db = GENOME_PREPARATION.out.sars_cov2_snpeff_db
        ch_sars_cov2_snpeff_config = GENOME_PREPARATION.out.sars_cov2_snpeff_config

        // RSV-A
        ch_rsv_a_fasta = GENOME_PREPARATION.out.rsv_a_fasta
        ch_rsv_a_fai = GENOME_PREPARATION.out.rsv_a_fai
        ch_rsv_a_chrom_sizes = GENOME_PREPARATION.out.rsv_a_chrom_sizes
        ch_rsv_a_gff = GENOME_PREPARATION.out.rsv_a_gff
        ch_rsv_a_snpeff_db = GENOME_PREPARATION.out.rsv_a_snpeff_db
        ch_rsv_a_snpeff_config = GENOME_PREPARATION.out.rsv_a_snpeff_config
        
        // RSV-B
        ch_rsv_b_fasta = GENOME_PREPARATION.out.rsv_b_fasta
        ch_rsv_b_fai = GENOME_PREPARATION.out.rsv_b_fai
        ch_rsv_b_chrom_sizes = GENOME_PREPARATION.out.rsv_b_chrom_sizes
        ch_rsv_b_gff = GENOME_PREPARATION.out.rsv_b_gff
        ch_rsv_b_snpeff_db = GENOME_PREPARATION.out.rsv_b_snpeff_db
        ch_rsv_b_snpeff_config = GENOME_PREPARATION.out.rsv_b_snpeff_config
        ch_rsv_combined_fasta = GENOME_PREPARATION.out.rsv_combined_fasta
        
        // H1N1 - Full genome
        ch_h1n1_fasta = GENOME_PREPARATION.out.h1n1_fasta
        ch_h1n1_fai = GENOME_PREPARATION.out.h1n1_fai
        ch_h1n1_chrom_sizes = GENOME_PREPARATION.out.h1n1_chrom_sizes
        
        // H1N1 - Segments
        ch_h1n1_segment_fasta = GENOME_PREPARATION.out.h1n1_segment_fasta
        ch_h1n1_segment_gff = GENOME_PREPARATION.out.h1n1_segment_gff
        ch_h1n1_segment_fai = GENOME_PREPARATION.out.h1n1_segment_fai
        ch_h1n1_segment_chrom_sizes = GENOME_PREPARATION.out.h1n1_segment_chrom_sizes
        ch_h1n1_segment_snpeff_db = GENOME_PREPARATION.out.h1n1_segment_snpeff_db
        ch_h1n1_segment_snpeff_config = GENOME_PREPARATION.out.h1n1_segment_snpeff_config
        
        // H3N2 - Full genome
        ch_h3n2_fasta = GENOME_PREPARATION.out.h3n2_fasta
        ch_h3n2_fai = GENOME_PREPARATION.out.h3n2_fai
        ch_h3n2_chrom_sizes = GENOME_PREPARATION.out.h3n2_chrom_sizes

        // H3N2 - Segments
        ch_h3n2_segment_fasta = GENOME_PREPARATION.out.h3n2_segment_fasta
        ch_h3n2_segment_gff = GENOME_PREPARATION.out.h3n2_segment_gff
        ch_h3n2_segment_fai = GENOME_PREPARATION.out.h3n2_segment_fai
        ch_h3n2_segment_chrom_sizes = GENOME_PREPARATION.out.h3n2_segment_chrom_sizes
        ch_h3n2_segment_snpeff_db = GENOME_PREPARATION.out.h3n2_segment_snpeff_db
        ch_h3n2_segment_snpeff_config = GENOME_PREPARATION.out.h3n2_segment_snpeff_config
        
        // H5N1 - Full genome
        ch_h5n1_fasta = GENOME_PREPARATION.out.h5n1_fasta
        ch_h5n1_fai = GENOME_PREPARATION.out.h5n1_fai
        ch_h5n1_chrom_sizes = GENOME_PREPARATION.out.h5n1_chrom_sizes

        // H5N1 - Segments
        ch_h5n1_segment_fasta = GENOME_PREPARATION.out.h5n1_segment_fasta
        ch_h5n1_segment_gff = GENOME_PREPARATION.out.h5n1_segment_gff
        ch_h5n1_segment_fai = GENOME_PREPARATION.out.h5n1_segment_fai
        ch_h5n1_segment_chrom_sizes = GENOME_PREPARATION.out.h5n1_segment_chrom_sizes
        ch_h5n1_segment_snpeff_db = GENOME_PREPARATION.out.h5n1_segment_snpeff_db
        ch_h5n1_segment_snpeff_config = GENOME_PREPARATION.out.h5n1_segment_snpeff_config

        ch_versions = ch_versions.mix(GENOME_PREPARATION.out.versions)
        }

    // Combine all segment references into single channels
    ch_all_segment_fasta = ch_h1n1_segment_fasta
        .mix(ch_h3n2_segment_fasta)
        .mix(ch_h5n1_segment_fasta)

    ch_all_segment_fai = ch_h1n1_segment_fai
        .mix(ch_h3n2_segment_fai)
        .mix(ch_h5n1_segment_fai)

    ch_all_segment_gff = ch_h1n1_segment_gff
        .mix(ch_h3n2_segment_gff)
        .mix(ch_h5n1_segment_gff)
    
    ch_all_segment_chrom_sizes = ch_h1n1_segment_chrom_sizes
        .mix(ch_h3n2_segment_chrom_sizes)
        .mix(ch_h5n1_segment_chrom_sizes)
    
    ch_all_segment_snpeff_db = ch_h1n1_segment_snpeff_db
        .mix(ch_h3n2_segment_snpeff_db)
        .mix(ch_h5n1_segment_snpeff_db)
    
    ch_all_segment_snpeff_config = ch_h1n1_segment_snpeff_config
        .mix(ch_h3n2_segment_snpeff_config)
        .mix(ch_h5n1_segment_snpeff_config)

    //
    // SUBWORKFLOW: Process reads (QC, trim)
    //
    if (!params.skip_read_processing) {
        READ_PROCESSING(ch_samplesheet)
        ch_versions = ch_versions.mix(READ_PROCESSING.out.versions)
    }
    ch_reads_qc = READ_PROCESSING.out.reads
    ch_reads_oligo = READ_PROCESSING.out.reads


    //
    // SUBWORKFLOW: Control oligo QC
    //
    if (!params.skip_genome_preparation) {
        panel_control_db = GENOME_PREPARATION.out.panel_control_db_fasta.map { fasta -> [[id: 'panel_control_db'], fasta] }
    }
    else {
        log.error("Genome preparation step was skipped. Panel control database FASTA is required for control oligo QC.")
        System.exit(1)
    }
    if (!params.skip_control_oligo_qc) {
        CONTROL_OLIGO_QC(ch_reads_oligo, panel_control_db)
        ch_versions = ch_versions.mix(CONTROL_OLIGO_QC.out.versions)
    }

    //
    // SUBWORKFLOW: Taxonomy classification and target extraction
    //
    if (!params.skip_taxonomy_classification) {
        TAXONOMY_CLASSIFICATION(ch_reads_qc)
        ch_versions = ch_versions.mix(TAXONOMY_CLASSIFICATION.out.versions)
        ch_kraken2_multiqc = TAXONOMY_CLASSIFICATION.out.kraken2_multiqc
    }

    ch_all_extracted_reads = TAXONOMY_CLASSIFICATION.out.sars_cov2_reads
        .mix(TAXONOMY_CLASSIFICATION.out.rsv_reads)
        .mix(TAXONOMY_CLASSIFICATION.out.flu_reads)

    //
    // SUBWORKFLOW: Influenza serotyping
    //
    if (!params.skip_influenza_serotyping) {
        INFLUENZA_SEROTYPING(
            ch_all_extracted_reads
        )
        ch_versions = ch_versions.mix(INFLUENZA_SEROTYPING.out.versions)
    }
    influenza_serotype_reads = INFLUENZA_SEROTYPING.out.serotype_reads
    ch_all_extracted_reads = ch_all_extracted_reads.mix(influenza_serotype_reads)
    
    //
    // SUBWORKFLOW: Alignment processing
    //
    if (!params.skip_alignment) {
        ALIGNMENT(
            ch_all_extracted_reads,
            ch_sars_cov2_fasta,
            ch_rsv_combined_fasta,
            ch_rsv_a_fasta,
            ch_rsv_b_fasta,
            ch_h1n1_fasta,
            ch_h3n2_fasta,
            ch_h5n1_fasta
        )
        ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

        //
        // Filter channels to get samples that passed minimum mapped reads threshold
        //
        ch_fail_mapping_multiqc = Channel.empty()
        ch_bam = Channel.empty()
        ch_bai = Channel.empty()


        ALIGNMENT.out.minimap2_flagstat_multiqc
            .map { meta, flagstat -> [meta] + WorkflowWasteflow.getFlagstatMappedReads(flagstat, params) }
            .set { ch_mapped_reads }

        ALIGNMENT.out.bam
            .join(ch_mapped_reads, by: [0])
            .map { meta, ofile, mapped, pass ->
                if (pass) {
                    [meta, ofile]
                }
            }
            .set { ch_bam }

        ALIGNMENT.out.bai
            .join(ch_mapped_reads, by: [0])
            .map { meta, ofile, mapped, pass ->
                if (pass) {
                    [meta, ofile]
                }
            }
            .set { ch_bai }

        ch_mapped_reads
            .branch { meta, mapped, pass ->
                pass: pass
                pass_mapped_reads[meta.id] = mapped
                return ["${meta.id}\t${mapped}"]
                fail: !pass
                fail_mapped_reads[meta.id] = mapped
                return ["${meta.id}\t${mapped}"]
            }
            .set { ch_pass_fail_mapped }

        ch_pass_fail_mapped.fail
            .collect()
            .map { tsv_data ->
                def header = ['Sample', 'Mapped reads']
                WorkflowCommons.multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_mapping_multiqc }
    }
    
    
    // Define influenza types and their BED files
    def influenza_bed_files = [
        'H1N1': params.genomes['H1N1'].bed,
        'H3N2': params.genomes['H3N2'].bed,
        'H5N1': params.genomes['H5N1'].bed
    ]

    // Parse all segment BED files
    ch_all_segments = Channel.empty()

    influenza_bed_files.each { virus, bed_path ->
        def ch_temp = Channel
            .fromPath(bed_path)
            .splitCsv(sep: '\t', header: false)
            .map { row -> 
                [
                    virus: virus,
                    accession: row[0],
                    segment: row[3]
                ]
            }
        
        ch_all_segments = ch_all_segments.mix(ch_temp)
    }

    
    // Filter influenza samples and split by segments
    ch_bam
        .filter { meta, bam -> meta.genome in ['H1N1', 'H3N2', 'H5N1'] }
        .join(ch_bai, by: [0])
        .combine(ch_all_segments)
        .filter { meta, bam, bai, segment_info ->
            meta.genome == segment_info.virus
        }
        .map { meta, bam, bai, segment_info ->
            def new_meta = meta.clone() + [
                segment: segment_info.segment,
                segment_accession: segment_info.accession
            ]
            [new_meta, bam, bai, segment_info.accession, segment_info.segment]
        }
        .set { ch_influenza_bam_for_split }
    

    SPLIT_BAM_BY_SEGMENT(ch_influenza_bam_for_split)
    ch_influenza_split_bam = SPLIT_BAM_BY_SEGMENT.out.bam
    
    
    //
    // MODULE: Reheader segment BAMs to match segment-specific references
    //
    
    
    // Prepare channel with BAM, BAI, and matching segment FASTA/FAI
    ch_influenza_split_bam
        .map { meta, bam -> 
            def key = "${meta.genome}_${meta.segment}"
            [key, meta, bam]
        }
        .join(
            SPLIT_BAM_BY_SEGMENT.out.bai
                .map { meta, bai -> 
                    def key = "${meta.genome}_${meta.segment}"
                    [key, meta, bai]
                },
            by: [0, 1]
        )
        .combine(
            ch_all_segment_fasta
                .map { ref_meta, fasta ->
                    def key = "${ref_meta.virus}_${ref_meta.segment}"
                    [key, fasta]
                },
            by: 0
        )
        .combine(
            ch_all_segment_fai
                .map { ref_meta, fai ->
                    def key = "${ref_meta.virus}_${ref_meta.segment}"
                    [key, fai]
                },
            by: 0
        )
        .map { key, meta, bam, bai, fasta, fai ->
            [meta, bam, bai, fasta, fai]
        }
        .set { ch_segments_for_reheader }

    REHEADER_SEGMENT_BAM(ch_segments_for_reheader)    
    ch_versions = ch_versions.mix(REHEADER_SEGMENT_BAM.out.versions)

    // Combine non-influenza BAMs with reheadered influenza segment BAMs
    ch_bam_reheader = ch_bam
        .filter { meta, bam -> 
            // Keep only non-influenza samples (SARS-CoV-2, RSV-A, RSV-B)
            !(meta.genome in ['H1N1', 'H3N2', 'H5N1'])
        }
        .mix(REHEADER_SEGMENT_BAM.out.bam)

    // Similarly for BAI
    ch_bai_reheader = ch_bai
        .filter { meta, bai -> 
            !(meta.genome in ['H1N1', 'H3N2', 'H5N1'])
        }
        .mix(REHEADER_SEGMENT_BAM.out.bai)

    // For QC, use the same combined channels
    ch_bam_for_qc = ch_bam_reheader
    ch_bai_for_qc = ch_bai_reheader


    // Run BAM QC with segment-specific references
    if (!params.skip_alignment && !params.skip_alignment_stats) {
        BAM_QC_METRICS(
            ch_bam_for_qc,
            ch_bai_for_qc,
            ch_sars_cov2_fasta,
            ch_sars_cov2_fai,
            ch_rsv_a_fasta,
            ch_rsv_a_fai,
            ch_rsv_b_fasta,
            ch_rsv_b_fai,
            ch_all_segment_fasta,
            ch_all_segment_fai
        )
        ch_versions = ch_versions.mix(BAM_QC_METRICS.out.versions)

        if (!params.skip_markduplicates) {
            // Update the reheadered channels with mark duplicates output
            ch_bam_reheader = BAM_QC_METRICS.out.bam
            ch_bai_reheader = BAM_QC_METRICS.out.bai
        }
    }
    
    //
    // SUBWORKFLOW: Call variants with FreeBayes + iVar + LoFreq
    //
    if (!params.skip_variants) {
        ch_variant_bam = ch_bam_reheader
        ch_vcf = Channel.empty()
        ch_tbi = Channel.empty()
        ch_ivar_counts_multiqc = Channel.empty()
        ch_bcftools_stats_multiqc = Channel.empty()
        ch_snpsift_txt = Channel.empty()
        ch_snpeff_multiqc = Channel.empty()

        VARIANT_CALLING(
            ch_variant_bam,
            ch_sars_cov2_fasta,
            ch_sars_cov2_fai,
            ch_sars_cov2_chrom_sizes,
            ch_sars_cov2_gff,
            ch_sars_cov2_snpeff_db,
            ch_sars_cov2_snpeff_config,
            ch_rsv_a_fasta,
            ch_rsv_a_fai,
            ch_rsv_a_chrom_sizes,
            ch_rsv_a_gff,
            ch_rsv_a_snpeff_db,
            ch_rsv_a_snpeff_config,
            ch_rsv_b_fasta,
            ch_rsv_b_fai,
            ch_rsv_b_chrom_sizes,
            ch_rsv_b_gff,
            ch_rsv_b_snpeff_db,
            ch_rsv_b_snpeff_config,
        )

        // Combine all variant calling outputs for SARS-CoV-2, RSV-A, RSV-B, and influenza types
        ch_vcf = VARIANT_CALLING.out.sars_cov_2_vcf
            .mix(VARIANT_CALLING.out.rsv_a_vcf)
            .mix(VARIANT_CALLING.out.rsv_b_vcf)
            

        ch_tbi = VARIANT_CALLING.out.sars_cov_2_tbi
            .mix(VARIANT_CALLING.out.rsv_a_tbi)
            .mix(VARIANT_CALLING.out.rsv_b_tbi)
            

        ch_ivar_counts_multiqc = VARIANT_CALLING.out.sars_cov_2_ivar_counts_multiqc
            .mix(VARIANT_CALLING.out.rsv_a_ivar_counts_multiqc)
            .mix(VARIANT_CALLING.out.rsv_b_ivar_counts_multiqc)
           

        ch_bcftools_stats_multiqc = VARIANT_CALLING.out.sars_cov_2_bcftools_stats_multiqc
            .mix(VARIANT_CALLING.out.rsv_a_bcftools_stats_multiqc)
            .mix(VARIANT_CALLING.out.rsv_b_bcftools_stats_multiqc)
            

        ch_snpeff_multiqc = VARIANT_CALLING.out.sars_cov_2_snpeff_multiqc
            .mix(VARIANT_CALLING.out.rsv_a_snpeff_multiqc)
            .mix(VARIANT_CALLING.out.rsv_b_snpeff_multiqc)
            

        ch_snpsift_txt = VARIANT_CALLING.out.sars_cov_2_snpsift_txt
            .mix(VARIANT_CALLING.out.rsv_a_snpsift_txt)
            .mix(VARIANT_CALLING.out.rsv_b_snpsift_txt)
            

        ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)
    }
    //
    // SUBWORKFLOW: Freyja variant analysis
    //
    if (!params.skip_variants && !params.skip_freyja) {
        FREYJA_ANALYSIS(
            ch_bam,
            ch_sars_cov2_fasta,
            ch_rsv_a_fasta,
            ch_rsv_b_fasta,
        )
        ch_freyja_organized = FREYJA_ANALYSIS.out.freyja_organized
        ch_versions = ch_versions.mix(FREYJA_ANALYSIS.out.versions)
    }


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'wasteflow2_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        summary_params = paramsSummaryMap(
            workflow,
            parameters_schema: "nextflow_schema.json"
        )
        ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description
            ? file(params.multiqc_methods_description, checkIfExists: true)
            : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description = Channel.value(
            methodsDescriptionText(ch_multiqc_custom_methods_description)
        )

        ch_multiqc_logo = params.multiqc_logo
            ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
            : Channel.empty()

        ch_multiqc_files = ch_multiqc_files.mix(
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
        )
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: false,
            )
        )

        // FastQC files,
        ch_multiqc_files = ch_multiqc_files.mix(READ_PROCESSING.out.fastqc_raw_zip.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(READ_PROCESSING.out.trim_json.collect { it[1] }.ifEmpty([]))

        // Kraken files
        ch_multiqc_files = ch_multiqc_files.mix(ch_kraken2_multiqc.collect { it[1] }.ifEmpty([]))

        // Add to MultiQC files
        ch_multiqc_files = ch_multiqc_files.mix(ch_freyja_organized.ifEmpty([]))

        ch_multiqc_config = file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? file(params.multiqc_config) : []
        MULTIQC(
            ch_multiqc_files.collect(),
            ch_multiqc_config,
            ch_multiqc_custom_config,
            ch_multiqc_logo.toList(),
            [],
            [],
        )
        multiqc_report = MULTIQC.out.report.toList()
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    emit:
    multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
