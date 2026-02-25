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
include { VARIANTS_QC             } from '../subworkflows/local/variants_qc'

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


    if (params.metadata) {
        // Parse metadata filter string into a map
        def filter_map = [:]
        if (params.metadata_filter) {
            params.metadata_filter.split(/\s+/).each { pair ->
                def parts = pair.split('=', 2)
                if (parts.size() == 2) {
                    filter_map[parts[0].trim()] = parts[1].trim()
                }
            }
        }
        
        // Load and parse metadata with filtering
        def metadata_map = Utils.loadMetadataWithFilters(
            file(params.metadata),
            filter_map,
            log
        )
        
        log.info "Filtering samples based on metadata file: ${params.metadata}"
        
        if (filter_map.size() > 0) {
            def filter_list = filter_map.collect { k, v -> "${k}=${v}" }
            log.info "Active filters: ${filter_list.join(', ')}"
        } else {
            log.info "No filters applied - including all samples from metadata that match sample IDs"
        }
        log.info "Found ${metadata_map.size()} samples matching filter criteria"
        log.info "All samples will be processed through READ_PROCESSING and CONTROL_OLIGO_QC"
        log.info "Filtering will be applied starting from TAXONOMY_CLASSIFICATION"
        
        // Create a closure that captures metadata_map
        ch_filter_closure = { meta, fastqs ->
            def normalized_id = meta.id
                .replaceAll(/_S\d+$/, '')
                .replaceAll('-', ' ')
            
            def should_process = metadata_map.containsKey(normalized_id)
            
            if (!should_process) {
                log.debug "Filtering out sample: ${meta.id} (normalized: ${normalized_id})"
            }
            
            return should_process
        }
        
    } else {
        log.info "No metadata file provided - processing all samples through entire pipeline"
        ch_filter_closure = null
    }

      
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
    
    // SARS-CoV-2
    // Remove the .map() - the channel is already [meta, db]
    ch_sars_cov2_snpeff_db_meta = GENOME_PREPARATION.out.sars_cov2_snpeff_db
    ch_sars_cov2_snpeff_config_meta = GENOME_PREPARATION.out.sars_cov2_snpeff_config    
    ch_sars_cov2_fasta_meta = GENOME_PREPARATION.out.sars_cov2_fasta
        .map { fasta -> [[genome: 'MN908947.3'], fasta] }  // This one still needs wrapping

    // Same for RSV-A
    ch_rsv_a_snpeff_db_meta = GENOME_PREPARATION.out.rsv_a_snpeff_db
    ch_rsv_a_snpeff_config_meta = GENOME_PREPARATION.out.rsv_a_snpeff_config
    ch_rsv_a_fasta_meta = GENOME_PREPARATION.out.rsv_a_fasta
        .map { fasta -> [[genome: 'PP109421.1'], fasta] }

    // Same for RSV-B
    ch_rsv_b_snpeff_db_meta = GENOME_PREPARATION.out.rsv_b_snpeff_db
    ch_rsv_b_snpeff_config_meta = GENOME_PREPARATION.out.rsv_b_snpeff_config
    ch_rsv_b_fasta_meta = GENOME_PREPARATION.out.rsv_b_fasta
        .map { fasta -> [[genome: 'OP975389.1'], fasta] }
    
    //
    // Step 2: Normalize segment SnpEff DBs (change 'virus' to 'genome' in meta)
    //
    
    // The segment channels have meta: [virus: 'H1N1', segment: 'PB2', id: 'H1N1_PB2']
    // We need to normalize to: [genome: 'H1N1', segment: 'PB2']
    
    ch_all_segment_snpeff_db_normalized = ch_all_segment_snpeff_db
        .map { meta, db -> 
            [[genome: meta.virus, segment: meta.segment], db]
        }
    
    ch_all_segment_snpeff_config_normalized = ch_all_segment_snpeff_config
        .map { meta, config -> 
            [[genome: meta.virus, segment: meta.segment], config]
        }
    
    ch_all_segment_fasta_normalized = ch_all_segment_fasta
        .map { meta, fasta ->
            [[genome: meta.virus, segment: meta.segment], fasta]
        }
    
    //
    // Step 3: Mix all SnpEff DBs into unified channels
    //
    
    ch_all_snpeff_db_unified = ch_all_segment_snpeff_db_normalized
        .mix(ch_sars_cov2_snpeff_db_meta)
        .mix(ch_rsv_a_snpeff_db_meta)
        .mix(ch_rsv_b_snpeff_db_meta)
    
    ch_all_snpeff_config_unified = ch_all_segment_snpeff_config_normalized
        .mix(ch_sars_cov2_snpeff_config_meta)
        .mix(ch_rsv_a_snpeff_config_meta)
        .mix(ch_rsv_b_snpeff_config_meta)
    
    ch_all_fasta_unified = ch_all_segment_fasta_normalized
        .mix(ch_sars_cov2_fasta_meta)
        .mix(ch_rsv_a_fasta_meta)
        .mix(ch_rsv_b_fasta_meta)

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
    if (params.metadata && ch_filter_closure != null) {
        ch_reads_for_taxonomy = ch_reads_qc.filter(ch_filter_closure)
    } else {
        ch_reads_for_taxonomy = ch_reads_qc
    }
    if (!params.skip_taxonomy_classification) {
        TAXONOMY_CLASSIFICATION(ch_reads_for_taxonomy)
        ch_versions = ch_versions.mix(TAXONOMY_CLASSIFICATION.out.versions)
        ch_kraken2_multiqc = TAXONOMY_CLASSIFICATION.out.kraken2_multiqc
    }

    ch_all_extracted_reads = TAXONOMY_CLASSIFICATION.out.sars_cov2_reads
        .mix(TAXONOMY_CLASSIFICATION.out.rsv_reads)
        .mix(TAXONOMY_CLASSIFICATION.out.flu_reads)

    //
    // SUBWORKFLOW: Influenza serotyping
    //
    
    influenza_serotype_reads = Channel.empty()
    if (!params.skip_influenza_serotyping) {
        INFLUENZA_SEROTYPING(
            ch_all_extracted_reads
        )
        influenza_serotype_reads = INFLUENZA_SEROTYPING.out.serotype_reads
        ch_versions = ch_versions.mix(INFLUENZA_SEROTYPING.out.versions)
    }
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
    
    // Sort all channels deterministically for reproducible matching
    ch_all_segment_fasta
        .toSortedList { a, b -> 
            def keyA = "${a[0].virus}_${a[0].segment}"
            def keyB = "${b[0].virus}_${b[0].segment}"
            keyA <=> keyB
        }
        .flatMap { it }
        .set { ch_sorted_segment_fasta }
    
    ch_all_segment_fai
        .toSortedList { a, b -> 
            def keyA = "${a[0].virus}_${a[0].segment}"
            def keyB = "${b[0].virus}_${b[0].segment}"
            keyA <=> keyB
        }
        .flatMap { it }
        .set { ch_sorted_segment_fai }
    
    // Sort split BAM outputs deterministically
    ch_influenza_split_bam
        .toSortedList { a, b ->
            def keyA = "${a[0].id}_${a[0].genome}_${a[0].segment}"
            def keyB = "${b[0].id}_${b[0].genome}_${b[0].segment}"
            keyA <=> keyB
        }
        .flatMap { it }
        .set { ch_sorted_split_bam }
    
    SPLIT_BAM_BY_SEGMENT.out.bai
        .toSortedList { a, b ->
            def keyA = "${a[0].id}_${a[0].genome}_${a[0].segment}"
            def keyB = "${b[0].id}_${b[0].genome}_${b[0].segment}"
            keyA <=> keyB
        }
        .flatMap { it }
        .set { ch_sorted_split_bai }
    
    // Create deterministic reference pairs [ref_meta, fasta, fai]
    ch_sorted_segment_fasta
        .combine(ch_sorted_segment_fai)
        .filter { fasta_meta, fasta, fai_meta, fai ->
            fasta_meta.virus == fai_meta.virus && fasta_meta.segment == fai_meta.segment
        }
        .map { fasta_meta, fasta, fai_meta, fai ->
            [fasta_meta, fasta, fai]
        }
        .set { ch_segment_refs }
    
    // Match BAM with BAI
    ch_sorted_split_bam
        .combine(ch_sorted_split_bai)
        .filter { bam_meta, bam, bai_meta, bai ->
            bam_meta.id == bai_meta.id && 
            bam_meta.genome == bai_meta.genome && 
            bam_meta.segment == bai_meta.segment
        }
        .map { bam_meta, bam, bai_meta, bai ->
            [bam_meta, bam, bai]
        }
        .set { ch_bam_bai_pairs }
    
    // Match with segment references
    ch_bam_bai_pairs
        .combine(ch_segment_refs)
        .filter { bam_meta, bam, bai, ref_meta, fasta, fai ->
            bam_meta.genome == ref_meta.virus && bam_meta.segment == ref_meta.segment
        }
        .map { bam_meta, bam, bai, ref_meta, fasta, fai ->
            // Extract segment accession as a simple value
            def seg_acc = bam_meta.segment_accession ?: bam_meta.segment
            
            // Return tuple with all inputs
            tuple(bam_meta, bam, bai, fasta, fai, seg_acc)
        }
        .multiMap { meta, bam, bai, fasta, fai, seg_acc ->
            tuple_input: tuple(meta, bam, bai, fasta, fai)
            seg_acc: seg_acc  // This will now be a clean scalar value
        }
        .set { ch_reheader_inputs }

        
    REHEADER_SEGMENT_BAM(
        ch_reheader_inputs.tuple_input,
        ch_reheader_inputs.seg_acc
    )

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
    ch_variant_bam = ch_bam_reheader
    ch_vcf = Channel.empty()
    ch_tbi = Channel.empty()
    ch_ivar_counts_multiqc = Channel.empty()
    ch_bcftools_stats_multiqc = Channel.empty()
    if (!params.skip_variants) {
        VARIANT_CALLING(
            ch_variant_bam,
            ch_sars_cov2_fasta,
            ch_sars_cov2_fai,
            ch_sars_cov2_chrom_sizes,
            ch_sars_cov2_gff,
            ch_rsv_a_fasta,
            ch_rsv_a_fai,
            ch_rsv_a_chrom_sizes,
            ch_rsv_a_gff,
            ch_rsv_b_fasta,
            ch_rsv_b_fai,
            ch_rsv_b_chrom_sizes,
            ch_rsv_b_gff,
            ch_all_segment_fasta,
            ch_all_segment_fai,
            ch_all_segment_chrom_sizes,
            ch_all_segment_gff,
        )            

        ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)
    }

    ch_vcf = VARIANT_CALLING.out.vcf
    ch_tbi = VARIANT_CALLING.out.tbi
    ch_ivar_counts_multiqc = VARIANT_CALLING.out.ivar_counts_multiqc
    ch_bcftools_stats_multiqc = VARIANT_CALLING.out.bcftools_stats

    VARIANTS_QC(
        ch_vcf,
        ch_all_snpeff_db_unified,         
        ch_all_snpeff_config_unified,     
        ch_all_fasta_unified              
    )
    snpeff_annotated_vcf = VARIANTS_QC.out.snpeff_vcf
    snpeff_annotated_tbi = VARIANTS_QC.out.snpeff_tbi
    snpeff_annotated_stats = VARIANTS_QC.out.snpeff_stats
    snpeff_annotated_csv = VARIANTS_QC.out.snpeff_csv
    snpeff_annotated_txt = VARIANTS_QC.out.snpeff_txt
    snpeff_annotated_html = VARIANTS_QC.out.snpeff_html
    snpsift_annotated_txt = VARIANTS_QC.out.snpsift_txt
    ch_versions = ch_versions.mix(VARIANTS_QC.out.versions)

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
