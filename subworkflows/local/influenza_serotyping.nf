#!/usr/bin/env nextflow

include { MINIMAP2_ALIGN                         } from '../../modules/nf-core/minimap2/align/main'
include { SEQKIT_GREP                            } from '../../modules/nf-core/seqkit/grep/main'
include { SEQKIT_STATS as SEQKIT_STATS_SEROTYPES } from '../../modules/nf-core/seqkit/stats/main'
include { ASSIGN_SEROTYPES                       } from '../../modules/local/assign_serotypes.nf'


workflow INFLUENZA_SEROTYPING {
    take:
    ch_all_extracted_reads // channel: [meta, reads] - flu samples from branched reads

    main:
    ch_versions = Channel.empty()
    ch_serotype_assignments = Channel.empty()
    ch_serotype_reads = Channel.empty()
    ch_assignment_plot = Channel.empty()
    ch_read_summary = Channel.empty()

    // Create channels for database files with validation
    ch_flu_db_fasta = Channel.fromPath(params.flu_db_fasta, checkIfExists: true)
        .ifEmpty { error("Flu database FASTA file not found: ${params.flu_db_fasta}") }

    ch_flu_db_info = Channel.fromPath(params.flu_db_info, checkIfExists: true)
        .ifEmpty { error("Flu database info file not found: ${params.flu_db_info}") }

    // Prepare flu samples from input channel
    ch_flu_reads = ch_all_extracted_reads
        .filter { meta, _reads ->
            meta.taxid == '11308'
        }
        .map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.genome = 'flu_db'
            [new_meta, reads]
        }

    // Prepare database for MINIMAP2
    ch_flu_db_fasta_meta = ch_flu_db_fasta.map { fasta -> [[id: 'flu_db'], fasta] }

    // Run MINIMAP2 alignment
    bam_format = 'true'
    bam_format_ext = 'bai'
    cigar_paf_format = 'false'
    cigar_bam = 'false'
    MINIMAP2_ALIGN(
        ch_flu_reads,
        ch_flu_db_fasta_meta.first(),
        bam_format,
        bam_format_ext,
        cigar_paf_format,
        cigar_bam,
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    ASSIGN_SEROTYPES(
        MINIMAP2_ALIGN.out.bam,
        ch_flu_db_info.first(),
    )
    
    ch_read_summary = ASSIGN_SEROTYPES.out.read_summary
    ch_serotype_assignments = ASSIGN_SEROTYPES.out.serotype_summary
    ch_assignment_plot = ASSIGN_SEROTYPES.out.assignment_plot
    ch_versions = ch_versions.mix(ASSIGN_SEROTYPES.out.versions)

    // ========================================================================
    // SEQKIT_GREP INTEGRATION - Extract reads matching assigned serotypes
    // ========================================================================
    
    // Step 1: Flatten the read lists (one tuple per serotype text file)
    read_lists_flattened = ASSIGN_SEROTYPES.out.read_lists
        .transpose()
        .filter { _meta, txt -> txt.size() > 0 }
    
    
    // Step 2: Add serotype information to metadata (keep original ID)
    read_lists_with_serotype = read_lists_flattened
        .map { meta, txt ->
            // Extract serotype from filename: sampleID_SEROTYPE.txt
            def prefix = meta.id + '_'
            def serotype = txt.baseName.replaceFirst(prefix, '')
            
            // Create new meta with serotype info (keep original ID)
            def new_meta = meta.clone()
            new_meta.serotype = serotype
            
            tuple(new_meta, txt)
        }
    
    // Step 3: Join with original reads using BOTH sample ID and ensure we preserve pairing
    // Key by sample ID only for the join
    read_lists_keyed = read_lists_with_serotype
        .map { meta, txt ->
            tuple(meta.id, meta, txt)
        }
    
    reads_keyed = ch_flu_reads
        .map { meta, reads ->
            tuple(meta.id, meta, reads)
        }
    
    reads_with_serotype = read_lists_keyed
        .combine(reads_keyed, by: 0)
        .map { _orig_id, meta_sero, txt, meta_reads, reads ->
            // Merge metadata, preserving single_end info from reads
            def final_meta = meta_sero + [
                single_end: meta_reads.single_end ?: false
            ]
            tuple(final_meta, txt, reads)
        }
    
    // Step 4: Split paired-end reads into separate R1/R2 invocations
    // Keep meta.id unchanged, add read_pair info for suffix only
    seqkit_inputs = reads_with_serotype
        .flatMap { meta, txt, reads ->
            if (meta.single_end) {
                // Single-end: one invocation
                [[meta, reads[0], txt]]
            } else {
                // Paired-end: add read_pair field but DON'T change meta.id
                def meta_r1 = meta.clone()
                meta_r1.read_pair = 'R1'
                
                def meta_r2 = meta.clone()
                meta_r2.read_pair = 'R2'
                
                [
                    [meta_r1, reads[0], txt],
                    [meta_r2, reads[1], txt]
                ]
            }
        }
    
    // Step 5: Call SEQKIT_GREP
    SEQKIT_GREP(
        seqkit_inputs.map { meta, read, _txt -> tuple(meta, read) },
        seqkit_inputs.map { _meta, _read, txt -> txt }
    )
    
    // Step 6: Group R1/R2 pairs back together by sample ID + serotype
    ch_serotype_reads = SEQKIT_GREP.out.filter
        .map { meta, fastq ->
            // Create grouping key from id + serotype
            def group_key = "${meta.id}_${meta.serotype}"
            tuple(group_key, meta, fastq)
        }
        .groupTuple(by: 0)
        .map { group_key, metas, fastqs ->
            // Take first meta (they should all have same id and serotype)
            def meta = metas[0].clone()
            // Remove read_pair from final meta if it exists
            meta.remove('read_pair')
            
            // Sort fastqs by R1/R2 if paired
            def sorted_fastqs = fastqs
            if (fastqs.size() == 2) {
                // Sort by the read_pair field from the metas
                def indexed = [metas, fastqs].transpose().sort { it[0].read_pair }
                sorted_fastqs = indexed.collect { it[1] }
                meta.single_end = false
            } else {
                meta.single_end = true
            }
            
            tuple(meta, sorted_fastqs)
        }
    
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)
    
    // Get stats of all extracted reads before branching
    ch_seqkit_sequences = ch_serotype_reads
        .map { meta, reads -> reads }
        .collect()
        .map { seqs -> [[id: "seqkit_stat_all_extracted"], seqs] }

    SEQKIT_STATS_SEROTYPES(ch_seqkit_sequences)
    ch_versions = ch_versions.mix(SEQKIT_STATS_SEROTYPES.out.versions)

    emit:
    read_summary         = ch_read_summary          // channel: [ val(meta), path(tsv) ]
    serotype_assignments = ch_serotype_assignments  // channel: [ val(meta), path(tsv) ]
    assignment_plot      = ch_assignment_plot       // channel: [ val(meta), path(pdf) ]
    serotype_reads       = ch_serotype_reads        // channel: [ val(meta), path(fastq) ]
    versions             = ch_versions              // channel: [ path(yml) ]
}