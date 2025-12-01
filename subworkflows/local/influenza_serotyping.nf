#!/usr/bin/env nextflow

include { MINIMAP2_ALIGN                         } from '../../modules/nf-core/minimap2/align/main'
include { PARSE_PAF_SEROTYPE                     } from '../../modules/local/parse_paf_serotype'
include { SEQKIT_GREP                            } from '../../modules/nf-core/seqkit/grep/main'

include { CLEAN_SEROTYPE_FILE                    } from '../../modules/local/clean_serotype_file'
include { SEQKIT_STATS as SEQKIT_STATS_SEROTYPES } from '../../modules/nf-core/seqkit/stats/main'


workflow INFLUENZA_SEROTYPING {
    take:
    ch_all_extracted_reads // channel: [meta, reads] - flu samples from branched reads

    main:
    ch_versions = Channel.empty()
    ch_serotype_assignments = Channel.empty()
    ch_serotype_reads = Channel.empty()

    if (!params.skip_influenza_serotyping) {
        // Create channels for database files with validation
        ch_flu_db_fasta = Channel.fromPath(params.flu_db_fasta, checkIfExists: true)
            .ifEmpty { error("Flu database FASTA file not found: ${params.flu_db_fasta}") }

        ch_flu_db_info = Channel.fromPath(params.flu_db_info, checkIfExists: true)
            .ifEmpty { error("Flu database info file not found: ${params.flu_db_info}") }

        // Prepare flu samples from input channel
        ch_flu_reads = ch_all_extracted_reads
            .filter { meta, reads ->
                def is_flu = meta.taxid == '11308'
                return is_flu
            }
            .map { meta, reads ->
                def new_meta = meta.clone()
                new_meta.genome = 'flu_db'
                [new_meta, reads]
            }
            .ifEmpty {
                log.warn("No flu samples found in input")
                Channel.empty()
            }

        // Prepare database for MINIMAP2
        ch_flu_db_fasta_meta = ch_flu_db_fasta.map { fasta -> [[id: 'flu_db'], fasta] }

        // Only proceed if we have flu reads
        if (ch_flu_reads) {
            // Run MINIMAP2 alignment
            MINIMAP2_ALIGN(
                ch_flu_reads,
                ch_flu_db_fasta_meta.first(),
                false,
                '',
                true,
                false,
            )
            ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

            // Parse PAF files and assign serotypes
            PARSE_PAF_SEROTYPE(
                MINIMAP2_ALIGN.out.paf,
                ch_flu_db_info.first(),
                params.score_threshold ?: 0.8,
            )
            ch_versions = ch_versions.mix(PARSE_PAF_SEROTYPE.out.versions)

            // Create input for serotype cleaning
            ch_serotype_input = PARSE_PAF_SEROTYPE.out.serotype_reads.flatMap { meta, files ->
                def file_list = files instanceof List ? files : [files]
                def entries = []

                file_list.each { file ->
                    def file_obj = file instanceof String ? new File(file) : file
                    def file_name = file_obj.name
                    def serotype = file_name.replaceAll("${meta.id}_", "").replaceAll(".txt", "")

                    if (file_obj.exists() && file_obj.size() > 0 && !serotype.trim().isEmpty()) {
                        def new_meta = meta.clone()
                        new_meta.serotype = serotype
                        entries.add([new_meta, file])
                    }
                }
                return entries
            }

            // Clean serotype files
            CLEAN_SEROTYPE_FILE(ch_serotype_input)

            // Prepare separate R1 and R2 channels for SEQKIT_GREP
            ch_seqkit_r1 = CLEAN_SEROTYPE_FILE.out.cleaned
                .map { meta, clean_file -> [meta.id, meta, clean_file] }
                .combine(
                    ch_flu_reads.map { meta, reads -> [meta.id, reads[0]] },
                    by: 0
                )
                .map { sample_id, meta, clean_file, r1_read ->
                    def r1_meta = meta.clone()
                    r1_meta.read_type = 'R1'
                    [r1_meta, [r1_read], clean_file]
                }

            ch_seqkit_r2 = CLEAN_SEROTYPE_FILE.out.cleaned
                .map { meta, clean_file -> [meta.id, meta, clean_file] }
                .combine(
                    ch_flu_reads.map { meta, reads -> [meta.id, reads[1]] },
                    by: 0
                )
                .map { sample_id, meta, clean_file, r2_read ->
                    def r2_meta = meta.clone()
                    r2_meta.read_type = 'R2'
                    [r2_meta, [r2_read], clean_file]
                }

            // Combine R1 and R2 channels
            ch_seqkit_input = ch_seqkit_r1.mix(ch_seqkit_r2)

            // Run SEQKIT_GREP (will run twice - once for R1, once for R2)
            SEQKIT_GREP(
                ch_seqkit_input.map { meta, reads, clean_file -> [meta, reads] },
                ch_seqkit_input.map { meta, reads, clean_file -> clean_file },
            )
            ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)

            // Combine R1 and R2 results back into pairs
            ch_serotype_reads = SEQKIT_GREP.out.filter
                .map { meta, filtered_reads ->
                    def key = "${meta.id}_${meta.serotype}"
                    [key, meta.read_type, filtered_reads, meta]
                }
                .groupTuple(by: 0)
                .map { key, read_types, files, metas ->
                    def meta = metas[0].clone()
                    meta.remove('read_type')
                    
                    def sorted_pairs = [read_types, files].transpose().sort { it[0] }
                    def sorted_files = sorted_pairs.collect { it[1] }
                    
                    [meta, sorted_files]
                }

            ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)
            ch_serotype_assignments = PARSE_PAF_SEROTYPE.out.summary
        }

        // Raw read statistics
        ch_seqkit_sequences = ch_serotype_reads
            .map { it[1] }
            .collect()
            .map { seqs -> [[id: "seqkit_stat_raw"], seqs] }

        SEQKIT_STATS_SEROTYPES(ch_seqkit_sequences)
        ch_versions = ch_versions.mix(SEQKIT_STATS_SEROTYPES.out.versions)
    }

    emit:
    serotype_assignments = ch_serotype_assignments // channel: [meta, serotype_summary] - serotype assignment summary
    serotype_reads       = ch_serotype_reads // channel: [meta, paired_files] - serotype read pairs
    paf_alignments       = MINIMAP2_ALIGN.out.paf // channel: [meta, paf] - PAF alignments
    bam_alignments       = MINIMAP2_ALIGN.out.bam // channel: [meta, bam] - BAM alignments
    versions             = ch_versions
}
