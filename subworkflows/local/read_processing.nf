#!/usr/bin/env nextflow

include { CAT_FASTQ                         } from '../../modules/nf-core/cat/fastq/main'
include { FASTQ_TRIM_FASTP_FASTQC           } from '../nf-core/fastq_trim_fastp_fastqc'
include { SEQKIT_STATS as SEQKIT_STATS_RAW  } from '../../modules/nf-core/seqkit/stats/main'
include { SEQKIT_STATS as SEQKIT_STATS_TRIM } from '../../modules/nf-core/seqkit/stats/main'

workflow READ_PROCESSING {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()
    ch_fail_reads_multiqc = Channel.empty()

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ(ch_samplesheet)
    ch_cat_fastq = CAT_FASTQ.out.reads
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    // Raw read statistics
    ch_seqkit_sequences = ch_cat_fastq
        .map { it[1] }
        .collect()
        .map { seqs -> [[id: "seqkit_stat_raw"], seqs] }

    SEQKIT_STATS_RAW(ch_seqkit_sequences)

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQ_TRIM_FASTP_FASTQC(
        ch_cat_fastq,
        [],
        params.save_trimmed_fail,
        params.discard_trimmed_pass,
        params.save_merged,
        false,
        false,
    )
    ch_variants_fastq = FASTQ_TRIM_FASTP_FASTQC.out.reads
    ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)

    // Trimmed read statistics
    ch_seqkit_sequences_trim = ch_variants_fastq
        .map { it[1] }
        .collect()
        .map { seqs -> [[id: "seqkit_stat_trim"], seqs] }
    SEQKIT_STATS_TRIM(ch_seqkit_sequences_trim)

    //
    // Filter empty FastQ files after adapter trimming
    //
    def fail_mapped_reads = [:]
    if (!params.skip_fastp) {
        ch_variants_fastq
            .join(FASTQ_TRIM_FASTP_FASTQC.out.trim_json)
            .map { meta, reads, json ->
                def pass = WorkflowWasteflow.getFastpReadsAfterFiltering(json) > 0
                [meta, reads, json, pass]
            }
            .set { ch_pass_fail_reads }

        ch_pass_fail_reads
            .map { meta, reads, json, pass ->
                if (pass) {
                    [meta, reads]
                }
            }
            .set { ch_variants_fastq }

        ch_pass_fail_reads
            .map { meta, reads, json, pass ->
                if (!pass) {
                    fail_mapped_reads[meta.id] = 0
                    def num_reads = WorkflowWasteflow.getFastpReadsBeforeFiltering(json)
                    return ["${meta.id}\t${num_reads}"]
                }
            }
            .collect()
            .map { tsv_data ->
                def header = ['Sample', 'Reads before trimming']
                WorkflowCommons.multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_reads_multiqc }
    }

    emit:
    reads              = ch_variants_fastq
    fail_reads_multiqc = ch_fail_reads_multiqc
    fastqc_raw_html    = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_html
    fastqc_raw_zip     = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_zip
    fastqc_trim_html   = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_trim_html
    fastqc_trim_zip    = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_trim_zip
    trim_json          = FASTQ_TRIM_FASTP_FASTQC.out.trim_json
    trim_html          = FASTQ_TRIM_FASTP_FASTQC.out.trim_html
    trim_log           = FASTQ_TRIM_FASTP_FASTQC.out.trim_log
    versions           = ch_versions
}
