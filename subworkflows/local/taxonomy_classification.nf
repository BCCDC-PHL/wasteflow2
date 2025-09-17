#!/usr/bin/env nextflow
// import the modules
include { KRAKEN2_BUILD                      } from '../../modules/local/kraken2_build'
include { KRAKEN2_KRAKEN2                    } from '../../modules/nf-core/kraken2/kraken2/main'
include { UNTAR as UNTAR_KRAKEN2_DB          } from '../../modules/nf-core/untar/main'
include { SEQKIT_STATS as SEQKIT_STATS_TAXQC } from '../../modules/nf-core/seqkit/stats/main'
// import the subworkflows
include { TAXONOMY_QC                        } from './taxonomy_qc'
include { EXTRACT_TARGETS                    } from './extract_targets'

workflow TAXONOMY_CLASSIFICATION {
    take:
    ch_reads // channel: [meta, reads]

    main:
    ch_versions = Channel.empty()
    ch_kraken2_multiqc = Channel.empty()
    ch_taxonomy_qc = ch_reads

    // Initialize individual channels
    ch_sars_cov2 = Channel.empty()
    ch_rsv = Channel.empty()
    ch_flu = Channel.empty()

    if (!params.skip_taxonomy_qc) {
        TAXONOMY_QC(ch_taxonomy_qc)
        ch_kraken2_multiqc = TAXONOMY_QC.out.kraken2_report
        ch_versions = ch_versions.mix(TAXONOMY_QC.out.versions)
    }

    if (!params.skip_target_extraction) {
        def target_taxids = params.kraken2_targets ? params.kraken2_targets.split(',').collect { it.trim() } : []
        ch_assignment = TAXONOMY_QC.out.kraken2_classified_reads_assignment
        ch_classified = TAXONOMY_QC.out.kraken2_classified_reads
        ch_report = TAXONOMY_QC.out.kraken2_report

        ch_joined = ch_assignment
            .join(ch_classified, by: 0)
            .join(ch_report, by: 0)
            .map { meta, assignment, classified, report ->
                tuple(meta, assignment, classified, report)
            }

        // Fan out each sample output to each taxid
        ch_target_inputs = ch_joined.flatMap { meta, assignment, classified, report ->
            target_taxids.collect { taxid ->
                if (taxid == null || taxid == '') {
                    log.warn("Found null or empty taxid for sample ${meta.id}")
                }
                tuple(
                    taxid,
                    tuple(meta, assignment),
                    tuple(meta, classified),
                    tuple(meta, report),
                )
            }
        }

        EXTRACT_TARGETS(ch_target_inputs)
        ch_versions = ch_versions.mix(EXTRACT_TARGETS.out.versions)

        // Get stats of all extracted reads before branching
        ch_seqkit_sequences = EXTRACT_TARGETS.out.extracted_reads
            .map { meta, reads -> reads }
            .collect()
            .map { seqs -> [[id: "seqkit_stat_all_extracted"], seqs] }

        SEQKIT_STATS_TAXQC(ch_seqkit_sequences)
        ch_versions = ch_versions.mix(SEQKIT_STATS_TAXQC.out.versions)

        // Branch extracted reads by pathogen and assign to individual channels
        EXTRACT_TARGETS.out.extracted_reads
            .branch { meta, reads ->
                sars_cov2: meta.taxid == '2697049'
                return tuple(meta, reads)
                rsv: meta.taxid == '3049954'
                return tuple(meta, reads)
                flu: meta.taxid == '11308'
                return tuple(meta, reads)
            }
            .set { branched_channels }

        // Assign branched channels to individual variables
        ch_sars_cov2 = branched_channels.sars_cov2
        ch_rsv = branched_channels.rsv
        ch_flu = branched_channels.flu
    }

    emit:
    sars_cov2_reads = ch_sars_cov2 // channel: [meta, reads] - SARS-CoV-2 reads
    rsv_reads       = ch_rsv // channel: [meta, reads] - RSV reads
    flu_reads       = ch_flu // channel: [meta, reads] - Influenza reads
    kraken2_multiqc = ch_kraken2_multiqc // channel: [meta, kraken2_report] - Kraken2 multiqc report
    versions        = ch_versions // channel: [meta, version] - workflow versions
}
