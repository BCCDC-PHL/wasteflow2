#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KRAKENTOOLS_EXTRACTKRAKENREADS } from '../../modules/nf-core/krakentools/extractkrakenreads/main'

workflow EXTRACT_TARGETS {
    take:
    input_data

    main:
    ch_versions = Channel.empty()

    ch_mapped_inputs = input_data.map { taxid, assignment, classified, report ->
        def (meta_a, path_assignment) = assignment
        def (meta_c, path_classified) = classified
        def (meta_r, path_report) = report

        // Add taxid to ALL metadata maps
        def new_meta_a = meta_a.clone()
        new_meta_a.taxid = taxid

        def new_meta_c = meta_c.clone()
        new_meta_c.taxid = taxid

        def new_meta_r = meta_r.clone()
        new_meta_r.taxid = taxid


        tuple(
            taxid,
            tuple(new_meta_a, path_assignment),
            tuple(new_meta_c, path_classified),
            tuple(new_meta_r, path_report),
        )
    }

    taxid_ch = ch_mapped_inputs.map { it[0] }
    assignment_ch = ch_mapped_inputs.map { it[1] }
    classified_ch = ch_mapped_inputs.map { it[2] }
    report_ch = ch_mapped_inputs.map { it[3] }

    // Call the krakentools extractkrakenreads module
    KRAKENTOOLS_EXTRACTKRAKENREADS(
        taxid_ch,
        assignment_ch,
        classified_ch,
        report_ch,
    )
    ch_versions = ch_versions.mix(KRAKENTOOLS_EXTRACTKRAKENREADS.out.versions)

    emit:
    extracted_reads = KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_kraken2_reads
    versions        = KRAKENTOOLS_EXTRACTKRAKENREADS.out.versions
}
