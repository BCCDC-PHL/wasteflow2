//
// Importing the modules required for the sub-workflow
//

include { KRAKEN2_BUILD                  } from '../../modules/nf-core/kraken2/build/main'
include { KRAKEN2_KRAKEN2                } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_EXTRACTKRAKENREADS } from '../../modules/nf-core/krakentools/extractkrakenreads/main'
include { KRAKENTOOLS_COMBINEKREPORTS    } from '../../modules/nf-core/krakentools/combinekreports/main'
include { KRAKENTOOLS_KREPORT2KRONA      } from '../../modules/nf-core/krakentools/kreport2krona/main'
include { BRACKEN_BUILD                  } from '../../modules/nf-core/bracken/build/main'
include { BRACKEN_BRACKEN                } from '../../modules/nf-core/bracken/bracken/main'
include { KRONA_KTUPDATETAXONOMY         } from '../../modules/nf-core/krona/ktupdatetaxonomy/main'
include { KRONA_KTIMPORTTAXONOMY         } from '../../modules/nf-core/krona/ktimporttaxonomy/main'
include { KRONA_KTIMPORTTEXT             } from '../../modules/nf-core/krona/ktimporttext/main'
include { UNTAR as UNTAR_KRAKEN2_DB      } from '../../modules/nf-core/untar/main'


workflow TAXONOMY_QC {
    take:
    reads

    main:
    ch_versions = Channel.empty()
    kraken2_report = Channel.empty()
    bracken_report = Channel.empty()

    //
    // Validate database paths
    //

    if (!params.skip_kraken2) {
        if (params.kraken2_db == null || !Utils.fileExists(params.kraken2_db)) {
            log.error("Path to Kraken2 database is not valid")
            exit(1)
        }
    }

    // Taxonomic classification
    save_output_fastqs = params.kraken2_variants_host_filter
    save_reads_assignment = params.kraken2_assembly_host_filter

    if (!params.skip_kraken2) {
        KRAKEN2_KRAKEN2(
            reads,
            params.kraken2_db,
            save_output_fastqs,
            save_reads_assignment,
        )
        kraken2_report = KRAKEN2_KRAKEN2.out.report
        ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)


        KRAKENTOOLS_KREPORT2KRONA(
            kraken2_report
        )
        ch_versions = ch_versions.mix(KRAKENTOOLS_KREPORT2KRONA.out.versions)
        ch_krona_txt = KRAKENTOOLS_KREPORT2KRONA.out.txt

        // If combining reports, combine here
        if (params.combine_kraken2_reports) {
            ch_combined = ch_krona_txt
                .map { it[1] }
                .collect()
                .map { files ->
                    // Join the files with single quotes and space
                    // String joinedFiles = files.collect { "'$it'" }.join(' ')
                    // if single file then make it [files] otherwise just files
                    [[id: 'combined_krona'], files instanceof List ? files : [files]]
                }
            ch_krona_txt = ch_combined
        }

        // Create Krona HTML
        KRONA_KTIMPORTTEXT(
            ch_krona_txt
        )
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions)
    }

    // Bracken abundance estimation
    if (!params.skip_bracken) {
        BRACKEN_BRACKEN(
            KRAKEN2_KRAKEN2.out.report,
            params.kraken2_db,
        )
        bracken_report = BRACKEN_BRACKEN.out.reports
        ch_versions = ch_versions.mix(BRACKEN_BRACKEN.out.versions)
    }

    emit:
    kraken2_report
    kraken2_classified_reads            = KRAKEN2_KRAKEN2.out.classified_reads_fastq
    kraken2_unclassified_reads          = KRAKEN2_KRAKEN2.out.unclassified_reads_fastq
    kraken2_classified_reads_assignment = KRAKEN2_KRAKEN2.out.classified_reads_assignment
    bracken_report
    versions                            = ch_versions // channel: [ versions.yml ]
}
