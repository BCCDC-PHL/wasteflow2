//
// This file holds several functions specific to the workflow/wasteflow.nf in the BCCDC-PHL/wasteflow2 pipeline
//
import nextflow.Nextflow
import groovy.json.JsonSlurper

class WorkflowWasteflow {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {
        WorkflowCommons.genomeExistsError(params, log)

        // Generic parameter validation
        

        if (!params.fasta) {
            Nextflow.error("Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file.")
        }

        if (!params.skip_kraken2 && !params.kraken2_db) {
            if (!params.kraken2_db_name) {
                Nextflow.error("Please specify a valid name to build Kraken2 database for host e.g. '--kraken2_db_name human'.")
            }
        }

    }

    //
    // Print warning if genome fasta has more than one sequence
    //
    public static void isMultiFasta(fasta_file, log) {
        def count = 0
        def line  = null
        fasta_file.withReader { reader ->
            while (line = reader.readLine()) {
                if (line.contains('>')) {
                    count++
                    if (count > 1) {
                        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                            "  This pipeline does not officially support multi-fasta genome files!\n\n" +
                            "  The parameters and processes are tailored for viral genome analysis.\n" +
                            "  Please amend the '--fasta' parameter.\n" +
                            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                        break
                    }
                }
            }
        }
    }

    //
    // Function that parses and returns the number of mapped reasds from flagstat files
    //
    public static ArrayList getFlagstatMappedReads(flagstat_file, params) {
        def mapped_reads = 0
        flagstat_file.eachLine { line ->
            if (line.contains(' mapped (')) {
                mapped_reads = line.tokenize().first().toInteger()
            }
        }

        def pass = false
        def logname = flagstat_file.getBaseName() - 'flagstat'
        if (mapped_reads > params.min_mapped_reads.toInteger()) {
            pass = true
        }
        return [ mapped_reads, pass ]
    }


    //
    // Function that parses fastp json output file to get total number of reads after trimming
    //
    public static Integer getFastpReadsAfterFiltering(json_file) {
        def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
        return json['after_filtering']['total_reads'].toInteger()
    }

    //
    // Function that parses fastp json output file to get total number of reads before trimming
    //
    public static Integer getFastpReadsBeforeFiltering(json_file) {
        def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
        return json['before_filtering']['total_reads'].toInteger()
    }
}
