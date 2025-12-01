//
// Uncompress and prepare reference genome files
//
/*
include { GUNZIP as GUNZIP_FASTA        } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF          } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRIMER_BED   } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRIMER_FASTA } from '../../modules/nf-core/gunzip/main'
include { UNTAR as UNTAR_BOWTIE2_INDEX  } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_NEXTCLADE_DB   } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_BLAST_DB       } from '../../modules/nf-core/untar/main'
include { BOWTIE2_BUILD                 } from '../../modules/nf-core/bowtie2/build/main'

include { BLAST_MAKEBLASTDB             } from '../../modules/nf-core/blast/makeblastdb/main'
include { BEDTOOLS_GETFASTA             } from '../../modules/nf-core/bedtools/getfasta/main'
include { CUSTOM_GETCHROMSIZES          } from '../../modules/nf-core/custom/getchromsizes/main'
include { CUSTOM_GETCHROMSIZES as CUSTOM_GETCHROMSIZES_SEGMENTS } from '../../modules/nf-core/custom/getchromsizes/main'

include { NEXTCLADE_DATASETGET          } from '../../modules/nf-core/nextclade/datasetget/main'
include { COLLAPSE_PRIMERS              } from '../../modules/local/collapse_primers'
include { SNPEFF_BUILD                  } from '../../modules/local/snpeff_build'
include { SNPEFF_BUILD as SNPEFF_BUILD_SEGMENT } from '../../modules/local/snpeff_build'

// New imports for segment handling
include { SPLIT_FASTA_BY_SEGMENT        } from '../../modules/local/split_fasta_by_segment'
include { SPLIT_GFF_BY_SEGMENT          } from '../../modules/local/split_gff_by_segment'



workflow PREPARE_GENOME {
    take:
    fasta
    gff
    primer_bed
    bowtie2_index
    nextclade_dataset
    nextclade_dataset_name
    nextclade_dataset_tag

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        GUNZIP_FASTA(
            [[:], fasta]
        )
        ch_fasta = GUNZIP_FASTA.out.gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    }
    else {
        ch_fasta = Channel.value(file(fasta))
    }

    //
    // Uncompress GFF annotation file
    //
    ch_gff = Channel.empty()
    if (gff) {
        if (gff.endsWith('.gz')) {
            GUNZIP_GFF(
                [[:], gff]
            )
            ch_gff = GUNZIP_GFF.out.gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        }
        else {
            ch_gff = Channel.value(file(gff))
        }
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES(
        ch_fasta.map { [[:], it] }
    )
    ch_fai = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)


    //
    // Prepare reference files required for Alignment
    //
    ch_bowtie2_index = Channel.empty()
    if (bowtie2_index) {
        if (bowtie2_index.endsWith('.tar.gz')) {
            UNTAR_BOWTIE2_INDEX(
                [[:], file(bowtie2_index)]
            )
            ch_bowtie2_index = UNTAR_BOWTIE2_INDEX.out.untar
            ch_versions = ch_versions.mix(UNTAR_BOWTIE2_INDEX.out.versions)
        }
        else {
            ch_bowtie2_index = [[:], file(bowtie2_index)]
        }
    }
    else {
        BOWTIE2_BUILD(
            ch_fasta.map { [[:], it] }
        )
        ch_bowtie2_index = BOWTIE2_BUILD.out.index
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }


    //
    // Prepare Nextclade dataset
    //
    ch_nextclade_db = Channel.empty()
    if (!params.skip_nextclade) {
        if (nextclade_dataset) {
            if (nextclade_dataset.endsWith('.tar.gz')) {
                UNTAR_NEXTCLADE_DB(
                    [[:], nextclade_dataset]
                )
                ch_nextclade_db = UNTAR_NEXTCLADE_DB.out.untar.map { it[1] }
                ch_versions = ch_versions.mix(UNTAR_NEXTCLADE_DB.out.versions)
            }
            else {
                ch_nextclade_db = Channel.value(file(nextclade_dataset))
            }
        }
        else if (nextclade_dataset_name) {
            NEXTCLADE_DATASETGET(
                nextclade_dataset_name,
                nextclade_dataset_tag,
            )
            ch_nextclade_db = NEXTCLADE_DATASETGET.out.dataset
            ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
        }
    }

    //
    // Make snpEff database
    //
    ch_snpeff_db = Channel.empty()
    ch_snpeff_config = Channel.empty()
    if (!params.skip_snpeff) {
        SNPEFF_BUILD(
            ch_fasta,
            ch_gff,
        )
        ch_snpeff_db = SNPEFF_BUILD.out.db
        ch_snpeff_config = SNPEFF_BUILD.out.config
        ch_versions = ch_versions.mix(SNPEFF_BUILD.out.versions)
    }

    emit:
    fasta         = ch_fasta // path: genome.fasta
    gff           = ch_gff // path: genome.gff
    fai           = ch_fai // path: genome.fai
    chrom_sizes   = ch_chrom_sizes // path: genome.sizes
    bowtie2_index = ch_bowtie2_index // channel: [ [:], bowtie2/index/ ]
    nextclade_db  = ch_nextclade_db // path: nextclade_db
    snpeff_db     = ch_snpeff_db // path: snpeff_db
    snpeff_config = ch_snpeff_config // path: snpeff.config
    versions      = ch_versions // channel: [ versions.yml ]
}
*/

include { GUNZIP as GUNZIP_FASTA        } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF          } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRIMER_BED   } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRIMER_FASTA } from '../../modules/nf-core/gunzip/main'
include { UNTAR as UNTAR_BOWTIE2_INDEX  } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_NEXTCLADE_DB   } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_BLAST_DB       } from '../../modules/nf-core/untar/main'
include { BOWTIE2_BUILD                 } from '../../modules/nf-core/bowtie2/build/main'

include { BLAST_MAKEBLASTDB             } from '../../modules/nf-core/blast/makeblastdb/main'
include { BEDTOOLS_GETFASTA             } from '../../modules/nf-core/bedtools/getfasta/main'
include { CUSTOM_GETCHROMSIZES          } from '../../modules/nf-core/custom/getchromsizes/main'
include { CUSTOM_GETCHROMSIZES as CUSTOM_GETCHROMSIZES_SEGMENTS } from '../../modules/nf-core/custom/getchromsizes/main'

include { NEXTCLADE_DATASETGET          } from '../../modules/nf-core/nextclade/datasetget/main'
include { COLLAPSE_PRIMERS              } from '../../modules/local/collapse_primers'
include { SNPEFF_BUILD                  } from '../../modules/local/snpeff_build'
include { SNPEFF_BUILD as SNPEFF_BUILD_SEGMENTS } from '../../modules/local/snpeff_build'

// Segment handling
include { SPLIT_FASTA_BY_SEGMENT        } from '../../modules/local/split_fasta_by_segment'

workflow PREPARE_GENOME {
    take:
    virus_name               // Virus name (e.g., 'H1N1', 'SARS-CoV-2', etc.)
    fasta
    gff
    primer_bed
    bowtie2_index
    nextclade_dataset
    nextclade_dataset_name
    nextclade_dataset_tag
    segments_bed             // BED file with segment info (null for non-segmented)
    segment_gffs             // Map of segment names to GFF files (null for non-segmented)

    main:

    ch_versions = Channel.empty()

    //
    // Check if this is a segmented genome that needs splitting
    //
    def is_segmented = segments_bed && virus_name

    //
    // Uncompress genome fasta file if required
    //
    ch_fasta_input = Channel.empty()
    if (fasta.endsWith('.gz')) {
        GUNZIP_FASTA(
            [[:], fasta]
        )
        ch_fasta_input = GUNZIP_FASTA.out.gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    }
    else {
        ch_fasta_input = Channel.value(file(fasta))
    }

    //
    // Uncompress GFF annotation file (for non-segmented genomes)
    //
    ch_gff_input = Channel.empty()
    if (gff && !is_segmented) {
        if (gff.endsWith('.gz')) {
            GUNZIP_GFF(
                [[:], gff]
            )
            ch_gff_input = GUNZIP_GFF.out.gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        }
        else {
            ch_gff_input = Channel.value(file(gff))
        }
    }

    //
    // Initialize output channels
    //
    ch_fasta_combined = Channel.empty()
    ch_fai_combined = Channel.empty()
    ch_chrom_sizes_combined = Channel.empty()
    ch_bowtie2_index_combined = Channel.empty()
    ch_gff_output = Channel.empty()
    ch_nextclade_db = Channel.empty()
    ch_snpeff_db = Channel.empty()
    ch_snpeff_config = Channel.empty()
    
    ch_fasta_segments = Channel.empty()
    ch_gff_segments = Channel.empty()
    ch_fai_segments = Channel.empty()
    ch_chrom_sizes_segments = Channel.empty()
    ch_snpeff_db_segments = Channel.empty()
    ch_snpeff_config_segments = Channel.empty()

    if (is_segmented) {
        //
        // SEGMENTED GENOME PROCESSING
        //
        
        // Step 1: Create chromosome sizes for COMBINED fasta (for alignment)
        CUSTOM_GETCHROMSIZES(
            ch_fasta_input.map { [[id: "${virus_name}_combined"], it] }
        )
        ch_fasta_combined = ch_fasta_input
        ch_fai_combined = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
        ch_chrom_sizes_combined = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
        ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

        // Step 2: Build Bowtie2 index for combined FASTA (for alignment)
        if (bowtie2_index) {
            if (bowtie2_index.endsWith('.tar.gz')) {
                UNTAR_BOWTIE2_INDEX(
                    [[:], file(bowtie2_index)]
                )
                ch_bowtie2_index_combined = UNTAR_BOWTIE2_INDEX.out.untar
                ch_versions = ch_versions.mix(UNTAR_BOWTIE2_INDEX.out.versions)
            }
            else {
                ch_bowtie2_index_combined = Channel.value([[:], file(bowtie2_index)])
            }
        }
        else {
            BOWTIE2_BUILD(
                ch_fasta_input.map { [[id: "${virus_name}_combined"], it] }
            )
            ch_bowtie2_index_combined = BOWTIE2_BUILD.out.index
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }

        // Step 3: Parse segment information and split FASTA
        ch_segments = Channel
            .fromPath(segments_bed)
            .splitCsv(sep: '\t', header: false)
            .map { row -> 
                [
                    virus: virus_name,
                    segment: row[3],
                    accession: row[0]
                ]
            }

        // Combine segments with the fasta file
        ch_fasta_for_split = ch_segments
            .combine(ch_fasta_input)
            .map { segment_info, fasta_file -> 
                [segment_info.virus, fasta_file, segment_info.accession, segment_info.segment] 
            }

        SPLIT_FASTA_BY_SEGMENT(ch_fasta_for_split)        
        ch_fasta_segments = SPLIT_FASTA_BY_SEGMENT.out.fasta
            .map { virus_val, segment_val, fasta_file ->
                [[virus: virus_val, segment: segment_val, id: "${virus_val}_${segment_val}"], fasta_file]
            }
        
        ch_versions = ch_versions.mix(SPLIT_FASTA_BY_SEGMENT.out.versions.first())

        // Step 4: Load segment GFF files from params map
        ch_gff_segments = Channel.empty()
        if (segment_gffs) {
            ch_gff_segments = ch_segments
                .map { segment_info ->
                    def gff_file = segment_gffs[segment_info.segment]
                    if (gff_file == null) {
                        error "No GFF file found for segment ${segment_info.segment} in virus ${virus_name}"
                    }
                    [[virus: virus_name, segment: segment_info.segment, id: "${virus_name}_${segment_info.segment}"], file(gff_file)]
                }
        }

        // Step 5: Create chromosome sizes for each segment
        CUSTOM_GETCHROMSIZES_SEGMENTS(
            ch_fasta_segments
        )
        ch_fai_segments = CUSTOM_GETCHROMSIZES_SEGMENTS.out.fai
        ch_chrom_sizes_segments = CUSTOM_GETCHROMSIZES_SEGMENTS.out.sizes
        ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES_SEGMENTS.out.versions)

        // Step 6: Build SnpEff databases for segments (reusing SNPEFF_BUILD)
        if (!params.skip_snpeff && segment_gffs) {
            ch_fasta_gff_segments = ch_fasta_segments
                .join(ch_gff_segments, by: [0])

            SNPEFF_BUILD_SEGMENTS(
                ch_fasta_gff_segments.map { meta, fasta_file, gff_file -> fasta_file },
                ch_fasta_gff_segments.map { meta, fasta_file, gff_file -> gff_file }
            )
            
            ch_snpeff_db_segments = SNPEFF_BUILD_SEGMENTS.out.db
            ch_snpeff_config_segments = SNPEFF_BUILD_SEGMENTS.out.config
            ch_versions = ch_versions.mix(SNPEFF_BUILD_SEGMENTS.out.versions.first())
        }

    } else {
        //
        // NON-SEGMENTED GENOME PROCESSING (SARS-CoV-2, RSV-A, RSV-B, etc.)
        //
        
        ch_fasta_to_process = ch_fasta_input.map { [[:], it] }
        ch_gff_to_process = ch_gff_input.map { [[:], it] }

        // Create chromosome sizes
        CUSTOM_GETCHROMSIZES(
            ch_fasta_to_process
        )
        ch_fasta_combined = ch_fasta_input
        ch_fai_combined = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
        ch_chrom_sizes_combined = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
        ch_gff_output = ch_gff_input
        ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

        // Build Bowtie2 index
        if (bowtie2_index) {
            if (bowtie2_index.endsWith('.tar.gz')) {
                UNTAR_BOWTIE2_INDEX(
                    [[:], file(bowtie2_index)]
                )
                ch_bowtie2_index_combined = UNTAR_BOWTIE2_INDEX.out.untar
                ch_versions = ch_versions.mix(UNTAR_BOWTIE2_INDEX.out.versions)
            }
            else {
                ch_bowtie2_index_combined = Channel.value([[:], file(bowtie2_index)])
            }
        }
        else {
            BOWTIE2_BUILD(
                ch_fasta_to_process
            )
            ch_bowtie2_index_combined = BOWTIE2_BUILD.out.index
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }

        // Prepare Nextclade dataset
        if (!params.skip_nextclade) {
            if (nextclade_dataset) {
                if (nextclade_dataset.endsWith('.tar.gz')) {
                    UNTAR_NEXTCLADE_DB(
                        [[:], nextclade_dataset]
                    )
                    ch_nextclade_db = UNTAR_NEXTCLADE_DB.out.untar.map { it[1] }
                    ch_versions = ch_versions.mix(UNTAR_NEXTCLADE_DB.out.versions)
                }
                else {
                    ch_nextclade_db = Channel.value(file(nextclade_dataset))
                }
            }
            else if (nextclade_dataset_name) {
                NEXTCLADE_DATASETGET(
                    nextclade_dataset_name,
                    nextclade_dataset_tag,
                )
                ch_nextclade_db = NEXTCLADE_DATASETGET.out.dataset
                ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
            }
        }

        // Make snpEff database
        if (!params.skip_snpeff && gff) {
            SNPEFF_BUILD(
                ch_fasta_to_process.map { meta, fasta_file -> fasta_file },
                ch_gff_to_process.map { meta, gff_file -> gff_file },
            )
            ch_snpeff_db = SNPEFF_BUILD.out.db
            ch_snpeff_config = SNPEFF_BUILD.out.config
            ch_versions = ch_versions.mix(SNPEFF_BUILD.out.versions)
        }
    }

    emit:
    // Standard outputs (always present)
    fasta         = ch_fasta_combined              // path: genome.fasta (combined for alignment)
    fai           = ch_fai_combined                // path: genome.fai
    chrom_sizes   = ch_chrom_sizes_combined        // path: genome.sizes
    bowtie2_index = ch_bowtie2_index_combined      // channel: [ [:], bowtie2/index/ ]
    gff           = ch_gff_output                  // path: genome.gff (only for non-segmented)
    nextclade_db  = ch_nextclade_db                // path: nextclade_db (only for non-segmented)
    snpeff_db     = ch_snpeff_db                   // path: snpeff_db (only for non-segmented)
    snpeff_config = ch_snpeff_config               // path: snpeff.config (only for non-segmented)
    
    // Segment-specific outputs (only for segmented genomes)
    segment_fasta        = ch_fasta_segments           // channel: [ [meta], segment.fasta ]
    segment_gff          = ch_gff_segments             // channel: [ [meta], segment.gff ]
    segment_fai          = ch_fai_segments             // channel: [ [meta], segment.fai ]
    segment_chrom_sizes  = ch_chrom_sizes_segments     // channel: [ [meta], segment.sizes ]
    segment_snpeff_db    = ch_snpeff_db_segments       // channel: [ path ]
    segment_snpeff_config = ch_snpeff_config_segments  // channel: [ path ]
    
    versions      = ch_versions                    // channel: [ versions.yml ]
}