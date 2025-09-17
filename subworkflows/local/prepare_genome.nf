//
// Uncompress and prepare reference genome files
//

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
include { NEXTCLADE_DATASETGET          } from '../../modules/nf-core/nextclade/datasetget/main'
include { COLLAPSE_PRIMERS              } from '../../modules/local/collapse_primers'
include { SNPEFF_BUILD                  } from '../../modules/local/snpeff_build'

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
