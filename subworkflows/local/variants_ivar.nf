//
// Variant calling with IVar, downstream processing and QC
//

include { IVAR_VARIANTS        } from '../../modules/nf-core/ivar/variants/main'
include { IVAR_VARIANTS_TO_VCF } from '../../modules/local/ivar_variants_to_vcf'
include { BCFTOOLS_SORT        } from '../../modules/nf-core/bcftools/sort/main'
include { VCF_TABIX_STATS      } from './vcf_tabix_stats'

workflow VARIANTS_IVAR {
    take:
    bam                 // channel: [ val(meta), [ bam ] ]
    fasta               // channel: /path/to/genome.fasta
    fai                 // channel: /path/to/genome.fai
    sizes               // channel: /path/to/genome.sizes
    gff                 // channel: /path/to/genome.gff
    ivar_multiqc_header // channel: /path/to/multiqc_header for ivar variants

    main:

    ch_versions = Channel.empty()

    //
    // Call variants
    //
    IVAR_VARIANTS(
        bam,
        fasta,
        fai,
        gff,
        params.save_mpileup,
    )
    ch_versions = ch_versions.mix(IVAR_VARIANTS.out.versions.first())


    IVAR_VARIANTS.out.tsv
    .branch { meta, tsv ->
        pass: WorkflowCommons.getNumLinesInFile(tsv) > 1
            return [meta, tsv]
        fail: true
            return [meta, tsv]
    }
    .set { ch_ivar_tsv_branched }

    // Use the passing samples for downstream processing
    def ch_ivar_tsv = ch_ivar_tsv_branched.pass

    // For segmented genomes, use the matched_fasta from meta
    // For non-segmented genomes, pair with the broadcast fasta
    ch_ivar_tsv
        .branch { meta, tsv ->
            segmented: meta.containsKey('matched_fasta')
                return [meta, tsv, meta.matched_fasta]
            full_genome: true
                return [meta, tsv]
        }
        .set { ch_tsv_branched }

    // For full genomes, add the broadcast fasta
    ch_tsv_branched.full_genome
        .combine(fasta.map { it })  
        .map { meta, tsv, fa ->
            [meta, tsv, fa]
        }
        .set { ch_full_genome_with_fasta }

    // Combine both branches
    ch_tsv_branched.segmented
        .mix(ch_full_genome_with_fasta)
        .map { meta, tsv, fa ->
            [meta, tsv, fa]
        }
        .set { ch_tsv_fasta_matched }

    //
    // Convert original iVar output to VCF, zip and index
    //
    IVAR_VARIANTS_TO_VCF(
        ch_tsv_fasta_matched.map { meta, tsv, fa -> [meta, tsv] },
        ch_tsv_fasta_matched.map { meta, tsv, fa -> fa },
        ivar_multiqc_header,
    )
    ch_versions = ch_versions.mix(IVAR_VARIANTS_TO_VCF.out.versions.first())

    BCFTOOLS_SORT(
        IVAR_VARIANTS_TO_VCF.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    VCF_TABIX_STATS(
        BCFTOOLS_SORT.out.vcf,
        [[:], []],
        [[:], []],
        [[:], []],
    )
    ch_versions = ch_versions.mix(VCF_TABIX_STATS.out.versions)

    emit:
    tsv          = ch_ivar_tsv // channel: [ val(meta), [ tsv ] ]
    vcf_orig     = IVAR_VARIANTS_TO_VCF.out.vcf // channel: [ val(meta), [ vcf ] ]
    log_out      = IVAR_VARIANTS_TO_VCF.out.log // channel: [ val(meta), [ log ] ]
    multiqc_tsv  = IVAR_VARIANTS_TO_VCF.out.tsv // channel: [ val(meta), [ tsv ] ]
    vcf          = BCFTOOLS_SORT.out.vcf // channel: [ val(meta), [ vcf ] ]
    tbi          = VCF_TABIX_STATS.out.tbi // channel: [ val(meta), [ tbi ] ]
    csi          = VCF_TABIX_STATS.out.csi // channel: [ val(meta), [ csi ] ]
    stats        = VCF_TABIX_STATS.out.stats // channel: [ val(meta), [ txt ] ]
    versions     = ch_versions // channel: [ versions.yml ]
}
