include { SNPEFF_ANN            } from '../../modules/local/snpeff_ann'
include { SNPSIFT_EXTRACTFIELDS } from '../../modules/local/snpsift_extractfields'
include { VCF_BGZIP_TABIX_STATS } from './vcf_bgzip_tabix_stats'

workflow SNPEFF_SNPSIFT {
    take:
    vcf_with_refs // channel: [ val(meta), vcf, snpeff_db, snpeff_config, fasta ]

    main:

    ch_versions = channel.empty()

  

    // Call SNPEFF_ANN with renamed files
    SNPEFF_ANN(vcf_with_refs)
    ch_versions = ch_versions.mix(SNPEFF_ANN.out.versions.first())

    // Compress, index, and stats
    VCF_BGZIP_TABIX_STATS(
        SNPEFF_ANN.out.vcf,
        [[:], []],
        [[:], []],
        [[:], []],
    )
    ch_versions = ch_versions.mix(VCF_BGZIP_TABIX_STATS.out.versions)

    // Extract fields with SnpSift
    SNPSIFT_EXTRACTFIELDS(
        VCF_BGZIP_TABIX_STATS.out.vcf
    )
    ch_versions = ch_versions.mix(SNPSIFT_EXTRACTFIELDS.out.versions.first())

    emit:
    csv         = SNPEFF_ANN.out.csv // channel: [ val(meta), [ csv ] ]
    txt         = SNPEFF_ANN.out.txt // channel: [ val(meta), [ txt ] ]
    html        = SNPEFF_ANN.out.html // channel: [ val(meta), [ html ] ]
    vcf         = VCF_BGZIP_TABIX_STATS.out.vcf // channel: [ val(meta), [ vcf.gz ] ]
    tbi         = VCF_BGZIP_TABIX_STATS.out.tbi // channel: [ val(meta), [ tbi ] ]
    csi         = VCF_BGZIP_TABIX_STATS.out.csi // channel: [ val(meta), [ csi ] ]
    stats       = VCF_BGZIP_TABIX_STATS.out.stats // channel: [ val(meta), [ txt ] ]
    snpsift_txt = SNPSIFT_EXTRACTFIELDS.out.txt // channel: [ val(meta), [ txt ] ]
    versions    = ch_versions // channel: [ versions.yml ]
}