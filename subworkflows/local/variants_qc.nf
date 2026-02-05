include { SNPEFF_SNPSIFT } from './snpeff_snpsift'

workflow VARIANTS_QC {
    take:
    vcf                // channel: [ val(meta), path(vcf) ]
    all_snpeff_db      // channel: [ val(meta), path(db) ] - from SNPEFF_BUILD
    all_snpeff_config  // channel: [ val(meta), path(config) ] - from SNPEFF_BUILD
    all_fasta          // channel: [ val(meta), path(fasta) ] - from segment processing

    main:

    ch_versions = channel.empty()

    //
    // Annotate variants
    //
    ch_snpeff_vcf = channel.empty()
    ch_snpeff_tbi = channel.empty()
    ch_snpeff_stats = channel.empty()
    ch_snpeff_csv = channel.empty()
    ch_snpeff_txt = channel.empty()
    ch_snpeff_html = channel.empty()
    ch_snpsift_txt = channel.empty()
    
    if (!params.skip_snpeff) {

        // Materialize all reference channels as Maps
        all_snpeff_db
            .map { meta, db ->
                def segment_key = (meta.segment == 'none' ? 'null' : meta.segment)
                def key = "${meta.genome}::${segment_key}"
                [key, db]
            }
            .set { ch_db_keyed }
        
        all_snpeff_config
            .map { meta, config ->
                def segment_key = (meta.segment == 'none' ? 'null' : meta.segment)
                def key = "${meta.genome}::${segment_key}"
                [key, config]
            }
            .set { ch_config_keyed }
        
        all_fasta
            .map { meta, fasta ->
                def segment_key = (meta.segment == 'none' ? 'null' : meta.segment)
                def key = "${meta.genome}::${segment_key}"
                [key, fasta]
            }
            .set { ch_fasta_keyed }
        
        // Prepare VCF channel
        vcf
            .map { meta, vcf_file ->
                def segment_key = (meta.segment == 'none' ? 'null' : meta.segment)
                def key = "${meta.genome}::${segment_key}"
                [key, meta, vcf_file]
            }
            .combine(ch_db_keyed, by: 0)
            .combine(ch_config_keyed, by: 0)
            .combine(ch_fasta_keyed, by: 0)
            .map { key, meta, vcf_file, db, config, fasta ->
                [meta, vcf_file, db, config, fasta]
            }
            .set { ch_vcf_with_all_refs }

        SNPEFF_SNPSIFT(ch_vcf_with_all_refs)
        
        // Collect outputs
        ch_snpeff_vcf = SNPEFF_SNPSIFT.out.vcf
        ch_snpeff_tbi = SNPEFF_SNPSIFT.out.tbi
        ch_snpeff_stats = SNPEFF_SNPSIFT.out.stats
        ch_snpeff_csv = SNPEFF_SNPSIFT.out.csv
        ch_snpeff_txt = SNPEFF_SNPSIFT.out.txt
        ch_snpeff_html = SNPEFF_SNPSIFT.out.html
        ch_snpsift_txt = SNPEFF_SNPSIFT.out.snpsift_txt
        ch_versions = ch_versions.mix(SNPEFF_SNPSIFT.out.versions)
    }

    emit:
    snpeff_vcf   = ch_snpeff_vcf   // channel: [ val(meta), path(vcf.gz) ]
    snpeff_tbi   = ch_snpeff_tbi   // channel: [ val(meta), path(tbi) ]
    snpeff_stats = ch_snpeff_stats // channel: [ val(meta), path(txt) ]
    snpeff_csv   = ch_snpeff_csv   // channel: [ val(meta), path(csv) ]
    snpeff_txt   = ch_snpeff_txt   // channel: [ val(meta), path(txt) ]
    snpeff_html  = ch_snpeff_html  // channel: [ val(meta), path(html) ]
    snpsift_txt  = ch_snpsift_txt  // channel: [ val(meta), path(txt) ]
    versions     = ch_versions     // channel: [ versions.yml ]
}