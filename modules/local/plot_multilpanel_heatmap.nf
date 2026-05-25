process PLOT_MULTIPANEL_COVERAGE_HEATMAP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a7/a7a2a97d06bed60d93000b0beef2fb15b80454425d764ef8408e3266dacf8046/data'
        : 'community.wave.seqera.io/library/matplotlib-base_seaborn-base_pip_matplotlib_pruned:14ed6c6e196df244'}"

    input:
    tuple val(meta), path(sars_coverage)
    tuple val(meta2), path(rsv_a_coverage)
    tuple val(meta3), path(rsv_b_coverage)
    tuple val(meta4), path(h1n1_coverage)
    tuple val(meta5), path(h3n2_coverage)
    tuple val(meta6), path(h5n1_coverage)
    tuple val(meta7), path(flu_b_vic_coverage)
    tuple val(meta8), path(measles_coverage)
    tuple val(meta), path(sars_bed)
    tuple val(meta2), path(rsv_a_bed)
    tuple val(meta3), path(rsv_b_bed)
    tuple val(meta4), path(h1n1_bed)
    tuple val(meta5), path(h3n2_bed)
    tuple val(meta6), path(h5n1_bed)
    tuple val(meta7), path(flu_b_vic_bed)
    tuple val(meta8), path(measles_bed)
    path metadata

    output:
    tuple val(meta), path("*_sarscov2_rsv_measles.png"), emit: heatmap_sars
    tuple val(meta4), path("*_influenza_A.png"), emit: heatmap_influenza_A
    tuple val(meta7), path("*_influenza_B.png"), emit: heatmap_influenza_B
    tuple val(meta), path("*_sarscov2_rsv_measles.tsv"), emit: tsv_sars
    tuple val(meta4), path("*_influenza_A.tsv"), emit: tsv_influenza_A
    tuple val(meta7), path("*_influenza_B.tsv"), emit: tsv_influenza_B

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metadata_arg = metadata && metadata.size() > 0 && !metadata.toString().contains('[]') ? "--metadata ${metadata}" : ""
    def sars_arg = sars_coverage ? "--sars ${sars_coverage}" : ""
    def rsv_a_arg = rsv_a_coverage ? "--rsvA ${rsv_a_coverage}" : ""
    def rsv_b_arg = rsv_b_coverage ? "--rsvB ${rsv_b_coverage}" : ""
    def h1n1_arg = h1n1_coverage ? "--h1n1 ${h1n1_coverage}" : ""
    def h3n2_arg = h3n2_coverage ? "--h3n2 ${h3n2_coverage}" : ""
    def h5n1_arg = h5n1_coverage ? "--h5n1 ${h5n1_coverage}" : ""
    def flu_b_vic_arg = flu_b_vic_coverage ? "--flu_b_vic ${flu_b_vic_coverage}" : ""
    def measles_arg = measles_coverage ? "--measles ${measles_coverage}" : ""
    def strict_arg = (params.strict_sample_names in [true, "true", "True"]) ? "--strict_sample_names" : ""
    

    """
    multi-panel_gene-level_heatmap.py \\
        ${sars_arg} \\
        ${rsv_a_arg} \\
        ${rsv_b_arg} \\
        ${h1n1_arg} \\
        ${h3n2_arg} \\
        ${h5n1_arg} \\
        ${flu_b_vic_arg} \\
        ${measles_arg} \\
        --sars_bed ${sars_bed} \\
        --rsvA_bed ${rsv_a_bed} \\
        --rsvB_bed ${rsv_b_bed} \\
        --h1n1_bed ${h1n1_bed} \\
        --h3n2_bed ${h3n2_bed} \\
        --h5n1_bed ${h5n1_bed} \\
        --flu_b_vic_bed ${flu_b_vic_bed} \\
        --measles_bed ${measles_bed} \\
        ${args} \\
        ${metadata_arg} \\
        ${strict_arg} \\
        --out_sars_rsv_measles ${prefix}_sarscov2_rsv_measles.png \\
        --out_influenza_a ${prefix}_influenza_A.png \\
        --out_influenza_b ${prefix}_influenza_B.png \\
        --out_sars_rsv_measles_tsv ${prefix}_sarscov2_rsv_measles.tsv \\
        --out_influenza_a_tsv ${prefix}_influenza_A.tsv \\
        --out_influenza_b_tsv ${prefix}_influenza_B.tsv 
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${params.run_id}"
    """
    touch ${prefix}_coverage_heatmap.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
    END_VERSIONS
    """
}
