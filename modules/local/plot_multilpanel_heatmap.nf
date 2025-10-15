process PLOT_MULTIPANEL_COVERAGE_HEATMAP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a7/a7a2a97d06bed60d93000b0beef2fb15b80454425d764ef8408e3266dacf8046/data':
        'community.wave.seqera.io/library/matplotlib-base_seaborn-base_pip_matplotlib_pruned:14ed6c6e196df244' }"

    input:
    tuple val(meta), path(sars_coverage) 
    tuple val(meta2), path(rsv_a_coverage)
    tuple val(meta3), path(rsv_b_coverage)
    tuple val(meta5), path(sars_bed)
    tuple val(meta6), path(rsv_a_bed)
    tuple val(meta7), path(rsv_b_bed)
    path(metadata)
    //val(batch_name)
    

    output:
    tuple val(meta), path("*.png"), emit: heatmap
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metadata_arg = (metadata && metadata.size() > 0 && !metadata.toString().contains('[]')) ? "--metadata ${metadata}" : ""
    def sars_arg = sars_coverage ? "--sars ${sars_coverage}" : ""
    def rsv_a_arg = rsv_a_coverage ? "--rsvA ${rsv_a_coverage}" : ""
    def rsv_b_arg = rsv_b_coverage ? "--rsvB ${rsv_b_coverage}" : ""

    
    """
    multi-panel_gene-level_heatmap.py \\
        ${sars_arg} \\
        ${rsv_a_arg} \\
        ${rsv_b_arg} \\
        --sars_bed ${sars_bed} \\
        --rsvA_bed ${rsv_a_bed} \\
        --rsvB_bed ${rsv_b_bed} \\
        ${args} \\
        ${metadata_arg} \\
        --out ${prefix}_coverage_heatmap.png 
        

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