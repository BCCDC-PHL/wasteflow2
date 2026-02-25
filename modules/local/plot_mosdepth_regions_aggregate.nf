process PLOT_MOSDEPTH_REGIONS_AGG {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a7/a7a2a97d06bed60d93000b0beef2fb15b80454425d764ef8408e3266dacf8046/data'
        : 'community.wave.seqera.io/library/matplotlib-base_seaborn-base_pip_matplotlib_pruned:14ed6c6e196df244'}"

    input:
    tuple val(meta), path(coverage_file)
    val pathogen

    output:
    tuple val(meta), path("*.png"), emit: boxplots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${params.run_id}"

    """
    plot_mosdepth_regions_compare.py \\
        --input ${coverage_file} \\
        --pathogen ${pathogen} \\
        --prefix ${prefix} 
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sars-cov-2_coverage_boxplot_EM.png
    touch ${prefix}_sars-cov-2_coverage_boxplot_research.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
