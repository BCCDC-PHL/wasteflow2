process PARSE_PAF_SEROTYPE {
    tag "${meta.id}"
    label 'process_high'

    conda "conda-forge::r-base=4.3.1 conda-forge::r-data.table=1.14.8 conda-forge::r-stringr=1.5.0 conda-forge::r-dplyr=1.1.3 conda-forge::r-ggplot2=3.4.4"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
        : 'community.wave.seqera.io/library/r-data.table_r-tidyverse:02853be7267d8cba'}"

    input:
    tuple val(meta), path(paf)
    path flu_db_info
    val score_threshold

    output:
    tuple val(meta), path("*_read_summary.tsv"), emit: summary
    tuple val(meta), path("*_read_serotype_assignment.pdf"), emit: plot
    tuple val(meta), path("*_*.txt"), emit: serotype_reads, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_pafs_influenza_A.R \\
        ${flu_db_info} \\
        ${paf} \\
        ${prefix} \\
        . \\
        ${score_threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ (.*//')
    END_VERSIONS
    """
}
