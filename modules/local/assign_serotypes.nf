process ASSIGN_SEROTYPES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pysam=0.22.0 conda-forge::polars=0.20.0 conda-forge::matplotlib=3.8.0 conda-forge::pyarrow"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b883fda47e52cc88fe8ffdc52fdf9d246ca4bf048d0ac7f043c89ce33a0d6cf9' :
        'community.wave.seqera.io/library/pysam_matplotlib_polars_pyarrow:b1a40fbc0eea3667' }"

    input:
    tuple val(meta), path(bam)
    path db_info

    output:
    tuple val(meta), path("${prefix}_per_read_summary.tsv")       , optional: true, emit: read_summary
    tuple val(meta), path("${prefix}_per_serotype_summary.tsv")   , optional: true, emit: serotype_summary
    tuple val(meta), path("${prefix}_read_serotype_assignment.pdf"), optional: true, emit: assignment_plot
    tuple val(meta), path("${prefix}_*.txt")                      , optional: true, emit: read_lists
    path "versions.yml"                                           , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def score_threshold = task.ext.score_thresh ?: 90
    """
    assign_serotypes.py \\
        --bam ${bam} \\
        --db-info ${db_info} \\
        --sample ${prefix} \\
        --outdir . \\
        --score-thresh ${score_threshold} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
        polars: \$(python -c "import polars; print(polars.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_per_read_summary.tsv
    touch ${prefix}_per_serotype_summary.tsv
    touch ${prefix}_read_serotype_assignment.pdf
    touch ${prefix}_ambiguous.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pysam: stub
        polars: stub
        matplotlib: stub
    END_VERSIONS
    """
}