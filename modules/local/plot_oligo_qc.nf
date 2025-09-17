process PLOT_OLIGO_QC {
    label 'process_single'

    conda "conda-forge::matplotlib=3.10.5 conda-forge::pandas=2.3.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5d/5d4923bc2b9c2d8c83468eddf1e32c0bfecffab8192e0e2ee69340fe267560cd/data'
        : 'community.wave.seqera.io/library/matplotlib_pandas:1c75679e093ee81f'}"

    input:
    path concatenated_counts

    output:
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Run visualization script
    plot_oligo_qc.py \
        --input_file ${concatenated_counts} \
        --control_ref CTRL_OLIGO \
        --out_prefix ${params.run_id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch viral_control_other_reads.png
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
