process SAMTOOLS_REHEADER_VIA_SAM {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::samtools=1.21"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0'
        : 'biocontainers/samtools:1.21--h50ea8bc_0'}"

    input:
    tuple val(meta), path(bam), path(header_file)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Convert BAM to SAM
    samtools view -h ${bam} > temp.sam
    
    # Extract only the alignment records (no header)
    samtools view ${bam} > temp_alignments.sam
    
    # Combine new header with alignment records
    cat ${header_file} temp_alignments.sam > temp_reheadered.sam
    
    # Convert back to BAM
    samtools view -bS temp_reheadered.sam > ${prefix}.bam
    
    # Clean up temporary files
    rm temp.sam temp_alignments.sam temp_reheadered.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
