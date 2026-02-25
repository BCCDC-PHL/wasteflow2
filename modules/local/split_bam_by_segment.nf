process SPLIT_BAM_BY_SEGMENT {
    tag "$meta.id - $segment_name"
    label 'process_low'

    conda "bioconda::samtools=1.19.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_1' :
        'biocontainers/samtools:1.19.2--h50ea8bc_1' }"

    input:
    tuple val(meta), path(bam), path(bai), val(accession), val(segment_name)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bam.bai"), emit: bai
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}.${meta.genome}.${segment_name}"
    """
    # Extract reads mapping to specific segment using accession ID
    # BAI file will be auto-detected if it follows naming convention
    samtools view -h -b ${bam} ${accession} > ${prefix}.bam
    
    # Index the output BAM
    samtools index ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}