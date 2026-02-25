process SPLIT_FASTA_BY_SEGMENT {
    tag "$virus - $segment_name"
    label 'process_low'

    conda "bioconda::samtools=1.19.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_1' :
        'biocontainers/samtools:1.19.2--h50ea8bc_1' }"

    input:
    tuple val(virus), path(fasta), val(accession), val(segment_name)

    output:
    tuple val(virus), val(segment_name), path("*.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = "${virus}.${segment_name}"
    """
    # Index FASTA
    samtools faidx ${fasta}
    
    # Extract specific segment
    samtools faidx ${fasta} ${accession} > ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}