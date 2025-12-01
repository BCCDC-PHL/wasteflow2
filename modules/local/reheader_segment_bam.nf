process REHEADER_SEGMENT_BAM {
    tag "$meta.id - $meta.segment"
    label 'process_medium'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"
    
    input:
    tuple val(meta), path(bam), path(bai), path(segment_fasta), path(segment_fai)
    
    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("${prefix}.bam.bai"), emit: bai
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.genome}.${meta.segment}"
    """
    # The split BAM has reads mapped to a specific segment, but the header
    # contains all 8 segments and reads reference the original sequence ID (e.g., ID 4 for NP)
    # We need to create a new BAM where reads map to the single-segment reference (ID 0)
    
    # Get the sequence name from the segment accession
    SEG_NAME=\$(awk '{print \$1}' ${segment_fai})
    
    # Convert BAM to SAM, extract reads, and rebuild with new reference
    # This remaps all sequence IDs to the new reference
    samtools view -h ${bam} | \\
        awk -v seg="\${SEG_NAME}" '
            BEGIN {OFS="\\t"}
            /^@/ {
                # Skip old @SQ headers
                if (\$1 == "@SQ") next
                # Keep other headers
                print
                next
            }
            # For alignment records, set reference to 0 (single segment reference)
            {
                \$3 = seg  # Reference name
                print
            }
        ' | \\
        samtools view -b -t ${segment_fai} - > ${prefix}.unsorted.bam
    
    # Sort and index
    samtools sort -o ${prefix}.bam ${prefix}.unsorted.bam
    samtools index ${prefix}.bam
    
    # Clean up
    rm ${prefix}.unsorted.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    
    stub:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.genome}.${meta.segment}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}