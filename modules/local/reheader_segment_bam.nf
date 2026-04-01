process REHEADER_SEGMENT_BAM {
    tag "$meta.id - $meta.segment"
    label 'process_medium'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"
    
    input:
    tuple val(meta), 
        path(bam), 
        path(bai), 
        path(segment_fasta),  // Remove stageAs
        path(segment_fai),    // Remove stageAs
        val(segment_accession)
    
    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("${prefix}.bam.bai"), emit: bai
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.genome}.${meta.segment}"

    // Use the explicit segment_accession parameter
    def seg_name = segment_accession

    """
    # ABSOLUTE DETERMINISTIC VERSION
    # All operations use explicit file names (no auto-generated temp files)

    # Step 1: Create header
    samtools view -H ${bam} | grep -v "^@SQ" > header.sam
    awk 'BEGIN {OFS="\\t"} {print "@SQ", "SN:" \$1, "LN:" \$2}' ${segment_fai} >> header.sam

    # Step 2: Extract and update alignments
    samtools view ${bam} | \\
        awk -v seg="${seg_name}" 'BEGIN {OFS="\\t"} {
            \$3 = seg
            print
        }' > alignments.sam

    # Step 3: Combine to unsorted BAM
    cat header.sam alignments.sam | \\
        samtools view -b -l 6 -t segment.fasta.fai -o unsorted.bam -

    # Step 4: Sort with MAXIMUM determinism
    # Key: Use current directory as temp location to avoid path issues
    samtools sort \\
        --threads 1 \\
        -l 6 \\
        -m 768M \\
        -o ${prefix}.bam \\
        unsorted.bam

    # Step 5: Index
    samtools index ${prefix}.bam

    # Clean up
    rm -f header.sam alignments.sam unsorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}