process SNPEFF_ANN {
    tag "${meta.id}.${meta.genome}.${meta.segment ?: 'nonsegmented'}"
    label 'process_medium'
    conda "bioconda::snpeff=5.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/snpeff:5.0--hdfd78af_1'
        : 'quay.io/biocontainers/snpeff:5.0--hdfd78af_1'}"

    input:
    tuple val(meta), 
          path(vcf), 
          path(snpeff_db),      
          path(snpeff_config),  
          path(fasta)

    output:
    tuple val(meta), path("*.snpeff.vcf")  , emit: vcf
    tuple val(meta), path("*.snpeff.csv")  , emit: csv
    tuple val(meta), path("*.snpeff.genes.txt")  , emit: txt
    tuple val(meta), path("*.snpeff.html") , emit: html
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Generate genome identifier - this matches the database name
    def genome_id = meta.containsKey('segment') ? "${meta.genome}.${meta.segment}" : "${meta.genome}"
    
    // Calculate available memory (80% of allocated memory for JVM heap)
    def avail_mem = 4
    if (!task.memory) {
        log.info('[SnpEff] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.')
    } else {
        avail_mem = (0.8 * task.memory.toGiga()).intValue()
    }
    
    """
    # Debug: Show what files are staged
    echo "=== Staged files ==="
    ls -la
    echo "=== Database directory ==="
    ls -la snpeff_db_${genome_id}/
    echo "=== Config file ==="
    cat snpeff_${genome_id}.config
    echo "===================="
    
    # Run SnpEff annotation
    snpEff ann \\
        -Xmx${avail_mem}g \\
        -config snpeff_${genome_id}.config \\
        -dataDir snpeff_db_${genome_id} \\
        -csvStats ${prefix}.snpeff.csv \\
        -stats ${prefix}.snpeff.txt \\
        -s ${prefix}.snpeff.html \\
        ${args} \\
        ${genome_id} \\
        ${vcf} \\
        > ${prefix}.snpeff.vcf

    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | sed 's/^.*SnpEff //; s/ .*\$//')
    END_VERSIONS
    """
}