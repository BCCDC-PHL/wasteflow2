process SNPEFF_BUILD {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::snpeff=5.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/snpeff:5.0--hdfd78af_1'
        : 'quay.io/biocontainers/snpeff:5.0--hdfd78af_1'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gff)

    output:
    tuple val(meta), path("snpeff_db_*"), emit: db
    tuple val(meta), path("snpeff_*.config"), emit: config
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    // Generate unique genome key from metadata
    def genome_key = meta.containsKey('segment') && meta.segment && meta.segment != 'none' ?  "${meta.genome}.${meta.segment}" : "${meta.genome}"
    
    // Detect GFF format
    def extension = gff.getExtension()
    def format = (extension == "gtf") ? "gtf22" : "gff3"

    // Calculate available memory
    def avail_mem = 4
    if (!task.memory) {
        log.info("[snpEff] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.")
    } else {
        avail_mem = task.memory.giga
    }
    
    """
    # Create directory structure with descriptive name
    mkdir -p snpeff_db_${genome_key}/genomes/
    cd snpeff_db_${genome_key}/genomes/
    ln -s ../../${fasta} ${genome_key}.fa
    
    cd ../../
    mkdir -p snpeff_db_${genome_key}/${genome_key}/
    cd snpeff_db_${genome_key}/${genome_key}/
    ln -s ../../${gff} genes.${extension}
    
    cd ../../
    
    # Create config with descriptive name
    cat > snpeff_${genome_key}.config << EOF
    data.dir = ./snpeff_db_${genome_key}
    ${genome_key}.genome : ${genome_key}
    EOF

    # Build database
    snpEff \\
        -Xmx${avail_mem}g \\
        build \\
        -config snpeff_${genome_key}.config \\
        -dataDir ./snpeff_db_${genome_key} \\
        -${format} \\
        ${args} \\
        -v \\
        ${genome_key}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}