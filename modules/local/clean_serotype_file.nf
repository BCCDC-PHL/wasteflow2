// Process to clean serotype files using your sed solution
process CLEAN_SEROTYPE_FILE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(serotype_file)

    output:
    tuple val(meta), path("${meta.id}_${meta.serotype}_clean.txt"), emit: cleaned

    script:
    """
    # Clean serotype file using awk - remove everything from first / or space onwards
    awk '{gsub(/[\\/\\ ].*/, ""); print}' ${serotype_file} > ${meta.id}_${meta.serotype}_clean.txt
    """
}
