process EXTRACT_HASHES {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(hto_matrix)

    output:
    tuple val(meta), path("*_hashes.txt"), emit: hashes

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix         = task.ext.prefix         ?: "${meta.id}"

    script:
    """
    zcat $hto_matrix | awk '{print \$2}' | paste -sd, > ${prefix}_hashes.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_hashes.txt
    """
}
