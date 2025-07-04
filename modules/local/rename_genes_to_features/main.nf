process RENAME_GENES_TO_FEATURES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.3':
        'biocontainers/coreutils:9.3' }"

    input:
    tuple val(meta), path(dir)

    output:
    tuple val(meta), path(prefix)

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -Lr ${dir} ${prefix}

    if [ -f "${prefix}/genes.tsv.gz" ]; then
        mv "${prefix}/genes.tsv.gz" "${prefix}/features.tsv.gz"
    fi
    """
}
