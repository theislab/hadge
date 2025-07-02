process RENAME_GENES_TO_FEATURES {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(dir)

    output:
    tuple val(meta), path(dir)

    script:
    """
    if [ -f "${dir}/genes.tsv.gz" ]; then
        mv "${dir}/genes.tsv.gz" "${dir}/features.tsv.gz"
    fi
    """
}
