process DROPLETUTILS_MTXCONVERT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c9f81df3cdd03c86a8133f74c0deb78719798c061895e4d9dd454f05e82ff93e/data'
        : 'community.wave.seqera.io/library/bioconductor-dropletutils:1.26.0--35a578ac06f1c531'}"

    input:
    tuple val(meta), path(input_mtx_dir)
    val write_csv

    output:
    tuple val(meta), path("*.csv"), emit: csv, optional: true
    tuple val(meta), path("*.h5"), emit: h5
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template("convert.R")

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ ${write_csv} == true ]]; then
        touch ${prefix}.csv
    fi
    touch ${prefix}.h5
    """
}
