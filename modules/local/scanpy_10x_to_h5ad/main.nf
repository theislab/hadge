process SCANPY_10X_TO_H5AD {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/168ecbbe27ccef766741ccbf937b0d2675be2e19b0565035e0719f1e9ea5ee95/data':
        'community.wave.seqera.io/library/scanpy:1.11.4--05fed4ff91d741c4' }"

    input:
    tuple val(meta), path(input_mtx_dir)

    output:
    tuple val(meta), path("*_hto_data.h5ad"), emit: h5ad
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'convert_to_h5ad.py'

    // stub:
    // prefix = task.ext.prefix ?: "${meta.id}"
    // """
    // touch ${prefix}_hto_data.h5ad
    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     python: $(python3 --version | cut -f 2 -d " ")
    //     scanpy: $(python3 -c "import scanpy; print(scanpy.__version__)")
    // END_VERSIONS
    // """
}
