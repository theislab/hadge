process HASH_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb8601e2171467026ea36c22328a15eb25025bbe686ff1a0ea04ab407c735aee/data':
        'community.wave.seqera.io/library/pegasusio_numpy_pandas_pyyaml_scanpy:e16c3756496aa20c' }"

    input:
    tuple val(meta),
        path(hto_matrix),
        path(htodemux_assignments), path (htodemux_classification),
        path(multiseq),
        path(bff),
        path(demuxem),
        path(gmmdemux_results), path(gmmdemux_config),
        path(hasheddrops_results), path(hasheddrops_id_to_hash),
        path(hashsolo)
    val(bff_methods)

    output:
    tuple val(meta), path("*_hashing_summary_assignment.csv")     , emit: assignment
    tuple val(meta), path("*_hashing_summary_classification.csv") , emit: classification
    tuple val(meta), path("*_hashing_overview_assignment.csv")    , emit: overview_assignment
    tuple val(meta), path("*_hashing_overview_classification.csv"), emit: overview_classification
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix         = task.ext.prefix         ?: "${meta.id}"
    hash_list      = "${meta.hto_names}".split(",")

    template 'hash_summary.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_hashing_summary_assignment.csv
    touch ${prefix}_hashing_summary_classification.csv
    touch ${prefix}_hashing_overview_assignment.csv
    touch ${prefix}_hashing_overview_classification.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        pandas: \$(python3 -c 'import pandas as pd; print(pd.__version__)')
        scanpy: \$(python3 -c 'import scanpy as sc; print(sc.__version__)')
        numpy: \$(python3 -c 'import numpy as np; print(np.__version__)')
        pegasusio: \$(python3 -c 'import pegasusio as io; print(io.__version__)')
    END_VERSIONS
    """
}
