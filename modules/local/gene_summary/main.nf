process GENE_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6d/6d63210b90bdadc321e15610f40c337ab08fa724719b7d4be0785944a86755fb/data':
        'community.wave.seqera.io/library/numpy_pandas_pyyaml_scanpy:d959777f7735763f' }"

    input:
    tuple val(meta),
        path(barcodes),
        path(vireo),
        path(demuxlet),
        path(freemuxlet),
        path(souporcell)

    output:
    tuple val(meta), path("*_genetic_summary_assignment.csv")     , emit: assignment
    tuple val(meta), path("*_genetic_summary_classification.csv") , emit: classification
    tuple val(meta), path("*_genetic_overview_assignment.csv")    , emit: overview_assignment
    tuple val(meta), path("*_genetic_overview_classification.csv"), emit: overview_classification
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    template 'gene_summary.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_genetic_summary_assignment.csv
    touch ${prefix}_genetic_summary_classification.csv
    touch ${prefix}_genetic_overview_assignment.csv
    touch ${prefix}_genetic_overview_classification.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        pandas: \$(python3 -c 'import pandas as pd; print(pd.__version__)')
        scanpy: \$(python3 -c 'import scanpy as sc; print(sc.__version__)')
        numpy: \$(python3 -c 'import numpy as np; print(np.__version__)')
    END_VERSIONS
    """
}
