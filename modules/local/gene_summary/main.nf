process GENE_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d863e56b5ce15b271e8c8666ec22217df5cfc57a9731cc23c7f92674dc7ab0c7/data':
        'community.wave.seqera.io/library/pegasusio_mudata_numpy_pandas_pruned:ecdbf7e42b2f3213' }"

    input:
    tuple val(meta),
        path(rna_matrix),
        path(hto_matrix),
        path(barcodes),
        path(vireo),
        path(demuxlet),
        path(freemuxlet)
    tuple val (generate_anndata), val(generate_mudata)


    output:
    tuple val(meta), path("*_genetic_summary_assignment.csv")    , emit: assignment    , optional: false
    tuple val(meta), path("*_genetic_summary_classification.csv"), emit: classification, optional: false
    tuple val(meta), path("*_genetic_summary.h5ad")              , emit: h5ad          , optional: true
    tuple val(meta), path("*_genetic_summary.h5mu")              , emit: h5mu          , optional: true
    path "versions.yml"                                          , emit: versions      , optional: false

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        pandas: \$(python3 -c 'import pandas as pd; print(pd.__version__)')
        scanpy: \$(python3 -c 'import scanpy as sc; print(sc.__version__)')
        numpy: \$(python3 -c 'import numpy as np; print(np.__version__)')
        mudata: \$(python3 -c 'import mudata as md; print(md.__version__)')
        pegasusio: \$(python3 -c 'import pegasusio as io; print(io.__version__)')
        yaml: \$(python3 -c 'import yaml; print(yaml.__version__)')
    END_VERSIONS
    """
}
