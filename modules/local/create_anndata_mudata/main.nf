process CREATE_ANNDATA_MUDATA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/31c261a4a1ed9c3b409457fe778a363fb941152f7307bfa76cb4c42d44235ddf/data':
        'community.wave.seqera.io/library/anndata_mudata_pandas_pyyaml_scanpy:e96a91e210372525' }"

    input:
    tuple val(meta),
        path(rna_matrix),
        path(hto_matrix),
        path(genetic_summary_assignment),
        path(genetic_summary_classification),
        path(hashing_summary_assignment),
        path(hashing_summary_classification)

    output:
    tuple val(meta), path("*_genetic.h5ad")            , emit: h5ad_genetic, optional: true
    tuple val(meta), path("*_hashing.h5ad")            , emit: h5ad_hashing, optional: true
    tuple val(meta), path("*_genetic_and_hashing.h5mu"), emit: h5mu        , optional: true
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    template 'create_anndata_mudata.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_genetic.h5ad
    touch ${prefix}_hashing.h5ad
    touch ${prefix}_genetic_and_hashing.h5mu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        pandas: \$(python3 -c 'import pandas as pd; print(pd.__version__)')
        scanpy: \$(python3 -c 'import scanpy as sc; print(sc.__version__)')
        mudata: \$(python3 -c 'import mudata as md; print(md.__version__)')
        anndata: \$(python3 -c 'import anndata as ad; print(ad.__version__)')
    END_VERSIONS
    """
}
