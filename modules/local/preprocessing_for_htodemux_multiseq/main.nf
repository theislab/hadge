process PREPROCESSING_FOR_HTODEMUX_MULTISEQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c9f81df3cdd03c86a8133f74c0deb78719798c061895e4d9dd454f05e82ff93e/data'
        : 'community.wave.seqera.io/library/r-seurat:5.3.0--eeb977835038859a'}"

    input:
    tuple val(meta), path(hto_matrix), path(rna_matrix)

    output:
    tuple val(meta), path("*_preprocessed.rds")        , emit: seurat_object
    tuple val(meta), path("*_params_preprocessing.csv"), emit: params
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // preprocessing parameters
    sel_method       = task.ext.sel_method       ?: "mean.var.plot"
    ndelim           = task.ext.ndelim           ?: "_"
    n_features       = task.ext.n_features       ?: "2000"
    assay            = task.ext.assay            ?: "HTO"
    margin           = task.ext.margin           ?: "2"
    norm_method      = task.ext.norm_method      ?: "CLR"
    preprocessOut    = task.ext.preprocessOut    ?: "preprocessed"
    gene_col         = task.ext.gene_col         ?: "2"

    // others
    prefix           = task.ext.prefix           ?: "${meta.id}"

    template 'pre_processing.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_preprocessed.rds
    touch ${prefix}_params_preprocessing.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
        r-argparse: \$(Rscript -e "library(argparse); cat(as.character(packageVersion('argparse')))")
    END_VERSIONS
    """
} 