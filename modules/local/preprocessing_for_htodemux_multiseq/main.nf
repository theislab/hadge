process PREPROCESSING_FOR_HTODEMUX_MULTISEQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b43d3afc47ad5e5d99bc97980d409e377f6b0595cb588b121076e67dca71d39/data':
        'community.wave.seqera.io/library/r-seurat_r-seuratobject:c1b3e7a7276bda09' }"


    input:
    tuple val(meta), path(rna_matrix), path(hto_matrix)

    output:
    tuple val(meta), path("*_preprocessed.rds")        , emit: seurat_object
    tuple val(meta), path("*_params_preprocessing.csv"), emit: params
    path "versions.yml"                                , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // preprocessing parameters
    sel_method  = task.ext.sel_method  ?: "mean.var.plot" // Selection method
    ndelim      = task.ext.ndelim      ?: "_"             // For the initial identity class for each cell, delimiter for the cell's column name
    n_features  = task.ext.n_features  ?: 2000            // Number of features to be used when finding variable features
    assay       = task.ext.assay       ?: "HTO"           // Assay name for hashing modality
    margin      = task.ext.margin      ?: 2               // Margin for normalisation
    norm_method = task.ext.norm_method ?: "CLR"           // Normalisation method
    gene_col    = task.ext.gene_col    ?: 2               // Specify which column of genes.tsv or features.tsv to use for gene names

    // others
    prefix      = task.ext.prefix      ?: "${meta.id}"

    template 'pre_processing.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_preprocessed.rds
    touch ${prefix}_params_preprocessing.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
