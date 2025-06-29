process PREPROCESSING_FOR_HTODEMUX_MULTISEQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-seurat:5.3.0--00f94834f5eea080':
        'community.wave.seqera.io/library/r-seurat:5.3.0--eeb977835038859a' }"


    input:
    tuple val(meta), path(rna_matrix), path(hto_matrix)

    output:
    tuple val(meta), path("*_preprocessed.rds")        , emit: seurat_object
    tuple val(meta), path("*_params_preprocessing.csv"), emit: params
    path "versions.yml"                                , emit: versions

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
    def prefix = task.ext.prefix ?: "${meta.id}"
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
