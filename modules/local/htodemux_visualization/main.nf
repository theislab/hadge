process HTODEMUX_VISUALIZATION {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-ggplot2_r-seurat:dac8c905972b98df':
        'community.wave.seqera.io/library/r-ggplot2_r-seurat:eefd54806320eae0' }"

    input:
    tuple val(meta), path(seurat_object), val(assay)

    output:
    tuple val(meta), path("*_ridge_htodemux.jpeg")         , emit: ridge_plot     , optional: true
    tuple val(meta), path("*_featureScatter_htodemux.jpeg"), emit: feature_scatter, optional: true
    tuple val(meta), path("*_violinPlot_htodemux.jpeg")    , emit: violin_plot    , optional: true
    tuple val(meta), path("*_tSNE_htodemux.jpeg")          , emit: tsne_plot      , optional: true
    tuple val(meta), path("*_heatMap_htodemux.jpeg")       , emit: heatmap_plot   , optional: true
    tuple val(meta), path("*_visual_params_htodemux.csv")  , emit: params
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Ridge Plot Parameters
    ridgePlot      = task.ext.ridgePlot      ?: true         // Generate ridge plot
    ridgeNCol      = task.ext.ridgeNCol      ?: 2            // Number of columns for ridge plot

    // Feature Scatter Plot Parameters
    featureScatter = task.ext.featureScatter ?: true         // Generate feature scatter plot
    scatterFeat1   = task.ext.scatterFeat1   ?: null         // Feature 1 for scatter plot
    scatterFeat2   = task.ext.scatterFeat2   ?: null         // Feature 2 for scatter plot

    // Violin Plot Parameters
    vlnPlot        = task.ext.vlnPlot        ?: true         // Generate violin plot
    vlnFeatures    = task.ext.vlnFeatures    ?: "nCount_RNA" // Features to plot (gene expression, metrics, PC scores, anything that can be retreived by FetchData)
    vlnLog         = task.ext.vlnLog         ?: true         // Plot the feature axis on log scale

    // TSNE Plot Parameters
    tSNE           = task.ext.tSNE           ?: true         // Generate a two dimensional tSNE embedding for HTOs
    tSNEIdents     = task.ext.tSNEIdents     ?: "Negative"   // What should we remove from the object (we have Singlet,Doublet and Negative)
    tSNEInvert     = task.ext.tSNEInvert     ?: true         // Invert tSNE selection
    tSNEVerbose    = task.ext.tSNEVerbose    ?: false        // Verbose tSNE
    tSNEApprox     = task.ext.tSNEApprox     ?: false        // Approximate tSNE
    tSNEDimMax     = task.ext.tSNEDimMax     ?: 2            // Max number of donors
    tSNEPerplexity = task.ext.tSNEPerplexity ?: 100          // Value for perplexity

    // Heatmap Parameters
    heatMap        = task.ext.heatMap        ?: true         // Generate heatmap
    heatMapNcells  = task.ext.heatMapNcells  ?: 500          // Number of cells for heatmap

    // Output Parameters
    prefix         = task.ext.prefix         ?: "${meta.id}"

    template 'htodemux_visualization.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ridge_htodemux.jpeg
    touch ${prefix}_featureScatter_htodemux.jpeg
    touch ${prefix}_violinPlot_htodemux.jpeg
    touch ${prefix}_tSNE_htodemux.jpeg
    touch ${prefix}_heatMap_htodemux.jpeg
    touch ${prefix}_visual_params_htodemux.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
    END_VERSIONS
    """
}
