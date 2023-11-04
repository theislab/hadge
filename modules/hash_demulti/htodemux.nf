#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process htodemux{
    publishDir "$projectDir/$params.outdir/${seurat_object.name.tokenize( '_' )[1]}/$params.mode/hash_demulti/htodemux", mode: 'copy'
    label 'small_mem'

    conda "conda-forge::r-base=4.1 conda-forge::r-seurat conda-forge::r-argparse"
    
    input:
        each seurat_object
        val assay
        val quantile
        val kfunc
        val nstarts
        val nsamples
        val seed
        val init
        val objectOutHTO
        val assignmentOutHTO
        
        //Ridge plot params
        val ridgePlot
        val ridgeNCol
        //Scatter features params
        val featureScatter
        val scatterFeat1
        val scatterFeat2
        //Violin plot params
        val vlnplot
        val vlnFeatures
        val vlnLog
        //tSNE
        val tsne
        val tsneIdents
        val tsneInvert
        val tsneVerbose
        val tsneApprox
        val tsneDimMax
        val tsnePerplexity
        //Heatmap
        val heatmap
        val heatmapNcells
        
    output:
        path "htodemux_${seurat_object.name.tokenize( '_' )[1]}"
        
    script:
        def sampleId = seurat_object.name.tokenize( '_' )[1]
        def init_val = init != 'None' ? " --init $init" : ''
        def vln_log = vlnLog != 'False' ?  "--vlnLog" : ''
        def invert = tsneInvert != 'False' ?  "--tSNEInvert" : ''
        def verbose = tsneVerbose != 'False' ?  "--tSNEVerbose" : ''
        def approx = tsneApprox != 'False' ?  "--tSNEApprox" : ''
        
        """
        mkdir htodemux_${sampleId}
        HTODemux.R --seuratObject $seurat_object --assay $assay --quantile $quantile --kfunc $kfunc \
                   --nstarts $nstarts --nsamples $nsamples --seed $seed $init_val --objectOutHTO $objectOutHTO \
                   --assignmentOutHTO $assignmentOutHTO --outputdir htodemux_${sampleId}
        HTODemux-visualisation.R --hashtagPath htodemux_${sampleId}/${objectOutHTO}.rds --assay $assay \
                   --ridgePlot $ridgePlot --ridgeNCol $ridgeNCol --featureScatter $featureScatter \
                   --scatterFeat1 $scatterFeat1 --scatterFeat2 $scatterFeat2 --vlnPlot $vlnplot \
                   --vlnFeatures $vlnFeatures $vln_log --tSNE $tsne --tSNEIdents $tsneIdents $invert $verbose $approx \
                   --tSNEDimMax $tsneDimMax --tSNEPerplexity $tsnePerplexity --heatMap $heatmap --heatMapNcells $heatmapNcells \
                   --outputdir htodemux_${sampleId}
        """

}

workflow htodemux_hashing{
    take:
        seurat_object
    main:
        quantile = params.quantile_htodemux
        assay = params.assay
        kfunc = params.kfunc
        nstarts = params.nstarts
        nsamples = params.nsamples
        seed = params.seed
        init = params.init
        objectOutHTO = params.objectOutHTO
        assignmentOutHTO = params.assignmentOutHTO
        
        ridgePlot = params.ridgePlot
        ridgeNCol = params.ridgeNCol
        featureScatter = params.featureScatter
        scatterFeat1 = params.scatterFeat1
        scatterFeat2 = params.scatterFeat2
        vlnplot = params.vlnplot
        vlnFeatures = params.vlnFeatures
        vlnLog = params.vlnLog
        
        tsne = params.tsne
        tsneIdents = params.tsneIdents
        tsneInvert = params.tsneInvert
        tsneVerbose = params.tsneVerbose
        tsneApprox = params.tsneApprox
        tsneDimMax = params.tsneDimMax
        tsnePerplexity = params.tsnePerplexity
        heatmap = params.heatmap
        heatmapNcells = params.heatmapNcells

        htodemux(seurat_object, assay, quantile, kfunc, nstarts, nsamples, seed, init, objectOutHTO, assignmentOutHTO, ridgePlot, ridgeNCol, featureScatter, scatterFeat1, scatterFeat2, vlnplot, vlnFeatures, vlnLog, tsne, tsneIdents, tsneInvert, tsneVerbose, tsneApprox, tsneDimMax, tsnePerplexity, heatmap, heatmapNcells)

    emit:
        htodemux.out.collect()
}
