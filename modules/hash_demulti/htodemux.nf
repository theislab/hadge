#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process htodemux{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti/htodemux", mode: 'copy'
    label 'seurat'
    input:
        each seurat_object
        each assay
        each quantile
        each kfunc
        each nstarts
        each nsamples
        each seed
        each init
        each objectOutHTO
        each assignmentOutHTO
        
        //Ridge plot params
        each ridgePlot
        each ridgeNCol
        //Scatter features params
        each featureScatter
        each scatterFeat1
        each scatterFeat2
        //Violin plot params
        each vlnplot
        each vlnFeatures
        each vlnLog
        //tSNE
        each tsne
        each tsneIdents
        each tsneInvert
        each tsneVerbose
        each tsneApprox
        each tsneDimMax
        each tsnePerplexity
        //Heatmap
        each heatmap
        each heatmapNcells
        
    output:
        path "htodemux_${task.index}"
        
    script:
        def init_val = init != 'NULL' ? " --init $init" : ''
        def vln_log = vlnLog != 'FALSE' ?  "--vlnLog" : ''
        def invert = tsneInvert != 'FALSE' ?  "--tSNEInvert" : ''
        def verbose = tsneVerbose != 'FALSE' ?  "--tSNEVerbose" : ''
        def approx = tsneApprox != 'FALSE' ?  "--tSNEApprox" : ''
        
        """
        mkdir htodemux_${task.index}
        HTODemux.R --seuratObject $seurat_object --assay $assay --quantile $quantile --kfunc $kfunc --nstarts $nstarts --nsamples $nsamples --seed $seed $init_val --objectOutHTO $objectOutHTO --assignmentOutHTO $assignmentOutHTO --outputdir htodemux_${task.index}
        HTODemux-visualisation.R --hashtagPath htodemux_${task.index}/${objectOutHTO}.rds --assay $assay --ridgePlot $ridgePlot --ridgeNCol $ridgeNCol --featureScatter $featureScatter --scatterFeat1 $scatterFeat1 --scatterFeat2 $scatterFeat2 --vlnPlot $vlnplot --vlnFeatures $vlnFeatures $vln_log --tSNE $tsne --tSNEIdents $tsneIdents $invert $verbose $approx --tSNEDimMax $tsneDimMax --tSNEPerplexity $tsnePerplexity --heatMap $heatmap --heatMapNcells $heatmapNcells --outputdir htodemux_${task.index}
        """

}

def split_input(input){
    if (input =~ /;/ ){
        Channel.from(input).map{ return it.tokenize(';')}.flatten()
    }
    else{
        Channel.from(input)
    }
}


workflow htodemux_hashing{
    take:
        seurat_object
    main:
        quantile = split_input(params.quantile_htodemux)
        assay = split_input(params.assay)
        kfunc = split_input(params.kfunc)
        nstarts = split_input(params.nstarts)
        nsamples = split_input(params.nsamples)
        seed = split_input(params.seed)
        init = split_input(params.init)
        objectOutHTO = split_input(params.objectOutHTO)
        assignmentOutHTO = split_input(params.assignmentOutHTO)
        
        ridgePlot = split_input(params.ridgePlot)
        ridgeNCol = split_input(params.ridgeNCol)
        featureScatter = split_input(params.featureScatter)
        scatterFeat1 = split_input(params.scatterFeat1)
        scatterFeat2 = split_input(params.scatterFeat2)
        vlnplot = split_input(params.vlnplot)
        vlnFeatures = split_input(params.vlnFeatures)
        vlnLog = split_input(params.vlnLog)
        
        tsne = split_input(params.tsne)
        tsneIdents = split_input(params.tsneIdents)
        tsneInvert = split_input(params.tsneInvert)
        tsneVerbose = split_input(params.tsneVerbose)
        tsneApprox = split_input(params.tsneApprox)
        tsneDimMax = split_input(params.tsneDimMax)
        tsnePerplexity = split_input(params.tsnePerplexity)
        heatmap = split_input(params.heatmap)
        heatmapNcells = split_input(params.heatmapNcells)

        htodemux(seurat_object, assay, quantile, kfunc, nstarts, nsamples, seed, init, objectOutHTO, assignmentOutHTO, ridgePlot, ridgeNCol, featureScatter, scatterFeat1, scatterFeat2, vlnplot, vlnFeatures, vlnLog, tsne, tsneIdents, tsneInvert, tsneVerbose, tsneApprox, tsneDimMax, tsnePerplexity, heatmap, heatmapNcells)

    emit:
        htodemux.out.collect()
}

workflow{
    htodemux_hashing()
}
