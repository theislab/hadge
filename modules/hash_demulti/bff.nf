#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bff{
    publishDir "$projectDir/$params.outdir/${seurat_object.name.tokenize( '_' )[1]}/$params.mode/hash_demulti/bff", mode:'copy'
    input:
        each seurat_object
        //Seurat process
        each assay
        each methods
        each methodsForConsensus
        each cellbarcodeWhitelist
        each metricsFile
        each doTSNE
        each doHeatmap
        each perCellSaturation
        each majorityConsensusThreshold
        each chemistry
        each callerDisagreementThreshold
        each rawFeatureMatrixH5
        each assignmentOutBff
        
    output:
        path "bff_${sampleId}"
        
        
    script:
        def h5_available = rawFeatureMatrixH5 != 'None' ? " --rawFeatureMatrixH5 ${rawFeatureMatrixH5}" : ''
        def sampleId = seurat_object.name.tokenize( '_' )[1]

        """
        mkdir bff_${sampleId}
        Rscript $baseDir/bin/bff.R --seuratObject $seurat_object --assay $assay --methods $methods --methodsForConsensus $methodsForConsensus \
        --cellbarcodeWhitelist $cellbarcodeWhitelist --cellbarcodeWhitelist $cellbarcodeWhitelist --metricsFile bff_${task.index}_$metricsFile \
        --doTSNE $doTSNE --doHeatmap $doHeatmap --perCellSaturation $perCellSaturation --majorityConsensusThreshold $majorityConsensusThreshold \
        --chemistry $chemistry --callerDisagreementThreshold $callerDisagreementThreshold $h5_available --outputdir bff_${sampleId} --assignmentOutBff $assignmentOutBff
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

workflow bff_hashing{
  take: 
        seurat_object
  main:
        assay = split_input(params.assay)
        methods = split_input(params.methods)
        methodsForConsensus = split_input(params.methodsForConsensus)
        cellbarcodeWhitelist = split_input(params.cellbarcodeWhitelist)
        metricsFile = split_input(params.metricsFile)
        doTSNE = split_input(params.doTSNE)
        doHeatmap = split_input(params.doHeatmap)
        perCellSaturation = split_input(params.perCellSaturation)
        majorityConsensusThreshold  = split_input(params.majorityConsensusThreshold)
        chemistry = split_input(params.chemistry)
        callerDisagreementThreshold = split_input(params.callerDisagreementThreshold)
        rawFeatureMatrixH5 = split_input(params.rawFeatureMatrixH5)
        assignmentOutBff = split_input(params.assignmentOutBff)
        

        bff(seurat_object,assay, methods, methodsForConsensus, metricsFile,cellbarcodeWhitelist,doTSNE,doHeatmap,perCellSaturation,majorityConsensusThreshold,chemistry,callerDisagreementThreshold,rawFeatureMatrixH5,assignmentOutBff)
  
  emit:
        bff.out.collect()
}


workflow{
    bff_hashing()

}
