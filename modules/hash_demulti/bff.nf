#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bff{
    publishDir "$projectDir/$params.outdir/${seurat_object.name.tokenize( '_' )[1]}/$params.mode/hash_demulti/bff", mode:'copy'
    label 'small_mem'
    input:
       
        tuple val(sampleId), path(hto_matrix, stageAs: 'hto_data'),
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
        each assignmentOutBff
        each preprocess_bff
        each barcodeWhitelist
        
    output:
        path "bff_${sampleId}"
        
        
    script:
        def sampleId = seurat_object.name.tokenize( '_' )[1]

        """
        mkdir bff_${sampleId}
        bff.R --fileHto hto_data --methods $methods --methodsForConsensus $methodsForConsensus \
        --cellbarcodeWhitelist $cellbarcodeWhitelist --metricsFile bff_${sampleId}_$metricsFile \
        --doTSNE $doTSNE --doHeatmap $doHeatmap --perCellSaturation $perCellSaturation --majorityConsensusThreshold $majorityConsensusThreshold \
        --chemistry $chemistry --callerDisagreementThreshold $callerDisagreementThreshold --outputdir bff_${sampleId} \
         --assignmentOutBff $assignmentOutBff --preprocess $preprocess_bff --barcodeWhitelist $barcodeWhitelist
        """

}


workflow bff_hashing{
  take: 
        hto_matrix
  main:
        methods = params.methods
        methodsForConsensus = params.methodsForConsensus
        cellbarcodeWhitelist = params.cellbarcodeWhitelist
        metricsFile = params.metricsFile
        doTSNE = params.doTSNE
        doHeatmap = params.doHeatmap
        perCellSaturation = params.perCellSaturation
        majorityConsensusThreshold  = params.majorityConsensusThreshold
        chemistry = params.chemistry
        callerDisagreementThreshold = params.callerDisagreementThreshold
        assignmentOutBff = params.assignmentOutBff
        preprocess_bff = params.preprocess_bff
        barcodeWhitelist = params.barcodeWhitelist

        bff(hto_matrix, methods, methodsForConsensus, cellbarcodeWhitelist, metricsFile,doTSNE,doHeatmap,perCellSaturation,majorityConsensusThreshold,chemistry,callerDisagreementThreshold,assignmentOutBff,preprocess_bff,barcodeWhitelist)
  
  emit:
        bff.out.collect()
}


workflow{
    bff_hashing()

}
