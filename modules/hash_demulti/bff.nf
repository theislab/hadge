#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bff{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti/bff", mode:'copy'
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
        
        
    output:
        path "bff_${task.index}"
        
    script:
        def generateGenderPlot = generate_gender_plot != 'None' ? " --generateGenderPlot ${generate_gender_plot}" : ''
        """
        mkdir bff_${task.index}
        bff.R --seuratObject $seurat_object --assay $assay --methods $methods --methodsForConsensus $methodsForConsensus \
        --cellbarcodeWhitelist $cellbarcodeWhitelist --cellbarcodeWhitelist $cellbarcodeWhitelist --metricsFile $metricsFile \
        --doTSNE $doTSNE --doHeatmap $doHeatmap --perCellSaturation $perCellSaturation --majorityConsensusThreshold $majorityConsensusThreshold \
        --chemistry $chemistry --callerDisagreementThreshold $callerDisagreementThreshold --outputdir bff_${task.index}
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


        bff(seurat_object,assay, methods, methodsForConsensus, metricsFile,cellbarcodeWhitelist,doTSNE,doHeatmap,perCellSaturation,majorityConsensusThreshold,chemistry,callerDisagreementThreshold)
  
  emit:
        bff.out.collect()
}


workflow{
    bff_hashing()

}
