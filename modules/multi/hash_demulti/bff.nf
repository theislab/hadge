#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process bff{
    publishDir "$params.outdir/$sampleId/$params.mode/hash_demulti/bff", mode:'copy'
    label 'small_mem'

    conda "$projectDir/conda/bff.yml"

    input:

        tuple val(sampleId), path(hto_matrix, stageAs: 'hto_data')
        val methods
        val methodsForConsensus
        val cellbarcodeWhitelist
        val metricsFile
        val doTSNE
        val doHeatmap
        val perCellSaturation
        val majorityConsensusThreshold
        val chemistry
        val callerDisagreementThreshold
        val assignmentOutBff
        val preprocess_bff
        val barcodeWhitelist

    output:
        path "bff_${sampleId}"

    script:

        def run_preprocess = preprocess_bff != 'False' ? " --preprocess_bff" : ''
        """
        mkdir bff_${sampleId}
        bff.R --fileHto hto_data --methods $methods --methodsForConsensus $methodsForConsensus \
        --cellbarcodeWhitelist $cellbarcodeWhitelist --metricsFile bff_${sampleId}_$metricsFile \
        --doTSNE $doTSNE --doHeatmap $doHeatmap --perCellSaturation $perCellSaturation --majorityConsensusThreshold $majorityConsensusThreshold \
        --chemistry $chemistry --callerDisagreementThreshold $callerDisagreementThreshold --outputdir bff_${sampleId} \
        --assignmentOutBff $assignmentOutBff ${run_preprocess} --barcodeWhitelist $barcodeWhitelist
        """
}

workflow bff_hashing {
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

        bff(hto_matrix, methods, methodsForConsensus, cellbarcodeWhitelist, metricsFile, doTSNE,
            doHeatmap, perCellSaturation, majorityConsensusThreshold, chemistry, callerDisagreementThreshold,
            assignmentOutBff, preprocess_bff, barcodeWhitelist)

  emit:
        bff.out.collect()
}

workflow {
    bff_hashing()
}
