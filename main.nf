#!/usr/bin/env nextflow

nextflow.enable.dsl=2
include { summary } from "$projectDir/modules/multi/gene_demultiplexing"
include { HADGE; SUMMARY } from "$projectDir/subworkflows/HADGE"

// Main entry point in the pipeline
workflow {
    HADGE()
}

// Enty point to only generate summary files
workflow SUMMARY_ONLY{
    SUMMARY()
}
