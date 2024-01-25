#!/usr/bin/env nextflow

nextflow.enable.dsl=2
include { run_multi } from "$projectDir/modules/multi_demultiplexing"
include {run_single} from "$projectDir/modules/single_demultiplexing"

workflow {
    // Here we decide if it is a single sample demultiplexing or multi input demutliplexing run.
    if (params.multi_input == null){
        // Single Mode
        run_single()
    }
    else{
        // Multi mode
        run_multi()
    }
}
