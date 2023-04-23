#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { gene_demultiplexing } from './modules/gene_demultiplexing'
include { hash_demultiplexing } from './modules/hash_demultiplexing'
include { donor_match } from './modules/donor_match'

process summary_all{
    publishDir "$projectDir/$params.outdir/$params.mode", mode: 'copy'
    input:
        path gene_demulti_result
        path hash_demulti_result
    output:
        path summary

    script:
        """
        summary.R --gene_demulti $gene_demulti_result --hash_demulti $hash_demulti_result
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

workflow{
    if (params.mode == "genetic"){
        gene_demultiplexing()
        if (params.match_donor == "True"){
            donor_match(gene_demultiplexing.out)
        }
    }
    else if (params.mode == "hashing"){
        hash_demultiplexing()
        if (params.match_donor == "True"){
            donor_match(hash_demultiplexing.out)
        }
    }
    else if (params.mode == "rescue"){
        hash_demultiplexing()
        gene_demultiplexing()
        gene_summary = gene_demultiplexing.out
        hash_summary = hash_demultiplexing.out
        summary_all(gene_summary, hash_summary)
        if (params.match_donor == "True"){
            donor_match(summary_all.out)
        }
    }
    else if (params.mode == "donor_match"){
        donor_match(params.demultiplexing_result)
    }
}
