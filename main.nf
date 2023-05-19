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

process generate_data{
    publishDir "$projectDir/$params.outdir/$params.mode/data_output", mode: 'copy'
    input:
        path assignment
        val generate_anndata
        val generate_mudata
        path rna_matrix
        path hto_matrix 
    output:
        path "adata_with_donor_matching.h5ad", optional: true
        path "mudata_with_donor_matching.h5mu", optional: true


    script:
        def generate_adata = ""
        def generate_mdata = ""

        if (generate_anndata == "True"){
            if(rna_matrix.name == "None"){
                error "Error: RNA count matrix is not given."
            }
            generate_adata = "--generate_anndata --read_rna_mtx $rna_matrix"
        }
        if (generate_mudata == "True"){
            if(rna_matrix.name == "None"){
                error "Error: RNA count matrix is not given."
            }
            if(hto_matrix.name == "None"){
                error "Error: HTO count matrix is not given."
            }
            generate_mdata = "--generate_mudata --read_rna_mtx $rna_matrix --read_hto_mtx $hto_matrix"
        }
        
        """
        generate_data.py --assignment $assignment $generate_adata $generate_mdata
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
        generate_data(donor_match.out, params.generate_anndata, params.generate_mudata, file(params.rna_matrix), file(params.hto_matrix))
    }
    else if (params.mode == "donor_match"){
        donor_match(params.demultiplexing_result)
        generate_data(donor_match.out, params.generate_anndata, params.generate_mudata, file(params.rna_matrix), file(params.hto_matrix))
    }
}
