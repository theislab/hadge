#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { run_hadge_multi} from './modules/multi_demultiplexing'
include { gene_demultiplexing } from './modules/gene_demultiplexing'
include { hash_demultiplexing } from './modules/hash_demultiplexing'
include { gene_demultiplexing as gene_demultiplexing_single } from './modules/single/gene_demultiplexing'
include { hash_demultiplexing as hash_demultiplexing_single } from './modules/single/hash_demultiplexing'
include { donor_match as donor_match_single } from './modules/single/donor_match'

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
        val rna_matrix
        val hto_matrix 
    output:
        path "adata_with_donor_matching.h5ad", optional: true
        path "mudata_with_donor_matching.h5mu", optional: true


    script:
        def generate_adata = ""
        def generate_mdata = ""

        if (generate_anndata == "True"){
            if(rna_matrix == "None"){
                error "Error: RNA count matrix is not given."
            }
            generate_adata = "--generate_anndata --read_rna_mtx $rna_matrix"
        }
        if (generate_mudata == "True"){
            if(rna_matrix == "None"){
                error "Error: RNA count matrix is not given."
            }
            if(hto_matrix == "None"){
                error "Error: HTO count matrix is not given."
            }
            generate_mdata = "--generate_mudata --read_rna_mtx $rna_matrix --read_hto_mtx $hto_matrix"
        }
        
        """
        generate_data.py --assignment $assignment $generate_adata $generate_mdata
        """
}

workflow run_hadge_single{
    if (params.mode == "genetic"){
        gene_demultiplexing_single()
        if (params.match_donor == "True"){
            donor_match_single(gene_demultiplexing_single.out)
        }
    }
    else if (params.mode == "hashing"){
        hash_demultiplexing_single(params.rna_matrix_raw, params.rna_matrix_filtered, params.hto_matrix_raw, params.hto_matrix_filtered)
        if (params.match_donor == "True"){
            donor_match_single(hash_demultiplexing_single.out)
        }
    }
    else if (params.mode == "rescue"){
        hash_demultiplexing_single(params.rna_matrix_raw, params.rna_matrix_filtered, params.hto_matrix_raw, params.hto_matrix_filtered)
        gene_demultiplexing_single()
        gene_summary = gene_demultiplexing_single.out
        hash_summary = hash_demultiplexing_single.out
        summary_all(gene_summary, hash_summary)
        if (params.match_donor == "True"){
            donor_match_single(summary_all.out)
        }
        generate_data(donor_match_single.out, params.generate_anndata, params.generate_mudata, 
                params.rna_matrix_filtered, params.hto_matrix_filtered)
    }
    else if (params.mode == "donor_match"){
        donor_match_single(params.demultiplexing_result)
        generate_data(donor_match_single.out, params.generate_anndata, params.generate_mudata, 
            params.rna_matrix_filtered, params.hto_matrix_filtered)
    }
}

workflow {
    if (params.multi_input == "None"){
        run_hadge_single()
    }
    else{
        run_hadge_multi()
    }

}