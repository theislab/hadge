#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { gene_demultiplexing } from "$projectDir/modules/single/gene_demultiplexing"
include { hash_demultiplexing } from "$projectDir/modules/single/hash_demultiplexing"
include { donor_match } from "$projectDir/modules/single/donor_match"


process generate_data{
    publishDir "$params.outdir/$params.mode/data_output", mode: 'copy'

    conda "pandas scanpy mudata"

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

process summary_all{
    publishDir "$params.outdir/$params.mode", mode: 'copy'
    label 'small_mem'

    conda "pandas scanpy mudata"

    input:
        path gene_demulti_result
        path hash_demulti_result
    output:
        path "summary"

    script:
        """
        summary.py --gene_demulti $gene_demulti_result --hash_demulti $hash_demulti_result
        """
}




workflow run_single{

    print("-----Running single sample-----")

    if (params.mode == "genetic"){

        // Performing genetic demultiplexing methodologies
        gene_demultiplexing()
        if (params.match_donor == "True"){
            donor_match(gene_demultiplexing.out)
        }
    }
    else if (params.mode == "hashing"){

        // Performing hashing demultplexing
        hash_demultiplexing(params.rna_matrix_raw, params.rna_matrix_filtered, params.hto_matrix_raw, params.hto_matrix_filtered)
        if (params.match_donor == "True"){
            donor_match(hash_demultiplexing.out)
        }
    }
    else if (params.mode == "rescue"){
        
        // Performing both hashing and genetic demultiplexing methods
        hash_demultiplexing(params.rna_matrix_raw, params.rna_matrix_filtered, params.hto_matrix_raw, params.hto_matrix_filtered)
        gene_demultiplexing()
        gene_summary = gene_demultiplexing.out
        hash_summary = hash_demultiplexing.out
        summary_all(gene_summary, hash_summary)

        if (params.match_donor == "True"){
            donor_match(summary_all.out)
            if (params.generate_anndata == "True" || params.generate_mudata == "True" ){
                generate_data(donor_match.out, params.generate_anndata, params.generate_mudata, 
                params.rna_matrix_filtered, params.hto_matrix_filtered)
            }
        }
    }
    else if (params.mode == "donor_match"){

        // Performing just donor matching
        donor_match(params.demultiplexing_result)
        if (params.generate_anndata == "True" || params.generate_mudata == "True" ){
            generate_data(donor_match.out, params.generate_anndata, params.generate_mudata, 
            params.rna_matrix_filtered, params.hto_matrix_filtered)
        }

    }
}