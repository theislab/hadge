#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { hash_demultiplexing } from './hash_demultiplexing'
include { gene_demultiplexing } from './gene_demultiplexing'
include { donor_match } from './donor_match'

process generate_data{
    publishDir "$projectDir/$params.outdir/$sampleId/$params.mode/data_output", mode: 'copy'
    label 'small_mem'

    conda "pandas scanpy mudata"

    input:
        tuple val(sampleId), val(hto_matrix), val(rna_matrix), path(assignment)
        val generate_anndata
        val generate_mudata
        
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
    publishDir "$projectDir/$params.outdir/$sampleId/$params.mode", mode: 'copy'
    label 'small_mem'
    
    conda "pandas scanpy mudata"

    input:
        tuple val(sampleId), path(gene_demulti_result), path(hash_demulti_result)
    output:
        tuple val(sampleId), path("summary")

    script:
        """
        summary.py --gene_demulti $gene_demulti_result --hash_demulti $hash_demulti_result
        """
}

workflow run_multi{
    if (params.mode == "genetic"){
        gene_demultiplexing()
        if (params.match_donor == "True"){
            gene_demultiplexing.out.view()
            Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, row.nsample, row.barcodes, "None", "None")}
                | join(gene_demultiplexing.out)
                | donor_match
        }
    }
    else if (params.mode == "hashing"){
        hash_demultiplexing()
        if (params.match_donor == "True"){
            Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row -> tuple(row.sampleId, row.nsample, row.barcodes, "None", "None")}
                | join(hash_demultiplexing.out)
                | donor_match
        }
    }
    else if (params.mode == "rescue"){
        hash_demultiplexing()
        gene_demultiplexing()
        gene_summary = gene_demultiplexing.out
        hash_summary = hash_demultiplexing.out
        input_summary_all = gene_summary.join(hash_summary)
        summary_all(input_summary_all)
        if (params.match_donor == "True"){
            Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row -> tuple(row.sampleId, row.nsample, row.barcodes, "None", "None")}
                | join(summary_all.out)
                | donor_match
        }
        if (params.generate_anndata == "True" || params.generate_mudata == "True" ){
        Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row -> tuple(row.sampleId, row.hto_matrix_filtered, row.rna_matrix_filtered)}
                | join(donor_match.out)
                | set {input_generate_data}
        generate_data(input_generate_data, params.generate_anndata, params.generate_mudata)
        }
    }
    else if (params.mode == "donor_match"){
        Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row -> tuple(row.sampleId, row.nsample, row.barcodes, row.celldata, row.vireo_parent_dir, row.demultiplexing_result)} \
                | donor_match
        if (params.generate_anndata == "True" || params.generate_mudata == "True" ){
            Channel.fromPath(params.multi_input) \
                    | splitCsv(header:true) \
                    | map { row -> tuple(row.sampleId, row.hto_matrix_filtered, row.rna_matrix_filtered)}
                    | join(donor_match.out)
                    | set {input_generate_data}
            generate_data(input_generate_data, params.generate_anndata, params.generate_mudata)
        }
    }
    
}
