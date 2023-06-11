#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { preprocessing_hashing as preprocessing_hashing_htodemux } from './hash_demulti/preprocess'
include { preprocessing_hashing as preprocessing_hashing_multiseq } from './hash_demulti/preprocess'
include { multiseq_hashing } from './hash_demulti/multiseq'
include { htodemux_hashing } from './hash_demulti/htodemux'
include { hash_solo_hashing } from './hash_demulti/hashsolo'
include { hashedDrops_hashing } from './hash_demulti/hashedDrops'
include { demuxem_hashing } from './hash_demulti/demuxem'

process summary{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti", mode: 'copy'
    label 'small_mem'
    label 'summary'
    input:
        val demuxem_result
        val hashsolo_result
        val htodemux_result
        val multiseq_result
        val hashedDrops_result
        val generate_anndata
        val generate_mudata
        path rna_matrix, stageAs: 'rna_data'
        path hto_matrix, stageAs: 'hto_data'
    
    output:
        path hash_summary

    script:
        def demuxem_files = ""
        def htodemux_files = ""
        def hashsolo_files = ""
        def multiseq_files = ""
        def hashedDrops_files = ""
        def generate_adata = ""
        def generate_mdata = ""
        
        if (demuxem_result != "no_result"){
            demuxem_files = "--demuxem ${demuxem_result.join(":")}"
        }
        if (hashsolo_result != "no_result"){
            hashsolo_files = "--hashsolo ${hashsolo_result.join(":")}"
        }
        if (htodemux_result != "no_result"){
            htodemux_files = "--htodemux ${htodemux_result.join(":")}"
        }
        if (multiseq_result != "no_result"){
            multiseq_files = "--multiseq ${multiseq_result.join(":")}"
        }
        if (hashedDrops_result != "no_result"){
            hashedDrops_files = "--hashedDrops ${hashedDrops_result.join(":")}"
        }
        if (generate_anndata == "True"){
            if(rna_matrix.name == "None"){
                error "Error: RNA count matrix is not given."
            }
            generate_adata = "--generate_anndata --read_rna_mtx rna_data"
        }
        if (generate_mudata == "True"){
            if(rna_matrix.name == "None"){
                error "Error: RNA count matrix is not given."
            }
            if(hto_matrix.name == "None"){
                error "Error: HTO count matrix is not given."
            }
            generate_mdata = "--generate_mudata --read_rna_mtx rna_data --read_hto_mtx hto_data"
        }
        
        """
        summary_hash.py $demuxem_files $htodemux_files $multiseq_files $hashedDrops_files $hashsolo_files $generate_adata $generate_mdata
        """
}


workflow hash_demultiplexing{
    take:
        rna_matrix_raw
        rna_matrix_filtered
        hto_matrix_raw
        hto_matrix_filtered
    main:
    
    if (params.htodemux == "True"){
        rna_matrix = params.rna_matrix_htodemux == "raw" ? rna_matrix_raw : rna_matrix_filtered
        hto_matrix = params.hto_matrix_htodemux == "raw" ? hto_matrix_raw : hto_matrix_filtered
        preprocessing_hashing_htodemux(hto_matrix, rna_matrix, params.hto_matrix_htodemux, params.rna_matrix_htodemux)
        htodemux_preprocess_out = preprocessing_hashing_htodemux.out
        htodemux_hashing(htodemux_preprocess_out)
        htodemux_out = htodemux_hashing.out
    }
    else{
        htodemux_out = channel.value("no_result")
    }
    
    if (params.multiseq == "True"){
        if (params.htodemux == "True" & params.hto_matrix_htodemux == params.hto_matrix_multiseq & 
            params.rna_matrix_htodemux == params.rna_matrix_multiseq){
            multiseq_preprocess_out = htodemux_preprocess_out
        }
        else{
            rna_matrix = params.rna_matrix_multiseq == "raw" ? rna_matrix_raw : rna_matrix_filtered
            hto_matrix = params.hto_matrix_multiseq == "raw" ? hto_matrix_raw : hto_matrix_filtered
            preprocessing_hashing_multiseq(hto_matrix, rna_matrix, params.hto_matrix_multiseq, params.rna_matrix_multiseq)
            multiseq_preprocess_out = preprocessing_hashing_multiseq.out
        }
        multiseq_hashing(multiseq_preprocess_out)
        multiseq_out = multiseq_hashing.out
    }
    else{
        multiseq_out = channel.value("no_result")
    }
    
    if (params.hashsolo == "True"){
        hashsolo_hto_input = params.hto_matrix_hashsolo = "raw" ? hto_matrix_raw : hto_matrix_filtered
        hashsolo_rna_input = params.rna_matrix_hashsolo == "False" ? channel.value("None") : 
                            (params.rna_matrix_hashsolo = "raw" ? rna_matrix_raw : rna_matrix_filtered)
        hash_solo_hashing(hashsolo_hto_input, hashsolo_rna_input)
        hashsolo_out = hash_solo_hashing.out
    }
    else{
        hashsolo_out = channel.value("no_result")
    }
    
    if (params.demuxem == "True"){
        demuxem_hto_input = params.hto_matrix_demuxem = "raw" ? hto_matrix_raw : hto_matrix_filtered
        demuxem_rna_input = params.rna_matrix_demuxem = "raw" ? rna_matrix_raw : rna_matrix_filtered
        demuxem_hashing(demuxem_hto_input, demuxem_rna_input)
        demuxem_out = demuxem_hashing.out
    }
    else{
        demuxem_out = channel.value("no_result")
    }
    
    if (params.hashedDrops == "True"){
        hashedDrops_hto_input = params.hto_matrix_hashedDrops = "raw" ? hto_matrix_raw : hto_matrix_filtered
        hashedDrops_hashing(hashedDrops_hto_input)
        hashedDrops_out = hashedDrops_hashing.out
    }
    else{
        hashedDrops_out = channel.value("no_result")
    }
    
    summary(demuxem_out, hashsolo_out, htodemux_out, multiseq_out, hashedDrops_out,
            params.generate_anndata, params.generate_mudata, rna_matrix_filtered, hto_matrix_filtered)
    emit:
        summary.out
}
