#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { preprocessing_hashing } from './hash_demulti/preprocess'
include { multiseq_hashing } from './hash_demulti/multiseq'
include { htodemux_hashing } from './hash_demulti/htodemux'
include { hash_solo_hashing } from './hash_demulti/hashsolo'
include { hashedDrops_hashing } from './hash_demulti/hashedDrops'
include { demuxem_hashing } from './hash_demulti/demuxem'
include { solo_hashing } from './hash_demulti/solo'

process summary{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti", mode: 'copy'
    label 'small_mem'
    
    input:
        val demuxem_result
        val hashsolo_result
        val htodemux_result
        val multiseq_result
        val hashedDrops_result
        val solo_result
        val generate_anndata
        val generate_mudata
        path rna_matrix
        path hto_matrix
    
    output:
        path hash_summary

    script:
        def demuxem_files = ""
        def htodemux_files = ""
        def hashsolo_files = ""
        def multiseq_files = ""
        def hashedDrops_files = ""
        def solo_files = ""
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
        if (solo_result != "no_result"){
            solo_files = "--solo ${solo_result.join(":")}"
        }
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
        summary_hash.py $demuxem_files $htodemux_files $multiseq_files $hashedDrops_files $hashsolo_files $solo_files $generate_adata $generate_mdata
        """
}



workflow hash_demultiplexing{
    main:
    if ((params.htodemux == "True" & params.htodemux_preprocess != "False")| \
       (params.multiseq == "True" & params.multiseq_preprocess != 'False')){
        preprocessing_hashing()
    }
    
    if (params.htodemux == "True"){
        rdsobj = params.htodemux_preprocess == 'True'? preprocessing_hashing.out: (params.htodemux_preprocess == 'False'? Channel.from(params.rdsObj_htodemux) : preprocessing_hashing.out.mix(Channel.from(params.rdsObj_htodemux)))
        htodemux_hashing(rdsobj)
        htodemux_out = htodemux_hashing.out
    }
    else{
        htodemux_out = channel.value("no_result")
    }
    
    if (params.multiseq == "True"){
        rdsobj = params.multiseq_preprocess == 'True'? preprocessing_hashing.out: (params.multiseq_preprocess == 'False'? Channel.from(params.rdsObj_multiseq) : preprocessing_hashing.out.mix(Channel.from(params.rdsObj_multiseq)))
        multiseq_hashing(rdsobj)
        multiseq_out = multiseq_hashing.out
    }
    else{
        multiseq_out = channel.value("no_result")
    }
    
    if (params.hashsolo == "True"){
        hash_solo_hashing()
        hashsolo_out = hash_solo_hashing.out
    }
    else{
        hashsolo_out = channel.value("no_result")
    }
    
    if (params.demuxem == "True"){
        demuxem_hashing()
        demuxem_out = demuxem_hashing.out
    }
    else{
        demuxem_out = channel.value("no_result")
    }
    
    if (params.hashedDrops == "True"){
        hashedDrops_hashing()
        hashedDrops_out = hashedDrops_hashing.out
    }
    else{
        hashedDrops_out = channel.value("no_result")
    }
    
    if (params.solo == "True"){
        solo_hashing()
        solo_out = solo_hashing.out
    }
    else{
        solo_out = channel.value("no_result")
    }
    
    summary(demuxem_out, hashsolo_out, htodemux_out, multiseq_out, hashedDrops_out, solo_out, 
            params.generate_anndata, params.generate_mudata,
            file(params.rna_matrix), file(params.hto_matrix))
    emit:
    summary.out
}
