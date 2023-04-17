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
    input:
        val demuxem_result
        val hashsolo_result
        val htodemux_result
        val multiseq_result
        val hashedDrops_result
        val solo_result
        val select
    
    output:
        path hash_summary

    script:
        def demuxem_files = ""
        def htodemux_files = ""
        def hashsolo_files = ""
        def multiseq_files = ""
        def hashedDrops_files = ""
        def solo_files = ""
        
        if (demuxem_result != "no_result"){
            demuxem_files = "--demuxem "
            for(r : demuxem_result) {
                demuxem_files  = demuxem_files + r + ":"
            }
        }
        if (hashsolo_result != "no_result"){
            hashsolo_files = "--hashsolo "
            for(r : hashsolo_result) {
                hashsolo_files = hashsolo_files + r + ":"
            }
        }
        if (htodemux_result != "no_result"){
            htodemux_files = "--htodemux "
            for(r : htodemux_result) {
                htodemux_files = htodemux_files + r + ":"
            }
        }
        if (multiseq_result != "no_result"){
            multiseq_files = "--multiseq "
            for(r : multiseq_result) {
                multiseq_files = multiseq_files + r + ":"
            }
        }
        if (hashedDrops_result != "no_result"){
            hashedDrops_files = "--hashedDrops "
            for(r : hashedDrops_result) {
                hashedDrops_files = hashedDrops_files + r + ":"
            }
        }
        if (solo_result != "no_result"){
            solo_files = "--solo "
            for(r : solo_result) {
                solo_files = solo_files + r + ":"
            }
        }
        
        """
        mkdir hash_summary && cd hash_summary
        summary_hash.R --select $select $demuxem_files $htodemux_files $multiseq_files $hashedDrops_files $hashsolo_files $solo_files
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
        rdsobj = params.multiseq_preprocess == 'True'? preprocessing_hashing.out: (params.multiseq_preprocess == 'False'? Channel.from(params.rdsObj_multi) : preprocessing_hashing.out.mix(Channel.from(params.rdsObj_multi)))
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
    
    summary(demuxem_out, hashsolo_out, htodemux_out, multiseq_out, hashedDrops_out, solo_out, params.select)
    emit:
    summary.out
}
