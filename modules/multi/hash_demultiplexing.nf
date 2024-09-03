#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { preprocessing_hashing as preprocessing_hashing_htodemux } from './hash_demulti/preprocess'
include { preprocessing_hashing as preprocessing_hashing_multiseq } from './hash_demulti/preprocess'
include { multiseq_hashing } from './hash_demulti/multiseq'
include { htodemux_hashing } from './hash_demulti/htodemux'
include { hash_solo_hashing } from './hash_demulti/hashsolo'
include { hashedDrops_hashing } from './hash_demulti/hashedDrops'
include { demuxem_hashing } from './hash_demulti/demuxem'
include { gmm_demux_hashing } from './hash_demulti/gmm_demux'
include { bff_hashing } from './hash_demulti/bff'

process summary{
    publishDir "$params.outdir/$sampleId/$params.mode/hash_demulti", mode: 'copy'
    label 'small_mem'
        
    conda "pandas scanpy mudata"

    input:
        tuple(val(sampleId), path(hto_matrix, stageAs: 'hto_data'), path(rna_matrix, stageAs: 'rna_data'), val(demuxem_result), val(hashedDrops_result), val(hashsolo_result), val(multiseq_result), val(htodemux_result), val(gmmDemux_result), val(bff_result))
        val generate_anndata
        val generate_mudata
        
    output:
        tuple val(sampleId), path("hash_summary")

    script:
        def htodemux_files = ""
        def hashsolo_files = ""
        def multiseq_files = ""
        def hashedDrops_files = ""
        def gmmDemux_files = ""
        def demuxem_files = ""
        def bff_files = ""
        def generate_adata = ""
        def generate_mdata = ""
        
        if (demuxem_result != "no_result"){
            demuxem_files = "--demuxem ${demuxem_result}"
        }
        if (hashedDrops_result != "no_result"){
            hashedDrops_files = "--hashedDrops ${hashedDrops_result}"
        }
        if (hashsolo_result != "no_result"){
            hashsolo_files = "--hashsolo ${hashsolo_result}"
        }
        if (multiseq_result != "no_result"){
            multiseq_files = "--multiseq ${multiseq_result}"
        }
        if (htodemux_result != "no_result"){
            htodemux_files = "--htodemux ${htodemux_result}"
        }
        if (gmmDemux_result != "no_result"){
            gmmDemux_files = "--gmm_demux ${gmmDemux_result}"
        }
        if (bff_result != "no_result"){
            bff_files = "--bff ${bff_result}"
        }
        
        if (generate_anndata == "True"){
            if(rna_matrix.name == "None"){
                error "Error: RNA count matrix is not given."
            }
            generate_adata = "--generate_anndata --read_rna_mtx rna_data"
            println "AnnData created"
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
            summary_hash.py $demuxem_files $htodemux_files $multiseq_files $hashedDrops_files $hashsolo_files $gmmDemux_files $bff_files $generate_adata $generate_mdata
        """
}


workflow hash_demultiplexing{
    take:
        input_channel
    main:
            
        if (params.htodemux == "True"){
            input_channel.splitCsv(header:true).map { row-> tuple(row.sampleId, params.hto_matrix_htodemux == "raw" ? row.hto_matrix_raw : row.hto_matrix_filtered,
                                        params.rna_matrix_htodemux == "raw" ? row.rna_matrix_raw : row.rna_matrix_filtered)}.set{input_list_preprocess_htodemux}

                    preprocessing_hashing_htodemux(input_list_preprocess_htodemux, params.hto_matrix_htodemux, params.rna_matrix_htodemux) 
                    htodemux_preprocess_out = preprocessing_hashing_htodemux.out
                    htodemux_hashing(htodemux_preprocess_out)
                    htodemux_out = htodemux_hashing.out
                    println htodemux_out
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
                input_channel.splitCsv(header:true).map { row-> tuple(row.sampleId, params.hto_matrix_multiseq == "raw" ? row.hto_matrix_raw : row.hto_matrix_filtered,
                                        params.rna_matrix_multiseq == "raw" ? row.rna_matrix_raw : row.rna_matrix_filtered)}.set {input_list_preprocess_multiseq}

                    preprocessing_hashing_multiseq(input_list_preprocess_multiseq, params.hto_matrix_multiseq, params.rna_matrix_multiseq) 
                    multiseq_preprocess_out = preprocessing_hashing_multiseq.out
            }
            multiseq_hashing(multiseq_preprocess_out)
            multiseq_out = multiseq_hashing.out
            println multiseq_out
        }
        else{
            multiseq_out = channel.value("no_result")
        }
        
        if (params.hashsolo == "True"){
            input_channel \
                    | splitCsv(header:true) \
                    | map { row-> tuple(row.sampleId, params.hto_matrix_hashsolo == "raw" ? row.hto_matrix_raw : row.hto_matrix_filtered,
                                        params.rna_matrix_hashsolo == "False" ? channel.value("None") : 
                                        (params.rna_matrix_hashsolo == "raw" ? row.rna_matrix_raw : row.rna_matrix_filtered)
                                        )}
                    | hash_solo_hashing
            hashsolo_out = hash_solo_hashing.out
        }
        else{
            hashsolo_out = channel.value("no_result")
        }
        
        if (params.demuxem == "True"){
            input_channel \
                    | splitCsv(header:true) \
                    | map { row-> tuple(row.sampleId, params.hto_matrix_demuxem == "raw" ? row.hto_matrix_raw : row.hto_matrix_filtered,
                                        params.rna_matrix_demuxem == "raw" ? row.rna_matrix_raw : row.rna_matrix_filtered)}
                    | demuxem_hashing
            demuxem_out = demuxem_hashing.out
        }
        else{
            demuxem_out = channel.value("no_result")
        }
        
        if (params.hashedDrops == "True"){
            input_channel \
                    | splitCsv(header:true) \
                    | map { row-> tuple(row.sampleId, params.hto_matrix_hashedDrops == "raw" ? row.hto_matrix_raw : row.hto_matrix_filtered )}
                    | hashedDrops_hashing
            hashedDrops_out = hashedDrops_hashing.out
        }
        else{ 
            hashedDrops_out = channel.value("no_result")
        }
        

        if (params.bff == "True"){
            input_channel \
                    | splitCsv(header:true) \
                    | map { row-> tuple(row.sampleId, params.hto_matrix_bff == "raw" ? row.hto_matrix_raw : row.hto_matrix_filtered )}
                    | bff_hashing
            bff_out = bff_hashing.out
        }
        else{
            bff_out = channel.value("no_result")
        }
        if (params.gmmDemux == "True"){
            input_channel \
                    | splitCsv(header:true) \
                    | map { row-> tuple(row.sampleId, 
                                        params.hto_matrix_gmm_demux == "raw" ? row.hto_matrix_raw : row.hto_matrix_filtered,
                                        params.hto_name_gmm )}
                    | gmm_demux_hashing
            gmmDemux_out = gmm_demux_hashing.out
        }
        else{
            gmmDemux_out = channel.value("no_result")
        }

        //////////
        //Summary
        //////////
        
        input_list_summary = input_channel.splitCsv(header:true).map { row-> tuple(row.sampleId, file(row.hto_matrix_filtered), file(row.rna_matrix_filtered))}

        htodemux_out_ch = htodemux_out.flatten().map{r1-> tuple(    "$r1".replaceAll(".*htodemux_",""), r1 )}
        multiseq_out_ch = multiseq_out.flatten().map{r1-> tuple(    "$r1".replaceAll(".*multiseq_",""), r1 )}
        hashsolo_out_ch = hashsolo_out.flatten().map{r1-> tuple(    "$r1".replaceAll(".*hashsolo_",""), r1 )}
        demuxem_out_ch = demuxem_out.flatten().map{r1-> tuple(    "$r1".replaceAll(".*demuxem_",""), r1 )}
        hashedDrops_out_ch = hashedDrops_out.flatten().map{r1-> tuple(    "$r1".replaceAll(".*hashedDrops_",""), r1 )}
        bff_out_ch = bff_out.flatten().map{r1-> tuple(    "$r1".replaceAll(".*bff_",""), r1 )}
        gmmDemux_out_ch = gmmDemux_out.flatten().map{r1-> tuple(    "$r1".replaceAll(".*gmmDemux",""), r1 )}
        
        summary_input = input_list_summary.join(demuxem_out_ch,by:0,remainder: true).join(hashedDrops_out_ch,by:0,remainder: true).join(hashsolo_out_ch,by:0,remainder: true).join(multiseq_out_ch,by:0,remainder: true).join(htodemux_out_ch,by:0,remainder: true).join(gmmDemux_out_ch,by:0,remainder: true).join(bff_out_ch,by:0,remainder: true)
        summary_input = summary_input.filter{ it[0] != 'no_result' }
        
        summary(summary_input,
                params.generate_anndata, params.generate_mudata)
                
    emit:
        summary.out
}
