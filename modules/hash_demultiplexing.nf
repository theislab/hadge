#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { preprocessing_hashing as preprocessing_hashing_htodemux } from './hash_demulti/preprocess'
include { preprocessing_hashing as preprocessing_hashing_multiseq } from './hash_demulti/preprocess'
include { multiseq_hashing } from './hash_demulti/multiseq'
include { htodemux_hashing } from './hash_demulti/htodemux'
include { hash_solo_hashing } from './hash_demulti/hashsolo'
include { hashedDrops_hashing } from './hash_demulti/hashedDrops'
include { demuxem_hashing } from './hash_demulti/demuxem'
include { gmm_demux_hashing } from './hash_demulti/gmm_demux'
include { bff_hashing } from './hash_demulti/bff'

process summary {
    publishDir "$projectDir/$params.outdir/$sampleId/$params.mode/hash_demulti", mode: 'copy'
    label 'small_mem'

    conda 'pandas scanpy mudata'

    input:
        tuple val(sampleId), path(hto_matrix, stageAs: 'hto_data'), path(rna_matrix, stageAs: 'rna_data')
        val demuxem_result
        val hashsolo_result
        val htodemux_result
        val multiseq_result
        val hashedDrops_result
        val bff_result
        val gmmDemux_result
        val generate_anndata
        val generate_mudata

    output:
        tuple val(sampleId), path('hash_summary')

    script:
        def htodemux_files = ''
        def hashsolo_files = ''
        def multiseq_files = ''
        def hashedDrops_files = ''
        def gmmDemux_files = ''
        def demuxem_files = ''
        def bff_files = ''
        def generate_adata = ''
        def generate_mdata = ''

    if (demuxem_result != 'no_result') {
            demuxem_res = demuxem_result.find { it.name.contains(sampleId) }
            demuxem_files = "--demuxem ${demuxem_res}"
    }
        if (hashsolo_result != 'no_result') {
            hashsolo_res = hashsolo_result.find { it.name.contains(sampleId) }
            hashsolo_files = "--hashsolo ${hashsolo_res}"
        }
        if (htodemux_result != 'no_result') {
            htodemux_res = htodemux_result.find { it.name.contains(sampleId) }
            htodemux_files = "--htodemux ${htodemux_res}"
        }
        if (multiseq_result != 'no_result') {
            multiseq_res = multiseq_result.find { it.name.contains(sampleId) }
            multiseq_files = "--multiseq ${multiseq_res}"
        }
        if (hashedDrops_result != 'no_result') {
            hashedDrops_res = hashedDrops_result.find { it.name.contains(sampleId) }
            hashedDrops_files = "--hashedDrops ${hashedDrops_res}"
        }
        if (gmmDemux_result != 'no_result') {
            gmmDemux_res = gmmDemux_result.find { it.name.contains(sampleId) }
            gmmDemux_files = "--gmm_demux ${gmmDemux_res}"
        }
        if (bff_result != 'no_result') {
            bff_res = bff_result.find { it.name.contains(sampleId) }
            bff_files = "--bff ${bff_res}"
        }
        if (generate_anndata == 'True') {
            if (rna_matrix.name == 'None') {
                error 'Error: RNA count matrix is not given.'
            }
            generate_adata = '--generate_anndata --read_rna_mtx rna_data'
        }
        if (generate_mudata == 'True') {
            if (rna_matrix.name == 'None') {
                error 'Error: RNA count matrix is not given.'
            }
            if (hto_matrix.name == 'None') {
                error 'Error: HTO count matrix is not given.'
            }
            generate_mdata = '--generate_mudata --read_rna_mtx rna_data --read_hto_mtx hto_data'
        }

        """
        summary_hash.py $demuxem_files $htodemux_files $multiseq_files $hashedDrops_files $hashsolo_files $gmmDemux_files $bff_files $generate_adata $generate_mdata --sampleId $sampleId
        """
}

workflow hash_demultiplexing {
    if (params.htodemux == 'True') {
        Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, params.hto_matrix_htodemux == 'raw' ? row.hto_matrix_raw : row.hto_matrix_filtered,
                                    params.rna_matrix_htodemux == 'raw' ? row.rna_matrix_raw : row.rna_matrix_filtered)}
                | set { input_list_preprocess_htodemux }
                preprocessing_hashing_htodemux(input_list_preprocess_htodemux, params.hto_matrix_htodemux, params.rna_matrix_htodemux)
                htodemux_preprocess_out = preprocessing_hashing_htodemux.out
                htodemux_hashing(htodemux_preprocess_out)
                htodemux_out = htodemux_hashing.out
    }
        else {
            htodemux_out = channel.value('no_result')
        }

    if (params.multiseq == 'True') {
        if (params.htodemux == 'True' & params.hto_matrix_htodemux == params.hto_matrix_multiseq &
            params.rna_matrix_htodemux == params.rna_matrix_multiseq) {
            multiseq_preprocess_out = htodemux_preprocess_out
            }
        else {
            Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, params.hto_matrix_multiseq == 'raw' ? row.hto_matrix_raw : row.hto_matrix_filtered,
                                    params.rna_matrix_multiseq == 'raw' ? row.rna_matrix_raw : row.rna_matrix_filtered)}
                | set { input_list_preprocess_multiseq }
                preprocessing_hashing_multiseq(input_list_preprocess_multiseq, params.hto_matrix_multiseq, params.rna_matrix_multiseq)
                multiseq_preprocess_out = preprocessing_hashing_multiseq.out
        }
        multiseq_hashing(multiseq_preprocess_out)
        multiseq_out = multiseq_hashing.out
    }
    else {
        multiseq_out = channel.value('no_result')
    }

    if (params.hashsolo == 'True') {
        Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, params.hto_matrix_hashsolo == 'raw' ? row.hto_matrix_raw : row.hto_matrix_filtered,
                                    params.rna_matrix_hashsolo == 'False' ? channel.value('None') :
                                    (params.rna_matrix_hashsolo == 'raw' ? row.rna_matrix_raw : row.rna_matrix_filtered)
                                    )}
                | hash_solo_hashing
        hashsolo_out = hash_solo_hashing.out
    }
    else {
        hashsolo_out = channel.value('no_result')
    }

    if (params.demuxem == 'True') {
        Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, params.hto_matrix_demuxem == 'raw' ? row.hto_matrix_raw : row.hto_matrix_filtered,
                                    params.rna_matrix_demuxem == 'raw' ? row.rna_matrix_raw : row.rna_matrix_filtered)}
                | demuxem_hashing
        demuxem_out = demuxem_hashing.out
    }
    else {
        demuxem_out = channel.value('no_result')
    }

    if (params.hashedDrops == 'True') {
        Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, params.hto_matrix_hashedDrops == 'raw' ? row.hto_matrix_raw : row.hto_matrix_filtered) }
                | hashedDrops_hashing
        hashedDrops_out = hashedDrops_hashing.out
    }
    else {
        hashedDrops_out = channel.value('no_result')
    }

    if (params.bff == 'True') {
        Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, params.hto_matrix_bff == 'raw' ? row.hto_matrix_raw : row.hto_matrix_filtered) }
                | bff_hashing
        bff_out = bff_hashing.out
        print('BFF path to output')
    }
    else {
        bff_out = channel.value('no_result')
    }
    if (params.gmmDemux == 'True') {
        Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, params.hto_matrix_gmm_demux == 'raw' ? row.hto_matrix_raw : row.hto_matrix_filtered, row.hto_name_gmm) }
                | gmm_demux_hashing
        gmmDemux_out = gmm_demux_hashing.out
    }
    else {
        gmmDemux_out = channel.value('no_result')
    }

    Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, file(row.hto_matrix_filtered), file(row.rna_matrix_filtered)) }
                | set { input_list_summary }
    summary(input_list_summary, demuxem_out, hashsolo_out, htodemux_out, multiseq_out, hashedDrops_out, bff_out, gmmDemux_out, params.generate_anndata, params.generate_mudata)

    emit:
        summary.out
}
