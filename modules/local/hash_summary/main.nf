process HASH_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-ggplot2_r-seurat:dac8c905972b98df':
        'community.wave.seqera.io/library/r-ggplot2_r-seurat:eefd54806320eae0' }"

    input:
    tuple val(meta), path(hto_matrix), path(rna_matrix), path(htodemux), path(multiseq), path(cellhashr), path(demuxem), path(gmmdemux), path(hasheddrops), path(hashsolo)
    val generate_anndata // boolean
    val generate_mudata // boolean


    output:
    tuple val(meta), path("*_hashing_assignment_summary.csv")    , emit: assignment
    tuple val(meta), path("*_hashing_classification_summary.csv"), emit: classification
    tuple val(meta), path("*_hashing_params_summary.json")       , emit: params
    script:
    prefix         = task.ext.prefix         ?: "${meta.id}"

    template 'hash_summary.py'
}
