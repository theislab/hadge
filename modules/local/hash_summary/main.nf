process HASH_SUMMARY {
    // TODO remove dubug for the pipeline
    debug true
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-ggplot2_r-seurat:dac8c905972b98df':
        'community.wave.seqera.io/library/pegasusio_anndata_mudata_numpy_pruned:9d13d0d12376624e' }"

    input:
    tuple val(meta), path(rna_matrix), path(hto_matrix), path(htodemux_assignments), path (htodemux_classification), path(multiseq), path(bff), path(demuxem), path(gmmdemux_results), path(gmmdemux_config), path(hasheddrops), path(hashsolo)
    val generate_anndata // boolean
    val generate_mudata // boolean


    output:
    tuple val(meta), path("*_hashing_summary_assignment.csv")    , emit: assignment
    tuple val(meta), path("*_hashing_summary_classification.csv"), emit: classification

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix         = task.ext.prefix         ?: "${meta.id}"

    template 'hash_summary.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_hashing_summary_assignment.csv
    touch ${prefix}_hashing_summary_classification.csv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
    END_VERSIONS
    """
}
