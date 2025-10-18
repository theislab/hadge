process DONOR_MATCH {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d863e56b5ce15b271e8c8666ec22217df5cfc57a9731cc23c7f92674dc7ab0c7/data':
        'community.wave.seqera.io/library/pegasusio_mudata_numpy_pandas_pruned:ecdbf7e42b2f3213' }"

    input:
        tuple val(meta), path(barcode_whitelist), val(cell_genotype), val(vireo_parent_dir), path(demultiplexing_result)
        val method1_name
        val method2_name
        val findVariants
        val variant_count
        val variant_pct


    //TODO do variants TRUE

    // TODO do without donor_mode then demultiplexing_result not necessary

    //TODO write donor_match with findVAriants false this will make the start easier
    // then also vireo_parent_dir, cell_genotype
    // write first part of the script and then lets see
    output:
    tuple val(meta), path("*_emptyDrops.png")         , emit: empty_drops_plot
    tuple val(meta), path("*_emptyDrops.csv")         , emit: empty_drops_csv
    tuple val(meta), path("*_emptyDrops.rds")         , emit: empty_drops_rds
    tuple val(meta), path("*_results_hasheddrops.csv"), emit: results
    tuple val(meta), path("*_id_to_hash.csv")         , emit: id_to_hash
    tuple val(meta), path("*_hasheddrops.rds")        , emit: rds
    tuple val(meta), path("*_plot_hasheddrops.png")   , emit: plot
    tuple val(meta), path("*_params_hasheddrops.csv") , emit: params
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template('donor_match.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_emptyDrops.png
    touch ${prefix}_emptyDrops.csv
    touch ${prefix}_emptyDrops.rds
    touch ${prefix}_results_hasheddrops.csv
    touch ${prefix}_id_to_hash.csv
    touch ${prefix}_hasheddrops.rds
    touch ${prefix}_plot_hasheddrops.png
    touch ${prefix}_params_hasheddrops.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
        dropletutils: \$(Rscript -e "library(DropletUtils); cat(as.character(packageVersion('DropletUtils')))")
    END_VERSIONS
    """
}

}
