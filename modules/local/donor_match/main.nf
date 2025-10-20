process DONOR_MATCH {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/45b060e69064c7a7894787b0cc29259bbb24357650d06b627912b56d3521899b/data':
        'community.wave.seqera.io/library/r-complexupset_r-data.table_r-pheatmap_r-r.utils_pruned:3bd8312041c22554' }"

    input:
        tuple val(meta), path(barcode_whitelist), path(demultiplexing_result), val(cell_genotype), val(vireo_parent_dir)
        val match_donor_method1
        val match_donor_method2
        val findVariants
        val variant_count
        val variant_pct

    //TODO do variants TRUE

    // TODO do without donor_mode then demultiplexing_result not necessary

    //TODO write donor_match with findVAriants false this will make the start easier
    // then also vireo_parent_dir, cell_genotype
    // write first part of the script and then lets see
    output:
    tuple val(meta), path("*_vs_*correlation_res.csv")               , emit: correlation_csv
    tuple val(meta), path("*_vs_*donor_match.csv")                   , emit: donor_match, optional: true
    tuple val(meta), path("*_vs_*concordance_heatmap.png")           , emit: concordance_heatmap, optional: true
    tuple val(meta), path("*_vs_*all_assignment_after_match.csv")    , emit: assignment_all_match, optional: true
    tuple val(meta), path("*_vs_*intersect_assignment_after_match.csv"), emit: assignment_intersect_match, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // stays like that if findVaraint is 0
    def cell_genotype_path = ''
    def vireo_parent_path = ''
    def ndonor = "${meta.nsample}"
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
