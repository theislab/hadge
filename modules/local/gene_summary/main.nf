process GENE_SUMMARY {
    tag "${meta.id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    // container "<a-small-python-pandas-scanpy image>"
    input:
    tuple val(meta),
        path(rna_matrix),
        path(hto_matrix),
        path(vireo),
        path(demuxlet),
        path(freemuxlet)
    tuple val (generate_anndata), val(generate_mudata)


    output:
    tuple val(meta), path("*_genetic_summary_assignment.csv")    , emit: assignment    , optional: false
    tuple val(meta), path("*_genetic_summary_classification.csv"), emit: classification, optional: false
    tuple val(meta), path("*_genetic_summary.h5ad")              , emit: h5ad          , optional: true
    tuple val(meta), path("*_genetic_summary.h5mu")              , emit: h5mu          , optional: true
    path "versions.yml"                                          , emit: versions      , optional: false

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    template 'gene_summary.py'

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
