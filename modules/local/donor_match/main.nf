process DONOR_MATCH {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/45b060e69064c7a7894787b0cc29259bbb24357650d06b627912b56d3521899b/data':
        'community.wave.seqera.io/library/r-complexupset_r-data.table_r-pheatmap_r-r.utils_pruned:3bd8312041c22554' }"

    //TODO findVariant = true not implemented
    input:
        tuple val(meta), path(barcode_whitelist), path(demultiplexing_result), val(cell_genotype), val(vireo_parent_dir)
        val match_donor_method1
        val match_donor_method2
        val findVariants
        val variant_count
        val variant_pct

    output:
    // best method combination for rescue/donor_match mode (has to be optional because runs with only genetic or hashing won't return this output)
    tuple val(meta), path("*_best_donor_match.csv")                      , emit: best_donor_match                     , optional:true
    tuple val(meta), path("*_best_all_assignment_after_match.csv")       , emit: best_assignment_after_match          , optional:true
    tuple val(meta), path("*_best_intersect_assignment_after_match.csv") , emit: best_intersect_assignment_after_match, optional:true
    tuple val(meta), path("*_score_record.csv")                          , emit: score_record                         , optional:true

    // comparison between deconvolution methods
    tuple val(meta), path("*/*_vs_*all_assignment_after_match.csv")      , emit: assignment_after_match               , optional: false
    tuple val(meta), path("*/*_vs_*intersect_assignment_after_match.csv"), emit: assignment_intersect_match           , optional: false
    tuple val(meta), path("*/*_vs_*correlation_res.csv")                 , emit: correlation                          , optional: false
    tuple val(meta), path("*/*_vs_*donor_match.csv")                     , emit: donor_match                          , optional: false
    tuple val(meta), path("*/*_vs_*concordance_heatmap.png")             , emit: concordance_heatmap                  , optional: false
    path "versions.yml"                                                  , emit: versions                             , optional: false

    when:
    task.ext.when == null || task.ext.when

    script:
    // TODO for findVariant = true (not used by findVariant = false)
    def cell_genotype_path = ''
    def vireo_parent_path = ''
    def ndonor = "${meta.n_sample}"
    template('donor_match.R')

    stub:
    //TODO is method1_vs_method2 correct for stub (number of new directories depends on the input data)?
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p method1_vs_method2

    touch ${prefix}_best_donor_match.csv
    touch ${prefix}_best_all_assignment_after_match.csv
    touch ${prefix}_best_intersect_assignment_after_match.csv
    touch ${prefix}_score_record.csv
    touch method1_vs_method2/${prefix}_method1_vs_method2_all_assignment_after_match.csv
    touch method1_vs_method2/${prefix}_method1_vs_method2_intersect_assignment_after_match.csv
    touch method1_vs_method2/${prefix}_method1_vs_method2_correlation_res.csv
    touch method1_vs_method2/${prefix}_method1_vs_method2_donor_match.csv
    touch method1_vs_method2/${prefix}_method1_vs_method2_concordance_heatmap.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(paste(R.version[['major']], R.version[['minor']], sep='.'))")
        r-complexupset: \$(Rscript -e "library(ComplexUpset); cat(as.character(packageVersion('ComplexUpset')))")
        r-data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')))")
        r-pheatmap: \$(Rscript -e "library(pheatmap); cat(as.character(packageVersion('pheatmap')))")
        r-tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        r-vcfr: \$(Rscript -e "library(vcfR); cat(as.character(packageVersion('vcfR')))")
    END_VERSIONS
    """
}
