process DONOR_MATCH {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d9138b380ca73daad0b5ad74a10b46324ca4f676efdf199f9dc9cb9145a4590c/data':
        'community.wave.seqera.io/library/r-data.table_r-pheatmap_r-tidyverse:ac2dbc33f827dbb9' }"

    input:
        tuple val(meta), path(demultiplexing_result)
        val match_donor_method1
        val match_donor_method2

    output:
    // donor match will return the optional outputs in rescue and donor_match mode but not in genetic or hashing mode
    // best combination of a genetic and hashing deconvolution method
    tuple val(meta), path("*_best_donor_match.csv")                      , emit: best_donor_match                     , optional:true
    tuple val(meta), path("*_best_all_assignment_after_match.csv")       , emit: best_all_assignment_after_match      , optional:true
    tuple val(meta), path("*_best_intersect_assignment_after_match.csv") , emit: best_intersect_assignment_after_match, optional:true
    tuple val(meta), path("*_score_record.csv")                          , emit: score_record                         , optional:true

    // comparison between deconvolution methods
    tuple val(meta), path("*/*_vs_*all_assignment_after_match.csv")      , emit: assignment_after_match               , optional:true
    tuple val(meta), path("*/*_vs_*intersect_assignment_after_match.csv"), emit: assignment_intersect_match           , optional:true
    tuple val(meta), path("*/*_vs_*correlation_res.csv")                 , emit: correlation
    tuple val(meta), path("*/*_vs_*donor_match.csv")                     , emit: donor_match
    tuple val(meta), path("*/*_vs_*concordance_heatmap.png")             , emit: concordance_heatmap
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('donor_match.R')

    stub:
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
        r-data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')))")
        r-pheatmap: \$(Rscript -e "library(pheatmap); cat(as.character(packageVersion('pheatmap')))")
        r-tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS
    """
}
