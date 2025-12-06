process FIND_VARIANTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/45b060e69064c7a7894787b0cc29259bbb24357650d06b627912b56d3521899b/data':
        'community.wave.seqera.io/library/r-complexupset_r-data.table_r-pheatmap_r-r.utils_pruned:3bd8312041c22554' }"

    input:
    tuple val(meta), path(best_intersect_assignment_after_match), path(cell_genotype), path(variants_vireo), path(demultiplexing_result)
    val variant_count
    val variant_pct

    output:
    tuple val(meta), path("*/*_matched_gt.csv")                       , emit: matched_gt
    tuple val(meta), path("*/*_unmatched_gt.csv")                     , emit: unmatched_gt
    tuple val(meta), path("*/*_informative_variants.csv")             , emit: informative_variants
    tuple val(meta), path("*_all_representative_variant_df.csv")      , emit: all_representative_variant_df
    tuple val(meta), path("*_donor_specific_variants_upset.png")      , emit: donor_specific_variants_upset
    tuple val(meta), path("*_donor_match_representative_variants.csv"), emit: donor_match_representative_variants
    tuple val(meta), path("*_vireo_representative_variants.csv")      , emit: vireo_representative_variants, optional: true
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('find_variants.R')

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
        r-complexupset: \$(Rscript -e "library(ComplexUpset); cat(as.character(packageVersion('ComplexUpset')))")
        r-data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')))")
        r-pheatmap: \$(Rscript -e "library(pheatmap); cat(as.character(packageVersion('pheatmap')))")
        r-tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        r-vcfr: \$(Rscript -e "library(vcfR); cat(as.character(packageVersion('vcfR')))")
    END_VERSIONS
    """
}
