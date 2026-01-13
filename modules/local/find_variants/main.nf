process FIND_VARIANTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/75/755e81f7b523df9db3a6f574fdb876ddfc1e1faf1e912260b99bca773f6dba2d/data':
        'community.wave.seqera.io/library/r-complexupset_r-data.table_r-tidyverse_r-vcfr:87602a1274fab432' }"

    input:
    tuple val(meta), path(best_intersect_assignment_after_match), path(cell_genotype), path(variants_vireo), path(demultiplexing_result)
    val variant_count
    val variant_pct

    output:
    tuple val(meta), path("*/*_matched_gt.csv")                 , emit: matched_gt
    tuple val(meta), path("*/*_unmatched_gt.csv")               , emit: unmatched_gt
    tuple val(meta), path("*/*_informative_variants.csv")       , emit: informative_variants
    tuple val(meta), path("*_all_representative_variants.csv")  , emit: all_representative_variants
    tuple val(meta), path("*_donor_specific_variants_upset.png"), emit: donor_specific_variants_upset
    tuple val(meta), path("*_donor_specific_variants.csv")      , emit: donor_specific_variants
    tuple val(meta), path("*_vireo_variants.csv")               , emit: vireo_variants, optional: true
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('find_variants.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p hto1

    touch ${prefix}_all_representative_variant.csv
    touch ${prefix}_donor_specific_variants_upset.png
    touch ${prefix}_donor_specific_representative_variants.csv
    touch ${prefix}_vireo_representative_variants.csv
    touch hto1/${prefix}_hto1_matched_gt.csv
    touch hto1/${prefix}_hto1_unmatched_gt.csv
    touch hto1/${prefix}_hto1_informative_variants.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(paste(R.version[['major']], R.version[['minor']], sep='.'))")
        r-complexupset: \$(Rscript -e "library(ComplexUpset); cat(as.character(packageVersion('ComplexUpset')))")
        r-data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')))")
        r-tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        r-vcfr: \$(Rscript -e "library(vcfR); cat(as.character(packageVersion('vcfR')))")
    END_VERSIONS
    """
}
