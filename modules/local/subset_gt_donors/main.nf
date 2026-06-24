process SUBSET_GT_DONORS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4e/4eb863d427b6e327abab4be2112d17760f86947eeaf3a214e1450bc680f14a49/data':
        'community.wave.seqera.io/library/bcftools:1.22--a51ee80717c2467e' }"

    input:
    tuple val(meta), path(subset_variants), val(output_basename), path(gt_donors_vcf), path(donor_match)

    output:
    tuple val(meta), path("*_${output_basename}.vcf.gz"), emit: donor_subset_vcf
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools sort $gt_donors_vcf -Oz -o ${prefix}_sorted.vcf.gz
    bcftools index ${prefix}_sorted.vcf.gz
    bcftools filter ${prefix}_sorted.vcf.gz -R $subset_variants -Oz -o ${prefix}_filtered.vcf.gz
    bcftools reheader ${prefix}_filtered.vcf.gz --samples $donor_match -o ${prefix}_${output_basename}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${output_basename}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
