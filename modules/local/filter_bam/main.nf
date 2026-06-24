process FILTER_BAM {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f4/f46d4b6a720d442979b57330a319a25863231b2c70c80348f7b1d1d7d422b1f6/data'
        : 'community.wave.seqera.io/library/bcftools_bedtools_samtools:f1acc4ec7fbdba9e'}"

    input:
    tuple val(meta), path(bam), path(barcodes)
    path vcf

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path 'versions.yml', emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    bcftools sort ${vcf} -Oz -o sorted.vcf.gz
    filter_bam_file_for_popscle_dsc_pileup.sh ${bam} ${barcodes} sorted.vcf.gz ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
