process POPSCLE_DSCPILEUP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/popscle:0.1beta--h2c78cec_0'
        : 'biocontainers/popscle:0.1beta--h2c78cec_0'}"

    input:
    tuple val(meta), path(bam), path(vcf)

    output:
    tuple val(meta), path("${prefix}"), emit: directory
    tuple val(meta), path("${prefix}/*.cel.gz"), emit: cel
    tuple val(meta), path("${prefix}/*.plp.gz"), emit: plp
    tuple val(meta), path("${prefix}/*.var.gz"), emit: var
    tuple val(meta), path("${prefix}/*.umi.gz"), emit: umi
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1'
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    """
    mkdir -p "${prefix}"
    popscle dsc-pileup \\
        --sam ${bam} \\
        --vcf ${vcf} \\
        --out ${prefix}/${prefix} \\
        ${args} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popscle dsc-pileup: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1'
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.cel.gz
    touch ${prefix}/${prefix}.var.gz
    touch ${prefix}/${prefix}.plp.gz
    touch ${prefix}/${prefix}.umi.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popscle dsc-pileup: ${VERSION}
    END_VERSIONS
    """
}
