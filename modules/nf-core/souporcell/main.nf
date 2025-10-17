process SOUPORCELL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9a/9a69c552c52aa5b3636a7a596f9406b2ec3e165809ccd58a012b9ea285ba6ecd/data' :
        'community.wave.seqera.io/library/souporcell_gxx:f648658dde2cdd53' }"

    input:
    tuple val(meta), path(bam), path(barcodes)
    tuple val(meta2), path(fasta)
    val(clusters)

    output:
    tuple val(meta), path("*/*.vcf"), emit: vcf
    tuple val(meta), path("*/*clusters.tsv"), emit: tsv_result
    tuple val(meta), path("*/*clusters_tmp.tsv"), emit: tsv_tmp
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    """
    mkdir -p temp
    export TMPDIR=./temp
    souporcell_pipeline.py \\
        -i $bam \\
        -b $barcodes \\
        -f $fasta \\
        -t $task.cpus \\
        -o $prefix \\
        -k $clusters \\
        $args


    # Build a single-line, YAML-safe version string
    souporcell_version=`souporcell_pipeline.py --version 2>&1 || echo 'N/A'`
    [ -z "\$souporcell_version" ] && souporcell_version=`souporcell --version 2>&1 || echo 'N/A'`
    # squash newlines -> space; replace double quotes with single quotes
    souporcell_version=`echo "\$souporcell_version" | tr '\\n' ' ' | sed 's/\"/'"'"'/g'`

    # Write versions.yml safely (Groovy expands ${task.process}; bash expands the rest)
    printf '"%s":\\n    souporcell: "%s"\\n' "${task.process}" "\$souporcell_version" > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    touch ${prefix}/cluster_genotypes.vcf
    touch ${prefix}/clusters.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        souporcell: \$(souporcell --version)
    END_VERSIONS
    """
}
