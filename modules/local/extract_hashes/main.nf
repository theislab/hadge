process EXTRACT_HASHES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/74/749b3cf99e0a33f46d2b49ab60d6e408ce467476c0ec86da1775e73ac4b5ba7b/data'
        : 'community.wave.seqera.io/library/coreutils_gawk_gzip:3d2dde6df78e314a'}"

    input:
    tuple val(meta), path(hto_dir)

    output:
    tuple val(meta), path("*_hashes.txt"), emit: hashes

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix         = task.ext.prefix         ?: "${meta.id}"
    """
    gunzip -c ${hto_dir}/features.tsv.gz | awk '{print \$2}' | paste -sd, - > ${prefix}_hashes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(gzip --version 2>&1 | head -n1 | sed 's/^.*gzip //; s/ .*\$//')
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
        coreutils: \$( cat --version | sed '1!d; s/.* //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_hashes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(gzip --version 2>&1 | head -n1 | sed 's/^.*gzip //; s/ .*\$//')
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
        coreutils: \$( cat --version | sed '1!d; s/.* //' )
    END_VERSIONS
    """
}
