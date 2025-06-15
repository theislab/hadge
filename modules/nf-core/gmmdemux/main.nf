process GMMDEMUX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gmm-demux:0.2.2.3--pyh7cba7a3_0'
        : 'biocontainers/gmm-demux:0.2.2.3--pyh7cba7a3_0'}"

    input:
    tuple val(meta), path(hto_matrix), val(hto_names), val(estimated_cells)
    val full_report
    val summary_report
    path classification_report
    path cell_list

    output:
    tuple val(meta), path("barcodes.tsv.gz"), emit: barcodes
    tuple val(meta), path("matrix.mtx.gz"), emit: matrix
    tuple val(meta), path("features.tsv.gz"), emit: features
    tuple val(meta), path("GMM_*.csv"), emit: classification_report
    tuple val(meta), path("GMM_*.config"), emit: config_report
    tuple val(meta), path("summary_report_*.txt"), emit: summary_report, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Docs: https://gmm-demux.readthedocs.io/en/latest/usage.html#command-line-tools
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def skip_arg = classification_report ? "--skip ${classification_report}" : ""
    def examine_arg = cell_list ? "--examine ${cell_list}" : ""
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    def VERSION = '0.2.2.3'
    def report_arg = full_report ? "-f ." : "-s ."
    def estimated_cells_arg = estimated_cells ? "--summary ${estimated_cells}" : ""
    def summary_report_arg = summary_report ? "-r ${prefix}_summary_report.txt" : ""
    """
    if [[ ${summary_report} == true ]]; then
        cat /dev/null > ${prefix}_summary_report.txt
    fi

    GMM-demux ${hto_matrix} ${hto_names} \\
        ${report_arg} \\
        ${summary_report_arg} \\
        ${estimated_cells_arg} \\
        ${skip_arg} \\
        ${examine_arg} \\
        -o . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GMM-Demux: ${VERSION}
    END_VERSIONS
    """

    stub:
    def VERSION = '0.2.2.3'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > barcodes.tsv.gz
    echo "" | gzip > features.tsv.gz
    echo "" | gzip > matrix.mtx.gz
    touch GMM_full.config
    touch GMM_full.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GMM-Demux: ${VERSION}
    END_VERSIONS
    """
}
