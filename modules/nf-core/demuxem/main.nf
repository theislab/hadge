process DEMUXEM {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0d/0d3f96aaa8437bfa1654570e1d2b84749f1ac14d68f97978acc19b3757af7f55/data'
        : 'community.wave.seqera.io/library/demuxem:0.1.7.post1--5ac55376ad7cb80e'}"

    input:
    tuple val(meta), path(input_raw_gene_bc_matrices_h5), path(input_hto_csv_file)
    val generate_gender_plot
    val genome
    val generate_diagnostic_plots

    output:
    tuple val(meta), path("*_demux.zarr.zip"), emit: zarr
    tuple val(meta), path("*.out.demuxEM.zarr.zip"), emit: out_zarr
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def generateGenderPlot = generate_gender_plot ? "--generate-gender-plot ${generate_gender_plot}" : ""
    def genome_file = genome ? "--genome ${genome}" : ""
    def diagnostic_plots = generate_diagnostic_plots ? "--generate-diagnostic-plots ${generate_diagnostic_plots}" : ""
    """
    demuxEM ${input_raw_gene_bc_matrices_h5} ${input_hto_csv_file} ${prefix} \\
        -p $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        echo \$(demuxEM --version  2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.out.demuxEM.zarr.zip
    touch ${prefix}_demux.zarr.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        echo \$(demuxEM --version  2>&1)
    END_VERSIONS
    """
}
