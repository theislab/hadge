process GENE_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // container "<a-small-python-pandas-scanpy image>"

    input:
    tuple val(meta), path(barcodes),
         path(vireo_donor_ids),
         path(vireo_summary),
         path(demuxlet_result),
         path(freemuxlet_result),
         path(souporcell_tsv),
         path(rna_matrix)                    // pass [] if you don’t want h5ad

    output:
    tuple val(meta), path("*_genetic_summary_assignment.csv"), emit: assignment
    tuple val(meta), path("*_genetic_summary_classification.csv"), emit: classification
    path "versions.yml", emit: versions
    // (optionally) tuple val(meta), path("*_genetic_summary.h5ad"), emit: h5ad, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${projectDir}/modules/local/genetic_summary/templates/genetic_summary.py

    cat > versions.yml <<EOF
    "${task.process}":
        python: "\$(python3 -V 2>&1)"
        pandas: "\$(python3 -c 'import pandas,sys;print(pandas.__version__)' 2>/dev/null || echo N/A)"
    EOF
    """
}
