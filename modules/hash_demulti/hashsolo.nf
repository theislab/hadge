#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
process hash_solo {
    publishDir "$projectDir/$params.outdir/$sampleId/$params.mode/hash_demulti/hashsolo", mode:'copy'
    label 'small_mem'

    conda "$projectDir/conda/hashsolo_py.yml"

    input:
        tuple val(sampleId), path(hto_data, stageAs: "hto_data_${params.hto_matrix_hashsolo}"), path(rna_data, stageAs: "rna_data_${params.rna_matrix_hashsolo}")
        val priors_negative
        val priors_singlet
        val priors_doublet
        val pre_existing_clusters
        val use_rna_data
        val number_of_noise_barcodes
        val assignmentOutHashSolo
        val plotOutHashSolo

    output:
        path "hashsolo_${sampleId}"

    script:
        def noise_barcodes = number_of_noise_barcodes != 'None' ? "--number_of_noise_barcodes $number_of_noise_barcodes" : ''
        def existing_clusters = pre_existing_clusters != 'None' ? "--pre_existing_clusters $pre_existing_clusters" : ''
        def clustering_data = use_rna_data != 'False' ? "--clustering_data rna_data_$params.rna_matrix_hashsolo" : ''
        """
        mkdir hashsolo_${sampleId}
        hashsolo.py --hto_data hto_data_${params.hto_matrix_hashsolo} --priors $priors_negative $priors_singlet $priors_doublet \
                    $existing_clusters $clustering_data $noise_barcodes \
                    --assignmentOutHashSolo $assignmentOutHashSolo \
                    --plotOutHashSolo $plotOutHashSolo --outputdir hashsolo_${sampleId}
        """
}

workflow hash_solo_hashing {
    take:
        input_list
    main:
        use_rna_data = params.use_rna_data
        priors_negative = params.priors_negative
        priors_singlet = params.priors_singlet
        priors_doublet = params.priors_doublet
        pre_existing_clusters = params.pre_existing_clusters
        number_of_noise_barcodes = params.number_of_noise_barcodes
        assignmentOutHashSolo = params.assignmentOutHashSolo
        plotOutHashSolo = params.plotOutHashSolo

        hash_solo(input_list, priors_negative, priors_singlet, priors_doublet, pre_existing_clusters, use_rna_data, number_of_noise_barcodes, assignmentOutHashSolo, plotOutHashSolo)
    emit:
        hash_solo.out.collect()
}
