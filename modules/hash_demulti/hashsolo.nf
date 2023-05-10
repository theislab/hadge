#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process hash_solo{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti/hashsolo", mode:'copy'
    label 'small_mem'
    
    input:
        each hto_data
        each priors_negative
        each priors_singlet
        each priors_doublet
        each pre_existing_clusters
        each rna_data
        each number_of_noise_barcodes
        each assignmentOutHashSolo
        each plotOutHashSolo
    output:
        path "hashsolo_${task.index}"

    script:
        def noise_barcodes = number_of_noise_barcodes != 'None' ? "--number_of_noise_barcodes $number_of_noise_barcodes" : ''
        def existing_clusters = pre_existing_clusters != 'None' ? "--pre_existing_clusters $pre_existing_clusters" : ''
        def clustering_data = rna_data != 'None' ? "--clustering_data $rna_data" : ''
        """
        mkdir hashsolo_${task.index}
        hashsolo.py --hto_data $hto_data --priors $priors_negative $priors_singlet $priors_doublet \
                    $existing_clusters $clustering_data $noise_barcodes \
                    --assignmentOutHashSolo $assignmentOutHashSolo \
                    --plotOutHashSolo $plotOutHashSolo --outputdir hashsolo_${task.index}
        """

}

def split_input(input){
    if (input =~ /;/ ){
        Channel.from(input).map{ return it.tokenize(';')}.flatten()
    }
    else{
        Channel.from(input)
    }
}

workflow hash_solo_hashing{
    main:
        hto_data = split_input(params.hto_matrix_hashsolo)
        priors_negative = split_input(params.priors_negative)
        priors_singlet = split_input(params.priors_singlet)
        priors_doublet = split_input(params.priors_doublet)
        pre_existing_clusters = split_input(params.pre_existing_clusters)
        rna_data = split_input(params.rna_matrix_hashsolo)
        number_of_noise_barcodes = split_input(params.number_of_noise_barcodes)
        assignmentOutHashSolo = split_input(params.assignmentOutHashSolo)
        plotOutHashSolo = split_input(params.plotOutHashSolo)

        hash_solo(hto_data, priors_negative, priors_singlet, priors_doublet, pre_existing_clusters, rna_data, number_of_noise_barcodes, assignmentOutHashSolo, plotOutHashSolo)
    emit:
        hash_solo.out.collect()
}


workflow{
    hash_solo_hashing()

}
