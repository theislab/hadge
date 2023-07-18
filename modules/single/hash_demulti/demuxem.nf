#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process demuxem{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti/demuxem", mode:'copy'
    label 'small_mem'
    
    input:
        path raw_rna_matrix_dir, stageAs: "rna_data_${params.rna_matrix_demuxem}"
        path raw_hto_matrix_dir, stageAs: "hto_data_${params.hto_matrix_demuxem}"
        val threads
        each alpha
        each alpha_noise
        each tol
        each min_num_genes
        each min_num_umis
        each min_signal
        each random_state
        each generate_gender_plot
        val objectOutDemuxem
        each filter_demuxem
    output:
        path "demuxem_${task.index}"
        
    script:
        def generateGenderPlot = generate_gender_plot != 'None' ? " --generateGenderPlot ${generate_gender_plot}" : ''
        """
        mkdir demuxem_${task.index}
        demuxem.py --rna_matrix_dir rna_data_${params.rna_matrix_demuxem} --hto_matrix_dir hto_data_${params.hto_matrix_demuxem} \
            --randomState $random_state --min_signal $min_signal --tol $tol \
            --min_num_genes $min_num_genes --min_num_umis $min_num_umis --alpha $alpha --alpha_noise $alpha_noise \
            --n_threads $threads $generateGenderPlot --objectOutDemuxem $objectOutDemuxem --outputdir demuxem_${task.index} \
            --filter_demuxem $filter_demuxem
        
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

workflow demuxem_hashing{
    take:
        hto_matrix
        rna_matrix
    main:
        threads = params.threads_demuxem
        alpha = split_input(params.alpha_demuxem)
        alpha_noise = split_input(params.alpha_noise)
        min_num_genes = split_input(params.min_num_genes)
        min_num_umis = split_input(params.min_num_umis)
        min_signal = split_input(params.min_signal)
        tol = split_input(params.tol)
        random_state = split_input(params.random_state)
        generate_gender_plot = split_input(params.generate_gender_plot)
        objectOutDemuxem = params.objectOutDemuxem
        filter_demuxem = split_input(params.filter_demuxem)

        demuxem(rna_matrix, hto_matrix, threads, alpha, alpha_noise, tol, min_num_genes, min_num_umis, 
                min_signal, random_state, generate_gender_plot, objectOutDemuxem,filter_demuxem)
  
  emit:
        demuxem.out.collect()
}