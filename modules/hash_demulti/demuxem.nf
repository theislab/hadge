#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process demuxem{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti/demuxem", mode:'copy'
    input:
        each raw_rna_matrix_dir
        each raw_hto_matrix_dir
        each threads
        each alpha
        each alpha_noise
        each tol
        each min_num_genes
        each min_num_umis
        each min_signal
        each random_state
        each generate_gender_plot
        each objectOutDemuxem
    output:
        path "demuxem_${task.index}"
        
    script:
        def generateGenderPlot = generate_gender_plot != 'None' ? " --generateGenderPlot ${generate_gender_plot}" : ''
        """
        mkdir demuxem_${task.index}
        demuxem.py --rna_matrix_dir $raw_rna_matrix_dir --hto_matrix_dir $raw_hto_matrix_dir --randomState $random_state --min_signal $min_signal --min_num_genes $min_num_genes --min_num_umis $min_num_umis --alpha $alpha --alpha_noise $alpha_noise --tol $tol --n_threads $threads $generateGenderPlot --objectOutDemuxem $objectOutDemuxem --outputdir demuxem_${task.index}
        
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
  main:
        raw_rna_matrix_dir = split_input(params.rna_matrix_demuxem)
        raw_hto_matrix_dir = split_input(params.hto_matrix_demuxem)
        threads = split_input(params.threads_demuxem)
        alpha = split_input(params.alpha_demuxem)
        alpha_noise = split_input(params.alpha_noise)
        min_num_genes = split_input(params.min_num_genes)
        min_num_umis = split_input(params.min_num_umis)
        min_signal = split_input(params.min_signal)
        tol = split_input(params.tol)
        random_state = split_input(params.random_state)
        generate_gender_plot = split_input(params.generate_gender_plot)
        objectOutDemuxem = split_input(params.objectOutDemuxem)

        demuxem(raw_rna_matrix_dir, raw_hto_matrix_dir, threads, alpha, alpha_noise, tol, min_num_genes, min_num_umis, min_signal, random_state, generate_gender_plot, objectOutDemuxem)
  
  emit:
        demuxem.out.collect()
}


workflow{
    demuxem_hashing()


}
