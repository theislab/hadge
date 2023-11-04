#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process demuxem{
    publishDir "$projectDir/$params.outdir/$sampleId/$params.mode/hash_demulti/demuxem", mode:'copy'
    label 'small_mem'

    conda "python=3.9 bioconda::pegasuspy pandas<2.0.0 demuxEM"
    
    input:
        tuple val(sampleId), path(raw_hto_matrix_dir, stageAs: "hto_data_${params.hto_matrix_demuxem}"), 
                             path(raw_rna_matrix_dir, stageAs: "rna_data_${params.rna_matrix_demuxem}")
        val threads
        val alpha
        val alpha_noise
        val tol
        val min_num_genes
        val min_num_umis
        val min_signal
        val random_state
        val generate_gender_plot
        val objectOutDemuxem
        val filter_demuxem
    output:
        path "demuxem_${sampleId}"
        
    script:
        def generateGenderPlot = generate_gender_plot != "None" ? " --generateGenderPlot ${generate_gender_plot}" : ''
        """
        mkdir demuxem_${sampleId}
        demuxem.py --rna_matrix_dir rna_data_${params.rna_matrix_demuxem} --hto_matrix_dir hto_data_${params.hto_matrix_demuxem} \
            --randomState $random_state --min_signal $min_signal --tol $tol \
            --min_num_genes $min_num_genes --min_num_umis $min_num_umis --alpha $alpha --alpha_noise $alpha_noise \
            --n_threads $threads $generateGenderPlot --objectOutDemuxem $objectOutDemuxem --outputdir demuxem_${sampleId} --filter_demuxem $filter_demuxem
        
        """

}

workflow demuxem_hashing{
    take:
        input_list
    main:
        threads = params.threads_demuxem
        alpha = params.alpha_demuxem
        alpha_noise = params.alpha_noise
        min_num_genes = params.min_num_genes
        min_num_umis = params.min_num_umis
        min_signal = params.min_signal
        tol = params.tol
        random_state = params.random_state
        generate_gender_plot = params.generate_gender_plot
        objectOutDemuxem = params.objectOutDemuxem
        filter_demuxem = params.filter_demuxem

        demuxem(input_list, threads, alpha, alpha_noise, tol, min_num_genes, min_num_umis, 
                min_signal, random_state, generate_gender_plot, objectOutDemuxem, filter_demuxem)
  
  emit:
        demuxem.out.collect()
}