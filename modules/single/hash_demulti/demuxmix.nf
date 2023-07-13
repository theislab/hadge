#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process demuxmix{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti/demuxmix", mode:'copy'
    label 'small_mem'
    input:
        path hto_matrix, stageAs: 'hto_data'
        path umi_matrix, stageAs: 'rna_data'
        val hto_raw_or_filtered
        val rna_raw_or_filtered
        each rna_available
        each assay
        val ndelim
        each model
        each alpha_demuxmix
        each beta_demuxmix
        each tol_demuxmix
        each maxIter_demuxmix
        each k_hto
        each k_rna
        each correctTails
        each assignmentOutDemuxmix
        
    output:
        path "demuxmix_${task.index}"
        
    script:
        
        """
        mkdir demuxmix_${task.index}
        demuxmix.R --fileUmi rna_data --fileHto hto_data --rna_available $rna_available --assay $assay --ndelim $ndelim --model $model --alpha_demuxmix $alpha_demuxmix \
            --beta_demuxmix $beta_demuxmix --tol_demuxmix $tol_demuxmix --maxIter_demuxmix $maxIter_demuxmix --correctTails $correctTails \
            --k_hto $k_hto  --k_rna $k_rna --outputdir demuxmix_${task.index} --assignmentOutDemuxmix $assignmentOutDemuxmix 
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

workflow demuxmix_hashing{
  take:
        hto_matrix
        rna_matrix
        hto_raw_or_filtered
        rna_raw_or_filtered
        rna_available
  main:
        assay = split_input(params.assay)
        ndelim = params.ndelim
        model = split_input(params.model)
        alpha_demuxmix =  split_input(params.alpha_demuxmix)
        beta_demuxmix = split_input(params.beta_demuxmix)
        tol_demuxmix = split_input(params.tol_demuxmix)
        maxIter_demuxmix = split_input(params.maxIter_demuxmix)
        k_hto = split_input(params.k_hto)
        k_rna = split_input(params.k_rna)
        correctTails = split_input(params.correctTails)
        assignmentOutDemuxmix = split_input(params.assignmentOutDemuxmix) 
        

        demuxmix(hto_matrix,rna_matrix,hto_raw_or_filtered,rna_raw_or_filtered,rna_available, assay,ndelim,model, alpha_demuxmix, beta_demuxmix, tol_demuxmix, maxIter_demuxmix, k_hto, k_rna,correctTails,assignmentOutDemuxmix )
  
  emit:
        demuxmix.out.collect()
}

workflow{
    demuxmix_hashing()
}
