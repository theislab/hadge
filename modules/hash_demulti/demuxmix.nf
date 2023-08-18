#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process demuxmix{
    publishDir "$projectDir/$params.outdir/${seurat_object.name.tokenize( '_' )[1]}/$params.mode/hash_demulti/demuxmix", mode:'copy'
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
        each gene_col
        
    output:
        path "demuxmix_${sampleId}"
        
    script:
        def sampleId = seurat_object.name.tokenize( '_' )[1]

        """
        mkdir demuxmix_${sampleId}

        demuxmix.R --fileUmi rna_data --fileHto hto_data --rna_available $rna_available --assay $assay --ndelim $ndelim --model $model --alpha_demuxmix $alpha_demuxmix \
            --beta_demuxmix $beta_demuxmix --tol_demuxmix $tol_demuxmix --maxIter_demuxmix $maxIter_demuxmix --correctTails $correctTails\
            --k_hto $k_hto  --k_rna $k_rna --outputdir demuxmix_${sampleId} --assignmentOutDemuxmix $assignmentOutDemuxmix --gene_col $gene_col
        """

}

workflow demuxmix_hashing{
  take:
      input_list
  main:
        assay = params.assay
        ndelim = params.ndelim
        model = params.model
        alpha_demuxmix =  params.alpha_demuxmix
        beta_demuxmix = params.beta_demuxmix
        tol_demuxmix = params.tol_demuxmix
        maxIter_demuxmix = params.maxIter_demuxmix
        k_hto = params.k_hto
        k_rna = params.k_rna
        correctTails = params.correctTails
        assignmentOutDemuxmix = params.assignmentOutDemuxmix 
        gene_col = params.gene_col
        

        demuxmix(input_list, assay,ndelim,model,alpha_demuxmix, beta_demuxmix, tol_demuxmix, maxIter_demuxmix, k_hto, k_rna,correctTails,assignmentOutDemuxmix, gene_col)
  
  emit:
        demuxmix.out.collect()
}

workflow{
    demuxmix_hashing()
}
