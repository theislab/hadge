#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process demuxmix{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti/demuxmix", mode:'copy'
    input:
        //shares pre-process Seurat
        each seurat_object
        //Same assay as Seurat
        each rna_available
        each assay
        each model
        each alpha_demuxmix
        each beta_demuxmix
        each tol_demuxmix
        each maxIter_demuxmix
        each k_hto
        each k_rna
        
    output:
        path "demuxmix_${task.index}"
        
    script:
        
        """
        mkdir demuxmix_${task.index}

        Rscript $baseDir/bin/demuxmix.R --seuratObject $seurat_object --rna_available $rna_available --assay $assay --model $model --alpha_demuxmix $alpha_demuxmix \
            --beta_demuxmix $beta_demuxmix --tol_demuxmix $tol_demuxmix --maxIter_demuxmix $maxIter_demuxmix \
            --k_hto $k_hto  --k_rna $k_rna --outputdir demuxmix_${task.index}
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
        seurat_object
  main:
        rna_available = split_input(params.rna_available)
        assay = split_input(params.assay)
        model = split_input(params.model)
        alpha_demuxmix =  split_input(params.alpha_demuxmix)
        beta_demuxmix = split_input(params.beta_demuxmix)
        tol_demuxmix = split_input(params.tol_demuxmix)
        maxIter_demuxmix = split_input(params.maxIter_demuxmix)
        k_hto = split_input(params.k_hto)
        k_rna = split_input(params.k_rna) 
        

        demuxmix(seurat_object,rna_available, assay,model, alpha_demuxmix, beta_demuxmix, tol_demuxmix, maxIter_demuxmix, k_hto, k_rna )
  
  emit:
        demuxmix.out.collect()
}

workflow{
    demuxmix_hashing()
}
