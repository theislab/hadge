#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process demuxmix{
    publishDir "$projectDir/$params.outdir/${seurat_object.name.tokenize( '_' )[1]}/$params.mode/hash_demulti/demuxmix", mode:'copy'
    label 'small_mem'
    
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
        each assignmentOutDemuxmix
        
    output:
        path "demuxmix_${sampleId}"
        
    script:
        def sampleId = seurat_object.name.tokenize( '_' )[1]

        """
        mkdir demuxmix_${sampleId}

        Rscript $baseDir/bin/demuxmix.R --seuratObject $seurat_object --rna_available $rna_available --assay $assay --model $model --alpha_demuxmix $alpha_demuxmix \
            --beta_demuxmix $beta_demuxmix --tol_demuxmix $tol_demuxmix --maxIter_demuxmix $maxIter_demuxmix \
            --k_hto $k_hto  --k_rna $k_rna --outputdir demuxmix_${sampleId} --assignmentOutDemuxmix $assignmentOutDemuxmix
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
        assignmentOutDemuxmix = split_input(params.assignmentOutDemuxmix) 
        

        demuxmix(seurat_object,rna_available, assay,model, alpha_demuxmix, beta_demuxmix, tol_demuxmix, maxIter_demuxmix, k_hto, k_rna,assignmentOutDemuxmix )
  
  emit:
        demuxmix.out.collect()
}

workflow{
    demuxmix_hashing()
}
