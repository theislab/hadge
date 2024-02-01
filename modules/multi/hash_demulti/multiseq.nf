#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process multi_seq{
    publishDir "$params.outdir/${seurat_object.name.tokenize( '_' )[1]}/$params.mode/hash_demulti/multiseq", mode: 'copy'
    label 'small_mem'
    
    conda "conda-forge::r-seurat conda-forge::r-argparse"

    input:
        each seurat_object
        val quantile
        val autoThresh
        val maxiter
        val qrangeFrom
        val qrangeTo
        val qrangeBy
        val verbose
        val assay
        val objectOutMulti
        val assignmentOutMulti

    output:
        path "multiseq_${seurat_object.name.tokenize( '_' )[1]}"

    script:
        def sampleId = seurat_object.name.tokenize( '_' )[1]
        def autoThr = autoThresh != 'False' ? " --autoThresh" : ''
        def verb = verbose != 'False' ? " --verbose" : ''
                
        """
        mkdir multiseq_${sampleId}
        MultiSeq.R --seuratObjectPath $seurat_object  --assay $assay --quantile $quantile $autoThr \
                   --maxiter $maxiter --qrangeFrom $qrangeFrom --qrangeTo $qrangeTo --qrangeBy $qrangeBy \
                   $verb --objectOutMulti $objectOutMulti --assignmentOutMulti $assignmentOutMulti --outputdir multiseq_${sampleId}
        """
}

workflow multiseq_hashing{
   take:
        seurat_object
   main:
        quantile = params.quantile_multi
        autoThresh = params.autoThresh
        maxIter = params.maxiter
        qrangeFrom = params.qrangeFrom
        qrangeTo = params.qrangeTo
        qrangeBy = params.qrangeBy
        verbose = params.verbose_multiseq
        assay = params.assay
        objectOutMulti = params.objectOutMulti
        assignmentOutMulti = params.assignmentOutMulti
        multi_seq(seurat_object, quantile, autoThresh, maxIter, qrangeFrom, qrangeTo, qrangeBy, verbose, assay, objectOutMulti, assignmentOutMulti)
   emit:
        multi_seq.out.collect()
}
