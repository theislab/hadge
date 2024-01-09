#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process demuxmix {
    publishDir "$projectDir/$params.outdir/$sampleId/$params.mode/hash_demulti/demuxmix", mode:'copy'
    label 'small_mem'

    conda "$projectDir/conda/demuxmix.yml"

    input:
        tuple val(sampleId), path(raw_hto_matrix_dir, stageAs: "hto_data_${params.hto_matrix_demuxmix}"),
                             path(raw_rna_matrix_dir, stageAs: "rna_data_${params.rna_matrix_demuxmix}")
        val hto_raw_or_filtered
        val rna_raw_or_filtered
        val rna_available
        val assay
        val ndelim
        val model
        val alpha_demuxmix
        val beta_demuxmix
        val tol_demuxmix
        val maxIter_demuxmix
        val k_hto
        val k_rna
        val correctTails
        val assignmentOutDemuxmix
        val gene_col

    output:
        path "demuxmix_${sampleId}"

    script:
        """
        mkdir demuxmix_${sampleId}

        demuxmix.R --fileUmi rna_data_${params.rna_matrix_demuxmix} --fileHto hto_data_${params.hto_matrix_demuxmix} --rna_available $rna_available --assay $assay \
                    --ndelim $ndelim --model $model --alpha_demuxmix $alpha_demuxmix --beta_demuxmix $beta_demuxmix --tol_demuxmix $tol_demuxmix \
                    --maxIter_demuxmix $maxIter_demuxmix --correctTails $correctTails --k_hto $k_hto --k_rna $k_rna --outputdir demuxmix_${sampleId} \
                    --assignmentOutDemuxmix $assignmentOutDemuxmix --gene_col $gene_col
        """
}

workflow demuxmix_hashing {
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

        demuxmix(input_list, assay, ndelim, model, alpha_demuxmix, beta_demuxmix, tol_demuxmix, maxIter_demuxmix, k_hto, k_rna, correctTails, assignmentOutDemuxmix, gene_col)

  emit:
        demuxmix.out.collect()
}

workflow {
    demuxmix_hashing()
}
