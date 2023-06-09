process preprocess{
    publishDir "$projectDir/$params.outdir/$sampleId/$params.mode/hash_demulti/preprocess", mode:'copy'
    label 'seurat'
    label 'small_mem'
    
    input:
        tuple val(sampleId), path(hto_matrix, stageAs: 'hto_data'), path(umi_matrix, stageAs: 'rna_data')
        val hto_raw_or_filtered
        val rna_raw_or_filtered
        val ndelim
        val selection_method
        val number_features
        val assay
        val margin
        val normalisation_method
        each rna_available
        each demuxmix_preprocess
        val preprocess_out
    output:
        path "preprocess_${sampleId}_hto_${hto_raw_or_filtered}_rna_${rna_raw_or_filtered}"

    script:
    """
        mkdir preprocess_${sampleId}_hto_${hto_raw_or_filtered}_rna_${rna_raw_or_filtered}
        pre_processing.R --fileUmi rna_data --fileHto hto_data --ndelim $ndelim \
                         --selectMethod $selection_method --numberFeatures $number_features --assay $assay \
                         --margin $margin --normalisationMethod $normalisation_method --rna_available $rna_available --demuxmix $demuxmix_preprocess --OutputFile $preprocess_out \
                         --outputdir preprocess_${sampleId}_hto_${hto_raw_or_filtered}_rna_${rna_raw_or_filtered}
    """


}

workflow preprocessing_hashing{
    take:
        input_list
        hto_raw_or_filtered
        rna_raw_or_filtered
    main:
        sel_method = params.sel_method
        ndelim = params.ndelim
        n_features = params.n_features
        assay = params.assay
        margin = params.margin
        norm_method = params.norm_method
        rna_available = split_input(params.rna_available)
        demuxmix_preprocess_mode = split_input(params.demuxmix_preprocess_mode)
        out_file = params.preprocessOut
        preprocess(input_list, hto_raw_or_filtered, rna_raw_or_filtered, ndelim, sel_method, n_features, assay, margin, norm_method,rna_available,demuxmix_preprocess_mode,out_file)
    emit:
        preprocess.out
}
  