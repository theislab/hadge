process preprocess{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti/preprocess", mode:'copy'
    label 'seurat'
    input:
        each rdsObject
        each umi_matrix
        each hto_matrix
        each ndelim
        each selection_method
        each number_features
        each assay
        each margin
        each normalisation_method
        each preprocess_out
    output:
        path "preprocess_${task.index}"

    script:
        def rds = rdsObject != "FALSE" ? "--rdsObject" : ""
    """
        mkdir preprocess_${task.index}
        pre_processing.R $rds --fileUmi $umi_matrix --fileHto $hto_matrix --ndelim $ndelim \
                                  --selectMethod $selection_method --numberFeatures $number_features --assay $assay \
                                  --margin $margin --normalisationMethod $normalisation_method --OutputFile $preprocess_out \
                                  --outputdir preprocess_${task.index}
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

workflow preprocessing_hashing{
    main:
        rdsObject = split_input(params.rdsObject)
        umi_matrix = split_input(params.umi_matrix_preprocess)
        hto_matrix =  split_input(params.hto_matrix_preprocess)
        sel_method = split_input(params.sel_method)
        ndelim = split_input(params.ndelim)
        n_features = split_input(params.n_features)
        assay = split_input(params.assay)
        margin = split_input(params.margin)
        norm_method = split_input(params.norm_method)
        out_file = split_input(params.preprocessOut)
        preprocess(rdsObject, umi_matrix, hto_matrix, ndelim, sel_method, n_features, assay, margin, norm_method, out_file)
    emit:
        preprocess.out.collect()
}
  
workflow{
    preprocessing_hashing()
}
