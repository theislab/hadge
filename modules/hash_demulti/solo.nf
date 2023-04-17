#!/usr/bin/env nextflow

nextflow.enable.dsl=2
process solo{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti/solo", mode:'copy'
    input:
        each rna_matrix_dir
        each max_epochs
        each lr
        each train_size
        each validation_size
        each batch_size
        each early_stopping
        each early_stopping_patience
        each early_stopping_min_delta
        each soft
        each include_simulated_doublets
        each assignmentOutSolo
        
    output:
        path "solo_${task.index}"
        
    script:

        """
        mkdir solo_${task.index}
        solo_demul.py --rna_matrix_dir $rna_matrix_dir --max_epochs $max_epochs --lr $lr --train_size $train_size \
                      --validation_size $validation_size --batch_size $batch_size --early_stopping $early_stopping \
                      --early_stopping_patience $early_stopping_patience --early_stopping_min_delta $early_stopping_min_delta \
                      --soft $soft --include_simulated_doublets $include_simulated_doublets \
                      --assignmentOutSolo $assignmentOutSolo --outputdir solo_${task.index}
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

workflow solo_hashing{
    main:
        rna_matrix_dir = split_input(params.rna_matrix_solo)
        max_epochs = split_input(params.max_epochs)
        lr = split_input(params.lr)
        train_size = split_input(params.train_size)
        validation_size = split_input(params.validation_size)
        batch_size = split_input(params.batch_size)
        early_stopping = split_input(params.early_stopping)
        early_stopping_patience = split_input(params.early_stopping_patience)
        early_stopping_min_delta = split_input(params.early_stopping_min_delta)
        soft = split_input(params.soft)
        include_simulated_doublets = split_input(params.include_simulated_doublets)
        assignmentOutSolo = split_input(params.assignmentOutSolo)
        solo(rna_matrix_dir, max_epochs, lr, train_size, validation_size, batch_size, early_stopping, early_stopping_patience, early_stopping_min_delta, soft, include_simulated_doublets, assignmentOutSolo)
        
    emit:
        solo.out.collect()
}



workflow{
    solo_hashing()
}
