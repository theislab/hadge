#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process gmm_demux{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti/gmm_demux", mode:'copy'
    label 'small_mem'
    input:
        path path_hto
        //HTO names as string separated by commas
        each hto_name_gmm
        //5 cases are available for the tool it all depends on the vars given
        //mode 2
        //need estimate number of cells in the single cell assay
        //obligatory
        each summary
        //need to be combined with summary to get a report as file
        each report_gmm 
        //mode 4
        // write csv or tsv - type of input
        each mode_GMM
        //case 5
        each extract 
        //float between 0 and 1
        each threshold_gmm
        each ambiguous
        
        
    
    output:
        path "gmm_demux_${task.index}"
        
    script:
        def extract_droplets = extract != 'None' ? " -x ${extract}" : ''
        def ambiguous_droplets = extract != 'None' ? " --ambiguous ${ambiguous}" : ''

        if(mode_GMM=="csv"){
            """
            mkdir gmm_demux_${task.index}
            touch gmm_demux_${task.index}_$report_gmm
            
            GMM-demux -c $path_hto $hto_name_gmm -u $summary --report gmm_demux_${task.index}_$report_gmm --full gmm_demux_${task.index} $extract_droplets -t $threshold_gmm
            gmm_demux_params.py --path_hto $path_hto --hto_name_gmm $hto_name_gmm --summary $summary --report gmm_demux_${task.index}_$report_gmm   --mode $mode_GMM  $extract_droplets --threshold_gmm $threshold_gmm $ambiguous_droplets  --outputdir gmm_demux_${task.index}
            
            """
        }else {
            """
            mkdir gmm_demux_${task.index}
            touch gmm_demux_${task.index}_$report_gmm
            
            GMM-demux $path_hto $hto_name_gmm -u $summary -r gmm_demux_${task.index}_$report_gmm --full gmm_demux_${task.index} -o gmm_demux_${task.index} $extract_droplets -t $threshold_gmm
            gmm_demux_params.py --path_hto $path_hto --hto_name_gmm $hto_name_gmm --summary $summary --report gmm_demux_${task.index}_$report_gmm --mode $mode_GMM $extract_droplets  --threshold_gmm $threshold_gmm $ambiguous_droplets --outputdir gmm_demux_${task.index}
            
            """
        }


}

def split_input(input){
    if (input =~ /;/ ){
        Channel.from(input).map{ return it.tokenize(';')}.flatten()
    }
    else{
        Channel.from(input)
    }
}

workflow gmm_demux_hashing{
    take:   
        hto_matrix
    main:
            hto_name_gmm = split_input(params.hto_name_gmm)
            summary = split_input(params.summary)
            report_gmm = split_input(params.report_gmm)
            mode = split_input(params.mode_GMM)
            extract = split_input(params.extract)
            threshold_gmm = split_input(params.threshold_gmm)
            ambiguous = split_input(params.ambiguous)

            gmm_demux(hto_matrix,hto_name_gmm,summary,report_gmm,mode,extract,threshold_gmm,ambiguous)
  
  emit:
        gmm_demux.out.collect()
}


workflow{
    gmm_demux_hashing()

}
