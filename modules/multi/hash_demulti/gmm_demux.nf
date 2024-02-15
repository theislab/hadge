#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process gmm_demux{
    publishDir "$params.outdir/$sampleId/$params.mode/hash_demulti/gmm_demux", mode:'copy'
    label 'small_mem'

    container 'mari3ga/gmm_demux:v1'
    
    input:
        tuple val(sampleId), path(filtered_hto_matrix_dir), val(hto_name_gmm)
        //HTO names as string separated by commas
        //val hto_name_gmm
        //mode 2
        //need estimate number of cells in the single cell assay
        //obligatory
        val summary
        //need to be combined with summary to get a report as file
        val report_gmm 
        //mode 4
        // write csv or tsv - type of input
        val mode_GMM
        //case 5
        val extract 
        //float between 0 and 1
        val threshold_gmm
        val ambiguous
        
        
    
    output:
        path "gmm_demux_${sampleId}"
        
    script:
        def extract_droplets = extract != 'None' ? " -x ${extract}" : ''
        def ambiguous_droplets = extract != 'None' ? " --ambiguous ${ambiguous}" : ''

        if(mode_GMM=="csv"){
            """
            mkdir gmm_demux_${sampleId}
            touch gmm_demux_${sampleId}_$report_gmm
            
            GMM-demux -c $filtered_hto_matrix_dir $hto_name_gmm -u $summary --report gmm_demux_${sampleId}_$report_gmm --full gmm_demux_${sampleId} $extract_droplets -t $threshold_gmm
            gmm_demux_params.py --path_hto $filtered_hto_matrix_dir --hto_name_gmm $hto_name_gmm --summary $summary --report gmm_demux_${sampleId}_$report_gmm   --mode $mode_GMM  $extract_droplets --threshold_gmm $threshold_gmm $ambiguous_droplets  --outputdir gmm_demux_${sampleId}
            
            """
        }else {
            """
            mkdir gmm_demux_${sampleId}
            touch gmm_demux_${sampleId}_$report_gmm

            print $filtered_hto_matrix_dir
            
            GMM-demux $filtered_hto_matrix_dir $hto_name_gmm -u $summary -r gmm_demux_${sampleId}_$report_gmm --full gmm_demux_${sampleId} -o gmm_demux_${sampleId} $extract_droplets -t $threshold_gmm
            gmm_demux_params.py --path_hto $filtered_hto_matrix_dir --hto_name_gmm $hto_name_gmm --summary $summary --report gmm_demux_${sampleId}_$report_gmm --mode $mode_GMM $extract_droplets  --threshold_gmm $threshold_gmm $ambiguous_droplets --outputdir gmm_demux_${sampleId}
            
            """
        }


}


workflow gmm_demux_hashing{
take: 
        input_list

  main:
        summary = params.summary
        report_gmm = params.report_gmm
        mode = params.mode_GMM
        extract = params.extract
        threshold_gmm = params.threshold_gmm
        ambiguous = params.ambiguous

        gmm_demux(input_list,summary,report_gmm,mode,extract,threshold_gmm,ambiguous)
  
  emit:
        gmm_demux.out.collect()
}


workflow{
    gmm_demux_hashing()

}
