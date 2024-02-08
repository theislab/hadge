process create_single_chanel_input{
    
    tag "${sampleId}"
    conda "$projectDir/conda/scsplit.yml"

    input:
            val sample_name
            val hto_matrix_raw
            val hto_matrix_filtered
            val rna_matrix_raw
            val rna_matrix_filtered
            val bam
            val bai
            val barcodes
            val fasta
            val fasta_index
            val nsample
            val cell_data
            val vcf_mixed
            val vcf_donor
            val vireo_parent_dir
            val demultiplexing_result
    output:
        path('hadge_single_input.csv'), emit: input_channel

    script:
    """
        echo "sampleId,rna_matrix_raw,rna_matrix_filtered,hto_matrix_raw,hto_matrix_filtered,bam,bam_index,barcodes,nsample,cell_data,vcf_mixed,vcf_donor,vireo_parent_dir,demultiplexing_result" > hadge_single_input.csv
        echo "${sample_name},${rna_matrix_raw},${rna_matrix_filtered},${hto_matrix_raw},${hto_matrix_filtered},${bam},${bai},${barcodes},${nsample},${cell_data},${vcf_mixed},${vcf_donor},${vireo_parent_dir},${demultiplexing_result}" >> hadge_single_input.csv
    """
   
}


