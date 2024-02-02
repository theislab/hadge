#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process scSplit{
    publishDir "$params.outdir/$sampleId/$params.mode/gene_demulti/scSplit", mode: 'copy'
    label 'big_mem'
    tag "${sampleId}"
    conda "$projectDir/conda/scsplit.yml"

    input:
        tuple val(sampleId), path(bam), path(bai), path(vcf), path(barcode), val(num), val(vcf_known)// file function not available, use val instead
        val tag
        val com
        val ref
        val alt
        val sub
        val ems
        val dbl
        val sample_geno
        val scsplit_out
       
       
    output:
        path "scsplit_${sampleId}"
   
    script:
        
        def common_data = com != 'None' ? "--com $com" : ''
        def common_data_name = com != 'None' ? com : 'Common variants are not given.'
        def vcf_known_data = vcf_known != 'None' ? "--vcf ${vcf_known}" : ''
        def vcf_known_data_name = vcf_known != 'None' ? vcf_known : 'Known variants are not given.'
        def sub_yesno = num == 0 ? "$sub": 'no_sub'
        def vcf_data = "-v $vcf"
        def bam_data = "-i $bam"
        def barcode_data = "-b $barcode"
        def tag_data = "--tag $tag"
        def num_data = "-n $num"
        def sub_data = num == 0 ? "--sub $sub": ''
        def ems_data = "--ems $ems"
        def dbl_data = dbl != 'None' ? "--dbl $dbl" : ''
        def out = "scsplit_${sampleId}/${scsplit_out}"
        
        """
        git clone https://github.com/jon-xu/scSplit
        mkdir scsplit_${sampleId}
        mkdir $out
        touch scsplit_${sampleId}/params.csv
        echo -e "Argument,Value \n vcf,$vcf \n bam,$bam \n barcode,$barcode \n common_data,${common_data_name} \n num,${num} \n sub,${sub_yesno} \n ems,${ems} \n dbl,${dbl} \n vcf_known_data,${vcf_known_data_name}" >> scsplit_${sampleId}/params.csv
        
        python scSplit/scSplit count ${vcf_data} ${bam_data} ${barcode_data} ${common_data} -r $ref -a $alt --out $out
        python scSplit/scSplit run -r $out/$ref -a $out/$alt ${num_data} ${sub_data} ${ems_data} ${dbl_data} ${vcf_known_data} --out $out
        
        if [[ "$sample_geno" != "False" ]]
        then
            python scSplit/scSplit genotype -r $out/$ref -a $out/$alt -p $out/scSplit_P_s_c.csv --out $out
        fi
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

workflow demultiplex_scSplit{
    take:
        input_list
    main:
        tag_scsplit = params.tag_group
        com_scsplit = params.common_variants_scSplit
        ref_scsplit = params.refscSplit
        alt_scsplit = params.altscSplit
        sub_scsplit = params.subscSplit
        ems_scsplit = params.emsscSplit
        dbl_scsplit = params.dblscSplit
        sample_geno = params.sample_geno
        scsplit_out = params.scsplit_out
        scSplit(input_list, tag_scsplit, com_scsplit, ref_scsplit, alt_scsplit, sub_scsplit, ems_scsplit, dbl_scsplit, 
            sample_geno, scsplit_out)
    emit:
        scSplit.out.collect()
    
}

