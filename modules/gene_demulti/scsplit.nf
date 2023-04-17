#!/usr/bin/env nextflow
'''
scSplit_loc = Channel.fromPath(params.scSplit_loc)
def vcf_known_data_name = vcf_known != 'False' ? file(vcf_known).baseName : 'known_variants_not_given'

'''

nextflow.enable.dsl=2

process scSplit{
    publishDir "$projectDir/$params.outdir/$params.mode/gene_demulti/scSplit", mode: 'copy'
    
    input:
        each vcf
        each bam
        each bai
        each barcode
        each tag
        each com
        each ref
        each alt
        each num
        each sub
        each ems
        each dbl
        each vcf_known
        each sample_geno
        each scsplit_out
        each scsplit_loc
       
       
    output:
        path "scsplit_${task.index}"
   
    script:
        
        def common_data = com != 'False' ? "--com $com" : ''
        def common_data_name = com != 'False' ? com : 'common_variants_not_given'
        def vcf_known_data = vcf_known != 'False' ? "--vcf ${vcf_known}" : ''
        def vcf_known_data_name = vcf_known != 'False' ? vcf_known : 'known_variants_not_given'
        def sub_yesno = num == 0 ? "$sub": 'no_sub'
        
        def vcf_data = "-v $vcf"
        def bam_data = "-i $bam"
        def barcode_data = "-b $barcode"
        def tag_data = "--tag $tag"
        def num_data = "-n $num"
        def sub_data = num == 0 ? "--sub $sub": ''
        def ems_data = "--ems $ems"
        def dbl_data = dbl != 'False' ? "--dbl $dbl" : ''
        def out = "scsplit_${task.index}/${scsplit_out}"
        
        """
        mkdir scsplit_${task.index}
        mkdir $out
        touch scsplit_${task.index}/params.csv
        echo -e "Argument,Value \n vcf,$vcf \n bam,$bam \n barcode,$barcode \n common_data, ${common_data_name} \n num, ${num} \n sub, ${sub_yesno} \n ems, ${ems} \n dbl, ${dbl} \n vcf_known_data, ${vcf_known_data_name}" >> scsplit_${task.index}/params.csv
        
        python $scsplit_loc/scSplit count ${vcf_data} ${bam_data} ${barcode_data} ${common_data} -r $ref -a $alt --out $out
        python $scsplit_loc/scSplit run -r $out/$ref -a $out/$alt ${num_data} ${sub_data} ${ems_data} ${dbl_data} ${vcf_known_data} --out $out
        
        if [[ "$sample_geno" != "False" ]]
        then
            python $scsplit_loc/scSplit genotype -r $out/$ref -a $out/$alt -p $out/scSplit_P_s_c.csv --out $out
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
        bam_scsplit
        vcf_scsplit
        bai_scsplit
        bar_scsplit
    main:
        tag_scsplit = split_input(params.tag_group)
        com_scsplit = split_input(params.common_variants_scSplit)
        ref_scsplit = split_input(params.refscSplit)
        alt_scsplit = split_input(params.altscSplit)
        num_scsplit = split_input(params.nsample)
        sub_scsplit = split_input(params.subscSplit)
        ems_scsplit = split_input(params.emsscSplit)
        dbl_scsplit = split_input(params.dblscSplit)
        vcf_known_scsplit = split_input(params.vcf_donor)
        sample_geno = split_input(params.sample_geno)
        scsplit_out = split_input(params.scsplit_out)
        scsplit_loc = split_input(params.scsplit_loc)
        scSplit(vcf_scsplit, bam_scsplit,bai_scsplit, bar_scsplit, tag_scsplit, com_scsplit, ref_scsplit, alt_scsplit, num_scsplit, sub_scsplit, ems_scsplit, dbl_scsplit, vcf_known_scsplit,sample_geno, scsplit_out, scsplit_loc)
    emit:
        scSplit.out.collect()
    
}

workflow{
    demultiplex_scSplit(channel.fromPath(params.bam), channel.fromPath(params.vcf_mixed), channel.fromPath(params.bai))
}
