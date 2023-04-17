#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process demuxlet {
    publishDir "$projectDir/$params.outdir/$params.mode/gene_demulti/demuxlet", mode: 'copy'
    input:
        each sam
        each tag_group
        each tag_UMI
        each vcf_ref
        each sm
        each sm_list
        each sam_verbose
        each vcf_verbose
        each skip_umi

        each cap_BQ
        each min_BQ
        each min_MQ
        each min_TD
        each excl_flag

        each group_list
        each min_total
        each min_uniq
        each min_umi
        each min_snp

        each plp
        each vcf_donor
        each field
        each geno_error_offset
        each geno_error_coeff
        each r2_info
        each min_mac
        each min_callrate
        
        each alpha
        each doublet_prior
        each demuxlet_out
        

    output:
        path "demuxlet_${task.index}"

    script:
        def samfile = plp == 'False' ? "--sam $sam" : ''
        def samfile_name = plp == 'False' ? sam.baseName : 'sam_file_not_used'
        def taggroup = tag_group != 'False' ? "--tag-group ${tag_group}" : ''
        def tagUMI = tag_UMI != 'False' ? "--tag-UMI ${tag_UMI}" : ''
        def vcfref = plp == 'True' ? "--vcf ${vcf_ref}" : ""
        def vcfref_name = plp == 'True' ? file(vcf_ref).baseName : "vcf_ref_file_not_used"
        def smlist = sm != 'False' ? "--sm $sm" : ''
        def sm_list_file = sm_list != 'False' ? "--sm-list ${sm_list}" : ''
        def sm_list_file_name = sm_list != 'False' ? file(sm_list).baseName : "sm_list_file_not_given"
        def samverbose = "--sam-verbose ${sam_verbose}"
        def vcfverbose = "--vcf-verbose ${vcf_verbose}"
        def skipumi = skip_umi != "False" ? "--skip-umi" : ""
        def capBQ = "--cap-BQ ${cap_BQ}"
        def minBQ = "--min-BQ ${min_BQ}"
        def minMQ = "--min-MQ ${min_MQ}"
        def minTD = "--min-TD ${min_TD}"
        def exclflag = "--excl-flag ${excl_flag}"
        def grouplist = group_list != 'False' ? "--group-list ${group_list}" : ''
        def grouplist_name = group_list != 'False' ? file(group_list).baseName : 'group_list_not_given'
        def mintotal = "--min-total ${min_total}"
        def minumi = "--min-umi ${min_umi}"
        def minuniq = "--min-uniq ${min_uniq}"
        def minsnp = "--min-snp ${min_snp}"
        def plpfile = (plp != 'True' && plp != 'False') ? "--plp $plp" : ''
        def plpfile_name = plp == 'True' ? "perform_plp" : ( plp != 'False' ? file(plp).baseName : "plp_file_not_used")
        def vcfdonor = "--vcf ${vcf_donor}"
        def vcfdonor_name = file(vcf_donor).baseName
        def fieldinfo = "--field $field"
        def genoerror_off = "--geno-error-offset ${geno_error_offset}"
        def genoerror_cof = "--geno-error-coeff ${geno_error_coeff}"
        def r2info = "--r2-info ${r2_info}"
        def minmac = "--min-mac ${min_mac}"
        def mincallrate = "--min-callrate ${min_callrate}"
        def alpha_value = alpha.replaceAll(/,/, " --alpha ")
        def doubletprior = "--doublet-prior ${doublet_prior}"
      
        """
        mkdir demuxlet_${task.index}
        mkdir demuxlet_${task.index}/plp
        touch demuxlet_${task.index}/params.csv
        barcode_num=\$(wc -l < "${group_list}")
        echo -e "Argument,Value \n samfile,${samfile_name} \n tag_group,${tag_group} \n tag_UMI, ${tag_UMI} \n vcf_ref, ${vcfref_name} \n sm, ${sm} \n sm_list_file, ${sm_list_file_name} \n sam_verbose, ${sam_verbose} \n vcf_verbose, ${vcf_verbose} \n skip_umi, ${skip_umi} \n cap_BQ, ${cap_BQ} \n min_BQ, ${min_BQ} \n min_MQ, ${min_MQ} \n min_TD, ${min_TD} \n excl_flag, ${excl_flag} \n grouplist, ${grouplist_name}_\${barcode_num} \n min_total, ${min_total} \n min_uniq, ${min_uniq} \n min_umi, ${min_umi} \n min_snp, ${min_snp} \n plpfile, ${plpfile_name} \n vcf_donor, ${vcfdonor_name} \n field, ${field} \n geno_error_offset, ${geno_error_offset} \n geno_error_coeff, ${geno_error_coeff} \n r2_info, ${r2_info} \n min_mac, ${min_mac} \n min_callrate, ${min_callrate} \n alpha, ${alpha} \n doublet_prior, ${doublet_prior}" >> demuxlet_${task.index}/params.csv
        if [[ "$plp" != "True" ]]
        then
            popscle demuxlet $samfile $taggroup $tagUMI $vcfdonor $fieldinfo ${genoerror_off} ${genoerror_cof} $r2info $minmac $mincallrate $smlist ${sm_list_file} --alpha ${alpha_value} $doubletprior $samverbose $vcfverbose $capBQ $minBQ $minMQ $minTD $exclflag $grouplist $mintotal $minumi $minsnp $plpfile --out demuxlet_${task.index}/${demuxlet_out}
        else
            popscle dsc-pileup $samfile ${taggroup} ${tagUMI} $vcfref ${smlist} ${sm_list_file} ${samverbose} ${vcfverbose} ${skipumi} ${capBQ} ${minBQ} ${minMQ} ${minTD} ${exclflag} ${grouplist} ${mintotal} ${minsnp} --out demuxlet_${task.index}/plp/${demuxlet_out}
            popscle demuxlet $taggroup $tagUMI --plp demuxlet_${task.index}/plp/${demuxlet_out} $vcfdonor $fieldinfo ${genoerror_off} ${genoerror_cof} $r2info $minmac $mincallrate $smlist ${sm_list_file} --alpha ${alpha_value} $doubletprior $samverbose $vcfverbose $capBQ $minBQ $minMQ $minTD $exclflag $grouplist $mintotal $minumi $minsnp $minuniq --out demuxlet_${task.index}/${demuxlet_out}
            
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

workflow demultiplex_demuxlet{
    take:
        sam
        group_list
    main:
        tag_group = split_input(params.tag_group)
        tag_UMI = split_input(params.tag_UMI)
        vcfref = split_input(params.vcf_ref)
        sm = split_input(params.sm)
        sm_list = split_input(params.sm_list)
        sam_verbose = split_input(params.sam_verbose)
        vcf_verbose = split_input(params.vcf_verbose)
        skip_umi = split_input(params.skip_umi)
        cap_BQ = split_input(params.cap_BQ)
        min_BQ = split_input(params.min_BQ)
        min_MQ = split_input(params.min_MQ)
        min_TD = split_input(params.min_TD)
        excl_flag = split_input(params.excl_flag)
        min_total = split_input(params.min_total)
        min_umi = split_input(params.min_umi)
        min_uniq = split_input(params.min_uniq)
        min_snp = split_input(params.min_snp)

        plp = split_input(params.plp)
        vcfdonor = split_input(params.vcf_donor)
        field = split_input(params.field)
        geno_error_offset = split_input(params.geno_error_offset)
        geno_error_coeff = split_input(params.geno_error_coeff)
        r2_info= split_input(params.r2_info)
        min_mac = split_input(params.min_mac)
        min_callrate = split_input(params.min_callrate)
        alpha = split_input(params.alpha)
        doublet_prior = split_input(params.doublet_prior)
        demuxlet_out = split_input(params.demuxlet_out)

        demuxlet(sam, tag_group, tag_UMI, vcfref, sm, sm_list, sam_verbose, vcf_verbose, skip_umi, cap_BQ, min_BQ, min_MQ, min_TD, excl_flag, group_list, min_total, min_uniq, min_umi, min_snp, plp, vcfdonor, field, geno_error_offset, geno_error_coeff, r2_info, min_mac, min_callrate, alpha, doublet_prior, demuxlet_out)
        
    emit:
        demuxlet.out.collect()
}

workflow{
    demultiplex_demuxlet(channel.fromPath(params.bam))
}
