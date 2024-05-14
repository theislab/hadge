#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process subset_bam_and_sort_vcf_based_on_reference{
    label 'small_mem'
    conda "-c conda-forge -c bioconda samtools=19.2 bedtools bcftools=1.19"
    tag "${sampleId}"
    input:
        tuple val(sampleId), path(sam), path(sam_index), path(barcodes), val(vcf)

    output:
        tuple val(sampleId), path("${sampleId}_dmx__filtered_bam_file.bam"), path("${sampleId}_dmx__filtered_bam_file.bam.csi"), path(barcodes), path("${sampleId}_dmx__samples.sorted_as_in_bam.vcf"), emit: input
    when:
         vcf !='None'
    script:
        """
            filter_bam_file_for_popscle_dsc_pileup.sh ${sam} ${barcodes} ${vcf} ${sampleId}_dmx__filtered_bam_file.bam
            sort_vcf_same_as_bam.sh ${sampleId}_dmx__filtered_bam_file.bam ${vcf} > ${sampleId}_dmx__samples.sorted_as_in_bam.vcf        
        """
}

process demuxlet {
    publishDir "$params.outdir/$sampleId/$params.mode/gene_demulti/demuxlet", mode: 'copy'
    label 'small_mem'
    tag "${sampleId}"
    conda "bioconda::popscle"

    input:
        tuple val(sampleId), path(sam), path(sam_index), path(group_list), val(vcf_donor)
        val tag_group
        val tag_UMI
        val sm
        val sm_list
        val sam_verbose
        val vcf_verbose
        val skip_umi

        val cap_BQ
        val min_BQ
        val min_MQ
        val min_TD
        val excl_flag

        val min_total
        val min_uniq
        val min_umi
        val min_snp

        val plp
        val field
        val geno_error_offset
        val geno_error_coeff
        val r2_info
        val min_mac
        val min_callrate
        
        val alpha
        val doublet_prior
        val demuxlet_out
        

    output:
        path "demuxlet_${sampleId}"

    when:
        vcf_donor !='None'

    script:
        def samfile = "--sam $sam"
        def taggroup = tag_group != 'None' ? "--tag-group ${tag_group}" : ''
        def tagUMI = tag_UMI != 'None' ? "--tag-UMI ${tag_UMI}" : ''
        def vcfref = plp == 'True' ? "--vcf ${vcf_donor}" : ""
        def vcfref_name = plp == 'True' ? vcf_donor : "No VCF Ref is used because plp is not performed."
        def smlist = sm != 'None' ? "--sm $sm" : ''
        def sm_list_file = sm_list != 'None' ? "--sm-list ${sm_list}" : ''
        def sm_list_file_name = sm_list != 'None' ? file(sm_list).baseName : "No sm list file is given"
        def samverbose = "--sam-verbose ${sam_verbose}"
        def vcfverbose = "--vcf-verbose ${vcf_verbose}"
        def skipumi = skip_umi != "False" ? "--skip-umi" : ""
        def capBQ = "--cap-BQ ${cap_BQ}"
        def minBQ = "--min-BQ ${min_BQ}"
        def minMQ = "--min-MQ ${min_MQ}"
        def minTD = "--min-TD ${min_TD}"
        def exclflag = "--excl-flag ${excl_flag}"
        def grouplist = "--group-list ${group_list}"
        def mintotal = "--min-total ${min_total}"
        def minumi = "--min-umi ${min_umi}"
        def minuniq = "--min-uniq ${min_uniq}"
        def minsnp = "--min-snp ${min_snp}"
        def plp_name = plp == 'True' ? "plp performed" : "plp not performed"
        def vcfdonor = "--vcf ${vcf_donor}"
        def fieldinfo = "--field $field"
        def genoerror_off = "--geno-error-offset ${geno_error_offset}"
        def genoerror_cof = "--geno-error-coeff ${geno_error_coeff}"
        def r2info = "--r2-info ${r2_info}"
        def minmac = "--min-mac ${min_mac}"
        def mincallrate = "--min-callrate ${min_callrate}"
        def alpha_value = alpha.replaceAll(/,/, " --alpha ")
        def doubletprior = "--doublet-prior ${doublet_prior}"
      
        """
        mkdir demuxlet_${sampleId}
        touch demuxlet_${sampleId}/params.csv
        echo -e "Argument,Value \n samfile,${sam} \n tag_group,${tag_group} \n tag_UMI,${tag_UMI} \n vcf_ref,${vcfref_name} \n sm,${sm} \n sm_list_file,${sm_list_file_name} \n sam_verbose,${sam_verbose} \n vcf_verbose,${vcf_verbose} \n skip_umi,${skip_umi} \n cap_BQ,${cap_BQ} \n min_BQ,${min_BQ} \n min_MQ,${min_MQ} \n min_TD,${min_TD} \n excl_flag,${excl_flag} \n grouplist,${group_list} \n min_total,${min_total} \n min_uniq,${min_uniq} \n min_umi,${min_umi} \n min_snp,${min_snp} \n plp,${plp_name} \n vcf_donor,${vcf_donor} \n field,${field} \n geno_error_offset,${geno_error_offset} \n geno_error_coeff,${geno_error_coeff} \n r2_info,${r2_info} \n min_mac,${min_mac} \n min_callrate,${min_callrate} \n alpha,${alpha} \n doublet_prior,${doublet_prior}" >> demuxlet_${sampleId}/params.csv

        if [[ "$plp" == "False" ]]
        then
            popscle demuxlet $samfile $taggroup $tagUMI $vcfdonor $fieldinfo ${genoerror_off} ${genoerror_cof} $r2info $minmac \
            $mincallrate $smlist ${sm_list_file} --alpha ${alpha_value} $doubletprior $samverbose $vcfverbose $capBQ $minBQ \
            $minMQ $minTD $exclflag $grouplist $mintotal $minumi $minsnp --out demuxlet_${sampleId}/${demuxlet_out}
        else
            mkdir demuxlet_${sampleId}/plp
            popscle dsc-pileup $samfile ${taggroup} ${tagUMI} $vcfref ${smlist} ${sm_list_file} ${samverbose} ${vcfverbose} \
            ${skipumi} ${capBQ} ${minBQ} ${minMQ} ${minTD} ${exclflag} ${grouplist} ${mintotal} ${minsnp} \
            --out demuxlet_${sampleId}/plp/${demuxlet_out}
            popscle demuxlet $taggroup $tagUMI --plp demuxlet_${sampleId}/plp/${demuxlet_out} $vcfdonor $fieldinfo \
            ${genoerror_off} ${genoerror_cof} $r2info $minmac $mincallrate $smlist ${sm_list_file} --alpha ${alpha_value} \
            $doubletprior $samverbose $vcfverbose $capBQ $minBQ $minMQ $minTD $exclflag $grouplist $mintotal $minumi $minsnp \
            $minuniq --out demuxlet_${sampleId}/${demuxlet_out}
            
        fi
        """
        
}

workflow demultiplex_demuxlet{
    take:
        input_list
    main:

    
        tag_group = params.tag_group
        tag_UMI = params.tag_UMI
        sm = params.sm
        sm_list = params.sm_list
        sam_verbose = params.sam_verbose
        vcf_verbose = params.vcf_verbose
        skip_umi = params.skip_umi
        cap_BQ = params.cap_BQ
        min_BQ = params.min_BQ
        min_MQ = params.min_MQ
        min_TD = params.min_TD
        excl_flag = params.excl_flag
        min_total = params.min_total
        min_umi = params.min_umi
        min_uniq = params.min_uniq
        min_snp = params.min_snp
        plp = params.plp
        field = params.field
        geno_error_offset = params.geno_error_offset
        geno_error_coeff = params.geno_error_coeff
        r2_info= params.r2_info
        min_mac = params.min_mac
        min_callrate = params.min_callrate
        alpha = params.alpha
        doublet_prior = params.doublet_prior
        demuxlet_out = params.demuxlet_out


        subset_bam_and_sort_vcf_based_on_reference(input_list)
        input_list = subset_bam_and_sort_vcf_based_on_reference.out.input
        demuxlet(input_list, tag_group, tag_UMI, sm, sm_list, sam_verbose, vcf_verbose, skip_umi, cap_BQ, min_BQ, 
            min_MQ, min_TD, excl_flag, min_total, min_uniq, min_umi, min_snp, plp, field, geno_error_offset, geno_error_coeff, 
            r2_info, min_mac, min_callrate, alpha, doublet_prior, demuxlet_out)
        
    emit:
        demuxlet.out.collect()
}
