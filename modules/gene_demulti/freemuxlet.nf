#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process freemuxlet {
    publishDir "$projectDir/$params.outdir/$params.mode/gene_demulti/freemuxlet", mode: 'copy'
    input:
        each sam
        each vcf
        each tag_group
        each tag_UMI
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
        each init_cluster
        each nsample
        each aux_files
        each verbose
        each doublet_prior
        each bf_thres
        each frac_init_clust
        each iter_init
        each keep_init_missing
        each freemuxlet_out

    output:
        path "freemuxlet_${task.index}"

    script:
        def samfile = plp == 'True' ? "--sam $sam" : ''
        def samfile_name = plp == 'True' ? sam.baseName : 'sam_file_not_used'
        def taggroup = tag_group != 'False' ? "--tag-group ${tag_group}" : ''
        def tagUMI = tag_UMI != 'False' ? "--tag-UMI ${tag_UMI}" : ''
        def vcffile = plp == 'True' ? "--vcf $vcf" : ''
        def vcffile_name = plp == 'True' ? file(vcf).baseName : "vcf_file_not_used"
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
        def minuniq = "--min-uniq ${min_uniq}"
        def minumi = "--min-umi ${min_umi}"
        def minsnp = "--min-snp ${min_snp}"
        def plpfile_name = plp != 'True' ? file(plp).baseName : "plp_file_not_given"
        def initcluster = init_cluster != 'False' ? "--init-cluster ${init_cluster}" : ''
        def n_sample = "--nsample $nsample"
        def auxfiles = aux_files != 'False' ? "--aux-files" : ''
        def verbose_info = "--verbose $verbose"
        def doubletprior = "--doublet-prior ${doublet_prior}"
        def bfthres = "--bf-thres ${bf_thres}"
        def frac_init_cluster = "--frac-init-clust ${frac_init_clust}"
        def iterinit = "--iter-init ${iter_init}"
        def keepinit_missing = keep_init_missing != "False" ? "--keep-init-missing" : ''
        
        """
        mkdir freemuxlet_${task.index}
        mkdir freemuxlet_${task.index}/plp
        touch freemuxlet_${task.index}/params.csv
        barcode_num=\$(wc -l < "${group_list}")
        echo -e "Argument,Value \n samfile,${samfile_name} \n tag_group,${tag_group} \n tag_UMI, ${tag_UMI} \n vcf_file, ${vcffile_name} \n sm, ${sm} \n sm_list_file, ${sm_list_file_name} \n sam_verbose, ${sam_verbose} \n vcf_verbose, ${vcf_verbose} \n skip_umi, ${skip_umi} \n cap_BQ, ${cap_BQ} \n min_BQ, ${min_BQ} \n min_MQ, ${min_MQ} \n min_TD, ${min_TD} \n excl_flag, ${excl_flag} \n grouplist, ${grouplist_name}_\${barcode_num} \n min_total, ${min_total} \n min_uniq, ${min_uniq} \n min_umi, ${min_umi} \n min_snp, ${min_snp} \n plpfile, ${plpfile_name} \n init_cluster, ${init_cluster} \n nsample, ${nsample} \n aux_files, ${aux_files} \n verbose, ${verbose} \n doublet_prior, ${doublet_prior} \n bf_thres, ${bf_thres} \n frac_init_clust, ${frac_init_clust} \n iter_init, ${iter_init} \n keep_init_missing, ${keep_init_missing}" >> freemuxlet_${task.index}/params.csv
        if [[ "$plp" != "True" ]]
        then
            popscle freemuxlet --plp ${plp} ${initcluster} ${n_sample} ${auxfiles} ${verbose_info} ${doubletprior} ${bfthres} ${frac_init_cluster} ${iterinit} ${keepinit_missing} ${capBQ} ${minBQ} ${grouplist} ${mintotal} ${minumi} ${minsnp} --out freemuxlet_${task.index}/${freemuxlet_out}
        else
            popscle dsc-pileup $samfile ${taggroup} ${tagUMI} $vcffile ${smlist} ${sm_list_file} ${samverbose} ${vcfverbose} ${skipumi} ${capBQ} ${minBQ} ${minMQ} ${minTD} ${exclflag} ${grouplist} ${mintotal} ${minuniq} ${minsnp} --out freemuxlet_${task.index}/plp/${freemuxlet_out}
            popscle freemuxlet --plp freemuxlet_${task.index}/plp/${freemuxlet_out} --out freemuxlet_${task.index}/${freemuxlet_out} ${initcluster} ${n_sample} ${auxfiles} ${verbose_info} ${doubletprior} ${bfthres} ${frac_init_cluster} ${iterinit} ${keepinit_missing} ${capBQ} ${minBQ} ${grouplist} ${mintotal} ${minumi} ${minsnp}
            
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


workflow demultiplex_freemuxlet{
    take:
        sam
        group_list
    main:
        vcf = split_input(params.vcf_ref)
        tag_group = split_input(params.tag_group)
        tag_UMI = split_input(params.tag_UMI)
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
        min_uniq = split_input(params.min_uniq)
        min_umi = split_input(params.min_umi)
        min_snp = split_input(params.min_snp)
        
        plp = split_input(params.plp_freemuxlet)

        init_cluster = split_input(params.init_cluster)
        nsample = split_input(params.nsample)
        aux_files = split_input(params.aux_files)
        verbose = split_input(params.verbose)
        doublet_prior = split_input(params.doublet_prior)
        bf_thres = split_input(params.bf_thres)
        frac_init_clust = split_input(params.frac_init_clust)
        iter_init = split_input(params.iter_init)
        keep_init_missing = split_input(params.keep_init_missing)
        freemuxlet_out = split_input(params.freemuxlet_out)

        freemuxlet(sam, vcf, tag_group, tag_UMI, sm, sm_list, sam_verbose, vcf_verbose, skip_umi, cap_BQ, min_BQ, min_MQ, min_TD, excl_flag, group_list, min_total, min_uniq, min_umi, min_snp, plp, init_cluster, nsample, aux_files, verbose, doublet_prior, bf_thres, frac_init_clust, iter_init, keep_init_missing, freemuxlet_out)
    
    emit:
        freemuxlet.out.collect()
}

workflow{
        demultiplex_freemuxlet(channel.fromPath(params.bam))
}
