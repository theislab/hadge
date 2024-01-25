#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process freemuxlet {
    publishDir "$projectDir/$params.outdir/$params.mode/gene_demulti/freemuxlet", mode: 'copy'
    label 'small_mem'

    conda "bioconda::popscle bioconda::samtools bioconda::bedtools bioconda::bcftools=1.9"

    input:
        each sam
        path(barcodes)
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
        def samfile = "--sam filtered_bam_file.bam"
        def taggroup = tag_group != 'None' ? "--tag-group ${tag_group}" : ''
        def tagUMI = tag_UMI != 'None' ? "--tag-UMI ${tag_UMI}" : ''
        def vcffile = "--vcf samples.sorted_as_in_bam.vcf"
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
        def minuniq = "--min-uniq ${min_uniq}"
        def minumi = "--min-umi ${min_umi}"
        def minsnp = "--min-snp ${min_snp}"
        def initcluster = init_cluster != 'None' ? "--init-cluster ${init_cluster}" : ''
        def n_sample = "--nsample $nsample"
        def auxfiles = aux_files != 'False' ? "--aux-files" : ''
        def verbose_info = "--verbose $verbose"
        def doubletprior = "--doublet-prior ${doublet_prior}"
        def bfthres = "--bf-thres ${bf_thres}"
        def frac_init_cluster = "--frac-init-clust ${frac_init_clust}"
        def iterinit = "--iter-init ${iter_init}"
        def keepinit_missing = keep_init_missing != "False" ? "--keep-init-missing" : ''
        
        """
        echo 'test5'
        filter_bam_file_for_popscle_dsc_pileup.sh ${sam} ${barcodes} ${vcf} filtered_bam_file.bam
        sort_vcf_same_as_bam.sh filtered_bam_file.bam ${vcf} > samples.sorted_as_in_bam.vcf
        
        mkdir freemuxlet_${task.index}
        mkdir freemuxlet_${task.index}/plp
        touch freemuxlet_${task.index}/params.csv
        echo -e "Argument,Value \n samfile,filtered_bam_file.bam \n tag_group,${tag_group} \n tag_UMI,${tag_UMI} \n vcf_file,samples.sorted_as_in_bam.vcf \n sm,${sm} \n sm_list_file,${sm_list_file_name} \n sam_verbose,${sam_verbose} \n vcf_verbose,${vcf_verbose} \n skip_umi,${skip_umi} \n cap_BQ,${cap_BQ} \n min_BQ,${min_BQ} \n min_MQ,${min_MQ} \n min_TD,${min_TD} \n excl_flag,${excl_flag} \n grouplist,${group_list} \n min_total,${min_total} \n min_uniq,${min_uniq} \n min_umi,${min_umi} \n min_snp,${min_snp} \n init_cluster,${init_cluster} \n nsample,${nsample} \n aux_files,${aux_files} \n verbose,${verbose} \n doublet_prior,${doublet_prior} \n bf_thres,${bf_thres} \n frac_init_clust,${frac_init_clust} \n iter_init,${iter_init} \n keep_init_missing,${keep_init_missing}" >> freemuxlet_${task.index}/params.csv
        
        popscle dsc-pileup $samfile ${taggroup} ${tagUMI} $vcffile ${smlist} ${sm_list_file} ${samverbose} \
            ${vcfverbose} ${skipumi} ${capBQ} ${minBQ} ${minMQ} ${minTD} ${exclflag} ${grouplist} ${mintotal} ${minuniq} \
            ${minsnp} --out freemuxlet_${task.index}/plp/${freemuxlet_out}
        popscle freemuxlet --plp freemuxlet_${task.index}/plp/${freemuxlet_out} --out freemuxlet_${task.index}/${freemuxlet_out} \
            ${initcluster} ${n_sample} ${auxfiles} ${verbose_info} ${doubletprior} ${bfthres} ${frac_init_cluster} ${iterinit} \
            ${keepinit_missing} ${capBQ} ${minBQ} ${grouplist} ${mintotal} ${minumi} ${minsnp}
            
        """
        
}


def split_input(input){
    if (input =~ /;/ ){
        Channel.from(input).map{return it.tokenize(';')}.flatten()
    }
    else{
        Channel.from(input)
    }
}


workflow demultiplex_freemuxlet{
    take:
        sam
    main:
        group_list = split_input(params.barcodes)
        vcf = split_input(params.common_variants_freemuxlet)
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
        init_cluster = split_input(params.init_cluster)
        nsample = split_input(params.nsample)
        aux_files = split_input(params.aux_files)
        verbose = split_input(params.verbose)
        doublet_prior = split_input(params.doublet_prior)
        bf_thres = split_input(params.bf_thres)
        frac_init_clust = split_input(params.frac_init_clust)
        iter_init = split_input(params.iter_init)
        keep_init_missing = split_input(params.keep_init_missing)
        freemuxlet_out = params.freemuxlet_out

        freemuxlet(sam, params.barcodes, vcf, tag_group, tag_UMI, sm, sm_list, sam_verbose, vcf_verbose, skip_umi, cap_BQ,
            min_BQ, min_MQ, min_TD, excl_flag, group_list, min_total, min_uniq, min_umi, min_snp, init_cluster, nsample, 
            aux_files, verbose, doublet_prior, bf_thres, frac_init_clust, iter_init, keep_init_missing, freemuxlet_out)
    
    emit:
        freemuxlet.out.collect()
}
