#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { data_preprocess } from './gene_demulti/samtools'
include { filter_variant } from './gene_demulti/bcftools'
include { variant_cellSNP } from './gene_demulti/cellsnp'
include { variant_freebayes } from './gene_demulti/freebayes'
include { demultiplex_demuxlet } from './gene_demulti/demuxlet'
include { demultiplex_freemuxlet } from './gene_demulti/freemuxlet'
include { demultiplex_scSplit } from './gene_demulti/scsplit'
include { demultiplex_souporcell } from './gene_demulti/souporcell'
include { demultiplex_vireo } from './gene_demulti/vireo'

def split_input(input){
    if (input =~ /;/ ){
        Channel.from(input).map{ return it.tokenize(';')}.flatten()
    }
    else{
        Channel.from(input)
    }
}


process subset_bam_to_comon_variants{
    
    label 'small_mem'
    conda "bioconda::samtools=1.19.2 bedtools bcftools=1.19"
    tag "${sampleId}"
    input:
        tuple val(sampleId), path(sam), path(sam_index), path(barcodes)
        path vcf

    output:
        tuple val(sampleId), path("${sampleId}__filtered_bam_file.bam"), path("${sampleId}__filtered_bam_file.bam.csi"), emit: input
    
    script:
        """ 
        
            bcftools sort ${vcf} -Oz -o sorted.vcf.gz
            filter_bam_file_for_popscle_dsc_pileup.sh ${sam} ${barcodes} sorted.vcf.gz ${sampleId}__filtered_bam_file.bam      
        """

}

process summary{
    publishDir "$projectDir/$params.outdir/$sampleId/$params.mode/gene_demulti", mode: 'copy'
    label 'small_mem'
    tag "${sampleId}"
    conda "pandas scanpy mudata"

    input:
        tuple val(sampleId), path(hto_matrix, stageAs: 'hto_data'), path(rna_matrix, stageAs: 'rna_data')
        val demuxlet_result
        val freemuxlet_result
        val vireo_result
        val souporcell_result
        val scsplit_result
        val generate_anndata
        val generate_mudata
  
    output:
        tuple val(sampleId), path("genetic_summary")

    script:
        def demuxlet_files = ""
        def freemuxlet_files = ""
        def vireo_files = ""
        def souporcell_files = ""
        def scsplit_files = ""
        def generate_adata = ""
        def generate_mdata = ""
        
        if (demuxlet_result != "no_result"){
            demuxlet_res = demuxlet_result.find{it.name.contains(sampleId)}
            demuxlet_files = "--demuxlet ${demuxlet_res}"
        }
        if (freemuxlet_result != "no_result"){
            freemuxlet_res = freemuxlet_result.find{it.name.contains(sampleId)}
            freemuxlet_files = "--freemuxlet ${freemuxlet_res}"
        }
        if (vireo_result != "no_result"){
            vireo_res = vireo_result.find{it.name.contains(sampleId)}
            vireo_files = "--vireo ${vireo_res}"
        }
        if (souporcell_result != "no_result"){
            souporcell_res = souporcell_result.find{it.name.contains(sampleId)}
            souporcell_files = "--souporcell ${souporcell_res}"
        }
        if (scsplit_result != "no_result"){
            scsplit_res = scsplit_result.find{it.name.contains(sampleId)}
            scsplit_files =  "--scsplit ${scsplit_res}"
        }
        if (generate_anndata == "True"){
            if(rna_matrix.name == "None"){
                error "Error: RNA count matrix is not given."
            }
            generate_adata = "--generate_anndata --read_rna_mtx rna_data"
        }
        if (generate_mudata == "True"){
            if(rna_matrix.name == "None"){
                error "Error: RNA count matrix is not given."
            }
            if(hto_matrix.name == "None"){
                error "Error: HTO count matrix is not given."
            }
            generate_mdata = "--generate_mudata --read_rna_mtx rna_data --read_hto_mtx hto_data"
        }
        """
        summary_gene.py $demuxlet_files $vireo_files $souporcell_files $scsplit_files $freemuxlet_files $generate_adata $generate_mdata
        """
}



workflow gene_demultiplexing {

    if ((params.demuxlet == "True" & params.demuxlet_preprocess == "True")      | \
       (params.freemuxlet == "True" & params.freemuxlet_preprocess == "True")   | \
       (params.scSplit == "True" & params.scSplit_preprocess  == "True")        | \
       (params.vireo == "True" & params.vireo_preprocess == "True")             | \
       (params.souporcell == "True" & params.souporcell_preprocess == "True"))  {

                Channel.fromPath(params.multi_input) \
                        | splitCsv(header:true) \
                        | map { row-> tuple(row.sampleId, row.bam)}
                        | data_preprocess
                qc_bam = data_preprocess.out.map{ it -> tuple( it.name.tokenize( '_' ).last(), it + "/sorted.bam", it + "/sorted.bam.bai") }
        }else{
            qc_bam = Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, row.bam, row.bam_index)}



        }

    input_param_cellsnp = Channel.fromPath(params.multi_input) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId, row.barcodes) }
    qc_bam_new = qc_bam.join(input_param_cellsnp)

    qc_bam = subset_bam_to_comon_variants(qc_bam_new,params.common_variants_freemuxlet)


    //////////
    //FreeBayes/ scSplit
    //////////
    if (params.scSplit == "True" & params.scSplit_variant == 'True'){

        freebayes_region = Channel.from(1..22, "X","Y").flatten()
        if (params.region != "None"){
            freebayes_region = split_input(params.region)
        }

        variant_freebayes(qc_bam, freebayes_region)
        filter_variant(variant_freebayes.out)
        freebayes_vcf = filter_variant.out.map{ it -> tuple(it[0], it[1] + "/filtered_sorted_total_chroms.vcf")}  
    }

    if (params.scSplit == "True"){


        input_bam_scsplit = qc_bam

        if (params.scSplit_variant == 'True'){
            input_vcf_scsplit = freebayes_vcf
        }
        else{

            input_vcf_scsplit = Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, row.vcf_mixed)}
        }

        input_param_scsplit = Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, row.barcodes, row.nsample, row.vcf_donor)}
        
        input_list_scsplit = input_bam_scsplit.join(input_vcf_scsplit)
        input_list_scsplit = input_list_scsplit.join(input_param_scsplit)
        demultiplex_scSplit(input_list_scsplit)
        scSplit_out = demultiplex_scSplit.out
    }
    else{
        scSplit_out = channel.value("no_result")
    }


    //////////
    //CellSNP/Vireo
    //////////
    if (params.vireo == "True" & params.vireo_variant == 'True'){

        variant_cellSNP(qc_bam_new)
        cellsnp_vcf = variant_cellSNP.out.out1.map{ it -> tuple( it.name.tokenize( '_' ).last(), it + "/*/cellSNP.cells.vcf") }

    }

    if (params.vireo == "True"){

        if (params.vireo_variant == 'True'){
            input_vcf_vireo = variant_cellSNP.out.cellsnp_input
        }
        else{
            input_vcf_vireo = Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, row.celldata)}
        }
        input_param_vireo = Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, row.nsample, row.vcf_donor)}
        
        input_list_vireo = input_vcf_vireo.join(input_param_vireo)
        demultiplex_vireo(input_list_vireo)
        vireo_out = demultiplex_vireo.out
    }
    else{
        vireo_out = channel.value("no_result")
    }


    //////////
    // Demuxlet/Freemuxlet
    // demuxlet (with genotypes) or freemuxlet (without genotypes)
    //////////

    if (params.demuxlet == "True"){

        input_bam_demuxlet = qc_bam

        input_param_demuxlet = Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, row.barcodes, row.vcf_donor)}
        input_list_demuxlet = input_bam_demuxlet.join(input_param_demuxlet)
        demultiplex_demuxlet(input_list_demuxlet)
        demuxlet_out = demultiplex_demuxlet.out
    }
    else{
        demuxlet_out = channel.value("no_result")
    }


    //////////
    //Freemuxlet
    //////////

    if (params.freemuxlet == "True"){

        input_bam_freemuxlet = qc_bam
        
        input_param_freemuxlet = Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, row.barcodes, row.nsample)}

        input_list_freemuxlet = input_bam_freemuxlet.join(input_param_freemuxlet)

        demultiplex_freemuxlet(input_list_freemuxlet)
        freemuxlet_out = demultiplex_freemuxlet.out
    }
    else{
        freemuxlet_out = channel.value("no_result")
    }


    //////////
    //Souporcell
    //////////

     if (params.souporcell == "True"){
        
        input_bam_souporcell = qc_bam

        input_param_souporcell = Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, row.barcodes, row.nsample, row.vcf_donor)}
                
        input_list_souporcell = input_bam_souporcell.join(input_param_souporcell)
        demultiplex_souporcell(input_list_souporcell)
        souporcell_out = demultiplex_souporcell.out
    }
    else{
        souporcell_out = channel.value("no_result")
    }

    //////////
    //Summary
    //////////
    
    Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, file(row.hto_matrix_filtered), file(row.rna_matrix_filtered))}
                | set {input_list_summary}

    summary(input_list_summary, demuxlet_out, freemuxlet_out, vireo_out, souporcell_out, scSplit_out, 
            params.generate_anndata, params.generate_mudata)
    
    
    emit:
        summary.out 
}

