#!/usr/bin/env nextflow
nextflow.enable.dsl=2
'''
BAM Index should in the format BAMFILE.bam.bai

'''
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

process summary{
    publishDir "$projectDir/$params.outdir/$params.mode/gene_demulti", mode: 'copy'
    input:
        val demuxlet_result
        val freemuxlet_result
        val vireo_result
        val souporcell_result
        val scsplit_result
    output:
        path gene_summary

    script:
        def demuxlet_files = ""
        def freemuxlet_files = ""
        def vireo_files = ""
        def souporcell_files = ""
        def scsplit_files = ""
        
        if (demuxlet_result != "no_result"){
            demuxlet_files = "--demuxlet "
            for(r : demuxlet_result) {
                demuxlet_files  = demuxlet_files + r + ":"
            }
        }
        if (freemuxlet_result != "no_result"){
            freemuxlet_files = "--freemuxlet "
            for(r : freemuxlet_result) {
                freemuxlet_files = freemuxlet_files + r + ":"
            }
        }
        if (vireo_result != "no_result"){
            vireo_files = "--vireo "
            for(r : vireo_result) {
                vireo_files = vireo_files + r + ":"
            }
        }
        if (souporcell_result != "no_result"){
            souporcell_files = "--souporcell "
            for(r : souporcell_result) {
                souporcell_files = souporcell_files + r + ":"
            }
        }
        if (scsplit_result != "no_result"){
            scsplit_files = "--scsplit "
            for(r : scsplit_result) {
                scsplit_files = scsplit_files + r + ":"
            }
        }
        
        """
        mkdir gene_summary && cd gene_summary
        summary_gene.R $demuxlet_files $vireo_files $souporcell_files $scsplit_files $freemuxlet_files
        """
}


workflow gene_demultiplexing {
    take:
    recover_barcodes
    main:
    input_bam = params.cellranger == 'True'? Channel.value(params.cellranger_dir).map{ return it + "/outs/possorted_genome_bam.bam"} : Channel.fromPath(params.bam)
    input_bai = params.cellranger == 'True'? Channel.value(params.cellranger_dir).map{ return it + "/outs/possorted_genome_bam.bam.bai"}: Channel.fromPath(params.bam).map{ return it + ".bai"}

    if ((params.demuxlet == "True" & params.demuxlet_preprocess != 'False')| \
       (params.freemuxlet == "True" & params.freemuxlet_preprocess != 'False')| \
       (params.scSplit == "True" & params.scSplit_preprocess != 'False') | \
       (params.souporcell == "True" & params.souporcell_preprocess != 'False')){
        data_preprocess(input_bam)
        qc_bam = data_preprocess.out.map{ return it + "/sorted.bam"}
        qc_bam_bai = data_preprocess.out.map{ return it + "/sorted.bam.bai"}
    }

    if (params.vireo == "True" &  params.vireo_variant != 'False'){
        if (params.vireo_variant == 'cellSNP' | params.vireo_variant == 'Both'){
            if(params.vireo_preprocess != 'False'){
                variant_cellSNP(qc_bam, recover_barcodes)
            }
            else{
                variant_cellSNP(input_bam, recover_barcodes)
            }
            cellsnp_vcf = variant_cellSNP.out.map{ return it + "/*/cellSNP.cells.vcf"}
        }
        else{
            if (params.scSplit == "True" & params.scSplit_variant != 'False'){
                if (params.scSplit_variant == 'cellSNP' | params.scSplit_variant == 'Both'){
                    if(params.scSplit_preprocess != 'False'){
                        variant_cellSNP(qc_bam, recover_barcodes)
                    }
                    else{
                        variant_cellSNP(input_bam, recover_barcodes)
                    }
                    cellsnp_vcf = variant_cellSNP.out.map{ return it + "/*/cellSNP.cells.vcf"}
                }
            }
        }
       
        if (params.vireo_variant == "freebayes" | params.vireo_variant == 'Both'){
            freebayes_region = Channel.from(1..22, "X","Y").flatten()
            if (params.region != "False"){
                freebayes_region = Channel.value(params.region)
            }
            if(params.vireo_preprocess != 'False'){
                variant_freebayes(qc_bam, qc_bam_bai, freebayes_region)
            }
            else{
                variant_freebayes(input_bam, input_bai, freebayes_region)
            }
            filter_variant(variant_freebayes.out, "True","True")
            freebayes_vcf = filter_variant.out.map{ return it + "/filtered_sorted_total_chroms.vcf"}

        }
        else{
            if (params.scSplit == "True" & params.scSplit_variant != 'False'){
                if (params.scSplit_variant == 'freebayes' | params.scSplit_variant == 'Both'){
                    freebayes_region = Channel.from(1..22, "X","Y").flatten()
                    if (params.region != "False"){
                        freebayes_region = split_input(params.region)
                    }
                    if(params.scSplit_preprocess != 'False'){
                        variant_freebayes(qc_bam, qc_bam_bai, freebayes_region)
                    }
                    else{
                        variant_freebayes(input_bam, input_bai, freebayes_region)
                    }
                    filter_variant(variant_freebayes.out, "True","True")
                    freebayes_vcf = filter_variant.out.map{ return it + "/filtered_sorted_total_chroms.vcf"}
                }
            }
        }
    }
    else{
        if (params.scSplit == "True" & params.scSplit_variant != 'False'){
            if (params.scSplit_variant == 'cellSNP' | params.scSplit_variant == 'Both'){
                if(params.scSplit_preprocess != 'False'){
                    variant_cellSNP(qc_bam, recover_barcodes)
                }
                else{
                    variant_cellSNP(input_bam, recover_barcodes)
                }
                cellsnp_vcf = variant_cellSNP.out.map{ return it + "/*/cellSNP.cells.vcf"}
            }
            if (params.scSplit_variant == "freebayes" | params.scSplit_variant == 'Both'){
                freebayes_region = Channel.from(1..22, "X","Y","MT").flatten()
                if (params.region != "False"){
                    freebayes_region = split_input(params.region)
                }
                if(params.scSplit_preprocess != 'False'){
                    variant_freebayes(qc_bam, qc_bam_bai, freebayes_region)
                }
                else{
                    variant_freebayes(input_bam, input_bai, freebayes_region)
                }
                filter_variant(variant_freebayes.out, "True","True")
                freebayes_vcf = filter_variant.out.map{ return it + "/filtered_sorted_total_chroms.vcf"}
            }
        }
    }
  
    
    if (params.demuxlet == "True"){
        bam = params.demuxlet_preprocess == 'True'? qc_bam: (params.demuxlet_preprocess == 'False'? input_bam : qc_bam.mix(input_bam))
        demultiplex_demuxlet(bam, recover_barcodes)
        demuxlet_out = demultiplex_demuxlet.out
    }
    else{
        demuxlet_out = channel.value("no_result")
    }
    
    
    if (params.freemuxlet == "True"){
        bam = params.freemuxlet_preprocess == 'True'? qc_bam: (params.freemuxlet_preprocess == 'False'? input_bam : qc_bam.mix(input_bam))
        demultiplex_freemuxlet(bam, recover_barcodes)
        freemuxlet_out = demultiplex_freemuxlet.out
    }
    else{
        freemuxlet_out = channel.value("no_result")
    }

    
    if (params.vireo == "True"){
        vcf = params.vireo_variant == 'False'? Channel.fromPath(params.celldata): \
            (params.vireo_variant == 'freebayes'? freebayes_vcf: \
            (params.vireo_variant == 'cellSNP'? variant_cellSNP.out.map{ return it + "/*/cellSNP.cells.vcf"} : \
            (params.vireo_variant == 'Both'? variant_cellSNP.out.map{ return it + "/*/cellSNP.cells.vcf"}.mix(freebayes_vcf) : \
            (params.vireo_variant == 'noCellSNP'? freebayes_vcf.mix(Channel.fromPath(params.celldata)) : \
            (params.vireo_variant == 'noFreebayes'? variant_cellSNP.out.map{ return it + "/*/cellSNP.cells.vcf"}.mix(Channel.fromPath(params.celldata)) : \
            variant_cellSNP.out.map{ return it + "/*/cellSNP.cells.vcf"}.mix(freebayes_vcf, Channel.fromPath(params.celldata)))))))
        demultiplex_vireo(vcf)
        vireo_out = demultiplex_vireo.out
    }
    else{
        vireo_out = channel.value("no_result")
    }

   
    if (params.scSplit == "True"){
        bam = params.scSplit_preprocess == 'True'? qc_bam: (params.scSplit_preprocess == 'False'? input_bam : qc_bam.mix(input_bam))
        bai = params.scSplit_preprocess == 'True'? qc_bam_bai: (params.scSplit_preprocess == 'False'? input_bai : qc_bam_bai.mix(input_bai))
        vcf = params.scSplit_variant == 'False'? Channel.fromPath(params.vcf_mixed): \
            (params.scSplit_variant == 'freebayes'? filter_variant.out.map{ return it + "/filtered_sorted_total_chroms.vcf"} : \
            (params.scSplit_variant == 'cellSNP'? variant_cellSNP.out.map{ return it + "/*/cellSNP.cells.vcf"} : \
            (params.scSplit_variant == 'Both'? variant_cellSNP.out.map{ return it + "/*/cellSNP.cells.vcf"}.mix(freebayes_vcf) : \
            (params.scSplit_variant == 'noCellSNP'? freebayes_vcf.mix(Channel.fromPath(params.vcf_mixed)) : \
            (params.scSplit_variant == 'noFreebayes'? variant_cellSNP.out.map{ return it + "/*/cellSNP.cells.vcf"}.mix(Channel.fromPath(params.vcf_mixed)) : \
            variant_cellSNP.out.map{ return it + "/*/cellSNP.cells.vcf"}.mix(freebayes_vcf, Channel.fromPath(params.vcf_mixed)))))))
        demultiplex_scSplit(bam, vcf, bai, recover_barcodes)
        scSplit_out = demultiplex_scSplit.out
    }
    else{
        scSplit_out = channel.value("no_result")
    }

    
    if (params.souporcell == "True"){
        bam = params.souporcell_preprocess == 'True'? qc_bam: (params.souporcell_preprocess == 'False'? input_bam : qc_bam.mix(input_bam))
        demultiplex_souporcell(bam, recover_barcodes)
        souporcell_out = demultiplex_souporcell.out
    }
    else{
        souporcell_out = channel.value("no_result")
    }

    summary(demuxlet_out, freemuxlet_out, vireo_out, souporcell_out, scSplit_out)
    emit:
    summary.out
}

