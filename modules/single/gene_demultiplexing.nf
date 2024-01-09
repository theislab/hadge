#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { data_preprocess } from './gene_demulti/samtools'
include { filter_variant } from './gene_demulti/bcftools'
include { variant_cellSNP } from './gene_demulti/cellsnp'
include { variant_freebayes } from './gene_demulti/freebayes'
include { demultiplex_demuxlet } from './gene_demulti/demuxlet'
include { demultiplex_freemuxlet } from './gene_demulti/freemuxlet'
include { demultiplex_scSplit } from './gene_demulti/scsplit'
include { demultiplex_souporcell } from './gene_demulti/souporcell'
include { demultiplex_vireo } from './gene_demulti/vireo'

def split_input(input) {
    if (input =~ /;/) {
        Channel.from(input).map { return it.tokenize(';') }.flatten()
    }
    else {
        Channel.from(input)
    }
}

process summary {
    publishDir "$projectDir/$params.outdir/$params.mode/gene_demulti", mode: 'copy'
    label 'small_mem'

    conda 'pandas scanpy mudata'

    input:
        val demuxlet_result
        val freemuxlet_result
        val vireo_result
        val souporcell_result
        val scsplit_result
        val generate_anndata
        val generate_mudata
        val rna_matrix
        val hto_matrix
    output:
        path genetic_summary

    script:
        def demuxlet_files = ''
        def freemuxlet_files = ''
        def vireo_files = ''
        def souporcell_files = ''
        def scsplit_files = ''
        def generate_adata = ''
        def generate_mdata = ''

        if (demuxlet_result != 'no_result') {
            demuxlet_files = "--demuxlet ${demuxlet_result.join(':')}"
        }
        if (freemuxlet_result != 'no_result') {
            freemuxlet_files = "--freemuxlet ${freemuxlet_result.join(':')}"
        }
        if (vireo_result != 'no_result') {
            vireo_files = "--vireo ${vireo_result.join(':')}"
        }
        if (souporcell_result != 'no_result') {
            souporcell_files = "--souporcell ${souporcell_result.join(':')}"
        }
        if (scsplit_result != 'no_result') {
            scsplit_files =  "--scsplit ${scsplit_result.join(':')}"
        }
        if (scsplit_result != 'no_result') {
            scsplit_files =  "--scsplit ${scsplit_result.join(':')}"
        }
        if (generate_anndata == 'True') {
            if (rna_matrix == 'None') {
                error 'Error: RNA count matrix is not given.'
            }
            generate_adata = "--generate_anndata --read_rna_mtx $rna_matrix"
        }
        if (generate_mudata == 'True') {
            if (rna_matrix == 'None') {
                error 'Error: RNA count matrix is not given.'
            }
            if (hto_matrix == 'None') {
                error 'Error: HTO count matrix is not given.'
            }
            generate_mdata = "--generate_mudata --read_rna_mtx $rna_matrix --read_hto_mtx $hto_matrix"
        }
        """
        summary_gene.py $demuxlet_files $vireo_files $souporcell_files $scsplit_files $freemuxlet_files $generate_adata $generate_mdata
        """
}

workflow gene_demultiplexing {
    main:
        input_bam = Channel.fromPath(params.bam)
        input_bai = Channel.fromPath(params.bai)

        if ((params.demuxlet == 'True' & params.demuxlet_preprocess != 'False') | \
        (params.freemuxlet == 'True' & params.freemuxlet_preprocess != 'False')| \
        (params.scSplit == 'True' & params.scSplit_preprocess != 'False') | \
        (params.vireo == 'True' & params.vireo_preprocess != 'False') | \
        (params.souporcell == 'True' & params.souporcell_preprocess != 'False')) {
            data_preprocess(input_bam)
            qc_bam = data_preprocess.out.map { return it + '/sorted.bam' }
            qc_bam_bai = data_preprocess.out.map { return it + '/sorted.bam.bai' }
        }

        if (params.vireo == 'True' & params.vireo_variant == 'True') {
            if (params.vireo_preprocess != 'False') {
                variant_cellSNP(qc_bam, qc_bam_bai)
            }
            else {
                variant_cellSNP(input_bam, input_bai)
            }
            cellsnp_vcf = variant_cellSNP.out.map { return it + '/*/cellSNP.cells.vcf' }
        }

        if (params.scSplit == 'True' & params.scSplit_variant == 'True') {
            freebayes_region = Channel.from(1..22, 'X', 'Y').flatten()
            if (params.region != 'None') {
                freebayes_region = split_input(params.region)
            }
            if (params.scSplit_preprocess != 'False') {
                variant_freebayes(qc_bam, qc_bam_bai, freebayes_region)
            }
            else {
                variant_freebayes(input_bam, input_bai, freebayes_region)
            }
            filter_variant(variant_freebayes.out)
            freebayes_vcf = filter_variant.out.map { return it + '/filtered_sorted_total_chroms.vcf' }
        }

        if (params.demuxlet == 'True') {
            bam = params.demuxlet_preprocess == 'True' ? qc_bam : input_bam //qc_bam.mix(input_bam))
            demultiplex_demuxlet(bam)
            demuxlet_out = demultiplex_demuxlet.out
        }
        else {
            demuxlet_out = channel.value('no_result')
        }

        if (params.freemuxlet == 'True') {
            bam = params.freemuxlet_preprocess == 'True' ? qc_bam : input_bam // qc_bam.mix(input_bam))
            demultiplex_freemuxlet(bam)
            freemuxlet_out = demultiplex_freemuxlet.out
        }
        else {
            freemuxlet_out = channel.value('no_result')
        }

        if (params.vireo == 'True') {
            vcf = params.vireo_variant != 'True' ? Channel.fromPath(params.celldata) :
                    variant_cellSNP.out.map { return it + '/*/cellSNP.cells.vcf' }
            //variant_cellSNP.out.map{ return it + "/*/cellSNP.cells.vcf"}.mix(Channel.fromPath(params.celldata)))
            demultiplex_vireo(vcf)
            vireo_out = demultiplex_vireo.out
        }
        else {
            vireo_out = channel.value('no_result')
        }

        if (params.scSplit == 'True') {
            bam = params.scSplit_preprocess == 'True' ? qc_bam : input_bam // qc_bam.mix(input_bam))
            bai = params.scSplit_preprocess == 'True' ? qc_bam_bai : input_bai // qc_bam_bai.mix(input_bai))
            vcf = params.scSplit_variant != 'True' ? Channel.fromPath(params.vcf_mixed) :
                    filter_variant.out.map { return it + '/filtered_sorted_total_chroms.vcf' }
            //freebayes_vcf.mix(Channel.fromPath(params.vcf_mixed)))
            demultiplex_scSplit(bam, vcf, bai)
            scSplit_out = demultiplex_scSplit.out
        }
        else {
            scSplit_out = channel.value('no_result')
        }

        if (params.souporcell == 'True') {
            bam = params.souporcell_preprocess == 'True' ? qc_bam : input_bam // qc_bam.mix(input_bam))
            demultiplex_souporcell(bam)
            souporcell_out = demultiplex_souporcell.out
        }
        else {
            souporcell_out = channel.value('no_result')
        }

        summary(demuxlet_out, freemuxlet_out, vireo_out, souporcell_out, scSplit_out,
                params.generate_anndata, params.generate_mudata,
                params.rna_matrix_filtered, params.hto_matrix_filtered)
    emit:
        summary.out
}
