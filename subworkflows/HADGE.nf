nextflow.enable.dsl=2

include { run_multi } from "$projectDir/modules/multi_demultiplexing"
include { summary } from "$projectDir/modules/multi/gene_demultiplexing"
include { donor_match } from "$projectDir/modules/multi/donor_match"
include {create_single_chanel_input}  from "$projectDir/modules/multi/preprocessing/preprocessing"


workflow HADGE {
    // Here we decide if it is a single sample demultiplexing or multi input demutliplexing run.
    if (params.multi_input == null){
        // Single Mode
        // Instead of running in single here we want to run a multi so that the chanels, workflows and processes are not redundant.
        create_single_chanel_input(
            params.sample_name,
            params.hto_matrix_raw,
            params.hto_matrix_filtered,
            params.rna_matrix_raw,
            params.rna_matrix_filtered,
            params.bam,
            params.bai,
            params.barcodes,
            params.fasta,
            params.fasta_index,
            params.nsamples_genetic,
            params.cell_data,
            params.vcf_mixed,
            params.vcf_donor,
            params.vireo_parent_dir,
            params.demultiplexing_result,
        )
        run_multi(create_single_chanel_input.out.input_channel)

    }
    else{
        // Multi mode
        input_channel = Channel.fromPath(params.multi_input)
        run_multi(input_channel)
    }
}



workflow SUMMARY{


    log.info('running summary only')

    Channel.fromPath(params.multi_input).splitCsv(header:true).map { row-> tuple(row.sampleId, file(row.hto_matrix_filtered), file(row.rna_matrix_filtered))}.set {input_list_summary}
    

    demuxlet_out = Channel.fromPath("${params.outdir}/*/genetic/gene_demulti/demuxlet/demuxlet_*", type: 'dir').collect().ifEmpty('no_result')
    freemuxlet_out= Channel.fromPath("${params.outdir}/*/genetic/gene_demulti/freemuxlet/freemuxlet_*", type: 'dir').collect().ifEmpty('no_result')
    vireo_out= Channel.fromPath("${params.outdir}/*/genetic/gene_demulti/vireo/vireo_*", type: 'dir').collect().ifEmpty('no_result')
    scSplit_out= Channel.fromPath("${params.outdir}/*/genetic/gene_demulti/scSplit/scsplit*", type: 'dir').collect().ifEmpty('no_result')
    souporcell_out= Channel.fromPath("${params.outdir}/*/genetic/gene_demulti/souporcell/souporcell_*", type: 'dir').collect().ifEmpty('no_result')

    demuxlet_out_ch = demuxlet_out.flatten().map{r1-> tuple(    "$r1".replaceAll(".*demuxlet_",""), r1 )}
    freemuxlet_out_ch = freemuxlet_out.flatten().map{r1-> tuple(    "$r1".replaceAll(".*freemuxlet_",""), r1 )}
    vireo_out_ch = vireo_out.flatten().map{r1-> tuple(    "$r1".replaceAll(".*vireo_",""), r1 )}
    scSplit_out_ch = scSplit_out.flatten().map{r1-> tuple(    "$r1".replaceAll(".*scsplit_",""), r1 )}
    souporcell_out_ch = souporcell_out.flatten().map{r1-> tuple(    "$r1".replaceAll(".*souporcell_",""), r1 )}

    summary_input = input_list_summary.join(souporcell_out_ch,by:0,remainder: true).join(scSplit_out_ch,by:0,remainder: true).join(vireo_out_ch,by:0,remainder: true).join(freemuxlet_out_ch,by:0,remainder: true).join(demuxlet_out_ch,by:0,remainder: true)
    summary_input = summary_input.filter{ it[0] != 'no_result' }

    summary(summary_input,
            params.generate_anndata, params.generate_mudata)

    // Channel.fromPath(params.multi_input) \
    //     | splitCsv(header:true) \
    //     | map { row-> tuple(row.sampleId, row.nsample, row.barcodes, "None", "None")}
    //     | join(summary.out)
    //     | donor_match

}