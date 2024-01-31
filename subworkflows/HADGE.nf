nextflow.enable.dsl=2

include { run_multi } from "$projectDir/modules/multi_demultiplexing"
include {run_single} from "$projectDir/modules/single_demultiplexing"
include { summary } from "$projectDir/modules/multi/gene_demultiplexing"
include { donor_match } from "$projectDir/modules/multi/donor_match"

workflow HADGE {
    // Here we decide if it is a single sample demultiplexing or multi input demutliplexing run.
    if (params.multi_input == null){
        // Single Mode
        run_single()
    }
    else{
        // Multi mode
        run_multi()
    }
}



workflow SUMMARY{

    Channel.fromPath(params.multi_input) \
                | splitCsv(header:true) \
                | map { row-> tuple(row.sampleId, file(row.hto_matrix_filtered), file(row.rna_matrix_filtered))}
                | set {input_list_summary}
    log.info('running summary only')

    // demuxlet_out = Channel.fromPath('/lustre/scratch125/humgen/teams/hgi/ip13/1k1k/hadge/hedge_multi/hadge_modified/lustre/scratch125/humgen/teams/hgi/ip13/1k1k/hadge/full_run/gt_absent/result/*/genetic/gene_demulti/demuxlet/demuxlet*', type: 'dir').collect()
    demuxlet_out = channel.value("no_result")
    freemuxlet_out= Channel.fromPath('/lustre/scratch125/humgen/teams/hgi/ip13/1k1k/hadge/hedge_multi/hadge_modified/lustre/scratch125/humgen/teams/hgi/ip13/1k1k/hadge/full_run/gt_absent/result/*/genetic/gene_demulti/freemuxlet/freemuxlet_*', type: 'dir').collect()
    vireo_out= Channel.fromPath('/lustre/scratch125/humgen/teams/hgi/ip13/1k1k/hadge/hedge_multi/hadge_modified/lustre/scratch125/humgen/teams/hgi/ip13/1k1k/hadge/full_run/gt_absent/result/*/genetic/gene_demulti/vireo/vireo_*', type: 'dir').collect()
    scSplit_out= Channel.fromPath('/lustre/scratch125/humgen/teams/hgi/ip13/1k1k/hadge/hedge_multi/hadge_modified/lustre/scratch125/humgen/teams/hgi/ip13/1k1k/hadge/full_run/gt_absent/result/*/genetic/gene_demulti/scSplit/scsplit*', type: 'dir').collect()
    souporcell_out= Channel.fromPath('/lustre/scratch125/humgen/teams/hgi/ip13/1k1k/hadge/hedge_multi/hadge_modified/lustre/scratch125/humgen/teams/hgi/ip13/1k1k/hadge/full_run/gt_absent/result/*/genetic/gene_demulti/souporcell/souporcell_*', type: 'dir').collect()
    // /lustre/scratch125/humgen/teams/hgi/ip13/1k1k/hadge/hedge_multi/hadge_modified/lustre/scratch125/humgen/teams/hgi/ip13/1k1k/hadge/full_run/gt_absent/result/pool1/genetic/gene_demulti/souporcell/souporcell_pool1
    // souporcell_out = channel.value("no_result")
    summary(input_list_summary, demuxlet_out, freemuxlet_out, vireo_out, souporcell_out, scSplit_out, 
            params.generate_anndata, params.generate_mudata)


    Channel.fromPath(params.multi_input) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId, row.nsample, row.barcodes, "None", "None")}
        | join(summary.out)
        | donor_match

}