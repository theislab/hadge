include { BAM_QC             } from '../bam_qc'
include { FILTER_BAM         } from '../../../modules/local/filter_bam'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index'
include { CELLSNP_MODEA      } from '../../../modules/nf-core/cellsnp/modea'
include { VIREO              } from '../../../modules/nf-core/vireo'
include { POPSCLE_DSCPILEUP  } from '../../../modules/nf-core/popscle/dscpileup'
include { POPSCLE_DEMUXLET   } from '../../../modules/nf-core/popscle/demuxlet'
include { POPSCLE_FREEMUXLET } from '../../../modules/nf-core/popscle/freemuxlet'
include { SOUPORCELL         } from '../../../modules/nf-core/souporcell/main'
include { GENE_SUMMARY       } from '../../../modules/local/gene_summary'

workflow GENETIC_DEMULTIPLEXING {
    take:
    ch_samplesheet  // channel: samplesheet read in from --input
    methods         // list of strings
    bam_qc          // boolean
    common_variants // file

    main:

    ch_versions = Channel.empty()

    ch_vireo_donor_id = Channel.empty()
    ch_vireo_summary = Channel.empty()
    ch_demuxlet_result = Channel.empty()
    ch_freemuxlet_result = Channel.empty()
    ch_souporcell = Channel.empty()

    if (bam_qc) {
        BAM_QC(ch_samplesheet.map { meta, bam, _barcodes, _vcf -> [meta, bam] })
        ch_versions = ch_versions.mix(BAM_QC.out.versions)

        ch_samplesheet = ch_samplesheet
            .join(BAM_QC.out.bam)
            .map { meta, _bam, barcodes, vcf, new_bam -> [meta, new_bam, barcodes, vcf] }
    }

    if (common_variants) {
        FILTER_BAM(
            ch_samplesheet.map { meta, bam, barcodes, _vcf ->
                [
                    meta,
                    bam,
                    barcodes,
                ]
            },
            common_variants,
        )
        ch_versions = ch_versions.mix(FILTER_BAM.out.versions)

        ch_samplesheet = ch_samplesheet
            .join(FILTER_BAM.out.bam)
            .map { meta, _bam, barcodes, vcf, new_bam -> [meta, new_bam, barcodes, vcf] }
    }


    if (methods.contains('vireo')) {
        SAMTOOLS_INDEX(ch_samplesheet.map { meta, bam, _barcodes, _vcf -> [meta, bam] })
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        CELLSNP_MODEA(
            ch_samplesheet.join(SAMTOOLS_INDEX.out.bai).map { meta, bam, barcodes, vcf, bai -> [meta, bam, bai, vcf, barcodes] }
        )
        ch_versions = ch_versions.mix(CELLSNP_MODEA.out.versions)

        VIREO(
            ch_samplesheet.join(CELLSNP_MODEA.out.cell).map { meta, _bam, _barcodes, vcf, cell -> [meta, cell, meta.n_samples, vcf, []] }
        )

        ch_vireo_donor_id = ch_vireo_donor_id.mix(VIREO.out.donor_ids)
        ch_vireo_summary = ch_vireo_summary.mix(VIREO.out.summary)
        ch_versions = ch_versions.mix(VIREO.out.versions)
    }
    if (methods.contains('demuxlet') || methods.contains('freemuxlet')) {
        ch_dscpileup = ch_samplesheet.map { meta, bam, _barcodes, vcf -> [meta, bam, vcf] }
        POPSCLE_DSCPILEUP(ch_dscpileup)
        ch_versions = ch_versions.mix(POPSCLE_DSCPILEUP.out.versions)

        if (methods.contains('demuxlet')) {
            ch_demuxlet = POPSCLE_DSCPILEUP.out.plp.join(ch_samplesheet).map { meta, plp, bam, _barcodes, vcf -> [meta, plp, bam, vcf] }
            POPSCLE_DEMUXLET(ch_demuxlet)
            ch_demuxlet_result = ch_demuxlet_result.mix(POPSCLE_DEMUXLET.out.demuxlet_result)
            ch_versions = ch_versions.mix(POPSCLE_DEMUXLET.out.versions)
        }
        if (methods.contains('freemuxlet')) {
            ch_freemuxlet = POPSCLE_DSCPILEUP.out.directory.join(ch_samplesheet).map { meta, plp_dir, _bam, _barcodes, _vcf -> [meta, plp_dir, meta.n_samples] }
            POPSCLE_FREEMUXLET(ch_freemuxlet)
            ch_freemuxlet_result = ch_freemuxlet_result.mix(POPSCLE_FREEMUXLET.out.result)
            ch_versions = ch_versions.mix(POPSCLE_FREEMUXLET.out.versions)
        }
    }

    if (methods.contains('souporcell')) {
        ch_soup_bam_barcodes = ch_samplesheet.map { meta, bam, barcodes, _vcf ->
        [ meta, bam, barcodes ]
        }

        ch_soup_fasta = ch_samplesheet.map { meta, _bam, _barcodes, _vcf ->
            [ meta, file(params.ref) ]
        }

        ch_soup_clusters = ch_samplesheet.map { meta, _bam, _barcodes, _vcf ->
            (params.souporcell_k ?: meta.n_samples)
        }

        SOUPORCELL(
            ch_soup_bam_barcodes,
            ch_soup_fasta,
            ch_soup_clusters
        )

        ch_souporcell = ch_souporcell.mix(SOUPORCELL.out.tsv)
        ch_versions = ch_versions.mix(SOUPORCELL.out.versions)
    }

    ch_summary = ch_samplesheet
        .join(ch_vireo_donor_id, remainder: true)
        .join(ch_vireo_summary, remainder: true)
        .join(ch_demuxlet_result, remainder: true)
        .join(ch_freemuxlet_result, remainder: true)
        .join(ch_souporcell, remainder: true)
        .map { tuple -> tuple.collect { it == null ? [] : it } }

    ch_summary.view()

    GENE_SUMMARY(
        ch_summary.map { meta, bam, barcodes, vcf,
                                      vireo_ids, vireo_sum,
                                      demuxlet_res, freemuxlet_res, souporcell_tsv ->
                        // If you also want .h5ad, pass a 10x RNA matrix dir here; else [].
                        def rna_mtx = []
                        [ meta, barcodes, vireo_ids, vireo_sum, demuxlet_res, freemuxlet_res, souporcell_tsv, rna_mtx ]
        }
    )

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
