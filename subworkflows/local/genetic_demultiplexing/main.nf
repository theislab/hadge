include { BAM_QC             } from '../bam_qc'
include { FILTER_BAM         } from '../../../modules/local/filter_bam'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index'
include { CELLSNP_MODEA      } from '../../../modules/nf-core/cellsnp/modea'
include { VIREO              } from '../../../modules/nf-core/vireo'
include { POPSCLE_DSCPILEUP  } from '../../../modules/nf-core/popscle/dscpileup'
include { POPSCLE_DEMUXLET   } from '../../../modules/nf-core/popscle/demuxlet'
include { POPSCLE_FREEMUXLET } from '../../../modules/nf-core/popscle/freemuxlet'
include { GENE_SUMMARY       } from '../../../modules/local/gene_summary'

workflow GENETIC_DEMULTIPLEXING {
    take:
    ch_samplesheet  // channel: samplesheet read in from --input
    methods         // list of strings
    bam_qc          // boolean
    common_variants // file

    main:

    ch_versions = Channel.empty()

    ch_vireo = Channel.empty()
    ch_demuxlet = Channel.empty()
    ch_freemuxlet = Channel.empty()
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
        ch_vireo = ch_vireo.mix(VIREO.out.donor_ids)
        ch_versions = ch_versions.mix(VIREO.out.versions)


    }
    if (methods.contains('demuxlet') || methods.contains('freemuxlet')) {
        ch_dscpileup = ch_samplesheet.map { meta, bam, _barcodes, vcf -> [meta, bam, vcf] }
        POPSCLE_DSCPILEUP(ch_dscpileup)
        ch_versions = ch_versions.mix(POPSCLE_DSCPILEUP.out.versions)

        if (methods.contains('demuxlet')) {
            ch_demuxlet_input = POPSCLE_DSCPILEUP.out.plp.join(ch_samplesheet).map { meta, plp, bam, _barcodes, vcf -> [meta, plp, bam, vcf] }
            POPSCLE_DEMUXLET(ch_demuxlet_input)
            ch_demuxlet = ch_demuxlet.mix(POPSCLE_DEMUXLET.out.demuxlet_result)
            ch_versions = ch_versions.mix(POPSCLE_DEMUXLET.out.versions)

        }
        if (methods.contains('freemuxlet')) {
            ch_freemuxlet_input = POPSCLE_DSCPILEUP.out.directory.join(ch_samplesheet).map { meta, plp_dir, _bam, _barcodes, _vcf -> [meta, plp_dir, meta.n_samples] }
            POPSCLE_FREEMUXLET(ch_freemuxlet_input)
            ch_freemuxlet = ch_freemuxlet.mix(POPSCLE_FREEMUXLET.out.result)
            ch_versions = ch_versions.mix(POPSCLE_FREEMUXLET.out.versions)
        }
    }

    if (methods.contains('souporcell')) {
        error("Souporcell not implemented")
    }

    ch_summary = ch_samplesheet
        .join(ch_vireo, remainder: true)
        .join(ch_demuxlet, remainder: true)
        .join(ch_freemuxlet, remainder: true)
        .map { tuple -> tuple.collect { it == null ? [] : it } }

    GENE_SUMMARY(ch_summary)
    ch_versions = ch_versions.mix(GENE_SUMMARY.out.versions)


    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
