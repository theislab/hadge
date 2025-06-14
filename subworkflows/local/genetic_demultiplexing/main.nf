include { BAM_QC             } from '../bam_qc'
include { FILTER_BAM         } from '../../../modules/local/filter_bam'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index'
include { CELLSNP_MODEA      } from '../../../modules/nf-core/cellsnp/modea'
include { VIREO              } from '../../../modules/nf-core/vireo'
include { POPSCLE_DSCPILEUP  } from '../../../modules/nf-core/popscle/dscpileup'
include { POPSCLE_DEMUXLET   } from '../../../modules/nf-core/popscle/demuxlet'
include { POPSCLE_FREEMUXLET } from '../../../modules/nf-core/popscle/freemuxlet'

workflow GENETIC_DEMULTIPLEXING {
    take:
    ch_samplesheet  // channel: samplesheet read in from --input
    methods         // list of strings
    bam_qc          // boolean
    common_variants // file

    main:

    ch_versions = Channel.empty()

    if (bam_qc) {
        BAM_QC(ch_samplesheet.map { meta, bam, _barcodes, _nsample, _vcf -> [meta, bam] })
        ch_versions = ch_versions.mix(BAM_QC.out.versions)

        ch_samplesheet = ch_samplesheet
            .join(BAM_QC.out.bam)
            .map { meta, _bam, barcodes, nsample, vcf, new_bam -> [meta, new_bam, barcodes, nsample, vcf] }
    }

    if (common_variants) {
        FILTER_BAM(
            ch_samplesheet.map { meta, bam, barcodes, _nsample, _vcf ->
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
            .map { meta, _bam, barcodes, nsample, vcf, new_bam -> [meta, new_bam, barcodes, nsample, vcf] }
    }

    if (methods.contains('demuxlet') || methods.contains('freemuxlet')) {
        POPSCLE_DSCPILEUP(ch_samplesheet.map { meta, bam, _barcodes, _nsample, vcf -> [meta, bam, vcf] })
        ch_versions = ch_versions.mix(POPSCLE_DSCPILEUP.out.versions)
    }

    if (methods.contains('vireo')) {
        SAMTOOLS_INDEX(ch_samplesheet.map { meta, bam, _barcodes, _nsample, _vcf -> [meta, bam] })
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        CELLSNP_MODEA(
            ch_samplesheet.join(SAMTOOLS_INDEX.out.bai).map { meta, bam, barcodes, _nsample, vcf, bai -> [meta, bam, bai, vcf, barcodes] }
        )
        ch_versions = ch_versions.mix(CELLSNP_MODEA.out.versions)

        VIREO(
            ch_samplesheet.join(CELLSNP_MODEA.out.cell).map { meta, _bam, _barcodes, nsample, vcf, cell -> [meta, cell, nsample, vcf, []] }
        )
        ch_versions = ch_versions.mix(VIREO.out.versions)
    }
    if (methods.contains('demuxlet')) {
        POPSCLE_DEMUXLET(
            POPSCLE_DSCPILEUP.out.plp.join(ch_samplesheet).map { meta, plp, bam, _barcodes, _nsample, vcf -> [meta, plp, bam, vcf] }
        )
        ch_versions = ch_versions.mix(POPSCLE_DEMUXLET.out.versions)
    }
    if (methods.contains('freemuxlet')) {
        POPSCLE_FREEMUXLET(
            POPSCLE_DSCPILEUP.out.directory.join(ch_samplesheet).map { meta, plp_dir, _bam, _barcodes, n_sample, _vcf -> [meta, plp_dir, n_sample] }
        )
        ch_versions = ch_versions.mix(POPSCLE_FREEMUXLET.out.versions)
    }
    if (methods.contains('souporcell')) {
        error("Souporcell not implemented")
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
