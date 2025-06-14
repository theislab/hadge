include { SAMTOOLS_INDEX    } from '../../../modules/nf-core/samtools/index'
include { POPSCLE_DSCPILEUP } from '../../../modules/nf-core/popscle/dscpileup'
include { POPSCLE_DEMUXLET  } from '../../../modules/nf-core/popscle/demuxlet'

workflow GENETIC_DEMULTIPLEXING {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    methods        // list of strings

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_INDEX(ch_samplesheet.map { meta, bam, _barcodes, _nsample, _vcf -> [meta, bam] })
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    if (methods.contains('demuxlet') || methods.contains('freemuxlet')) {
        POPSCLE_DSCPILEUP(ch_samplesheet.map { meta, bam, _barcodes, _nsample, vcf -> [meta, bam, vcf] })
        ch_versions = ch_versions.mix(POPSCLE_DSCPILEUP.out.versions)
    }

    if (methods.contains('vireo')) {
        error("Vireo not implemented")
    }
    if (methods.contains('demuxlet')) {
        POPSCLE_DEMUXLET(ch_samplesheet.map { meta, bam, _barcodes, _nsample, vcf -> [meta, [], bam, vcf] })
        ch_versions = ch_versions.mix(POPSCLE_DEMUXLET.out.versions)
    }
    if (methods.contains('freemuxlet')) {
        error("Freemuxlet not implemented")
    }
    if (methods.contains('souporcell')) {
        error("Souporcell not implemented")
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
