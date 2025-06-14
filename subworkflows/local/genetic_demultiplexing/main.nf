include { POPSCLE_DSCPILEUP } from '../../../modules/nf-core/popscle/dscpileup'
include { SAMTOOLS_INDEX    } from '../../../modules/nf-core/samtools/index'

workflow GENETIC_DEMULTIPLEXING {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    methods        // list of strings

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_INDEX(ch_samplesheet.map { meta, bam, _barcodes, _nsample, _vcf -> [meta, bam] })
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    if (methods.contains('demuxlet') || methods.contains('freemuxlet')) {
    }

    if (methods.contains('vireo')) {
        error("Vireo not implemented")
    }
    if (methods.contains('demuxlet')) {
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
