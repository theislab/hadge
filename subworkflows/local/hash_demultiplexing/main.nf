include { DROPLETUTILS_MTXCONVERT as MTXCONVERT_RNA } from '../../../modules/local/dropletutils/mtxconvert'
include { DROPLETUTILS_MTXCONVERT as MTXCONVERT_HTO } from '../../../modules/local/dropletutils/mtxconvert'
include { DEMUXEM                                   } from '../../../modules/nf-core/demuxem'

workflow HASH_DEMULTIPLEXING {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    methods        // list of strings

    main:

    ch_versions = Channel.empty()

    MTXCONVERT_RNA(ch_samplesheet.map { meta, rna, _hto -> [meta, rna] }, false)
    ch_versions = ch_versions.mix(MTXCONVERT_RNA.out.versions)

    MTXCONVERT_HTO(ch_samplesheet.map { meta, _rna, hto -> [meta, hto] }, true)
    ch_versions = ch_versions.mix(MTXCONVERT_HTO.out.versions)

    if (methods.contains('htodemux')) {
    }
    if (methods.contains('multiseq')) {
    }
    if (methods.contains('cellhashr')) {
    }
    if (methods.contains('demuxem')) {
        DEMUXEM(
            MTXCONVERT_RNA.out.h5.join(MTXCONVERT_HTO.out.csv),
            true,
            [],
            true,
        )
        ch_versions = ch_versions.mix(DEMUXEM.out.versions)
    }
    if (methods.contains('gmm-demux')) {
    }
    if (methods.contains('hasheddrops')) {
    }
    if (methods.contains('hashsolo')) {
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
