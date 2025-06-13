include { DEMUXEM } from '../../../modules/nf-core/demuxem'

workflow HASH_DEMULTIPLEXING {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    methods        // list of strings

    main:

    ch_samplesheet.view()

    ch_versions = Channel.empty()

    if (methods.contains('htodemux')) {
    }
    if (methods.contains('multiseq')) {
    }
    if (methods.contains('cellhashr')) {
    }
    if (methods.contains('demuxem')) {
        DEMUXEM(ch_samplesheet, "test", true, [], true)
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
