workflow HASH_DEMULTIPLEXING {
    take:
    methods // list of strings

    main:

    ch_versions = Channel.empty()

    if (methods.contains('htodemux')) {
    }
    if (methods.contains('multiseq')) {
    }
    if (methods.contains('cellhashr')) {
    }
    if (methods.contains('demuxem')) {
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
