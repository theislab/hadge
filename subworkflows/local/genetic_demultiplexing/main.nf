workflow GENETIC_DEMULTIPLEXING {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    methods        // list of strings

    main:

    ch_versions = Channel.empty()

    if (methods.contains('vireo')) {
    }
    if (methods.contains('demuxlet')) {
    }
    if (methods.contains('freemuxlet')) {
    }
    if (methods.contains('souporcell')) {
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
