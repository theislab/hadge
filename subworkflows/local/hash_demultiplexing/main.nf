workflow HASH_DEMULTIPLEXING {

    take:
    methods // list of strings

    main:

    ch_versions = Channel.empty()

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
