workflow DONOR_MATCHING {
    main:

    ch_versions = Channel.empty()

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
