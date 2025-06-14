include { SAMTOOLS_VIEW  } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index'
include { UMITOOLS_DEDUP } from '../../../modules/nf-core/umitools/dedup'
include { SAMTOOLS_SORT  } from '../../../modules/nf-core/samtools/sort'

workflow BAM_QC {
    take:
    ch_bam

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_VIEW(ch_bam.map { meta, bam -> [meta, bam, []] }, [[], []], [], 'bai')
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    SAMTOOLS_INDEX(SAMTOOLS_VIEW.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    UMITOOLS_DEDUP(SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_INDEX.out.bai), true)
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)

    SAMTOOLS_SORT(UMITOOLS_DEDUP.out.bam, [[], []])
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    emit:
    bam      = SAMTOOLS_SORT.out.bam
    versions = ch_versions // channel: [ versions.yml ]
}
