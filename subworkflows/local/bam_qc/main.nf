include { SAMTOOLS_VIEW  } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index'
include { UMITOOLS_DEDUP } from '../../../modules/nf-core/umitools/dedup'
include { SAMTOOLS_SORT  } from '../../../modules/nf-core/samtools/sort'

workflow BAM_QC {
    take:
    ch_bam

    main:

    SAMTOOLS_VIEW(
        ch_bam.map { meta, bam -> [meta, bam, []] },
        [[], [], []],
        [[], []],
        [[], []],
        'bai'
    )

    SAMTOOLS_INDEX(SAMTOOLS_VIEW.out.bam)

    UMITOOLS_DEDUP(SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_INDEX.out.bai), true)

    SAMTOOLS_SORT(UMITOOLS_DEDUP.out.bam, [[], [], []], '')

    emit:
    bam = SAMTOOLS_SORT.out.bam
}
