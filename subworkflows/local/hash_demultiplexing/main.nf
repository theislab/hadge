include { DROPLETUTILS_MTXCONVERT as MTXCONVERT_RNA } from '../../../modules/local/dropletutils/mtxconvert'
include { DROPLETUTILS_MTXCONVERT as MTXCONVERT_HTO } from '../../../modules/local/dropletutils/mtxconvert'
include { DEMUXEM                                   } from '../../../modules/nf-core/demuxem'
include { GMMDEMUX                                  } from '../../../modules/nf-core/gmmdemux'

workflow HASH_DEMULTIPLEXING {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    methods        // list of strings

    main:

    ch_versions = Channel.empty()

    if (methods.contains('htodemux')) {
        error("HtoDemux not implemented")
    }
    if (methods.contains('multiseq')) {
        error("MultiSeq not implemented")
    }
    if (methods.contains('cellhashr')) {
        error("CellHashR not implemented")
    }
    if (methods.contains('demuxem')) {
        ch_samplesheet.map { meta, rna, hto ->
            {
                if (!rna) {
                    error("RNA matrix not provided for sample ${meta.id}, but this is required for DemuxEM. Please check your input samplesheet.")
                }
                if (!hto) {
                    error("HTO matrix not provided for sample ${meta.id}, but this is required for DemuxEM. Please check your input samplesheet.")
                }
            }
        }

        MTXCONVERT_RNA(ch_samplesheet.map { meta, rna, _hto -> [meta, rna] }, false)
        ch_versions = ch_versions.mix(MTXCONVERT_RNA.out.versions)

        MTXCONVERT_HTO(ch_samplesheet.map { meta, _rna, hto -> [meta, hto] }, true)
        ch_versions = ch_versions.mix(MTXCONVERT_HTO.out.versions)

        DEMUXEM(
            MTXCONVERT_RNA.out.h5.join(MTXCONVERT_HTO.out.csv),
            params.demuxem_gender_genes,
            params.genome ?: [],
            true,
        )
        ch_versions = ch_versions.mix(DEMUXEM.out.versions)
    }
    if (methods.contains('gmm-demux')) {
        ch_gmmdemux = ch_samplesheet.map { meta, _rna, hto -> [meta, hto, "MS-11,MS-12", meta.n_cells] }

        ch_gmmdemux.map { meta, hto, hto_names, _estimated_cells ->
            {
                if (!hto) {
                    error("HTO matrix not provided for sample ${meta.id}, but this is required for GMM-Demux. Please check your input samplesheet.")
                }
                if (!hto_names) {
                    error("HTO names not provided for sample ${meta.id}, but this is required for GMM-Demux. Please check your input samplesheet.")
                }
            }
        }
        GMMDEMUX(
            ch_gmmdemux,
            true,
            true,
            [],
            [],
        )
        ch_versions = ch_versions.mix(GMMDEMUX.out.versions)
    }
    if (methods.contains('hasheddrops')) {
        error("HashedDrops not implemented")
    }
    if (methods.contains('hashsolo')) {
        error("HashSolo not implemented")
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
