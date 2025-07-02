include { UNTAR as UNTAR_RNA                                       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_HTO                                } from '../../../modules/nf-core/untar'
include { RENAME_GENES_TO_FEATURES as RENAME_GENES_TO_FEATURES_RNA } from '../../../modules/local/rename_genes_to_features'
include { RENAME_GENES_TO_FEATURES as RENAME_GENES_TO_FEATURES_HTO } from '../../../modules/local/rename_genes_to_features'
include { DROPLETUTILS_MTXCONVERT as MTXCONVERT_RNA } from '../../../modules/local/dropletutils/mtxconvert'
include { DROPLETUTILS_MTXCONVERT as MTXCONVERT_HTO } from '../../../modules/local/dropletutils/mtxconvert'
include { PREPROCESSING_FOR_HTODEMUX_MULTISEQ       } from '../../../modules/local/preprocessing_for_htodemux_multiseq'
include { HTODEMUX                                  } from '../../../modules/nf-core/htodemux'
include { HTODEMUX_VISUALIZATION                    } from '../../../modules/local/htodemux_visualization'
include { MULTISEQDEMUX                             } from '../../../modules/nf-core/multiseqdemux'
include { DEMUXEM                                   } from '../../../modules/nf-core/demuxem'
include { GMMDEMUX                                  } from '../../../modules/nf-core/gmmdemux'
include { HASHEDDROPS                               } from '../../../modules/nf-core/hasheddrops'

workflow HASH_DEMULTIPLEXING {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    methods        // list of strings

    main:

    ch_versions = Channel.empty()

    ch_samplesheet.map { meta, rna, hto ->
        {
            if (!rna) {
                error("RNA matrix not provided for sample ${meta.id}, but this is required for hash demultiplexing. Please check your input samplesheet.")
            }
            if (!hto) {
                error("HTO matrix not provided for sample ${meta.id}, but this is required for hash demultiplexing. Please check your input samplesheet.")
            }
        }
    }

    ch_rna = ch_samplesheet.map { meta, rna, _hto -> 
        // add _rna to the id to prevent input file name collision of preprocessing and hasheddrops (both modules take two matrices as input)
        def new_meta = meta.clone()
        new_meta.id = "${meta.id}_rna"
        [new_meta, rna]
    }
    .branch { _meta, rna ->
        tar: rna.endsWith('.tar.gz')
        directory: true
    }

    ch_hto = ch_samplesheet.map { meta, _rna, hto -> 
        // add _hto to the id to prevent input file name collision of preprocessing and hasheddrops (both modules take two matrices as input)
        def new_meta = meta.clone()
        new_meta.id = "${meta.id}_hto"
        [new_meta, hto]
    }
    .branch { _meta, hto ->
        tar: hto.endsWith('.tar.gz')
        directory: true
    }

    UNTAR_RNA(ch_rna.tar)
    ch_versions = ch_versions.mix(UNTAR_RNA.out.versions)

    UNTAR_HTO(ch_hto.tar)
    ch_versions = ch_versions.mix(UNTAR_HTO.out.versions)


    // remove the changes to meta.id
    ch_rna = UNTAR_RNA.out.untar.map { meta, rna -> 
        def inital_id = meta.id.split("_")[0]
        [meta + [id: inital_id], rna]
    }
    ch_hto = UNTAR_HTO.out.untar.map { meta, hto -> 
        def inital_id = meta.id.split("_")[0]
        [meta + [id: inital_id], hto]
    }

    // rename genes.tsv to features.tsv to avoid Seurat 5.3 file missing error
    // ch_rna = RENAME_GENES_TO_FEATURES_RNA(ch_rna)
    // ch_hto = RENAME_GENES_TO_FEATURES_HTO(ch_hto)

    UNTAR_RNA.out.untar.view { "RNA untar output: ${it}" }
    UNTAR_HTO.out.untar.view { "HTO untar output: ${it}" }

    ch_samplesheet = ch_samplesheet.map { meta, _rna, _hto -> [meta] }.join(ch_rna).join(ch_hto)

    // ch_samplesheet.view { "Samplesheet input: ${it}" }

    if (methods.contains('htodemux') || methods.contains('multiseq')) {
        PREPROCESSING_FOR_HTODEMUX_MULTISEQ(
            ch_samplesheet
        )
        ch_versions = ch_versions.mix(PREPROCESSING_FOR_HTODEMUX_MULTISEQ.out.versions)

        if (methods.contains('htodemux')) {
            HTODEMUX(
                PREPROCESSING_FOR_HTODEMUX_MULTISEQ.out.seurat_object.map { meta, seurat_object -> [meta, seurat_object, "HTO"] }
            )
            ch_versions = ch_versions.mix(HTODEMUX.out.versions)

            HTODEMUX_VISUALIZATION(
                HTODEMUX.out.rds.map { meta, seurat_object -> [meta, seurat_object, "HTO"] }
            )
            ch_versions = ch_versions.mix(HTODEMUX_VISUALIZATION.out.versions)
        }
        if (methods.contains('multiseq')) {
            MULTISEQDEMUX(
                PREPROCESSING_FOR_HTODEMUX_MULTISEQ.out.seurat_object.map { meta, seurat_object -> [meta, seurat_object, "HTO"] }
            )
            ch_versions = ch_versions.mix(MULTISEQDEMUX.out.versions)
        }
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
        HASHEDDROPS(
            ch_samplesheet.map { meta, rna, hto -> 
            [meta, hto, "FALSE", rna] }
        )
        ch_versions = ch_versions.mix(HASHEDDROPS.out.versions)
    }
    if (methods.contains('hashsolo')) {
        error("HashSolo not implemented")
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
