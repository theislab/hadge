# hadge

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

![Caption](docs/source/_static/images/pipeline.png)

Read the [documentation](https://hadge.readthedocs.io/en/latest/) and the associated [manuscript](https://www.biorxiv.org/content/10.1101/2023.07.23.550061v1)

## Quickstart

```bash
git clone https://github.com/theislab/hadge.git
sh hadge/test_data/download_data.sh
nextflow run main.nf -profile test
```

## Citation

hadge: a comprehensive pipeline for donor deconvolution in single cell

Fabiola Curion, Xichen Wu, Lukas Heumos, Mariana Gonzales Andre, Lennard Halle, Melissa Grant-Peters, Charlotte Rich-Griffin, Hing-Yuen Yeung, Calliope A. Dendrou, Herbert B. Schiller, Fabian J. Theis

bioRxiv 2023.07.23.550061; doi: https://doi.org/10.1101/2023.07.23.550061
