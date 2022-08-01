
# ScNaST
ScNaST (Single-cell Nanopore Spatial Transcriptomics)

This repository contains *(i)* a pipeline for the annotation of genomes using long read transcriptomics data ([assembly](workflow/assembly/)), *(ii)* a pipeline to generate single-cell feature-count matrices at the isoform level ([sc-nanopore-utils](workflow/sc-nanopore-utils/)), and *(iii)* a set of tools and notebooks to facilitate the integration of long read data obtained with the Visium Spatial Gene Expression molecular profiling 10x Genomics protocol ([visiumTools](workflow/visiumtools/)).


## Supplementary resources 

This repository contains in addition supplementary material for

> Boileau E, Li X, Naarmann-de Vries IS, Becker C, Casper R, Altmüller J, Leuschner F and Dieterich C (2022)
> Full-Length Spatial Transcriptomics Reveals the Unexplored Isoform Diversity of the Myocardium Post-MI.
> Front. Genet. 13:912572. doi: [10.3389/fgene.2022.912572](https://www.frontiersin.org/articles/10.3389/fgene.2022.912572/full) .
> Research Topic: Advances in Single-Cell and Spatial Analyses Using Nanopore Sequencing


## How to use this repository

### Installation

This repository is not *installable* ( *i.e.* it is not a package ). We intend to eventually integrate these tools with the [scNapbar](https://github.com/dieterich-lab/single-cell-nanopore) workflow. Currently, both the [assembly](workflow/assembly/) and the [sc-nanopore-utils](workflow/sc-nanopore-utils/) quantification workflow can be run in the same environment as the one created for [scNapbar](https://github.com/dieterich-lab/single-cell-nanopore), but may
require additional dependencies (see details for each pipeline). To use the [visiumTools](workflow/visiumtools/), we recommend to install additional packages in a separate virtual environment.

Clone the pipeline toolset:

```bash
git clone https://github.com/dieterich-lab/ScNaST.git
```

### Dependencies

- See [scNapbar](https://github.com/dieterich-lab/single-cell-nanopore).
- Additional dependencies are currently handled via the environment modules, see [assembly](workflow/assembly/) and [sc-nanopore-utils](workflow/sc-nanopore-utils/).
- In addition, for the analysis of Visium Spatial Gene Expression (Nanopore), see [visiumTools](workflow/visiumtools/).


## Content

### workflow

The pipeline toolset: [assembly](workflow/assembly/), [sc-nanopore-utils](workflow/sc-nanopore-utils/), and [visiumTools](workflow/visiumtools/).
To run each workflow, descend into the respective directory and follow the instructions.

**Note:** All paths and options are specified in the `config.yaml` file for all tools/workflows! The actual location of this file is irrelevant; you just 
need to make sure it is referenced correctly when calling a script.

### paper

Analysis and plotting scripts, additional results used in the manuscript. Analysis notebooks are under [visiumTools](workflow/visiumtools/).

Additional data is available at

> Boileau, Etienne, Li, Xue, Naarmann-de Vries, Isabel, Becker, Christian, Casper, Ramona, Altmüller, Janine, Leuschner, Florian, & Dieterich, Christoph. (2022).
> Full-Length Spatial Transcriptomics Reveals the Unexplored Isoform Diversity of the Myocardium Post-MI [Data set]. 
> Zenodo. [https://doi.org/10.5281/zenodo.6546361](https://doi.org/10.5281/zenodo.6546361)

