
VisiumTools
===========

This directory contains tools (scripts and notebooks) to prepare and/or analyze long-read data from the [scNaST Nanopore barcode/UMI assignment and transcript isoform quantification](https://github.com/dieterich-lab/ScNaST/workflow/sc-nanopore-utils/) workflow and obtained with the Visium Spatial Gene Expression molecular profiling protocol (10x Genomics).

These tools are typically used jointly with the **scNaST** assembly and transcript isoform quantification workflow and use the same configuration file, see `config.yaml` (parent directory).


Getting started
===============


## Dependencies

The Python3 scripts can be run in the same environment as the one created for [scNapbar](https://github.com/dieterich-lab/single-cell-nanopore),
however additional dependencies are required, such as `PyTables`, and can be installed via `pip` or `conda`.

However, since the analysis notebooks depend on many more packages, we recommend to create a separate python virtual environment *e.g.*


```bash
python3 -m venv /path/to/new/virtual/environment
source /path/to/new/virtual/environment/bin/activate
# install required packages
pip install -r environment.txt
```


## Usage

These scripts are currently not *installable* ( *i.e.* they are not a package ), and **MUST** be called from this directory ( *visiumtools* ). 
In addition, you might have to *e.g.* `chmod +x fmt_nanopore_visium.py`. This is not a clean pythonic solution, but only a temporary workaround
as we intend to eventually integrate these with the [scNapbar](https://github.com/dieterich-lab/single-cell-nanopore) pipeline.

**utils**

| script                       | description |
| -----------------------------|-------------|
| fmt_nanopore_visium.py       | Format scNaST output to use with Scanpy read_visium or Seurat Load10X_Spatial |


**notebooks**

| notebook          | description |
| ------------------|-------------|






