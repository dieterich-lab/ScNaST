
Nanopore barcode/UMI assignment - Transcript isoform quantification
===================================================================

Starting from long read alignments and barcode/UMI assignments from [scNapbar](https://github.com/dieterich-lab/single-cell-nanopore), the
main script `get_iso_mtx.py` integrates several steps for single-cell transcript isoform quantification, including *(i)* splitting the alignment BAM file into multiple BAM files, one per cell barcode, based on the barcode assignments, *(ii)* converting to FASTQ, and mapping to the transcriptome with [minimap2](https://lh3.github.io/minimap2/minimap2.html), *(iii)* quantifying abundances with [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) in alignment-based mode, and finally *(iv)* converting Salmon's estimate of the number of reads mapping to each transcript (NumReads) to a CellRanger-like output (feature-barcode matrices) for downstrean analyses.

This workflow is typically used jointly with the *assembly* workflow and uses the same configuration file, see `config.yaml`.


Getting started
===============

## Input

- The input alignments in BAM format (output from [scNapbar](https://github.com/dieterich-lab/single-cell-nanopore))
- The assignment file `real.label` from [scNapbar](https://github.com/dieterich-lab/single-cell-nanopore)
- The transcript sequences, transcriptome index, and GTF annotations (output from assembly workflow)

## Output

- `feature_bc_matrix` (one directory up from the BAM file into a directory called "analysis"), containing `features.tsv.gz`, `barcodes.tsv.gz`, and `matrix.mtx.gz`

Most intermediate files ( split BAM files, transcriptome alignments, *etc.* ) are removed by default. To keep them use `--keep-all`. Salmon output is kept by default, and must be removed by hand.


## Dependencies

The pipeline is run in the same environment created for [scNapbar](https://github.com/dieterich-lab/single-cell-nanopore).
Additional dependencies ([Salmon](https://salmon.readthedocs.io/en/latest/salmon.html)) are currently handled via the environment modules.

**Note:** [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) > 1.5.0 is required to use the `--ont` flag (set by default) for long-reads quantification. The *EffectiveLength* of transcripts is not used, and this field in the output is replaced by a fixed value of 100.


## Usage

These scripts are not *installable*. You might have to *e.g.* `chmod +x get_iso_mtx.py`.

Edit `config.yml`, and 

```bash
./get_iso_mtx.py ../config.yaml
```

and/or `runjob` to submit using [Slurm](https://slurm.schedmd.com/documentation.html)

```bash
sbatch runjob
```

## Issues

1. Salmon overrides the link to the C math standard library and causes havoc when using environment modules.

See [Salmon modifies the link to the 64-bit x86 C math library](https://github.com/COMBINE-lab/salmon/issues/710). Temporary solution: use `conda` or delete salmon's copy of *libm.so.6*.

2. Salmon long read error model

For less than 10% of reads, we get 

```
[jointLog] [warning] (in update()) CIGAR string for read [] seems inconsistent. It implied an error rate greater than 1
```

It was previously recommended to disable the default error model when quantifying alignments from long reads (`--noErrorModel --noLengthCorrection`).
From v1.5.0, the error model enabled with the `--ont` flag is designed for the alignment of long reads. This flag also disables the length effect in the generative model, but it is unclear if we still should pass the `--noErrorModel` flag in addition to the `--ont` flag?

Why do we get this error for some reads only? Does this arise from the long read error model? See also this [issue](https://github.com/COMBINE-lab/salmon/issues/289).



