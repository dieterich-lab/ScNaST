
Assembly workflow for long read data 
====================================

This `snakemake` pipeline generates a GTF annotation file from existing long read alignments. 
It is used to process the output from the [scNapbar](https://github.com/dieterich-lab/single-cell-nanopore) workflow.
Samples (path to alignment in BAM format) and options are specified in a configuration file. 
For each sample, the alignment file is processed by [StringTie](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) in long read mode (with or without a a reference annotation file to guide the assembly process) to generate a GTF annotation. 
The annotations are merged into a non-redundant set of transcripts and compared with a reference annotation (if available) using [GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml).

If running all rules, single-exon transcripts are removed from the annotation, transcript sequences are extracted using [GffRead](http://ccb.jhu.edu/software/stringtie/gff.shtml), and a [minimap2](https://lh3.github.io/minimap2/minimap2.html) index is generated (for transcriptome mapping).


Currently no filtering is done on the BAM files prior to calling [StringTie](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual), but it is possible to specify custom thresholds ( *e.g* isoform abundance, minimum read coverage *etc.* ), see `config.yaml`. Options are not checked for consistency.


Getting started
===============

## Input

- The input alignments in BAM format (output from [scNapbar](https://github.com/dieterich-lab/single-cell-nanopore))
- The input genome in fasta format
- The existing annotations in GFF3 or GFF2 format

## Output

- `gffcompare/gffcmp.multi_exons.annotated.gtf` under the directory specified by `strg_base`
- mmi index (minimap2) under the directory specified by `idx_base` (with name given by `idx_name`)


## Dependencies

The pipeline is run in the same environment created for [scNapbar](https://github.com/dieterich-lab/single-cell-nanopore).
Additional dependencies [StringTie](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual), [GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml), [GffRead](http://ccb.jhu.edu/software/stringtie/gff.shtml) are currently handled via the environment modules.


## Usage

Edit `config.yml`, and 

```bash
snakemake --configfile ../config.yaml -j <num_cores> all --use-envmodules
```

or edit `config.json` and/or `runjob` to submit using [Slurm](https://slurm.schedmd.com/documentation.html)

```bash
sbatch runjob
```

