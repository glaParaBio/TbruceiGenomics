<!-- vim-markdown-toc GFM -->

* [Description](#description)
* [Setup](#setup)
* [Run](#run)

<!-- vim-markdown-toc -->

Description
===========

This workflow detects base modification in one or more samples sequenced with
Nanopore technology. The workhorse of the procedure is the
[tombo](https://github.com/nanoporetech/tombo) suite.

Setup
=====

Use [bioconda](https://bioconda.github.io/user/install.html) and [mamba](https://github.com/mamba-org/mamba) to
create a separate environment and install the project dependencies listed in
`requirements.txt`:

```
mamba create --yes -n base_mods
mamba activate base_mods
# See https://stackoverflow.com/questions/70958434/unexpected-python-paths-in-conda-environment/70961159
export PYTHONNOUSERSITE=1
mamba install -n base_mods --yes --file requirements.txt

# See https://www.biostars.org/p/9504758/ for the reason to deactivate/activate 
mamba deactivate
mamba activate base_mods
export PYTHONNOUSERSITE=1
```

Alternatively, install the dependencies in `requirements.txt` manually.

Run
===

```
mamba activate base_mods
export PYTHONNOUSERSITE=1

snakemake -p -n -j 10 --rerun-trigger mtime -C \
    sample_sheet=$PWD/sample_sheet.tsv \
    genomes=$PWD/genomes.tsv \
        -d output
```

[sample_sheet.tsv](sample_sheet.tsv) is a tabular, tab-separeted file with the information of libraries and files to analyse. Columns are:


* `library_id`: Name of the library

* `control_id`: The control library for this `library_id`

* `experiment_id`: Identifier of the nanopore experiement

* `barcode`: Barcode of this library

* `run_dir`: Full path to the Nonopore run directory

* `fast5_dir`: Path of fast5 directory relative to `run_dir`

* `sequencing_summary`: Path to the sequencing summary file relative to `run_dir`. This file is typically produced by guppy basecaller 

* `fastq_dir`: Path of fastq directory relative to `run_dir` 

* `genome`: Identifier of the reference genome (matching the corresponding column in `genome.tsv`)

[genomes.tsv](genomes.tsv) is a tabular, tab-separated file with the information of the reference genome(s). Columns are:

* `genome`: Identifier of the reference genome (matching the corresponding column in `sample_sheet.tsv`)

* `fasta`: Full path or URL to the reference genome
