<!-- vim-markdown-toc GFM -->

* [Description](#description)
* [Demultiplex Nanopore run](#demultiplex-nanopore-run)
* [Setup](#setup)
* [Run](#run)

<!-- vim-markdown-toc -->

Description
===========

This workflow detects base modification in one or more samples sequenced with
Nanopore technology. The workhorse of the procedure is the
[tombo](https://github.com/nanoporetech/tombo) suite.

Demultiplex Nanopore run
========================

We demultiplex the Nanopore runs to produce the sequencing summary table
indicating which read belongs to each barcode (*i.e.* sample). We also output
fast5 files and use these as input for tombo. (Note that tombo does not use
fastq files. The pipeline will split fast5 files to one read per file as
required by tombo). This is an example command for basecalling and
demultiplexing:

```
guppydir='/export/projects/III-data/wcmp_bioinformatics/db291g/applications/ont-guppy/ont-guppy-cpu_6.0.1_linux64'

cd /export/projects/III-data/wcmp_bioinformatics/db291g/data/20211218_marija_basej/20211213_1616_MN23371_FAQ95558_218cbaa2

nohup $guppydir/bin/guppy_basecaller \
    -c $guppydir/data/dna_r9.4.1_450bps_fast.cfg \
    -i . \
    -s ./guppy \
    --barcode_kits 'EXP-NBD104' \
    --num_callers 50 \
    --recursive \
    --compress_fastq \
    --fast5_out \
    --disable_qscore_filtering \
    --trim_barcode &
```

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
