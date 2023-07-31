<!-- vim-markdown-toc GFM -->

* [Description](#description)
* [Genome annotation](#genome-annotation)
* [Explore annotation](#explore-annotation)
    * [Set up](#set-up)
    * [Run](#run)

<!-- vim-markdown-toc -->

Description
===========

The procedure described here annotates a genome assembly and then characterises
the "intermediate" chromosomes. In particular, we test whether the intermediate
chromosomes are enriched for any GO term or PFAM domain. We also ask whether
VSG genes tend to be present more densely on the intermediate chromosomes.

* See if the *de novo* genes are more or less represented in the intermediates

* Read length distribution of genes or proteins

Genome annotation
=================

The genome assembly listed in
(annotation_sample_sheet.tsv)[annotation_sample_sheet.tsv] was annotated with
the annotation pipeline at
https://github.com/glaParaBio/genomeAnnotationPipeline to detect protein coding
genes. The output required for the steps below is
`{genome_id}/hmmer/augustus.hints.gff3` (or `{genome_id}/hmmer/merge.gff3`).

The sample sheet required by the pipeline is [annotation_sample_sheet.tsv](annotation_sample_sheet.tsv).

Explore annotation
==================

Set up
------

Use bioconda and mamba to create a separate environment and install the project dependencies listed in requirements.txt:

```
conda config --set channel_priority strict
conda create -n intermediate_ctgs
conda activate intermediate_ctgs
mamba install -n intermediate_ctgs --yes --file requirements.txt
```

Alternatively, install the dependencies in requirements.txt manually.

Run
---

```
snakemake -p -n -j 5 --rerun-trigger mtime \
    --config gff3=/path/to/canu_tb427_ont.contigs_pilon_4x/hmmer/augustus.hints.gff3 \
             intercontigs=$PWD/data/canu_tb427_ont.contigs_pilon_4x.intermediate_contigs.txt \
             fai=$PWD/data/canu_tb427_ont.contigs_pilon_4x.fasta.fai \
    -d output
```

`/path/to/canu_tb427_ont.contigs_pilon_4x/hmmer/augustus.hints.gff3` is the
output of the genome annotation described above. The file
[canu_tb427_ont.contigs_pilon_4x.intermediate_contigs.txt](data/canu_tb427_ont.contigs_pilon_4x.intermediate_contigs.txt)
lists the "intermendiate contigs", one per line.
[canu_tb427_ont.contigs_pilon_4x.fasta.fai](data/canu_tb427_ont.contigs_pilon_4x.fasta.fai). 
