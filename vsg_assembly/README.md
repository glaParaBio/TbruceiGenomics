<!-- vim-markdown-toc GFM -->

* [Description](#description)
* [Set up](#set-up)
* [Run](#run)

<!-- vim-markdown-toc -->

Description
===========

This pipeline identifies VSG in one or more assemblies of *T. brucei* genome
using known coding sequences of VSG. Known VSGs from
[TriTrypDB-47_TbruceiLister427_2018_AnnotatedCDSs.fasta](https://tritrypdb.org/common/downloads/release-47/TbruceiLister427_2018)
are selected by parsing the fasta header for appropriate keywords (i.e.
*"variant surface glycoprotein"*, see rule `extract_vsg` in the
[Snakefile](Snakefile)).

These reference sequences are blasted against the genome(s) to annotate. We
retain blast hits with 90% identity and where the coverage of the query VSG is
at least 95%; see rule `blast` for details of the blast command.

Finally, the output of blast is converted to bed format (rule `blastToBed`).
These are some example lines of the output:

```
#sseqid    | m.sstart | m.send | qseqid                                                                                     | N
contig_102 |     2216 |   3676 | Tb427_000556700.1                                                                          | 1
contig_102 |    10074 |  11633 | Tb427_000155900.1,Tb427_000556900.1,Tb427_000574900.1                                      | 3
contig_102 |    14180 |  15713 | Tb427_000078500.1,Tb427_000078600.1,Tb427_000155700.1,Tb427_000248600.1,Tb427_000574200.1  | 5
contig_102 |    17066 |  18472 | Tb427_000155600.1,Tb427_000574100.1                                                        | 2
contig_102 |    23382 |  24496 | Tb427_000558000.1,Tb427_000573900.1                                                        | 2
contig_102 |    27576 |  28962 | Tb427_000573700.1                                                                          | 1
contig_102 |    30218 |  31720 | Tb427_000077500.1,Tb427_000128300.1,Tb427_000154800.1,Tb427_000573600.1                    | 4
```

The first three columns are the coordinates of the detected VSG on the
assembly. In the 4th and 5th column there are, respectively, the reference VSG
genes with a blast hit and the count of these genes.

Set up
======

Use [bioconda](https://bioconda.github.io/user/install.html) and [mamba](https://github.com/mamba-org/mamba) to
create a separate environment and install the project dependencies listed in
`requirements.txt`:

```
mamba create --yes -n 20200508_marija_vsg_assembly
mamba activate 20200508_marija_vsg_assembly
mamba install --freeze-installed -n 20200508_marija_vsg_assembly --yes --file requirements.txt
```

Alternatively, install the dependencies in `requirements.txt` manually.

Run 
===

Edit `sample_sheet.tsv` with the details of the assemblies to annotate. This is
a tab-separated file with columns:

* `assembly`: Name of your choice for the assembly - must be unique within the
  sample sheet

* `fasta`: Full path to genome fasta file to annotate

Additional columns are ignored. It is strongly recommended to avoid names and
file paths containing spaces and special characters.

To execute the pipeline:

```
cd /path/to/Snakefile
snakemake -p -n -C ss="$PWD/sample_sheet.tsv" -j 10 -d output
```

`-d` is a working directory of your choice where the main output will go. With
`-n` option, execute in dry-run mode (and do nothing).

The main output files are `blast/*.vsg.mrg.bed`, one file per assembly.

If necessary, edit the [Snakefile](Snakefile) to tweak the parameters of the analysis.
