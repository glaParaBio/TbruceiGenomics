#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(ggbeeswarm))

parser <- ArgumentParser()
parser$add_argument('--gff', help='GFF file annotated with Pfam domains')
parser$add_argument('--intercontigs', help='File listing the intermediate contigs to test for enrichment')
parser$add_argument('--contigsize', help='File with contig names (1st column) and contig sizes (2nd column). E.g. the fai index')
parser$add_argument('--plotout', help='Write plot to this file')
parser$add_argument('--tblout', help='Write table of domain counts to this file')

xargs <- parser$parse_args()

# xargs <- list(gff='canu_tb427_ont.contigs_pilon_4x/hmmer/augustus.hints.gff3', 
#    intercontigs='/export/III-data/wcmp_bioinformatics/db291g/git_repos/wcip/20220909_marija_annotation/data/canu_tb427_ont.contigs_pilon_4x.intermediate_contigs.txt', 
#    contigsize='ref/canu_tb427_ont.contigs_pilon_4x.fasta.fai')

get_gff_attr_ <- function(x, key, split=NULL) {
    xattr <- strsplit(x, ';')[[1]]
    xkey <- sprintf('^%s=', key)
    value <- grep(xkey, xattr, value=TRUE)
    if(length(value) == 0) {
        return(NA)
    }
    stopifnot(length(value) == 1)
    value <- sub(xkey, '', value)
    if(!is.null(split)) {
        value <- strsplit(value, split)[[1]]
    }
    return(value)
}

get_gff_attr <- function(x, key, split=NULL) {
    sapply(x, get_gff_attr_, key=key, split=split, USE.NAMES=FALSE)
}

VSG_DOMS <- c('Trypanosomal VSG domain', 'Trypanosome variant surface glycoprotein (A-type)', 'Trypanosome variant surface glycoprotein C-terminal domain')

gff <- fread(cmd=sprintf('grep -v "^#" %s', xargs$gff))
inter <- fread(xargs$intercontigs, header=FALSE, col.name='contig')
contig_size <- fread(xargs$contigsize, header=FALSE, col.names=c('contig', 'size'), select=1:2)

pfam <- gff[V3 == 'protein_match']
pfam[, txid := get_gff_attr(V9, 'Parent')]
pfam[, signature_desc := get_gff_attr(V9, 'signature_desc')]
pfam[, is_vsg := signature_desc %in% VSG_DOMS]

mrna <- gff[V3 == 'mRNA', list(contig=V1, start=V4-1, end=V5, strand=V7, txid=get_gff_attr(V9, 'ID'), gene_id=get_gff_attr(V9, 'Parent'))]
# Only keep mrna with Pfam match
mrna <- mrna[txid %in% pfam$txid]
mrna[, is_vsg := txid %in% pfam[is_vsg == TRUE]$txid]

genes <- mrna[, list(is_vsg=sum(is_vsg) > 0), by=list(contig, gene_id)]

contig_smry <- genes[, list(n_genes=.N, n_vsg=sum(is_vsg)), by=contig]
contig_smry <- merge(contig_smry, contig_size, by='contig')
contig_smry[, genes_100kb := 100000 * n_genes/size]
contig_smry[, vsg_100kb := 100000 * n_vsg/size]
contig_smry[, is_intermediate := contig %in% inter$contig]

fit.vsg <- glm.nb(n_vsg ~ is_intermediate + offset(log(size)), data=contig_smry)
fit.genes <- glm.nb(n_genes ~ is_intermediate + offset(log(size)), data=contig_smry)
fit.vsg.geneoffset <- glm.nb(n_vsg ~ is_intermediate + offset(log(n_genes)), data=contig_smry)

dat <- melt(contig_smry, id.vars=c('contig', 'is_intermediate'), measure.vars=c('genes_100kb', 'vsg_100kb'), variable.name='genes', value.name='genes_100kb')
dat[, contig_type := ifelse(is_intermediate, 'Intermediate', 'Others')]
dat[, genes := ifelse(genes == 'genes_100kb', 'All genes', ifelse(genes == 'vsg_100kb', 'VSG', NA))]

gg <- ggplot(data=dat, aes(x=contig_type, y=genes_100kb)) +
    geom_point(data=dat[, list(genes_100kb=median(genes_100kb)), by=list(genes, contig_type)], pch='-', size=15, colour='red') +
    geom_quasirandom(size=0.5) +
    ylab('Genes per 100kb') +
    xlab('Contig type') +
    facet_wrap(~genes) +
    theme_light() +
    theme(strip.text=element_text(colour='black'))
ggsave(xargs$plotout, height=8, width=14, units='cm')

pfam[, signature_desc := get_gff_attr(V9, 'signature_desc')]
pfam[, name := get_gff_attr(V9, 'Name')]
pfam <- merge(pfam[, list(txid, signature_desc, name)], mrna[, list(txid, contig, gene_id)], by='txid')
pfam[, is_intermediate := contig %in% inter$contig]
gene_pfam <- unique(pfam[, list(signature_desc, name, contig, gene_id, is_intermediate)])

all_intermediate <- length(unique(gene_pfam[is_intermediate == TRUE]$gene_id))
all_nointermediate <- length(unique(gene_pfam[is_intermediate == FALSE]$gene_id))
pfam_cnt <- gene_pfam[, list(intrmd_dom=sum(is_intermediate), not_intrmd_dom=sum(is_intermediate == FALSE)), list(name, signature_desc)]
pfam_cnt[, intrmd_not_dom := all_intermediate - intrmd_dom]
pfam_cnt[, not_intrmd_not_dom := all_nointermediate - not_intrmd_dom]

ftlist <- function(a, b, c, d) {
    m <- matrix(c(a, b, c, d), nrow=2, byrow=TRUE)
    ft <- fisher.test(m, conf.int=FALSE)
    if(ft$estimate == 0) {
        ft$estimate <- 1/sum(m)
    }
    if(ft$estimate == Inf) {
        ft$estimate <- sum(m)
    }
    return(list(fisher.pval=ft$p.value, log2_odds_ratio=log2(ft$estimate)))
}

fttab <- pfam_cnt[, ftlist(intrmd_dom, not_intrmd_dom, intrmd_not_dom, not_intrmd_not_dom), by=list(name)]
pfam_cnt <- merge(pfam_cnt, fttab, by='name')
pfam_cnt <- pfam_cnt[order(fisher.pval)]
pfam_cnt[, signature_desc := URLdecode(signature_desc)]
fwrite(pfam_cnt, xargs$tblout, sep='\t')
