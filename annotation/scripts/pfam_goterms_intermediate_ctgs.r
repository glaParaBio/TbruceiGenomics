#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(topGO))

parser <- ArgumentParser()
parser$add_argument('--gff', help='GFF file annotated with GO terms')
parser$add_argument('--intercontigs', help='File listing the intermediate contigs to test for enrichment')
parser$add_argument('--tsvout', help='Write output table to this file')
parser$add_argument('--rdata-out', help='Write Rdata objects to this dir')

xargs <- parser$parse_args()

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

lod_enrichment_score <- function(sig_go, sig_nogo, nosig_go, nosig_nogo) {
    # Almost identical to:
    #lod <- log2((sig_go / sig_nogo) / (nosig_go / nosig_nogo))
    #return(lod)
    m <- matrix(c(sig_go, sig_nogo, nosig_go, nosig_nogo), nrow=2, byrow=TRUE)
    ft <- fisher.test(m+1, conf.int=FALSE)
    return(log2(ft$estimate))
}

# gff <- fread(sprintf('grep -v "^#" canu_tb427_ont.contigs_pilon_4x/hmmer/augustus.hints.gff3'))
# inter <- fread('/export/projects/III-data/wcmp_bioinformatics/db291g/git_repos/wcip/20220909_marija_annotation/data/canu_tb427_ont.contigs_pilon_4x.intermediate_contigs.txt', header=FALSE, col.name='contig')
gff <- fread(cmd=sprintf('grep -v "^#" %s', xargs$gff))
inter <- fread(xargs$intercontigs, header=FALSE, col.name='contig')

go <- gff[V3 == 'protein_match', list(contig=V1, txid=get_gff_attr(V9, 'Parent'), Name=get_gff_attr(V9, 'Name'), go_term=get_gff_attr(V9, 'Ontology_term'))][!is.na(go_term)]
go <- unique(go) # Because you can have a tx with copies of the same domain
stopifnot(nrow(go[, list(txid, Name)]) == nrow(unique(go[, list(txid, Name)])))
go <- go[, list(go_term=strsplit(go_term, ',')[[1]]), by=list(contig, txid, Name)]
go[, is_intermediate := contig %in% inter$contig]

gene2tx <- gff[V3 == 'mRNA', list(gene_id=get_gff_attr(V9, 'Parent'), txid=get_gff_attr(V9, 'ID'))]
stopifnot(go$txid %in% gene2tx$txid)
go <- merge(go, gene2tx, by='txid')
genes <- unique(go[, list(gene_id, is_intermediate)])

gene2GO <- list()
for(g in genes$gene_id) {
    gene2GO[[g]] <- unique(go[gene_id == g]$go_term)
}

allgenes <- as.numeric(genes$is_intermediate)
names(allgenes) <- genes$gene_id

gout <- list()
for(is_intermediate in c(0, 1)){
    for(alg in c('parentChild', 'classic')) {
        for(go_set in c('BP', 'MF', 'CC')) {
            topGOdata <- new("topGOdata",
                ontology = go_set,
                allGenes = allgenes, 
                geneSelectionFun= function(x) {return(x == is_intermediate)},
                nodeSize = 5,
                annot = annFUN.gene2GO, gene2GO= gene2GO)
            topGOresult <- runTest(topGOdata, algorithm=alg, statistic='fisher')

            go_category <- ifelse(go_set == 'BP', 'Biological Process', ifelse(go_set == 'MF', 'Molecular Function', ifelse(go_set == 'CC', 'Cellular Component', NA)))
            enriched_in <- ifelse(as.logical(is_intermediate) == TRUE, 'intermediate_ctgs', 'other_ctgs')

            outdir <- file.path(xargs$rdata_out, enriched_in, alg, gsub(' ', '_', go_category))
            system(sprintf('mkdir -p %s', outdir))
            save(topGOdata, is_intermediate, file=file.path(outdir, 'topGOdata.RData'))
            save(topGOresult, file=file.path(outdir, 'topGOresult.RData'))

            got <- data.table(GenTable(topGOdata, p.value= topGOresult, topNodes= length(topGOdata@graph@nodes), numChar= 1000))
            got[, go_category := go_category]
            got[, algorithm := sprintf('pvalue.%s', alg)]
            got[, enriched_in := enriched_in]
            
            if(is_intermediate == 1){
                allsig <- length(unique(go[is_intermediate == TRUE]$gene_id))
                allnosig <- length(unique(go[is_intermediate == FALSE]$gene_id))
            } else if (is_intermediate == 0) {
                allsig <- length(unique(go[is_intermediate == FALSE]$gene_id))
                allnosig <- length(unique(go[is_intermediate == TRUE]$gene_id))
            } else {
                stop(sprintf('Invalid value: %s', is_intermediate))
            }

            got[, lod_enrichment_score := lod_enrichment_score(sig_go = Significant,  
                                               sig_nogo = allsig - Significant,
                                               nosig_go = Annotated - Significant,
                                               nosig_nogo = allnosig - Annotated + Significant
                                               ), by=GO.ID]
            gout[[length(gout)+1]] <- got
        }
    }
}
gout <- rbindlist(gout)
gout[, p.value := as.numeric(p.value)]
setnames(gout, c('GO.ID', 'Term', 'Significant'), c('go_term', 'go_name', 'test_set'))

gout <- dcast(data=gout, go_term + go_name + Annotated + test_set + enriched_in + Expected + go_category + lod_enrichment_score ~ algorithm, value.var='p.value')[order(pvalue.parentChild)]

fwrite(gout, xargs$tsvout, sep='\t')
