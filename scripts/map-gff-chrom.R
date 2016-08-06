#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
})

read.chrom.mapping <- function(filename) {
    x <- read.table(filename, sep="\t", stringsAsFactors = FALSE)
    x <- setNames(x[[2]], x[[1]])
    x[x != ""]
}

map.seqlevels <- function(gr, new_seqlevels, keep_unmatched=FALSE) {
    new_seqlevels <- new_seqlevels[names(new_seqlevels) %in% seqlevels(gr)]
    if (!keep_unmatched) {
        gr <- keepSeqlevels(gr, names(new_seqlevels))
    }
    gr <- renameSeqlevels(gr, new_seqlevels)
    return(gr)
}
