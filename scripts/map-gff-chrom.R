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
    seqlevels_to_rename <- intersect(seqlevels(gr), names(new_seqlevels))
    new_seqlevels <- new_seqlevels[seqlevels_to_rename]
    if (keep_unmatched) {
        seqlevels_to_keep <- setdiff(seqlevels(gr), seqlevels_to_rename)
        new_seqlevels <- c(new_seqlevels, setNames(nm=seqlevels_to_keep))
    } else {
        gr <- keepSeqlevels(gr, seqlevels_to_rename)
    }
    gr <- renameSeqlevels(gr, new_seqlevels)
    return(gr)
}
