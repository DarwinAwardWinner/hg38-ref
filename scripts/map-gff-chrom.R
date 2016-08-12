#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(GenomeInfoDb)
    library(GenomicRanges)
    library(GenomicFeatures)
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

makeTxDbFromBiomart <- function (biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
    transcript_ids = NULL, circ_seqs = DEFAULT_CIRC_SEQS, filter = NULL,
    id_prefix = "ensembl_", host = "www.ensembl.org", port = 80,
    taxonomyId = NA, miRBaseBuild = NA, chrom.mapping=NULL)
{
    mart <- .useMart2(biomart = biomart, dataset = dataset, host = host,
        port = port)
    id_prefix <- .normarg_id_prefix(id_prefix)
    filter <- .add_tx_id_filter(filter, transcript_ids, id_prefix)
    valid_filter_names <- listFilters(mart, what = "name")
    invalid_filter_names <- setdiff(names(filter), valid_filter_names)
    if (length(invalid_filter_names) != 0L) {
        in1string <- paste0(invalid_filter_names, collapse = ", ")
        stop(wmsg("Invalid filter name(s): ", in1string, "\n\nPlease use the listFilters() function from the ",
            "biomaRt package to get valid filter names."))
    }
    is_full_dataset <- length(filter) == 0L
    recognized_attribs <- recognizedBiomartAttribs(id_prefix)
    transcripts <- .makeBiomartTranscripts(filter, mart, transcript_ids,
                                           recognized_attribs, id_prefix)
    transcripts_tx_id <- transcripts$tx_id
    names(transcripts_tx_id) <- transcripts$tx_name
    chrominfo <- .makeBiomartChrominfo(mart, extra_seqnames = transcripts$tx_chrom,
        circ_seqs = circ_seqs, host, port)
    if (!is_full_dataset) {
        keep_idx <- which(chrominfo[, "chrom"] %in% transcripts$tx_chrom)
        chrominfo <- S4Vectors:::extract_data_frame_rows(chrominfo,
            keep_idx)
    }
    splicings <- .makeBiomartSplicings(filter, mart, transcripts_tx_id,
        recognized_attribs, id_prefix = id_prefix)
    utr_anomaly <- splicings$utr_anomaly
    if (!is.null(utr_anomaly)) {
        invalid_tx <- unique(splicings[utr_anomaly != 0L, "tx_id"])
        if (length(invalid_tx) != 0L) {
            message("Drop transcripts with UTR anomalies (",
                length(invalid_tx), " transcripts) ... ", appendLF = FALSE)
            keep_idx1 <- !(transcripts$tx_id %in% invalid_tx)
            transcripts <- S4Vectors:::extract_data_frame_rows(transcripts,
                keep_idx1)
            transcripts_tx_id <- transcripts_tx_id[keep_idx1]
            keep_idx2 <- !(splicings$tx_id %in% invalid_tx)
            splicings <- S4Vectors:::extract_data_frame_rows(splicings,
                keep_idx2)
            message("OK")
        }
        splicings$utr_anomaly <- NULL
    }
    genes <- .makeBiomartGenes(filter, mart, transcripts_tx_id,
        recognized_attribs, id_prefix)
    metadata <- .prepareBiomartMetadata(mart, is_full_dataset,
        host, port, taxonomyId, miRBaseBuild)
    ## browser()
    if (!is.null(chrom.mapping)) {
        transcripts <- transcripts[transcripts$tx_chrom %in% names(chrom.mapping),]
        chrominfo <- chrominfo[chrominfo$chrom %in% names(chrom.mapping),]
        transcripts_tx_id <- transcripts$tx_id
        splicings <- splicings[splicings$tx_id %in% transcripts_tx_id,]
        genes <- genes[genes$tx_id %in% transcripts_tx_id,]
        levels(chrominfo$chrom) <- chrom.mapping[levels(chrominfo$chrom)]
        levels(transcripts$tx_chrom) <- chrom.mapping[levels(transcripts$tx_chrom)]
    }
    message("Make the TxDb object ... ", appendLF = FALSE)
    txdb <- makeTxDb(transcripts, splicings, genes = genes, chrominfo = chrominfo,
        metadata = metadata, reassign.ids = TRUE)
    message("OK")
    txdb
}
environment(makeTxDbFromBiomart) <- new.env(parent=environment(GenomicFeatures::makeTxDbFromBiomart))
