suppressMessages({
  library(magrittr)
  library(assertthat)
  library(biomaRt)
})

build.annot.from.orgdb <- function(orgdb, keytype=c("ENTREZID", "SYMBOL", "UNIGENE")) {
    library(AnnotationDbi)
    keytype <- match.arg(keytype)
    assert_that(keytype %in% keytypes(orgdb))
    ## Add gene annotations as metadata on the exons-by-gene list
    metacols <- c("ENTREZID", "SYMBOL", "GENENAME", "UNIGENE", "ENSEMBL", "REFSEQ") %>%
        setdiff(keytype) %>%
        intersect(keytypes(orgdb))
    allkeys <- keys(orgdb, keytype)
    gene.annot <- DataFrame(row.names=allkeys, setNames(list(allkeys), keytype))
    suppressMessages({
        gene.annot[metacols] <- lapply(metacols, mapIds, x=orgdb,
                                       keys=gene.annot[[keytype]],
                                       keytype=keytype, multiVals="CharacterList")
    })
    for (i in names(gene.annot)) {
        vals <- gene.annot[[i]]
        if (all(lengths(vals) <= 1)) {
            ## Single-value column: fill in empty values with NA, then
            ## unlist
            vals[lengths(vals) == 0] <- NA
            assert_that(all(lengths(vals) == 1))
            gene.annot[[i]] <- unlist(vals)
        } else {
            ## Multi-value column: Remove NA values
            vals %<>% endoapply(na.omit)
            gene.annot[[i]] <- vals
        }
    }

    annot.meta <- metadata(orgdb) %$%
        setNames(value, name) %>%
        .[names(.) %in%
          c("ORGANISM", "SPECIES", "EGSOURCEDATE", "EGSOURCENAME",
            "EGSOURCEURL", "TAXID", "ENSOURCEDATE", "ENSOURCENAME",
            "ENSOURCEURL", "UPSOURCENAME", "UPSOURCEURL",
            "UPSOURCEDATE")] %>% as.list
    annot.meta$KEYTYPE <- keytype
    annot.meta$PACKAGE <- orgdb$packageName
    metadata(gene.annot) %<>% c(annot.meta)

    return(gene.annot)
}

build.annot.from.biomart <-
    function(host="www.ensembl.org",
             keytype=c("ENSEMBL", "ENTREZID", "SYMBOL", "UNIGENE"),
             mart=useEnsembl(biomart="ensembl", host=host,
                             dataset="hsapiens_gene_ensembl"))
{

    ## Name: Biomart attr name;
    ## Value: corresponding data frame colname
    metacols <- c(ensembl_gene_id="ENSEMBL",
                  entrezgene="ENTREZID",
                  external_gene_name="SYMBOL",
                  wikigene_description="GENENAME",
                  unigene="UNIGENE",
                  refseq_mrna="REFSEQ")

    keytype <- match.arg(keytype)
    key_attr_name <- names(metacols)[match(keytype, metacols)]
    metacols %<>% .[names(.) != key_attr_name]

    gene.annot <- getBM(attributes=key_attr_name, mart=mart) %>%
        setNames(keytype) %>% set_rownames(.[[1]])

    for (i in names(metacols)) {
        bm <- getBM(attributes=c(key_attr_name, i),
                    mart=mart) %>%
            {List(split(as.character(.[[2]]), .[[1]]))} %>%
            .[! (. == "" | is.na(.)) ]
        if (all(lengths(bm) <= 1)) {
            ## Single-value column: fill in empty values with NA, then
            ## unlist
            bm[lengths(bm) == 0] <- NA
            assert_that(all(lengths(bm) == 1))
            bm <-  unlist(bm)
        } else {
            ## Multi-value column: fill in missing IDs with empty vector
            missing.ens <- setdiff(gene.annot$ENSEMBL, names(bm))
            missing.ids <- List(character(0)) %>% rep(length(missing.ens)) %>% setNames(missing.ens)
            bm %<>% c(missing.ids)
        }
        gene.annot[[metacols[i]]] <- bm[gene.annot$ENSEMBL]
    }

    rownames(gene.annot) <- gene.annot[[keytype]]

    annot.meta <- list(
        KEYTYPE=keytype,
        MARTHOST=mart@host,
        MART=mart@biomart,
        MARTDATASET=mart@dataset)
    metadata(gene.annot) %<>% c(annot.meta)

    return(gene.annot)
}
