suppressMessages({
  library(magrittr)
  library(assertthat)
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
        ids <- gene.annot[[i]]
        if (all(lengths(ids) <= 1)) {
            ids[lengths(ids) == 0] <- NA
            assert_that(all(lengths(ids) == 1))
            gene.annot[[i]] <- unlist(ids)
        }
    }
    return(gene.annot)
}

build.annot.from.biomart <-
    function(host="www.ensembl.org",
             keytype=c("ENSEMBL", "ENTREZID", "SYMBOL", "UNIGENE"),
             mart=useEnsembl(biomart="ensembl", host=host,
                             dataset="hsapiens_gene_ensembl"))
{
    library(biomaRt)

    ## Names: Biomart attr name;
    ## Values: corresponding data frame colname
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
            bm[lengths(bm) == 0] <- NA
            assert_that(all(lengths(bm) == 1))
            bm <-  unlist(bm)
        } else {
            missing.ens <- setdiff(gene.annot$ENSEMBL, names(bm))
            missing.ids <- List(character(0)) %>% rep(length(missing.ens)) %>% setNames(missing.ens)
            bm %<>% c(missing.ids)

        }
        gene.annot[[metacols[i]]] <- bm[gene.annot$ENSEMBL]
    }
    return(gene.annot)
}

## Defunct code for auto-determining ID type

## if (cmdopts$txdb_geneid_type == "auto") {
##     idtype.from.meta <- metadata(txdb) %>% as.data.frame %>%
##         filter(name == "Type of Gene ID") %$% value
##     ## Try to find key words in the gene ID type
##     idtype <- idtype.from.meta %>%
##         str_detect(fixed(names(known.gene.id.types), ignore_case = TRUE)) %>%
##         which %>% .[1] %>% na.omit %>% known.gene.id.types[.]
##     if (length(idtype) != 1) {
##         idtype <- tryCatch({
##             identify.ids(names(annot), db="org.Hs.eg.db", idtypes=known.gene.id.types)
##         }, error=function(...) {
##             NULL
##         })
##     }
## } else {
##     idtype <- known.gene.id.types[cmdopts$txdb_geneid_type]
## }
