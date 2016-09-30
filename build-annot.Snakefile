'''Rules for creating various annotation file formats'''

include: 'localutils.py'
include: 'tool_versions.py'

rule make_gencode_txdb:
    input: gff='gencode.v{release}.gff3',
    output: dbfile='TxDb.Hsapiens.gencode.hg38.v{release}.sqlite3',
    version: BIOC_VERSION
    run:
        import rpy2.rinterface
        rpy2.rinterface.set_writeconsole_warnerror(lambda x: sys.stderr.write(x))
        from rpy2.robjects import r
        from rpy2.robjects import globalenv as r_env
        from rpy2.rinterface import StrSexpVector
        txdb_meta = {
            'TaxID': '9606',
            'Data source': 'Gencode',
            'Gencode Release': wildcards.release,
            'Resource URL': 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/',
            'Type of Gene ID': 'Ensembl Gene ID',
        }
        txdb_meta_df = r['data.frame'](
            name=StrSexpVector(list(txdb_meta.keys())),
            value=StrSexpVector(list(txdb_meta.values())),
            stringsAsFactors=False,
        )
        r_env['txdb.meta'] = txdb_meta_df
        r_env['taxid'] = txdb_meta['TaxID']
        r_env['input.gff'] = input.gff
        r_env['output.dbfile'] = output.dbfile
        r('''
        suppressMessages({
            library(rtracklayer)
            library(GenomicFeatures)
            library(BSgenome.Hsapiens.UCSC.hg38)
            library(stringr)
            library(magrittr)
        })
        gff <- import(input.gff, format='GFF3')
        seqlevels(gff) <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
        seqinfo(gff) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
        # Remove the version number from the end of the IDs
        mcols(gff)[c("ID", "Parent", "gene_id", "transcript_id", "exon_id", "protein_id")] %<>%
            lapply(. %>% str_replace("\\\\.[0-9]+", ""))
        txdb <- makeTxDbFromGRanges(gff, taxonomyId=taxid, metadata=txdb.meta)
        saveDb(txdb, output.dbfile)
        ''')

rule make_ensembl_txdb:
    input: mapping='chrom_mapping_GRCh38_ensembl2UCSC.txt'
    output: dbfile='TxDb.Hsapiens.ensembl.hg38.v{release}.sqlite3'
    params: ensembl_host='e{release}.ensembl.org'
    version: BIOC_VERSION
    run:
        release = int(wildcards.release)
        if release < 76:
            raise ValueError('Only Ensembl releases 76 and higher are built on hg38')
        from rpy2.robjects import r
        from rpy2.robjects import globalenv as r_env
        import rpy2.rinterface
        rpy2.rinterface.set_writeconsole_warnerror(lambda x: sys.stderr.write(x))
        r_env['chrom.mapping.file'] = input.mapping
        r_env['ensembl.host'] = params.ensembl_host
        r_env['output.dbfile'] = output.dbfile
        r('''
        suppressMessages({
            library(biomaRt)
            library(GenomicFeatures)
            library(BSgenome.Hsapiens.UCSC.hg38)
        })
        source("scripts/map-gff-chrom.R")
        chrom.mapping <- read.chrom.mapping(chrom.mapping.file)
        txdb <- makeTxDbFromBiomart(
            biomart="ENSEMBL_MART_ENSEMBL",
            dataset="hsapiens_gene_ensembl",
            host=ensembl.host,
            chrom.mapping=chrom.mapping)
        saveDb(txdb, output.dbfile)
        ''')

rule make_ensembl_genemeta:
    output: datafile='genemeta.ensembl.{release}.RDS'
    params: ensembl_host='e{release}.ensembl.org'
    version: BIOC_VERSION
    run:
        from rpy2.robjects import r
        from rpy2.robjects import globalenv as r_env
        import rpy2.rinterface
        rpy2.rinterface.set_writeconsole_warnerror(lambda x: sys.stderr.write(x))
        r_env['output.datafile'] = output.datafile
        r_env['ensembl.host'] = params.ensembl_host
        r('''
        source("scripts/build-annot.R")
        annot <- build.annot.from.biomart(host=ensembl.host, keytype="ENSEMBL")
        saveRDS(annot, output.datafile)
        ''')

rule make_orgdb_genemeta:
    output: datafile='genemeta.org.{name}.db.RDS'
    params: orgdbname='org.{name}.db'
    # Are orgdb packages versioned independently of BioC?
    version: BIOC_VERSION
    run:
        from rpy2.robjects import r
        from rpy2.robjects import globalenv as r_env
        import rpy2.rinterface
        rpy2.rinterface.set_writeconsole_warnerror(lambda x: sys.stderr.write(x))
        r_env['output.datafile'] = output.datafile
        r_env['orgdbname'] = params.orgdbname
        r('''
        source("scripts/build-annot.R")
        library(orgdbname, character.only=TRUE)
        orgdb <- get(orgdbname)
        annot <- build.annot.from.orgdb(orgdb, keytype="ENTREZID")
        saveRDS(annot, output.datafile)
        ''')
