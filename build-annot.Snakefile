'''Rules for creating various annotation file formats'''

from localutils import *
from tool_versions import *

rule make_gencode_txdb:
    input: gff='gencode.v{release}.gff3'
    output: dbfile='TxDb.Hsapiens.Gencode.hg38.v{release}.sqlite3'
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
    output: dbfile='TxDb.Hsapiens.Ensembl.hg38.v{release}.sqlite3'
    run:
        release = int(wildcards.release)
        if release < 76:
            raise ValueError('Only Ensembl releases 76 and higher are built on hg38')
        host = 'e{release}.ensembl.org'.format(release=release)
        from rpy2.robjects import r
        from rpy2.robjects import globalenv as r_env
        import rpy2.rinterface
        rpy2.rinterface.set_writeconsole_warnerror(lambda x: sys.stderr.write(x))
        r_env['chrom.mapping.file'] = input.mapping
        r_env['ensembl.host'] = host
        r_env['output.dbfile'] = output.dbfile
        r('''
        suppressMessages({
            library(biomaRt)
            library(GenomicFeatures)
            library(BSgenome.Hsapiens.UCSC.hg38)
        })
        source("scripts/map-gff-chrom.R")
        txdb <- makeTxDbFromBiomart(
            biomart="ENSEMBL_MART_ENSEMBL",
            dataset="hsapiens_gene_ensembl",
            host=ensembl.host)
        chrom.mapping <- read.chrom.mapping(chrom.mapping.file)
        txdb <- map.seqlevels(txdb, chrom.mapping)
        saveDb(txdb, output.dbfile)
        ''')
