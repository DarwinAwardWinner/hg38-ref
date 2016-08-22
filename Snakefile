from localutils import *

from snakemake.utils import min_version
min_version("3.7.1")

include: 'download.Snakefile'
include: 'convert.Snakefile'
include: 'build-annot.Snakefile'
include: 'build-index.Snakefile'
include: 'rulegraph.Snakefile'

rule all:
    input: BT1=bt1_index_files('BT1_index_hg38.analysisSet', 'index', large=True),
           BT2=bt2_index_files('BT2_index_hg38.analysisSet', 'index', large=True),
           bwa=bwa_index_files('BWA_index_hg38.analysisSet', 'index'),
           BBMap=bbmap_index_files('BBMap_index_hg38.analysisSet'),
           HISAT2=hisat2_index_files('HISAT2_index_grch38_snp_tran', 'index', large=False),
           STAR_kg=star_index_files('STAR_index_hg38.analysisSet_knownGene'),
           tophat2_kg=tophat2_index_files('TH2_index_hg38.analysisSet_knownGene', 'index'),
           salmon_kg=salmon_index_files('Salmon_index_hg38.analysisSet_knownGene'),
           STAR_ens85=star_index_files('STAR_index_hg38.analysisSet_ensembl.85'),
           tophat2_ens85=tophat2_index_files('TH2_index_hg38.analysisSet_ensembl.85', 'index'),
           salmon_ens85=salmon_index_files('Salmon_index_hg38.analysisSet_ensembl.85'),
           ensembl_txdb='TxDb.Hsapiens.ensembl.hg38.v85.sqlite3',
