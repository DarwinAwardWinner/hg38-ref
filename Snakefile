from localutils import *

include: 'download.Snakefile'
include: 'convert.Snakefile'
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
           STAR_gen25=star_index_files('STAR_index_hg38.analysisSet_gencode.v25'),
           tophat2_gen25=tophat2_index_files('TH2_index_hg38.analysisSet_gencode.v25', 'index'),
           salmon_gen25=salmon_index_files('Salmon_index_hg38.analysisSet_gencode.v25'),
