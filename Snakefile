from localutils import *

include: "download.Snakefile"
include: "convert.Snakefile"
include: "build-index.Snakefile"

rule all:
    input: BBMap="BBMap_index_hg38.analysisSet/ref/genome/1/summary.txt",
           HISAT2="HISAT2_index_grch38_snp_tran/index.1.ht2",
           STAR_kg="STAR_index_hg38.analysisSet_knownGene/SA",
           salmon_kg='Salmon_index_hg38.analysisSet_knownGene/sa.bin',
           STAR_gen25="STAR_index_hg38.analysisSet_gencode.v25/SA",
           salmon_gen25='Salmon_index_hg38.analysisSet_gencode.v25/sa.bin',
