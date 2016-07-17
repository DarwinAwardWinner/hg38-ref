from localutils import *

include: "download.Snakefile"
include: "convert.Snakefile"
include: "build-index.Snakefile"

rule all:
    input: STAR_ref="STAR_index_hg38.analysisSet_knownGene/SA",
           BBMap_ref="BBMap_index_hg38.analysisSet/ref/genome/1/summary.txt",
