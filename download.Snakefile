'''Rules for downloading specific files from the internet.'''

rule get_hg38_analysisSet_twoBit:
    output: 'hg38.analysisSet.2bit'
    shell: 'curl -o {output:q} http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/analysisSet/hg38.analysisSet.2bit'

rule get_knownGene_gtf:
    output: 'knownGene.gtf'
    shell: 'genePredToGtf -addComments -utr hg38 knownGene {output:q}'
