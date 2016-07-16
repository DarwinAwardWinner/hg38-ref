rule get_hg38_analysisSet_twoBit:
    output: 'hg38.analysisSet.2bit'
    shell: 'curl -o {output:q} http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/analysisSet/hg38.analysisSet.2bit'

rule twoBit_to_fasta:
    input: '{basename}.2bit'
    output: '{basename}.fa'
    shell: 'twoBitToFa {input:q} {output:q}'

rule get_knownGene_gtf:
    output: 'knownGene.gtf'
    shell: 'genePredToGtf -addComments -utr hg38 knownGene {output:q}'

rule gtf_to_gff3:
    input: '{basename}.gtf'
    output: '{basename}.gff'
    shell: 'gffread -FOE {input:q} -o {output:q}'
