rule all:
    input: "STAR_index_hg38.analysisSet_knownGene/SA"

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

# TODO: Write a function to compute the input files for a STAR index.
# Also do so with BBMap, etc.
star_index_files = (
    'chrLength.txt',
    'chrNameLength.txt',
    'chrName.txt',
    'chrStart.txt',
    'Genome',
    'genomeParameters.txt',
    'SA',
    'SAindex',
    'sjdbInfo.txt',
    'sjdbList.out.tab',
)

rule build_star_index:
    input: genome_fa='{genome_build}.fa', transcriptome_gff='{transcriptome}.gff'
    output: expand('STAR_index_{{genome_build}}_{{transcriptome}}/{filename}', filename=star_index_files)
    threads: 16
    shell: '''
    mkdir -p STAR_index_{wildcards.genome_build:q}_{wildcards.transcriptome:q} && \
        STAR --runMode genomeGenerate \
        --genomeDir STAR_index_{wildcards.genome_build:q}_{wildcards.transcriptome:q} \
        --genomeFastaFiles {input.genome_fa:q} \
        --sjdbGTFfile {input.transcriptome_gff:q} \
        --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent \
        --sjdbOverhang 100 \
        --runThreadN {threads:q}
    '''
