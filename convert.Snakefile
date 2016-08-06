'''Generic rules for converting between file types.'''

from rpy2 import robjects

from localutils import *
from tool_versions import *

rule twoBit_to_fasta:
    input: '{basename}.2bit'
    output: '{basename}.fa'
    shell: 'twoBitToFa {input:q} {output:q}'

rule gtf_to_gff3:
    input: '{basename}.gtf'
    output: '{basename}.gff3'
    version: CUFFLINKS_VERSION
    shell: 'gffread -FOE {input:q} -o {output:q}'

rule index_fa:
    input: '{basename}.fa'
    output: '{basename}.fa.fai'
    version: SAMTOOLS_VERSION
    shell: 'samtools faidx {input:q}'

rule extract_transcript_seqs:
    input: genome_fa='{genome_build}.fa', transcriptome_gff='{transcriptome}.gff3',
           genome_fai='{genome_build}.fa.fai'
    output: '{genome_build}_{transcriptome}_transcripts.fa'
    version: CUFFLINKS_VERSION
    shell: 'gffread -w {output:q} -g {input.genome_fa:q} {input.transcriptome_gff:q}'

# rule make_gencode_txdb:
#     input: 'gencode.v{release}.gff3'
#     output: 'gencode.v{release}.txdb.rds'
#     run:
#         robjects.globalenv['input'] <- input
#         robjects.globalenv['output'] <- output
#         infile <-
#         R('''suppressWarnings({
#         library(rtracklayer)
#         library(GenomicFeatures)
#         gff <- import(input, format="gff3")
#         txdb <- make
#         })''')
