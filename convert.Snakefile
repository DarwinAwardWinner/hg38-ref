'''Generic rules for converting between file types.'''

from rpy2.robjects import r

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
    priority: 0
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

# Convert GENCODE annotation from GENCODE to UCSC
rule fix_gencode_annot_chrom_ids:
    input: gff='gencode.v{release}_raw.gff3', mapping='chrom_mapping_GRCh38_gencode2UCSC.txt'
    output: gff='gencode.v{release,\\d+}.gff3'
    version: BIOC_VERSION
    run:
        r['source']("scripts/map-gff-chrom.R")
        mapping = r['read.chrom.mapping'](input.mapping)
        gr = r['import'](input.gff, format="GFF3")
        fixed_gr = r['map.seqlevels'](gr, mapping, keep_unmatched=False)
        r['export'](fixed_gr, output.gff, format="GFF3")
