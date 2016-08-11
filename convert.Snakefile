'''Generic rules for converting between file types.'''

from atomicwrites import atomic_write

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
    run:
        mapping = read_chrom_mapping(input.mapping)
        with open(input.gff, "r") as infile, \
             atomic_write(output.gff, overwrite=True) as outfile:
            for line in infile:
                line = line.strip('\n')
                # Line indicating chromosome length: chr is 2nd field.
                if line.startswith('##sequence-region'):
                    fields = line.split()
                    try:
                        fields[1] = mapping[fields[1]]
                    except KeyError:
                        continue
                    line = ' '.join(fields)
                # Regular comment: pass through unchanged
                elif line.startswith('#'):
                    pass
                # Normal GFF line: chr is 1st field
                else:
                    fields = line.split('\t')
                    try:
                        fields[0] = mapping[fields[0]]
                    except KeyError:
                        continue
                    line = '\t'.join(fields)
                outfile.write(line + '\n')
