'''Generic rules for converting between file types.'''

rule twoBit_to_fasta:
    input: '{basename}.2bit'
    output: '{basename}.fa'
    shell: 'twoBitToFa {input:q} {output:q}'

rule gtf_to_gff3:
    input: '{basename}.gtf'
    output: '{basename}.gff'
    shell: 'gffread -FOE {input:q} -o {output:q}'

rule extract_transcript_seqs:
    input: genome_fa='{genome_build}.fa', transcriptome_gff='{transcriptome}.gff'
    output: '{genome_build}_{transcriptome}_transcripts.fa'
    shell: 'gffread -w {output:q} -g {input.genome_fa:q} {input.transcriptome_gff:q}'
