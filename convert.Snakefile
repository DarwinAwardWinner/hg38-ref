'''Generic rules for converting between file types.'''

rule twoBit_to_fasta:
    input: '{basename}.2bit'
    output: '{basename}.fa'
    shell: 'twoBitToFa {input:q} {output:q}'

rule gtf_to_gff3:
    input: '{basename}.gtf'
    output: '{basename}.gff'
    shell: 'gffread -FOE {input:q} -o {output:q}'
