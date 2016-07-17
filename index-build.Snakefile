import os.path

from localutils import *

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

# TODO: Place the Log file somewhere better than the root
rule build_star_index:
    input: genome_fa='{genome_build}.fa', transcriptome_gff='{transcriptome}.gff'
    output: expand('STAR_index_{{genome_build}}_{{transcriptome}}/{filename}', filename=star_index_files)
    params: outdir='STAR_index_{wildcards.genome_build}_{wildcards.transcriptome}'
    threads: 16
    run:
        ensure_empty_dir(params.outdir)
        shell('''
        STAR --runMode genomeGenerate \
            --genomeDir {params.outdir:q} \
            --genomeFastaFiles {input.genome_fa:q} \
            --sjdbGTFfile {input.transcriptome_gff:q} \
            --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent \
            --sjdbOverhang 100 \
            --runThreadN {threads:q}
        ''')

# TODO: Make it configurable
BBMAP=os.path.expanduser("~/opt/bbmap/bbmap.sh")

# File names are not really certain for bbmap, so just list the few
# that are. These will be used as indicator files. TODO: Maybe use
# dynamic files to include every file produced in the directory?
bbmap_index_files = (
    "ref/genome/1/info.txt",
    "ref/genome/1/summary.txt",
)

rule build_bbmap_index:
    input: genome_fa='{genome_build}.fa'
    output: expand('BBMap_index_{{genome_build}}/{filename}', filename=bbmap_index_files)
    params: outdir='BBMap_index_{wildcards.genome_build}'
    threads: 16
    run:
        ensure_empty_dir(params.outdir)
        shell('''
        {BBMAP:q} ref={input.genome_fa:q} \
            path={params.outdir:q} \
            t={threads}
        ''')
