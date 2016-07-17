'''Rules for building indices for various alignment tools.'''

import os.path

from localutils import *

_bt1_index_file_infixes = ('1', '2', '3', '4', 'rev.1', 'rev.2',)
def bt1_index_files(path, prefix='index', large=True):
    suffix = "ebwt"
    if large:
        suffix += "l"
    fnames = ('{prefix}.{infix}.{suffix}'.format(prefix=prefix, infix=x, suffix=suffix)
              for x in _bt1_index_file_infixes)
    return tuple(os.path.join(path, f) for f in fnames)
rule build_bowtie1_index:
    input: genome_fa='{genome_build}.fa'
    output: bt1_index_files('BT1_index_{genome_build}', 'index', large=True)
    params: outdir='BT1_index_{wildcards.genome_build}',
            basename='BT1_index_{wildcards.genome_build}/index'
    shell: '''
    mkdir -p {params.outdir:q} && \
        bowtie-build --large-index {input.genome_fa:q} {params.basename:q}
    '''

_bt2_index_file_infixes = ('1', '2', '3', '4', 'rev.1', 'rev.2',)
def bt2_index_files(path, prefix='index', large=True):
    suffix = "bt2"
    if large:
        suffix += "l"
    fnames = ('{prefix}.{infix}.{suffix}'.format(prefix=prefix, infix=x, suffix=suffix)
              for x in _bt2_index_file_infixes)
    return tuple(os.path.join(path, f) for f in fnames)

rule build_bowtie2_index:
    input: genome_fa='{genome_build}.fa'
    output: bt2_index_files('BT2_index_{genome_build}', 'index', large=True)
    params: outdir='BT2_index_{wildcards.genome_build}',
            basename='BT2_index_{wildcards.genome_build}/index'
    shell: '''
    mkdir -p {params.outdir:q} && \
        bowtie2-build --large-index {input.genome_fa:q} {params.basename:q}
    '''

bwa_index_files = (
    'index.amb',
    'index.ann',
    'index.bwt',
    'index.pac',
    'index.sa',
)
rule build_bwa_index:
    input: genome_fa='{genome_build}.fa'
    output: expand('BWA_index_{{genome_build}}/{filename}', filename=bwa_index_files)
    params: outdir='BWA_index_{wildcards.genome_build}',
            basename='BWA_index_{wildcards.genome_build}/index'
    shell: '''
    mkdir -p {params.outdir:q} && \
        bwa index -a bwtsw -p {params.basename:q} {input.genome_fa:q}
    '''

# TODO: Make it configurable
BBMAP=os.path.expanduser("~/opt/bbmap/bbmap.sh")

# File names are not all consistent for bbmap, so just list the few
# that are. These will be used as indicator files for the presence of
# the index.
_bbmap_index_files = (
    "ref/genome/1/info.txt",
    "ref/genome/1/summary.txt",
)
def bbmap_index_files(path):
    return tuple(os.path.join(path, f) for f in _bbmap_index_files)

rule build_bbmap_index:
    input: genome_fa='{genome_build}.fa'
    output: bbmap_index_files('BBMap_index_{genome_build}')
    params: outdir='BBMap_index_{wildcards.genome_build}'
    threads: 16
    run:
        ensure_empty_dir(params.outdir)
        shell('''
        {BBMAP:q} ref={input.genome_fa:q} \
            path={params.outdir:q} \
            t={threads}
        ''')

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
