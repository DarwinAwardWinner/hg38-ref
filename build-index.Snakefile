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

_bt2_index_file_infixes = ('1', '2', '3', '4', 'rev.1', 'rev.2',)
def bt2_index_files(path, prefix='index', large=True):
    suffix = "bt2"
    if large:
        suffix += "l"
    fnames = ('{prefix}.{infix}.{suffix}'.format(prefix=prefix, infix=x, suffix=suffix)
              for x in _bt2_index_file_infixes)
    return tuple(os.path.join(path, f) for f in fnames)

# These are in addition to the bowtie2 files for the same prefix
_tophat2_index_suffixes = ('.gff', '.fa', '.fa.tlst', '.fa.ver')
def tophat2_index_files(path, prefix='index', large=True):
    fnames = ('{prefix}{suffix}'.format(prefix=prefix, suffix=s)
              for s in _tophat2_index_suffixes)
    paths = tuple(os.path.join(path, f) for f in fnames)
    return paths + bt2_index_files(path, prefix, large)

_bwa_index_suffixes = (
    '.amb',
    '.ann',
    '.bwt',
    '.pac',
    '.sa',
)
def bwa_index_files(path, prefix='index'):
    fnames = ('{prefix}{suffix}'.format(prefix=prefix, suffix=s)
              for s in _bwa_index_suffixes)
    return tuple(os.path.join(path, f) for f in fnames)

# File names are not all consistent for bbmap, so just list the few
# that are. These will be used as indicator files for the presence of
# the index.
_bbmap_index_filenames = (
    "ref/genome/1/info.txt",
    "ref/genome/1/summary.txt",
)
def bbmap_index_files(path):
    return tuple(os.path.join(path, f) for f in _bbmap_index_filenames)

_star_index_filenames = (
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
def star_index_files(path):
    '''Return a tuple of all STAR index files for path.'''
    return tuple(os.path.join(path, f) for f in _star_index_filenames)

_salmon_index_filenames = (
    'hash.bin',
    'header.json',
    'kintervals.bin',
    'rsd.bin',
    'sa.bin',
    'txpInfo.bin',
    'versionInfo.json',
)
def salmon_index_files(path):
    return tuple(os.path.join(path, f) for f in _salmon_index_filenames)

_hisat2_index_file_infixes = tuple(map(str, range(1,9)))
def hisat2_index_files(path, prefix='index', large=False):
    suffix = "ht2"
    if large:
        suffix += "l"
    fnames = ('{prefix}.{infix}.{suffix}'.format(prefix=prefix, infix=x, suffix=suffix)
              for x in _hisat2_index_file_infixes)
    return tuple(os.path.join(path, f) for f in fnames)

rule all_indices:
    input: BT1=bt1_index_files('BT1_index_hg38.analysisSet', 'index', large=True),
           BT2=bt2_index_files('BT2_index_hg38.analysisSet', 'index', large=True),
           bwa=bwa_index_files('BWA_index_hg38.analysisSet', 'index'),
           bbmap=bbmap_index_files('BBMap_index_hg38.analysisSet'),
           STAR=star_index_files('STAR_index_hg38.analysisSet_knownGene'),
           TH2=tophat2_index_files('TH2_index_hg38.analysisSet_knownGene', 'index', large=True),
           HISAT2=hisat2_index_files('HISAT2_index_grch38_snp_tran', 'index', large=False),
           salmon=salmon_index_files('Salmon_index_hg38.analysisSet_knownGene'),

rule build_bowtie1_index:
    input: genome_fa='{genome_build}.fa'
    output: bt1_index_files('BT1_index_{genome_build}', 'index', large=True)
    params: outdir='BT1_index_{genome_build}',
            basename='BT1_index_{genome_build}/index'
    shell: '''
    mkdir -p {params.outdir:q} && \
        bowtie-build --large-index {input.genome_fa:q} {params.basename:q}
    '''

rule build_bowtie2_index:
    input: genome_fa='{genome_build}.fa'
    output: bt2_index_files('BT2_index_{genome_build}', 'index', large=True)
    params: outdir='BT2_index_{genome_build}',
            basename='BT2_index_{genome_build}/index'
    shell: '''
    mkdir -p {params.outdir:q} && \
        bowtie2-build --large-index {input.genome_fa:q} {params.basename:q}
    '''

rule build_bwa_index:
    input: genome_fa='{genome_build}.fa'
    output: bwa_index_files('BWA_index_{genome_build}', 'index')
    params: outdir='BWA_index_{genome_build}',
            basename='BWA_index_{genome_build}/index'
    shell: '''
    mkdir -p {params.outdir:q} && \
        bwa index -a bwtsw -p {params.basename:q} {input.genome_fa:q}
    '''

# TODO: Make it configurable
BBMAP=os.path.expanduser("~/opt/bbmap/bbmap.sh")

rule build_bbmap_index:
    input: genome_fa='{genome_build}.fa'
    output: bbmap_index_files('BBMap_index_{genome_build}')
    params: outdir='BBMap_index_{genome_build}'
    threads: 16
    run:
        ensure_empty_dir(params.outdir)
        shell('''
        mkdir -p {params.outdir:q} && \
          {BBMAP:q} ref={input.genome_fa:q} \
            path={params.outdir:q} \
            t={threads}
        ''')

# TODO: Place the Log file somewhere better than the root
rule build_star_index:
    input: genome_fa='{genome_build}.fa', transcriptome_gff='{transcriptome}.gff'
    output: star_index_files('STAR_index_{genome_build}_{transcriptome}')
    params: outdir='STAR_index_{genome_build}_{transcriptome}'
    threads: 16
    run:
        ensure_empty_dir(params.outdir)
        shell('''
        mkdir -p {params.outdir:q} && \
          STAR --runMode genomeGenerate \
            --genomeDir {params.outdir:q} \
            --genomeFastaFiles {input.genome_fa:q} \
            --sjdbGTFfile {input.transcriptome_gff:q} \
            --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent \
            --sjdbOverhang 100 \
            --runThreadN {threads:q}
        ''')

rule build_tophat2_index:
    input: bt2_idx=bt2_index_files('BT2_index_{genome_build}', 'index', large=True),
           transcriptome_gff='{transcriptome}.gff'
    output: tophat2_index_files('TH2_index_{genome_build}_{transcriptome}', 'index')
    params: outdir='TH2_index_{genome_build}_{transcriptome}',
            bt2_index_basename=os.path.join('BT2_index_{genome_build}', 'index'),
            basename='TH2_index_{genome_build}_{transcriptome}/index',
    shell: '''
    tophat2 --GTF {input.transcriptome_gff:q} \
      --transcriptome-index={wildcards.basename:q} \
      {params.bt2_index_basename:q}
    '''

# Unless you have 200GB of RAM, you want to use the pre-built index.
rule unpack_hisat2_index:
    input: 'grch38_snp_tran.tar.gz'
    output: hisat2_index_files('HISAT2_index_grch38_snp_tran','index')
    params: outdir='HISAT2_index_grch38_snp_tran'
    shell: '''
    mkdir -p {params.outdir:q} && \
      tar -C {params.outdir:q} --strip-components 1 \
        -s /^genome_snp_tran/index/ \
        -xzvf {input:q}
    '''

rule build_salmon_index:
    input: transcriptome_fa='{transcriptome_build}_transcripts.fa'
    output: salmon_index_files('Salmon_index_{transcriptome_build}')
    params: outdir='Salmon_index_{transcriptome_build}'
    shell: '''
    mkdir -p {params.outdir:q} && \
      salmon index --transcripts {input.transcriptome_fa:q} \
        --index {params.outdir:q} --threads {threads:q} \
        --perfectHash --type quasi
    '''
