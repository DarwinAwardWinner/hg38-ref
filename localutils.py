'''Utility functions for use in Snakemake workflows.'''

import os
import os.path
import shutil
import sys

from subprocess import check_output

# TODO: Make this configurable
BBMAP=os.path.expanduser("~/opt/bbmap/bbmap.sh")

def check_output_decoded(*args, **kwargs):
    return check_output(*args, **kwargs).decode(sys.getdefaultencoding())

def ensure_empty_dir(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)

def read_chrom_mapping(filename):
    mapping = {}
    with open(filename, "r") as infile:
        for line in infile:
            (from_id, to_id) = line.split("\t")
            to_id = to_id.strip()
            if to_id:
                mapping[from_id] = to_id
    return mapping

# Series of functions for listing the expected index files of various
# aligner tools

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
_tophat2_index_suffixes = ('.gff', '.fa', '.fa.tlst', '.ver')
def tophat2_index_files(path, prefix='index'):
    fnames = ('{prefix}{suffix}'.format(prefix=prefix, suffix=s)
              for s in _tophat2_index_suffixes)
    paths = tuple(os.path.join(path, f) for f in fnames)
    return paths + bt2_index_files(path, prefix, large=False)

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
    'indexing.log',
    'quasi_index.log',
    'rsd.bin',
    'sa.bin',
    'txpInfo.bin',
    'versionInfo.json',
)
# salmon index --perfectHash
_salmon_perfectHash_index_filenames = (
    'hash_info.bph',
    'hash_info.val',
    'header.json',
    'indexing.log',
    'quasi_index.log',
    'rsd.bin',
    'sa.bin',
    'txpInfo.bin',
    'versionInfo.json',
)

def salmon_index_files(path, perfectHash=False):
    if perfectHash:
        return tuple(os.path.join(path, f) for f in _salmon_perfectHash_index_filenames)
    else:
        return tuple(os.path.join(path, f) for f in _salmon_index_filenames)

_hisat2_index_file_infixes = tuple(map(str, range(1,9)))
def hisat2_index_files(path, prefix='index', large=False):
    suffix = "ht2"
    if large:
        suffix += "l"
    fnames = ('{prefix}.{infix}.{suffix}'.format(prefix=prefix, infix=x, suffix=suffix)
              for x in _hisat2_index_file_infixes)
    return tuple(os.path.join(path, f) for f in fnames)
