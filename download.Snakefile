'''Rules for downloading specific files from the internet.'''

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

HTTP = HTTPRemoteProvider()
FTP = FTPRemoteProvider()

rule get_hg38_analysisSet_twoBit:
    input: HTTP.remote('hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/analysisSet/hg38.analysisSet.2bit', insecure=True)
    output: 'hg38.analysisSet.2bit'
    shell: 'mv {input:q} {output:q}'

rule get_knownGene_gtf:
    output: 'knownGene.gtf'
    shell: 'genePredToGtf -addComments -utr hg38 knownGene {output:q}'

rule get_hisat2_index_tar:
    input: FTP.remote('ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_snp_tran.tar.gz')
    output: 'grch38_snp_tran.tar.gz'
    shell: 'mv {input:q} {output:q}'
