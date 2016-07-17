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

rule unpack_hisat2_index:
    input: 'grch38_snp_tran.tar.gz'
    output: expand('HISAT2_index_grch38_snp_tran/index.{num}.ht2', list(range(1, 9)))
    params: outdir='HISAT2_index_grch38_snp_tran',
            src_basename='genome_snp_tran',
            dest_basename='index',
    shell: '''
    mkdir -p {params.outdir:q} && \
      tar -C {params.outdir:q} --strip-components 1 \
        -s /^genome_snp_tran/index/ \
        -xzvf {input:q}
    '''
