'''Rules for downloading specific files from the internet.'''

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

HTTP = HTTPRemoteProvider()
FTP = FTPRemoteProvider()

rule get_hg38_analysisSet_twoBit:
    input: HTTP.remote('hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/analysisSet/hg38.analysisSet.2bit', insecure=True, static=True)
    output: 'hg38.analysisSet.2bit'
    shell: 'mv {input:q} {output:q}'

rule get_knownGene_gtf:
    output: 'knownGene.gtf'
    shell: 'genePredToGtf -addComments -utr hg38 knownGene {output:q}'

rule get_gencode_annotation_gff:
    input: FTP.remote('ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_{release}/gencode.v{release}.chr_patch_hapl_scaff.annotation.gff3.gz', static=True)
    output: 'gencode.v{release,\\d+}_raw.gff3'
    shell: 'zcat < {input:q} > {output:q}'

# http://uswest.ensembl.org/info/data/ftp/index.html
rule get_ensembl_annotation_gff:
    input: FTP.remote('ftp.ensembl.org/pub/release-{release}/gff3/homo_sapiens/Homo_sapiens.GRCh38.{release}.chr_patch_hapl_scaff.gff3.gz', static=True)
    output: 'ensembl.{release,\\d+}_raw.gff3'
    shell: 'zcat < {input:q} > {output:q}'

rule get_chrom_mapping:
    input: HTTP.remote('raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/{genome_build}_{from_ids}2{to_ids}.txt', static=True)
    output: 'chrom_mapping_{genome_build}_{from_ids}2{to_ids}.txt'
    shell: 'mv {input:q} {output:q}'

# Need to download the HISAT2 index instead of building it because
# building requires 200GB of RAM.
rule get_hisat2_index_tar:
    input: FTP.remote('ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_snp_tran.tar.gz', static=True)
    output: 'grch38_snp_tran.tar.gz'
    shell: 'mv {input:q} {output:q}'
