'''Rules for building indices for various alignment tools.'''

import os.path

from localutils import *
from tool_versions import *

rule all_indices:
    input: BT1=bt1_index_files('BT1_index_hg38.analysisSet', 'index', large=True),
           BT2=bt2_index_files('BT2_index_hg38.analysisSet', 'index', large=True),
           bwa=bwa_index_files('BWA_index_hg38.analysisSet', 'index'),
           bbmap=bbmap_index_files('BBMap_index_hg38.analysisSet'),
           STAR=star_index_files('STAR_index_hg38.analysisSet_knownGene'),
           TH2=tophat2_index_files('TH2_index_hg38.analysisSet_knownGene', 'index'),
           HISAT2=hisat2_index_files('HISAT2_index_grch38_snp_tran', 'index', large=False),
           salmon=salmon_index_files('Salmon_index_hg38.analysisSet_knownGene'),

rule build_bowtie1_index:
    input: genome_fa='{genome_build}.fa'
    output: bt1_index_files('BT1_index_{genome_build}', 'index', large=True)
    params: outdir='BT1_index_{genome_build}',
            basename='BT1_index_{genome_build}/index'
    version: BOWTIE1_VERSION
    run:
        ensure_empty_dir(params.outdir)
        shell('''bowtie-build --large-index {input.genome_fa:q} {params.basename:q}''')

rule build_bowtie2_index:
    input: genome_fa='{genome_build}.fa'
    output: bt2_index_files('BT2_index_{genome_build}', 'index', large=True)
    params: outdir='BT2_index_{genome_build}',
            basename='BT2_index_{genome_build}/index'
    version: BOWTIE2_VERSION
    run:
        ensure_empty_dir(params.outdir)
        shell('''bowtie2-build --large-index {input.genome_fa:q} {params.basename:q}''')

rule build_bwa_index:
    input: genome_fa='{genome_build}.fa'
    output: bwa_index_files('BWA_index_{genome_build}', 'index')
    params: outdir='BWA_index_{genome_build}',
            basename='BWA_index_{genome_build}/index'
    version: BWA_VERSION
    run:
        ensure_empty_dir(params.outdir)
        shell('''bwa index -a bwtsw -p {params.basename:q} {input.genome_fa:q}''')

rule build_bbmap_index:
    input: genome_fa='{genome_build}.fa'
    output: bbmap_index_files('BBMap_index_{genome_build}')
    params: outdir='BBMap_index_{genome_build}'
    version: BBMAP_VERSION
    threads: 16
    run:
        ensure_empty_dir(params.outdir)
        shell('''
        {BBMAP:q} ref={input.genome_fa:q} \
          path={params.outdir:q} \
          t={threads}
        ''')

# TODO: Place the Log file somewhere better than the root
rule build_star_index:
    input: genome_fa='{genome_build}.fa', transcriptome_gff='{transcriptome}.gff3'
    output: index_files=star_index_files('STAR_index_{genome_build}_{transcriptome}'),
            log_file='STAR_index_{genome_build}_{transcriptome}/Log.out'
    params: outdir='STAR_index_{genome_build}_{transcriptome}'
    threads: 16
    version: STAR_VERSION
    run:
        ensure_empty_dir(params.outdir)
        # This is to run STAR from the index's subdirectory to ensure
        # that the log file is produced there.
        basedir = os.path.abspath(os.path.curdir)
        shell('''
        cd {params.outdir:q} && \
          STAR --runMode genomeGenerate \
            --genomeDir ./ \
            --genomeFastaFiles {basedir:q}/{input.genome_fa:q} \
            --sjdbGTFfile {basedir:q}/{input.transcriptome_gff:q} \
            --sjdbGTFfeatureExon exon \
            --sjdbGTFtagExonParentTranscript Parent \
            --sjdbGTFtagExonParentGene gene_id \
            --sjdbOverhang 100 \
            --runThreadN {threads:q}
        ''')

rule build_tophat2_index:
    input: bt2_idx=bt2_index_files('BT2_index_{genome_build}', 'index', large=True),
           transcriptome_gff='{transcriptome}.gff3'
    output: tophat2_index_files('TH2_index_{genome_build}_{transcriptome}', 'index')
    params: outdir='TH2_index_{genome_build}_{transcriptome}',
            bt2_index_basename=os.path.join('BT2_index_{genome_build}', 'index'),
            basename='TH2_index_{genome_build}_{transcriptome}/index',
    # Can tophat do parallel building of a transcriptome index?
    threads: 1
    version: TOPHAT2_VERSION
    run:
        ensure_empty_dir(params.outdir)
        shell('''
        tophat2 --GTF {input.transcriptome_gff:q} \
          --output-dir={params.outdir:q} \
          --transcriptome-index={params.basename:q} \
          {params.bt2_index_basename:q}
        ''')

# Unless you have 200GB of RAM, you want to use the pre-built index.
rule unpack_hisat2_index:
    input: 'grch38_snp_tran.tar.gz'
    output: hisat2_index_files('HISAT2_index_grch38_snp_tran','index')
    params: outdir='HISAT2_index_grch38_snp_tran'
    # This is a downloaded index, so there is no associated hisat2 version
    run:
        ensure_empty_dir(params.outdir)
        shell('''
        tar -C {params.outdir:q} --strip-components 1 \
          --transform 's/genome_snp_tran/index/' \
          -xzvf {input:q}
        ''')

rule build_salmon_index:
    input: transcriptome_fa='{transcriptome_build}_transcripts.fa'
    output: salmon_index_files('Salmon_index_{transcriptome_build}')
    params: outdir='Salmon_index_{transcriptome_build}'
    version: SALMON_VERSION
    run:
        ensure_empty_dir(params.outdir)
        shell('''
        salmon index --transcripts {input.transcriptome_fa:q} \
          --index {params.outdir:q} --threads {threads:q} \
          --perfectHash --type quasi
        ''')

# Salmon genemap files, required to get gene-level quantification
# results
rule make_salmon_genemap_from_knownGene:
    output: 'Salmon_genemap_hg38.analysisSet_knownGene/genemap.txt'
    run:
        from rpy2.robjects import r
        from rpy2.robjects import globalenv as r_env
        import rpy2.rinterface
        rpy2.rinterface.set_writeconsole_warnerror(lambda x: sys.stderr.write(x))
        r_env['outfile'] = output[0]
        r('''
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)
        txg <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, 'gene')
        txg_df <- data.frame(TXID=unlist(txg)$tx_name, GeneID=rep(names(txg), lengths(txg)))
        write.table(txg_df, outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
        ''')

rule make_salmon_genemap_from_ensembl:
    input: txdb='TxDb.Hsapiens.ensembl.hg38.v{release}.sqlite3'
    output: 'Salmon_index_hg38.analysisSet_ensembl.{release}/genemap.txt'
    run:
        from rpy2.robjects import r
        from rpy2.robjects import globalenv as r_env
        import rpy2.rinterface
        rpy2.rinterface.set_writeconsole_warnerror(lambda x: sys.stderr.write(x))
        r_env['input.txdb'] = input.txdb
        r_env['outfile'] = output[0]
        r('''
        library(GenomicFeatures)
        txdb <- loadDb(input.txdb)
        txg <- transcriptsBy(txdb, 'gene')
        txg_df <- data.frame(TXID=unlist(txg)$tx_name, GeneID=rep(names(txg), lengths(txg)))
        write.table(txg_df, outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
        ''')
