#!/bin/bash

cd "$HOME/references/hg38"

JOBS=1
OPTS="-q new -l nodes=1:ppn=${JOBS},mem=30gb,walltime=48:00:00 -N TH:IDX:hg38 -j oe"

qsub <<EOF
#!/bin/bash
#PBS $OPTS
cd \$PBS_O_WORKDIR
module load bowtie2/2.2.9 tophat/2.1.0
mkdir -p hg38-tophat
tophat --GTF knownGene.gtf --transcriptome-index hg38-tophat/hg38 ./hg38-bt2/hg38 && \
    echo "Finished building tophat index"
EOF
