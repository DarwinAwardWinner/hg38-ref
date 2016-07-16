#!/bin/bash

cd "$HOME/references/hg38"

JOBS=1
OPTS="-q new -l nodes=1:ppn=${JOBS},mem=30gb,walltime=48:00:00 -N BT1:IDX:hg38 -j oe"

CHROMS="$(echo analysisSet/hg38.analysisSet.chroms/*.fa | perl -lape 's/\s+/,/g')"

qsub <<EOF
#!/bin/bash
#PBS $OPTS
cd \$PBS_O_WORKDIR
module load bowtie/1.1.2
mkdir -p hg38-bt1
bowtie-build "$CHROMS" hg38-bt1/hg38 && \
    echo "Finished building bowtie1 index"
EOF
