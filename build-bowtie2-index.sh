#!/bin/bash

cd "$HOME/references/hg38"

JOBS=1
OPTS="-q new -l nodes=1:ppn=${JOBS},mem=30gb,walltime=48:00:00 -N BT2:IDX:hg38 -j oe"

CHROMS="$(echo analysisSet/hg38.analysisSet.chroms/*.fa | perl -lape 's/\s+/,/g')"

qsub <<EOF
#!/bin/bash
#PBS $OPTS
cd \$PBS_O_WORKDIR
module load bowtie2/2.2.9
mkdir -p hg38-bt2
bowtie2-build "$CHROMS" hg38-bt2/hg38 && \
    echo "Finished building bowtie2 index"
EOF
