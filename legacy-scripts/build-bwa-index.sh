#!/bin/bash


cd "$HOME/references/hg38"

JOBS=1
OPTS="-q new -l nodes=1:ppn=${JOBS},mem=30gb,walltime=48:00:00 -N BWA:IDX:hg38 -j oe"

CHROMS="$(echo analysisSet/hg38.analysisSet.chroms/*.fa)"

qsub <<EOF
#!/bin/bash
#PBS $OPTS
cd \$PBS_O_WORKDIR
module load bwa
mkdir -p hg38-bwa
cat $CHROMS | \
    bwa index -a bwtsw -p hg38-bwa/hg38 /dev/stdin && \
    echo "Finished building bwa index"
EOF
