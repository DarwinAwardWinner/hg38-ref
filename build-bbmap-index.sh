#!/bin/bash


cd "$HOME/references/hg38"

JOBS=8
OPTS="-q new -l nodes=1:ppn=${JOBS},mem=30gb,walltime=48:00:00 -N BBMAP:IDX:hg38 -j oe"

CHROMS="$(echo analysisSet/hg38.analysisSet.chroms/*.fa)"

BBPATH="$HOME/opt/bbmap"

qsub <<EOF
#!/bin/bash
#PBS $OPTS
cd \$PBS_O_WORKDIR
mkdir -p hg38-bbmap
cat $CHROMS | \
    $BBPATH/bbmap.sh ref=stdin.fa path=hg38-bbmap t=$JOBS && \
    echo "Finished building BBMap index"
EOF
