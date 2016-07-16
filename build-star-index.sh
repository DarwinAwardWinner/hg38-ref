#!/bin/bash

cd "$HOME/references/hg38"

JOBS=8
OPTS="-q new -l nodes=1:ppn=${JOBS},mem=40gb,walltime=48:00:00 -N STAR:IDX:hg38"

qsub <<EOF
#!/bin/bash
#PBS $OPTS
cd \$PBS_O_WORKDIR
module load star
mkdir -p hg38-star
STAR --runMode genomeGenerate --genomeDir ./hg38-star \
  --genomeFastaFiles ./analysisSet/hg38.analysisSet.chroms/*.fa \
  --sjdbGTFfile ./knownGene.gff \
  --sjdbGTFfeatureExon exon \ --sjdbGTFtagExonParentTranscript Parent \
  --sjdbOverhang 100 \
  --runThreadN ${JOBS} \
    echo "Finished building STAR genome"
EOF
