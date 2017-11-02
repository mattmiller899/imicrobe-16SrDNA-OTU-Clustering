#!/bin/bash

module load singularity

INPUT_DIR=$1
OUTPUT_DIR=$2

echo "starting directory : `pwd`"
echo "`ls -l`"
echo "input directory    : ${INPUT_DIR}"
echo "output directory   : ${OUTPUT_DIR}"

export LAUNCHER_DIR="$HOME/src/launcher"
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_WORKDIR=${OUTPUT_DIR}
export LAUNCHER_RMI=SLURM

export LAUNCHER_JOB_FILE=`pwd`/launcher_jobfile_${SLURM_JOB_ID}
echo ${LAUNCHER_JOB_FILE}

xz --decompress imicrobe-16SrDNA-OTU-Clustering.img.xz

echo "`ls -l`"

time singularity run imicrobe-16SrDNA-OTU-Clustering.img \
  --input-dir ${INPUT_DIR} \
  --work-dir ${OUTPUT_DIR} \
  --core-count 2 \
  --cutadapt-min-length 100 \
  --pear-min-overlap 200\
  --pear-max-assembly-length 270 \
  --pear-min-assembly-length 220 \
  --vsearch-filter-maxee 1 \
  --vsearch-filter-trunclen 245 \
  --vsearch-derep-minuniquesize 3

#sleep 10
#export LAUNCHER_PPN=2

#$LAUNCHER_DIR/paramrun
echo "Ended launcher"
