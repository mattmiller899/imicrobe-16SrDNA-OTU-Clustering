#!/bin/bash

rm -rf mock_run3/pipeline_output
mkdir -p mock_run3/pipeline_output

time pipeline \
  --input-dir ~/host/project/imicrobe/apps/imicrobe-16SrDNA-OTU-Clustering/mock_run3/input \
  --work-dir ~/host/project/imicrobe/apps/imicrobe-16SrDNA-OTU-Clustering/mock_run3/pipeline_output \
  --core-count 1 \
  --cutadapt-min-length 100 \
  --pear-min-overlap 200\
  --pear-max-assembly-length 270 \
  --pear-min-assembly-length 220 \
  --vsearch-filter-maxee 1 \
  --vsearch-filter-trunclen 245
