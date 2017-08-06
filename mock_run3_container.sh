#!/bin/bash

#rm -rf mock_run3/output_container
mkdir mock_run3/output_container

time singularity run singularity/imicrobe-16SrDNA-OTU-Clustering.img \
  --input-dir ./mock_run3/input \
  --work-dir ./mock_run3/output_container \
  --core-count 2 \
  --cutadapt-min-length 100 \
  --pear-min-overlap 200\
  --pear-max-assembly-length 270 \
  --pear-min-assembly-length 220 \
  --vsearch-filter-maxee 1 \
  --vsearch-filter-trunclen 245 \
  --vsearch-derep-minuniquesize 3
