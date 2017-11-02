#!/bin/bash

#rm -rf mock_run3/output_pipeline
mkdir -p mock_run3/output_pipeline

time pipeline \
  --input-dir ~/host/project/imicrobe/apps/imicrobe-16SrDNA-OTU-Clustering/mock_run3/input \
  --work-dir ~/host/project/imicrobe/apps/imicrobe-16SrDNA-OTU-Clustering/mock_run3/output_pipeline \
  --core-count 2 \
  --cutadapt-min-length 100 \
  --pear-min-overlap 200\
  --pear-max-assembly-length 270 \
  --pear-min-assembly-length 220 \
  --uchime-ref-db-fp ~/host/project/silva/SILVA_128_SSURef_Nr99_tax_silva.fasta.gz \
  --vsearch-filter-maxee 1 \
  --vsearch-filter-trunclen 245 \
  --vsearch-derep-minuniquesize 3
