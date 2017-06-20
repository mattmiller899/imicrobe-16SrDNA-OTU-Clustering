#!/bin/bash

rm -rf mock_run3/output
mkdir -p mock_run3/output

time pipeline \
  --input-dir ~/host/project/imicrobe/apps/imicrobe-16SrDNA-OTU-Clustering/mock_run3/input \
  --output-dir ~/host/project/imicrobe/apps/imicrobe-16SrDNA-OTU-Clustering/mock_run3/output \
  --core-count 1
