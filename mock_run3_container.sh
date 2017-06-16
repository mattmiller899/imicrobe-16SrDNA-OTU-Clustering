#!/bin/bash

rm -rf mock_run3/output
mkdir mock_run3/output

time singularity run singularity/imicrobe-16SrDNA-OTU-Clustering.img \
  --input-dir ~/host/project/imicrobe/apps/imicrobe-16SrDNA-OTU-Clustering/mock_run3/input \
  --output-dir ~/host/project/imicrobe/apps/imicrobe-16SrDNA-OTU-Clustering/mock_run3/output \
  --core-count 1
