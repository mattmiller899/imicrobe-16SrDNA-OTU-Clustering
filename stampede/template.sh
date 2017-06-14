#!/bin/bash

echo "Started $(date)"

sh run.sh ${INPUT_DIR} ${FORWARD_PRIMER} ${REVERSE_PRIMER} ${MIN_OVERLAP}

echo "Ended $(date)"
exit 0
