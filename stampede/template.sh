#!/bin/bash

echo "Started $(date)"

sh run.sh ${INPUT_DIR} ${FORWARD_PRIMER} ${REVERSE_PRIMER}

echo "Ended $(date)"
exit 0
