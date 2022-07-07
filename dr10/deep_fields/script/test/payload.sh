#!/bin/bash

module load parallel

if [[ -z "${SLURM_NODEID}" ]]; then
    echo "need \$SLURM_NODEID set"
    exit
fi
if [[ -z "${SLURM_NNODES}" ]]; then
    echo "need \$SLURM_NNODES set"
    exit
fi

cat $1 |                                               \
awk -v NNODE="$SLURM_NNODES" -v NODEID="$SLURM_NODEID" \
'NR % NNODE == NODEID' |                               \
parallel --jobs 24 shifter --image docker:legacysurvey/legacypipe:DR10.0.0 ./runbrick-cosmos-test-100.sh {}
