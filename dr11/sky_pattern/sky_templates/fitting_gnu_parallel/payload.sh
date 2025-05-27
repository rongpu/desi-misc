#!/bin/bash

source /global/common/software/desi/desi_environment.sh 24.11

if [[ -z "${SLURM_NODEID}" ]]; then
    echo "need \$SLURM_NODEID set"
    exit
fi
if [[ -z "${SLURM_NNODES}" ]]; then
    echo "need \$SLURM_NNODES set"
    exit
fi

# ntasks=$(($(wc -l < tasks.txt)+1))

cat $1 |                                               \
awk -v NNODE="$SLURM_NNODES" -v NODEID="$SLURM_NODEID" \
'NR % NNODE == NODEID' |                               \
# parallel --jobs 1 ./task.sh $ntasks {} $SLURM_NNODES $SLURM_NODEID
parallel --jobs 1 python sky_template_fitting.py {}
