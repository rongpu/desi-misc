#!/bin/bash
#SBATCH --qos=premium
#SBATCH --nodes=10
#SBATCH --cpus-per-task=64
#SBATCH --constraint=haswell
#SBATCH --time=10:00:00
#SBATCH --ntasks=10

cd $HOME/jobs/blobmask_taskfarmer
export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=1

runcommands.sh blobmask_tasks_g.txt
# runcommands.sh blobmask_tasks_r.txt
# runcommands.sh blobmask_tasks_z.txt
