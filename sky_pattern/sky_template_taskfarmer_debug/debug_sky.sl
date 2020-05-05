#!/bin/bash
#SBATCH --qos=debug
#SBATCH --nodes=7
#SBATCH --cpus-per-task=64
#SBATCH --constraint=haswell
#SBATCH --time=00:30:00
#SBATCH --ntasks=7

cd $HOME/jobs/sky_template_taskfarmer_debug
export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=32

runcommands.sh sky_task.txt
