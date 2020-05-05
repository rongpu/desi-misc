#!/bin/bash
#SBATCH --qos=debug
#SBATCH --nodes=16
#SBATCH --cpus-per-task=64
#SBATCH --constraint=haswell
#SBATCH --time=00:30:00
#SBATCH --ntasks=16

cd $HOME/jobs/sky_template_taskfarmer
export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=1

runcommands.sh sky_task.txt
