#!/bin/bash
#SBATCH --qos=premium
#SBATCH --nodes=6
#SBATCH --cpus-per-task=64
#SBATCH --constraint=haswell
#SBATCH --time=04:30:00
#SBATCH --ntasks=6

cd $HOME/jobs/sky_template_taskfarmer
export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=1

runcommands.sh sky_task.txt
