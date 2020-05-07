#!/bin/bash
#SBATCH --qos=premium
#SBATCH --nodes=7
#SBATCH --cpus-per-task=64
#SBATCH --constraint=haswell
#SBATCH --time=10:00:00
#SBATCH --ntasks=7

cd $HOME/jobs/sky_template_taskfarmer_no_z
export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=1

runcommands.sh sky_task.txt
