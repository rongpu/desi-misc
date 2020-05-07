#!/bin/bash
#SBATCH --qos=debug
#SBATCH --nodes=4
#SBATCH --cpus-per-task=64
#SBATCH --constraint=haswell
#SBATCH --time=00:10:00
#SBATCH --ntasks=4

cd $HOME/temp/taskfarmer_debug
export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=1

runcommands.sh task.txt
