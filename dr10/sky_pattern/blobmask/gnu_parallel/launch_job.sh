#!/bin/bash
#SBATCH --qos=regular
#SBATCH --nodes=40
#SBATCH --ntasks-per-node 1
#SBATCH --constraint=haswell
#SBATCH --time=8:00:00
#SBATCH --account=m3592

srun --wait=0 payload.sh tasks.txt