#!/bin/bash
#SBATCH --qos=regular
#SBATCH --nodes=10
#SBATCH --ntasks-per-node 1
#SBATCH --constraint=cpu
#SBATCH --time=8:00:00
#SBATCH --account=m3592

srun --wait=0 slurm_payload.sh runbrick-commands.txt 4 5

