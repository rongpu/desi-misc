#!/bin/bash
#SBATCH --qos=debug
#SBATCH --nodes=8
#SBATCH --ntasks-per-node 1
#SBATCH --constraint=cpu
#SBATCH --time=0:30:00
#SBATCH --account=m3592

srun --wait=0 slurm_payload.sh runbrick-commands.txt 4 5

