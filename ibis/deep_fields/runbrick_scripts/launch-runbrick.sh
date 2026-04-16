#!/bin/bash
#SBATCH --qos=regular
#SBATCH --nodes=20
#SBATCH --ntasks-per-node 1
#SBATCH --constraint=cpu
#SBATCH --time=14:00:00
#SBATCH --account=m3592

srun --wait=0 payload-runbrick.sh runbrick-commands-remain.txt
