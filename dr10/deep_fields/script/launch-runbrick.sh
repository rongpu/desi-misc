#!/bin/bash
#SBATCH --qos=regular
#SBATCH --nodes=6
#SBATCH --ntasks-per-node 1
#SBATCH --constraint=cpu
#SBATCH --time=12:00:00
#SBATCH --account=desi

srun --wait=0 payload-runbrick.sh tasks-runbrick-cosmos.sh
