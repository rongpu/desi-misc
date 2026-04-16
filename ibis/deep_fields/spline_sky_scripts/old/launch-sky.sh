#!/bin/bash
#SBATCH --qos=regular
#SBATCH --nodes=4
#SBATCH --ntasks-per-node 1
#SBATCH --constraint=cpu
#SBATCH --time=4:00:00
#SBATCH --account=desi

srun --wait=0 payload-sky.sh splinesky_cosmos_commands.sh
