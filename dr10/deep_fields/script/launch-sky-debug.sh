#!/bin/bash
#SBATCH --qos=debug
#SBATCH --nodes=4
#SBATCH --ntasks-per-node 1
#SBATCH --constraint=cpu
#SBATCH --time=00:10:00
#SBATCH --account=desi

srun --wait=0 payload_sky.sh splinesky_cosmos_commands.sh
