#!/bin/bash
#SBATCH --qos=debug
#SBATCH --nodes=8
#SBATCH --ntasks-per-node 1
#SBATCH --constraint=cpu
#SBATCH --time=00:30:00
#SBATCH --account=m3592

srun --wait=0 payload-runbrick.sh runbrick-commands-edge_bricks_first.txt
