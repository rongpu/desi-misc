#!/bin/bash
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --constraint=haswell
#SBATCH --time=10:00:00
#SBATCH --array=0-5

srun python save_ccd_blob_mask_sbatch.py $SLURM_ARRAY_TASK_ID

