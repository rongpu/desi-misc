#!/bin/bash
#SBATCH --qos=premium
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --constraint=haswell
#SBATCH --time=8:00:00
#SBATCH --array=0-3

source /project/projectdirs/desi/software/desi_environment.sh 19.12
export PYTHONPATH=/global/homes/r/rongpu/modules/lib/python3.6/site-packages/fitsio-1.1.0-py3.6-linux-x86_64.egg

srun python python sky_template_fitting.py 6 $SLURM_ARRAY_TASK_ID
