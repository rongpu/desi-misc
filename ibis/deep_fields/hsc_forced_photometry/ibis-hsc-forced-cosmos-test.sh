#! /bin/bash

brick=$1

export COSMO=/dvs_ro/cfs/cdirs/cosmo

# survey-dir contains HSC CCD list and the HSC coadds for outlier rejection
export survey_dir=/dvs_ro/cfs/cdirs/cosmo/work/users/dstn/ODIN/2024-a/hsc-co-4

catdir=$COSMO/work/legacysurvey/ibis/reductions/ibis-dr1
outdir=$SCRATCH/ibis-tractor/ibis-dr1-hsc-wide-forced-test

export LARGEGALAXIES_CAT=$COSMO/work/legacysurvey/sga/2025/SGA2025-ellipse-v0.80-dr11-south.kd.fits

mkdir -p ${outdir}/logs-forced
echo Logging to ${outdir}/logs-forced/${brick}.log

# Don't add ~/.local/ to Python's sys.path
export PYTHONNOUSERSITE=1
# Force MKL single-threaded
# https://software.intel.com/en-us/articles/using-threaded-intel-mkl-in-multi-thread-application
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
# To avoid problems with MPI and Python multiprocessing
export MPICH_GNI_FORK_MODE=FULLCOPY
export KMP_AFFINITY=disabled

# # Config directory nonsense
export TMPCACHE=$(mktemp -d)
mkdir $TMPCACHE/cache
mkdir $TMPCACHE/config
# astropy
export XDG_CACHE_HOME=$TMPCACHE/cache
export XDG_CONFIG_HOME=$TMPCACHE/config
mkdir $XDG_CACHE_HOME/astropy
cp -r $HOME/.astropy/cache $XDG_CACHE_HOME/astropy
mkdir $XDG_CONFIG_HOME/astropy
cp -r $HOME/.astropy/config $XDG_CONFIG_HOME/astropy
# matplotlib
export MPLCONFIGDIR=$TMPCACHE/matplotlib
mkdir $MPLCONFIGDIR
cp -r $HOME/.config/matplotlib $MPLCONFIGDIR

# # cosmos:
# --bands g,r2,i2,z,y \
# # xmm-lss:
# --bands g,r,i,i2,z,y \

# export LEGACYPIPE_DIR=/src/legacypipe/py
# use the newer local copy (8f8f3f0) to implement the outliers_mask bug fix
export LEGACYPIPE_DIR=/global/cfs/cdirs/desicollab/users/rongpu/moregit/legacypipe/py
export PYTHONPATH="$LEGACYPIPE_DIR:$PYTHONPATH"

python -O $LEGACYPIPE_DIR/legacypipe/forced_photom_brickwise.py \
       --brick ${brick} \
       --bands g,r2,i2,z,y \
       --survey-dir ${survey_dir} \
       --catalog-dir ${catdir} \
       --outdir ${outdir} \
       --threads 32 \
       >> ${outdir}/logs-forced/${brick}.log 2>&1

# Save the return value from the python command -- otherwise we
# exit 0 because the rm succeeds!
status=$?

# /Config directory nonsense
rm -R $TMPCACHE

exit $status
