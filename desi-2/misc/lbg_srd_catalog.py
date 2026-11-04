# one-time setup (the README suggests DESI2_SRD_BASE as the parent directory of the repo)
export PYTHONPATH=/global/cfs/cdirs/desicollab/users/rongpu/moregit/desi-2-srd/code:$PYTHONPATH
export DESI2_SRD_BASE=/global/cfs/cdirs/desicollab/users/rongpu/moregit
export DESI2_SRD_DATA=$DESI2_SRD_BASE/desi-2-srd/data
export PYTHONPATH=$DESI2_SRD_BASE/desi-2-srd/code:$PYTHONPATH

# needs scikit-learn, fitsio, astropy, pandas, matplotlib
# use a conda/desi environment that has them (the README mentions: conda install scikit-learn; pip install fitsio)

mkdir -p $SCRATCH/lbg_rf_new
cd $DESI2_SRD_BASE/desi-2-srd/code/lbg_target_selection

python random_forest_lbg.py --target-density 900 --lsst-years 4 \
    --output-filename $SCRATCH/lbg_rf_new/toto.json \
    --output-target-catalog $SCRATCH/lbg_rf_new/lbg-target-catalog.fits


python random_forest_lbg_deepmag.py --target-density 900 --lsst-years 4 \
    --output-filename $SCRATCH/lbg_rf_new/toto-deepmag.json \



python random_forest_lbg.py --target-density 900 --lsst-years 4 \
    --input-catalog-filename /pscratch/sd/r/rongpu/tmp/COSMOS_11bands-SExtractor-Lephare-HSC_PDR3_dud-clean.fits \
    --output-filename $SCRATCH/lbg_rf_new/toto.json \
    --output-target-catalog $SCRATCH/lbg_rf_new/lbg-target-catalog.fits

