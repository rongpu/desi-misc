source /global/cfs/cdirs/desi/software/desi_environment.sh master
export REDUXDIR=/global/cfs/cdirs/desi/spectro/redux/daily
export OUTDIR=$SCRATCH/desi/minisv2/4_exp_coadd
mkdir -p $OUTDIR

################# desi_coadd_spectra commands: #################

time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051053/cframe-[brz]0-*.fits $REDUXDIR/exposures/20200219/00051054/cframe-[brz]0-*.fits $REDUXDIR/exposures/20200219/00051055/cframe-[brz]0-*.fits $REDUXDIR/exposures/20200219/00051056/cframe-[brz]0-*.fits  -o $OUTDIR/70003/coadd-0-4exp-subset-0.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051057/cframe-[brz]0-*.fits $REDUXDIR/exposures/20200219/00051058/cframe-[brz]0-*.fits $REDUXDIR/exposures/20200219/00051059/cframe-[brz]0-*.fits $REDUXDIR/exposures/20200219/00051060/cframe-[brz]0-*.fits  -o $OUTDIR/70003/coadd-0-4exp-subset-1.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051053/cframe-[brz]3-*.fits $REDUXDIR/exposures/20200219/00051054/cframe-[brz]3-*.fits $REDUXDIR/exposures/20200219/00051055/cframe-[brz]3-*.fits $REDUXDIR/exposures/20200219/00051056/cframe-[brz]3-*.fits  -o $OUTDIR/70003/coadd-3-4exp-subset-0.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051057/cframe-[brz]3-*.fits $REDUXDIR/exposures/20200219/00051058/cframe-[brz]3-*.fits $REDUXDIR/exposures/20200219/00051059/cframe-[brz]3-*.fits $REDUXDIR/exposures/20200219/00051060/cframe-[brz]3-*.fits  -o $OUTDIR/70003/coadd-3-4exp-subset-1.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051053/cframe-[brz]6-*.fits $REDUXDIR/exposures/20200219/00051054/cframe-[brz]6-*.fits $REDUXDIR/exposures/20200219/00051055/cframe-[brz]6-*.fits $REDUXDIR/exposures/20200219/00051056/cframe-[brz]6-*.fits  -o $OUTDIR/70003/coadd-6-4exp-subset-0.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051057/cframe-[brz]6-*.fits $REDUXDIR/exposures/20200219/00051058/cframe-[brz]6-*.fits $REDUXDIR/exposures/20200219/00051059/cframe-[brz]6-*.fits $REDUXDIR/exposures/20200219/00051060/cframe-[brz]6-*.fits  -o $OUTDIR/70003/coadd-6-4exp-subset-1.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051053/cframe-[brz]7-*.fits $REDUXDIR/exposures/20200219/00051054/cframe-[brz]7-*.fits $REDUXDIR/exposures/20200219/00051055/cframe-[brz]7-*.fits $REDUXDIR/exposures/20200219/00051056/cframe-[brz]7-*.fits  -o $OUTDIR/70003/coadd-7-4exp-subset-0.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051057/cframe-[brz]7-*.fits $REDUXDIR/exposures/20200219/00051058/cframe-[brz]7-*.fits $REDUXDIR/exposures/20200219/00051059/cframe-[brz]7-*.fits $REDUXDIR/exposures/20200219/00051060/cframe-[brz]7-*.fits  -o $OUTDIR/70003/coadd-7-4exp-subset-1.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051053/cframe-[brz]9-*.fits $REDUXDIR/exposures/20200219/00051054/cframe-[brz]9-*.fits $REDUXDIR/exposures/20200219/00051055/cframe-[brz]9-*.fits $REDUXDIR/exposures/20200219/00051056/cframe-[brz]9-*.fits  -o $OUTDIR/70003/coadd-9-4exp-subset-0.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051057/cframe-[brz]9-*.fits $REDUXDIR/exposures/20200219/00051058/cframe-[brz]9-*.fits $REDUXDIR/exposures/20200219/00051059/cframe-[brz]9-*.fits $REDUXDIR/exposures/20200219/00051060/cframe-[brz]9-*.fits  -o $OUTDIR/70003/coadd-9-4exp-subset-1.fits

################# redrock commands: #################

# salloc -N 1 -C haswell -q interactive -t 04:00:00

srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-0-4exp-subset-0.h5 -z $OUTDIR/70003/zbest-0-4exp-subset-0.fits $OUTDIR/70003/coadd-0-4exp-subset-0.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-0-4exp-subset-1.h5 -z $OUTDIR/70003/zbest-0-4exp-subset-1.fits $OUTDIR/70003/coadd-0-4exp-subset-1.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-3-4exp-subset-0.h5 -z $OUTDIR/70003/zbest-3-4exp-subset-0.fits $OUTDIR/70003/coadd-3-4exp-subset-0.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-3-4exp-subset-1.h5 -z $OUTDIR/70003/zbest-3-4exp-subset-1.fits $OUTDIR/70003/coadd-3-4exp-subset-1.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-6-4exp-subset-0.h5 -z $OUTDIR/70003/zbest-6-4exp-subset-0.fits $OUTDIR/70003/coadd-6-4exp-subset-0.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-6-4exp-subset-1.h5 -z $OUTDIR/70003/zbest-6-4exp-subset-1.fits $OUTDIR/70003/coadd-6-4exp-subset-1.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-7-4exp-subset-0.h5 -z $OUTDIR/70003/zbest-7-4exp-subset-0.fits $OUTDIR/70003/coadd-7-4exp-subset-0.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-7-4exp-subset-1.h5 -z $OUTDIR/70003/zbest-7-4exp-subset-1.fits $OUTDIR/70003/coadd-7-4exp-subset-1.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-9-4exp-subset-0.h5 -z $OUTDIR/70003/zbest-9-4exp-subset-0.fits $OUTDIR/70003/coadd-9-4exp-subset-0.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-9-4exp-subset-1.h5 -z $OUTDIR/70003/zbest-9-4exp-subset-1.fits $OUTDIR/70003/coadd-9-4exp-subset-1.fits
