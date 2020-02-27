source /global/cfs/cdirs/desi/software/desi_environment.sh master
export REDUXDIR=/global/cfs/cdirs/desi/spectro/redux/daily
export OUTDIR=$SCRATCH/desi/minisv2/single_exp_coadd
mkdir -p $OUTDIR

################# desi_coadd_spectra commands: #################

time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051053/cframe-[brz]0-*.fits -o $OUTDIR/70003/coadd-0-00051053.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051053/cframe-[brz]3-*.fits -o $OUTDIR/70003/coadd-3-00051053.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051053/cframe-[brz]6-*.fits -o $OUTDIR/70003/coadd-6-00051053.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051053/cframe-[brz]7-*.fits -o $OUTDIR/70003/coadd-7-00051053.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051053/cframe-[brz]9-*.fits -o $OUTDIR/70003/coadd-9-00051053.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051054/cframe-[brz]0-*.fits -o $OUTDIR/70003/coadd-0-00051054.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051054/cframe-[brz]3-*.fits -o $OUTDIR/70003/coadd-3-00051054.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051054/cframe-[brz]6-*.fits -o $OUTDIR/70003/coadd-6-00051054.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051054/cframe-[brz]7-*.fits -o $OUTDIR/70003/coadd-7-00051054.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051054/cframe-[brz]9-*.fits -o $OUTDIR/70003/coadd-9-00051054.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051055/cframe-[brz]0-*.fits -o $OUTDIR/70003/coadd-0-00051055.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051055/cframe-[brz]3-*.fits -o $OUTDIR/70003/coadd-3-00051055.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051055/cframe-[brz]6-*.fits -o $OUTDIR/70003/coadd-6-00051055.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051055/cframe-[brz]7-*.fits -o $OUTDIR/70003/coadd-7-00051055.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051055/cframe-[brz]9-*.fits -o $OUTDIR/70003/coadd-9-00051055.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051056/cframe-[brz]0-*.fits -o $OUTDIR/70003/coadd-0-00051056.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051056/cframe-[brz]3-*.fits -o $OUTDIR/70003/coadd-3-00051056.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051056/cframe-[brz]6-*.fits -o $OUTDIR/70003/coadd-6-00051056.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051056/cframe-[brz]7-*.fits -o $OUTDIR/70003/coadd-7-00051056.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051056/cframe-[brz]9-*.fits -o $OUTDIR/70003/coadd-9-00051056.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051057/cframe-[brz]0-*.fits -o $OUTDIR/70003/coadd-0-00051057.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051057/cframe-[brz]3-*.fits -o $OUTDIR/70003/coadd-3-00051057.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051057/cframe-[brz]6-*.fits -o $OUTDIR/70003/coadd-6-00051057.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051057/cframe-[brz]7-*.fits -o $OUTDIR/70003/coadd-7-00051057.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051057/cframe-[brz]9-*.fits -o $OUTDIR/70003/coadd-9-00051057.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051058/cframe-[brz]0-*.fits -o $OUTDIR/70003/coadd-0-00051058.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051058/cframe-[brz]3-*.fits -o $OUTDIR/70003/coadd-3-00051058.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051058/cframe-[brz]6-*.fits -o $OUTDIR/70003/coadd-6-00051058.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051058/cframe-[brz]7-*.fits -o $OUTDIR/70003/coadd-7-00051058.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051058/cframe-[brz]9-*.fits -o $OUTDIR/70003/coadd-9-00051058.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051059/cframe-[brz]0-*.fits -o $OUTDIR/70003/coadd-0-00051059.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051059/cframe-[brz]3-*.fits -o $OUTDIR/70003/coadd-3-00051059.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051059/cframe-[brz]6-*.fits -o $OUTDIR/70003/coadd-6-00051059.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051059/cframe-[brz]7-*.fits -o $OUTDIR/70003/coadd-7-00051059.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051059/cframe-[brz]9-*.fits -o $OUTDIR/70003/coadd-9-00051059.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051060/cframe-[brz]0-*.fits -o $OUTDIR/70003/coadd-0-00051060.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051060/cframe-[brz]3-*.fits -o $OUTDIR/70003/coadd-3-00051060.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051060/cframe-[brz]6-*.fits -o $OUTDIR/70003/coadd-6-00051060.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051060/cframe-[brz]7-*.fits -o $OUTDIR/70003/coadd-7-00051060.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051060/cframe-[brz]9-*.fits -o $OUTDIR/70003/coadd-9-00051060.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051072/cframe-[brz]0-*.fits -o $OUTDIR/70003/coadd-0-00051072.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051072/cframe-[brz]3-*.fits -o $OUTDIR/70003/coadd-3-00051072.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051072/cframe-[brz]6-*.fits -o $OUTDIR/70003/coadd-6-00051072.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051072/cframe-[brz]7-*.fits -o $OUTDIR/70003/coadd-7-00051072.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051072/cframe-[brz]9-*.fits -o $OUTDIR/70003/coadd-9-00051072.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051073/cframe-[brz]0-*.fits -o $OUTDIR/70003/coadd-0-00051073.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051073/cframe-[brz]3-*.fits -o $OUTDIR/70003/coadd-3-00051073.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051073/cframe-[brz]6-*.fits -o $OUTDIR/70003/coadd-6-00051073.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051073/cframe-[brz]7-*.fits -o $OUTDIR/70003/coadd-7-00051073.fits
time desi_coadd_spectra --coadd-cameras -i $REDUXDIR/exposures/20200219/00051073/cframe-[brz]9-*.fits -o $OUTDIR/70003/coadd-9-00051073.fits

################# redrock commands: #################

# salloc -N 1 -C haswell -q interactive -t 04:00:00

srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-0-00051053.h5 -z $OUTDIR/70003/zbest-0-00051053.fits $OUTDIR/70003/coadd-0-00051053.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-3-00051053.h5 -z $OUTDIR/70003/zbest-3-00051053.fits $OUTDIR/70003/coadd-3-00051053.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-6-00051053.h5 -z $OUTDIR/70003/zbest-6-00051053.fits $OUTDIR/70003/coadd-6-00051053.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-7-00051053.h5 -z $OUTDIR/70003/zbest-7-00051053.fits $OUTDIR/70003/coadd-7-00051053.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-9-00051053.h5 -z $OUTDIR/70003/zbest-9-00051053.fits $OUTDIR/70003/coadd-9-00051053.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-0-00051054.h5 -z $OUTDIR/70003/zbest-0-00051054.fits $OUTDIR/70003/coadd-0-00051054.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-3-00051054.h5 -z $OUTDIR/70003/zbest-3-00051054.fits $OUTDIR/70003/coadd-3-00051054.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-6-00051054.h5 -z $OUTDIR/70003/zbest-6-00051054.fits $OUTDIR/70003/coadd-6-00051054.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-7-00051054.h5 -z $OUTDIR/70003/zbest-7-00051054.fits $OUTDIR/70003/coadd-7-00051054.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-9-00051054.h5 -z $OUTDIR/70003/zbest-9-00051054.fits $OUTDIR/70003/coadd-9-00051054.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-0-00051055.h5 -z $OUTDIR/70003/zbest-0-00051055.fits $OUTDIR/70003/coadd-0-00051055.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-3-00051055.h5 -z $OUTDIR/70003/zbest-3-00051055.fits $OUTDIR/70003/coadd-3-00051055.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-6-00051055.h5 -z $OUTDIR/70003/zbest-6-00051055.fits $OUTDIR/70003/coadd-6-00051055.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-7-00051055.h5 -z $OUTDIR/70003/zbest-7-00051055.fits $OUTDIR/70003/coadd-7-00051055.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-9-00051055.h5 -z $OUTDIR/70003/zbest-9-00051055.fits $OUTDIR/70003/coadd-9-00051055.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-0-00051056.h5 -z $OUTDIR/70003/zbest-0-00051056.fits $OUTDIR/70003/coadd-0-00051056.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-3-00051056.h5 -z $OUTDIR/70003/zbest-3-00051056.fits $OUTDIR/70003/coadd-3-00051056.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-6-00051056.h5 -z $OUTDIR/70003/zbest-6-00051056.fits $OUTDIR/70003/coadd-6-00051056.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-7-00051056.h5 -z $OUTDIR/70003/zbest-7-00051056.fits $OUTDIR/70003/coadd-7-00051056.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-9-00051056.h5 -z $OUTDIR/70003/zbest-9-00051056.fits $OUTDIR/70003/coadd-9-00051056.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-0-00051057.h5 -z $OUTDIR/70003/zbest-0-00051057.fits $OUTDIR/70003/coadd-0-00051057.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-3-00051057.h5 -z $OUTDIR/70003/zbest-3-00051057.fits $OUTDIR/70003/coadd-3-00051057.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-6-00051057.h5 -z $OUTDIR/70003/zbest-6-00051057.fits $OUTDIR/70003/coadd-6-00051057.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-7-00051057.h5 -z $OUTDIR/70003/zbest-7-00051057.fits $OUTDIR/70003/coadd-7-00051057.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-9-00051057.h5 -z $OUTDIR/70003/zbest-9-00051057.fits $OUTDIR/70003/coadd-9-00051057.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-0-00051058.h5 -z $OUTDIR/70003/zbest-0-00051058.fits $OUTDIR/70003/coadd-0-00051058.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-3-00051058.h5 -z $OUTDIR/70003/zbest-3-00051058.fits $OUTDIR/70003/coadd-3-00051058.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-6-00051058.h5 -z $OUTDIR/70003/zbest-6-00051058.fits $OUTDIR/70003/coadd-6-00051058.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-7-00051058.h5 -z $OUTDIR/70003/zbest-7-00051058.fits $OUTDIR/70003/coadd-7-00051058.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-9-00051058.h5 -z $OUTDIR/70003/zbest-9-00051058.fits $OUTDIR/70003/coadd-9-00051058.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-0-00051059.h5 -z $OUTDIR/70003/zbest-0-00051059.fits $OUTDIR/70003/coadd-0-00051059.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-3-00051059.h5 -z $OUTDIR/70003/zbest-3-00051059.fits $OUTDIR/70003/coadd-3-00051059.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-6-00051059.h5 -z $OUTDIR/70003/zbest-6-00051059.fits $OUTDIR/70003/coadd-6-00051059.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-7-00051059.h5 -z $OUTDIR/70003/zbest-7-00051059.fits $OUTDIR/70003/coadd-7-00051059.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-9-00051059.h5 -z $OUTDIR/70003/zbest-9-00051059.fits $OUTDIR/70003/coadd-9-00051059.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-0-00051060.h5 -z $OUTDIR/70003/zbest-0-00051060.fits $OUTDIR/70003/coadd-0-00051060.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-3-00051060.h5 -z $OUTDIR/70003/zbest-3-00051060.fits $OUTDIR/70003/coadd-3-00051060.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-6-00051060.h5 -z $OUTDIR/70003/zbest-6-00051060.fits $OUTDIR/70003/coadd-6-00051060.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-7-00051060.h5 -z $OUTDIR/70003/zbest-7-00051060.fits $OUTDIR/70003/coadd-7-00051060.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-9-00051060.h5 -z $OUTDIR/70003/zbest-9-00051060.fits $OUTDIR/70003/coadd-9-00051060.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-0-00051072.h5 -z $OUTDIR/70003/zbest-0-00051072.fits $OUTDIR/70003/coadd-0-00051072.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-3-00051072.h5 -z $OUTDIR/70003/zbest-3-00051072.fits $OUTDIR/70003/coadd-3-00051072.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-6-00051072.h5 -z $OUTDIR/70003/zbest-6-00051072.fits $OUTDIR/70003/coadd-6-00051072.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-7-00051072.h5 -z $OUTDIR/70003/zbest-7-00051072.fits $OUTDIR/70003/coadd-7-00051072.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-9-00051072.h5 -z $OUTDIR/70003/zbest-9-00051072.fits $OUTDIR/70003/coadd-9-00051072.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-0-00051073.h5 -z $OUTDIR/70003/zbest-0-00051073.fits $OUTDIR/70003/coadd-0-00051073.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-3-00051073.h5 -z $OUTDIR/70003/zbest-3-00051073.fits $OUTDIR/70003/coadd-3-00051073.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-6-00051073.h5 -z $OUTDIR/70003/zbest-6-00051073.fits $OUTDIR/70003/coadd-6-00051073.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-7-00051073.h5 -z $OUTDIR/70003/zbest-7-00051073.fits $OUTDIR/70003/coadd-7-00051073.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi -o $OUTDIR/70003/redrock-9-00051073.h5 -z $OUTDIR/70003/zbest-9-00051073.fits $OUTDIR/70003/coadd-9-00051073.fits
