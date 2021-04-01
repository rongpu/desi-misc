################# redrock commands: #################

srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80607/coadd-00068663.fits -o $OUTDIR/80607/redrock-00068663.fits -z $OUTDIR/80607/zbest-00068663.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80607/coadd-00068664.fits -o $OUTDIR/80607/redrock-00068664.fits -z $OUTDIR/80607/zbest-00068664.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80607/coadd-00068666.fits -o $OUTDIR/80607/redrock-00068666.fits -z $OUTDIR/80607/zbest-00068666.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80607/coadd-00068844.fits -o $OUTDIR/80607/redrock-00068844.fits -z $OUTDIR/80607/zbest-00068844.fits

srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80609/coadd-00068336.fits -o $OUTDIR/80609/redrock-00068336.fits -z $OUTDIR/80609/zbest-00068336.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80609/coadd-00068337.fits -o $OUTDIR/80609/redrock-00068337.fits -z $OUTDIR/80609/zbest-00068337.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80609/coadd-00068338.fits -o $OUTDIR/80609/redrock-00068338.fits -z $OUTDIR/80609/zbest-00068338.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80609/coadd-00068339.fits -o $OUTDIR/80609/redrock-00068339.fits -z $OUTDIR/80609/zbest-00068339.fits
