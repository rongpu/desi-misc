################# redrock commands: #################

srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80607/coadd-00068028.fits -o $OUTDIR/80607/redrock-00068028.fits -z $OUTDIR/80607/zbest-00068028.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80607/coadd-00068662.fits -o $OUTDIR/80607/redrock-00068662.fits -z $OUTDIR/80607/zbest-00068662.fits

srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80609/coadd-00067781.fits -o $OUTDIR/80609/redrock-00067781.fits -z $OUTDIR/80609/zbest-00067781.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80609/coadd-00067783.fits -o $OUTDIR/80609/redrock-00067783.fits -z $OUTDIR/80609/zbest-00067783.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80609/coadd-00068063.fits -o $OUTDIR/80609/redrock-00068063.fits -z $OUTDIR/80609/zbest-00068063.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80609/coadd-00068064.fits -o $OUTDIR/80609/redrock-00068064.fits -z $OUTDIR/80609/zbest-00068064.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80609/coadd-00068065.fits -o $OUTDIR/80609/redrock-00068065.fits -z $OUTDIR/80609/zbest-00068065.fits
srun -N 1 -n 32 -c 2 rrdesi_mpi $OUTDIR/80609/coadd-00068334.fits -o $OUTDIR/80609/redrock-00068334.fits -z $OUTDIR/80609/zbest-00068334.fits
