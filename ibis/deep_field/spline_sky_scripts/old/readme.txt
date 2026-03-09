# launch-sky.sh is not used
cd /pscratch/sd/r/rongpu/ibis-dr1/splinesky
salloc -N 2 -C cpu -q interactive -t 04:00:00
srun --wait=0 --ntasks-per-node 1 payload-sky.sh runcalibs-sky-commands.txt ; exit
