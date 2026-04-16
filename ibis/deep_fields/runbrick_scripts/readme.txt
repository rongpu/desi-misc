# To run on interactive nodes
salloc -N 4 -C cpu -q interactive -t 04:00:00
srun --account=m3592 --wait=0 --ntasks-per-node 1 payload-runbrick.sh runbrick-commands-remain.txt
