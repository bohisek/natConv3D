#!/bin/bash
# sets home directory
DATADIR="/storage/plzen1/home/$LOGNAME/palabos-v1.5r1/examples/own/2018.03.01_HE"
cd $DATADIR

#cat $PBS_NODEFILE |uniq > nodes.txt
#export OMP_NUM_THREADS=$PBS_NUM_PPN
# --hostfile nodes.txt

module add openmpi-2.1.1-intel

mpirun -np 100 ./natConv3D 3.0151e+9 > results.out
