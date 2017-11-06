#!/bin/bash
#QSUB -q gr10260f
#QSUB -W 24:00
#QSUB -A p=20:t=1:c=1:m=3072M
#QSUB -rn
#QSUB -J mp-30311-born_effective_charge
#QSUB -e err.log
#QSUB -o std.log

mpirun ~/vasp544mpi
sleep 60