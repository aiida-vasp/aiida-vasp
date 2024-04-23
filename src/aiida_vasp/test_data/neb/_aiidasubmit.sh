#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -m n
#$ -N aiida-49859
#$ -V
#$ -o _scheduler-stdout.txt
#$ -e _scheduler-stderr.txt
#$ -pe mpi 2
#$ -l h_rt=00:30:00

'/home/bonan/appdir/vasp.5.4.4/bin/vasp_std'  > 'stdout' 2>&1
