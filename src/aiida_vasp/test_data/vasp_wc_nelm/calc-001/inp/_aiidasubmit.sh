#!/bin/bash
exec > _scheduler-stdout.txt
exec 2> _scheduler-stderr.txt


source /opt/intel/bin/compilervars.sh intel64

'mpirun' '-np' '2' '/home/bonan/appdir/VASP/bin/vasp_std_mkl'  > 'vasp_output' 2>&1
