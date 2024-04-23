#!/bin/bash
exec > _scheduler-stdout.txt
exec 2> _scheduler-stderr.txt


export AIIDA_PATH=/tmp/tmp9gxj6gz9/.aiida

'mpirun' '-np' '3' '/home/bonan/appdir/vtst-code/vasp.6.2.0/bin/vasp_std'  > 'stdout' 2>&1
