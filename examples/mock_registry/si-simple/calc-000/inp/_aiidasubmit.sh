#!/bin/bash
exec > _scheduler-stdout.txt
exec 2> _scheduler-stderr.txt


export VASP_MOCK_CODE_BASE=/home/bonan/aiida_envs/aiida-2.0/aiida-vasp/examples/mock_registry

'mpirun' '-np' '1' '/home/bonan/appdir/vtst-code/vasp.6.2.0/bin/vasp_std'  > 'vasp_output' 2>&1
