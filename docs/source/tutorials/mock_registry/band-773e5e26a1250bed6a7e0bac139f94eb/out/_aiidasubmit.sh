#!/bin/bash
exec > _scheduler-stdout.txt
exec 2> _scheduler-stderr.txt


'mpirun' '-np' '1' '/home/bonan/miniconda3/envs/aiida-2.4-dev/bin/mock-vasp'  > 'vasp_output' 2>&1
