#!/bin/bash
exec > _scheduler-stdout.txt
exec 2> _scheduler-stderr.txt


export AIIDA_PATH=/tmp/tmp4q2o_4kq/.aiida

'/home/bonan/appdir/vtst-code/vasp.6.2.0/bin/vasp_std'  > 'vasp_output' 2>&1
