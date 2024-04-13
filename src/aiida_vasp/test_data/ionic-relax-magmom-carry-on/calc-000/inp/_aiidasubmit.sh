#!/bin/bash
exec > _scheduler-stdout.txt
exec 2> _scheduler-stderr.txt


'/home/bonan/appdir/vtst-code/vasp.6.2.0/bin/vasp_std'  > 'vasp_output' 2>&1
