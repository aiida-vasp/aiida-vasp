#!/bin/bash
exec > _scheduler-stdout.txt
exec 2> _scheduler-stderr.txt
export MOCK_VASP_UPLOAD_PREFIX=mock_silicon_hybrid_no_relax


'/home/bonan/miniconda3/envs/aiida-2.4-dev/bin/mock-vasp'  > 'vasp_output' 2>&1
