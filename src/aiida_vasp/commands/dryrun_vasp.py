"""
Module for dry-running a VASP calculation
"""

import shutil
import subprocess as sb
import tempfile
import time
from pathlib import Path

import click
import yaml
from parsevasp.kpoints import Kpoints

# pylint:disable=too-many-branches,consider-using-with


@click.command('dryrun-vasp')
@click.option(
    '--input-dir',
    help='Where the VASP input is, default to the current working directory.',
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default='.',
    show_default=True,
)
@click.option('--vasp-exe', help='Executable for VASP', default='vasp_std', show_default=True)
@click.option(
    '--timeout',
    help='Timeout in seconds to terminate VASP',
    default=10,
    show_default=True,
)
@click.option('--work-dir', help='Working directory for running', show_default=True)
@click.option('--keep', help='Wether to the dryrun files', is_flag=True, show_default=True)
@click.option(
    '--force',
    help='Force the run even if the working directory exists.',
    is_flag=True,
    show_default=True,
)
def cmd_dryrun_vasp(input_dir, vasp_exe, timeout, work_dir, keep, force):
    """
    A simple tool to dryrun a VASP calculation. The calculation will be run for
    up to <timeout> seconds. The underlying VASP process will be terminated once it enters
    the main loop, which is signalled by the appearance of a `INWAV` keyword in the OUTCAR.
    """
    result = dryrun_vasp(
        input_dir=input_dir,
        vasp_exe=vasp_exe,
        timeout=timeout,
        work_dir=work_dir,
        keep=keep,
        force=force,
    )
    with open(Path(input_dir) / 'dryrun.yaml', 'w', encoding='utf-8') as fhandle:
        yaml.dump(result, fhandle, Dumper=yaml.SafeDumper)


def dryrun_vasp(input_dir, vasp_exe='vasp_std', timeout=10, work_dir=None, keep=False, force=False):
    """
    Perform a "dryrun" for a VASP calculation - get the number of kpoints, bands and
    estimated memory usage.
    """
    input_dir = Path(input_dir)
    if not work_dir:
        tmpdir = tempfile.mkdtemp()  # tmpdir is the one to remove when finished
        work_dir = Path(tmpdir) / 'vasp_dryrun'
    else:
        work_dir = Path(work_dir)
        if work_dir.resolve() == input_dir.resolve():
            raise ValueError('The working directory cannot be the input directory!')
        if work_dir.exists():
            if not force:
                raise FileExistsError(f'Working directory {work_dir} exists already! Please remove it first.')
            shutil.rmtree(work_dir)
        tmpdir = str(work_dir)
    shutil.copytree(str(input_dir), str(work_dir))

    # Add the DRYRUNCAR for triggering the dryrun interface
    (Path(work_dir) / 'DRYRUNCAR').write_text('LDRYRUN = .TRUE.\n')

    process = sb.Popen(vasp_exe, cwd=str(work_dir))
    launch_start = time.time()
    outcar = work_dir / 'OUTCAR'
    time.sleep(3.0)  # Sleep for 3 seconds to wait for VASP creating the file
    dryrun_finish = False
    try:
        while (time.time() - launch_start < timeout) and not dryrun_finish:
            with open(outcar, encoding='utf-8') as fhandle:
                for line in fhandle:
                    if 'INWAV' in line or 'Terminating' in line:
                        dryrun_finish = True
                        break
            # Stop if VASP is terminated or crashed
            if process.poll() is not None:
                break
            time.sleep(0.2)
    except Exception as error:
        raise error
    finally:
        # Once we are out side the loop, kill VASP process
        process.kill()
    result = parse_outcar(outcar)
    ibzkpt = parse_ibzkpt(work_dir / 'IBZKPT')
    result['kpoints_and_weights_ibzkpt'] = ibzkpt

    if not keep:
        shutil.rmtree(tmpdir)

    return result


def parse_ibzkpt(ibzkpt_path):
    """
    Parsing the IBZKPT file
    """

    kpoints = Kpoints(file_path=str(ibzkpt_path))
    tmp = kpoints.get_dict()['points']
    kpoints_and_weights = [elem[0].tolist() + [elem[1]] for elem in tmp]
    total_weight = sum(tmp[3] for tmp in kpoints_and_weights)

    # Normalise the kpoint weights
    normalised = []
    for entry in kpoints_and_weights:
        normalised.append(entry[:3] + [entry[3] / total_weight])

    return normalised


def parse_outcar(outcar_path):
    """
    Parse the header part of the OUTCAR

    Returns:
        A dictionary of the parsed information
    """
    output_dict = {
        'POTCARS': [],
    }
    with open(outcar_path, encoding='utf-8') as fhandle:
        lines = fhandle.readlines()
    for line_number, line in enumerate(lines):
        if 'POTCAR:' in line:
            content = line.split(maxsplit=1)[1].strip()
            if content not in output_dict['POTCARS']:
                output_dict['POTCARS'].append(content)
        elif 'NKPTS' in line:
            tokens = line.strip().split()
            output_dict['num_kpoints'] = int(tokens[tokens.index('NKPTS') + 2])
            output_dict['num_bands'] = int(tokens[-1])
        elif 'dimension x,y,z NGX =' in line:
            tokens = line.strip().split()
            output_dict['NGX'] = int(tokens[tokens.index('NGX') + 2])
            output_dict['NGY'] = int(tokens[tokens.index('NGY') + 2])
            output_dict['NGZ'] = int(tokens[tokens.index('NGZ') + 2])
        elif 'FFT grid for exact exchange' in line:
            tokens = lines[line_number + 1].replace(';', '').strip().split()
            output_dict['EX NGX'] = int(tokens[tokens.index('NGX') + 2])
            output_dict['EX NGY'] = int(tokens[tokens.index('NGY') + 2])
            output_dict['EX NGZ'] = int(tokens[tokens.index('NGZ') + 2])
        elif 'NPLWV' in line:
            try:
                output_dict['num_plane_waves'] = int(line.split()[-1])
            except ValueError:
                pass
        elif 'k-points in reciprocal lattice and weights:' in line:
            kblock = lines[line_number + 1 : line_number + 1 + output_dict['num_kpoints']]
            k_list = [[float(token) for token in subline.strip().split()] for subline in kblock]
            output_dict['kpoints_and_weights'] = k_list
        elif 'maximum and minimum number of plane-waves per node :' in line:
            output_dict['plane_waves_min_max'] = [float(token) for token in line.split()[-2:]]
        elif 'total amount of memory used by VASP MPI-rank0' in line:
            output_dict['max_ram_rank0'] = float(line.split()[-2])
            for subline in lines[line_number + 3 : line_number + 9]:
                tokens = subline.replace(':', '').split()
                output_dict['mem_' + tokens[0]] = float(tokens[-2])
    return output_dict
