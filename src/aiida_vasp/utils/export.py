import gzip
import os
import re
import shutil
from pathlib import Path

from aiida import orm
from aiida.repository import FileType

from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.parsers.content_parsers.potcar import MultiPotcarIo

from .aiida_utils import ensure_node_first_arg, ensure_node_kwargs


def export_vasp(process, folder, decompress=True, include_potcar=True):
    """
    Export a VASP calculation, works for both `VaspCalculation` or `VaspWorkChain`
    """

    # Dispatch export function based on process type
    if process.process_type.endswith('vasp'):
        export_vasp_calc(process, folder, decompress=decompress, include_potcar=include_potcar)
    elif process.process_type.endswith('neb'):
        export_neb(process, folder, decompress=decompress, include_potcar=include_potcar)
    elif process.process_type.endswith('relax'):
        export_relax(process, folder, decompress=decompress, include_potcar=include_potcar)
    else:
        raise TypeError(f'Unsupported process type: {process.process_type}')


@ensure_node_kwargs
@ensure_node_first_arg
def export_vasp_calc(node, folder, decompress=False, include_potcar=True):
    """
    Export a AiiDA VASP calculation

    :param node: A VaspCalculation node or VaspWorkChain node
    """
    from aiida.common.links import LinkType
    from aiida.orm import CalcJobNode, WorkChainNode

    folder = Path(folder)
    folder.mkdir(exist_ok=True)

    # Inputs
    retrieved = node.get_outgoing(link_label_filter='retrieved').one().node
    if isinstance(node, CalcJobNode):
        calcjob = node
    elif isinstance(node, WorkChainNode):
        # In this case the node is an workchain we export the
        # 'retrieved' output link and trace to its ancestor
        calcjob = retrieved.get_incoming(link_label_filter='retrieved', link_type=LinkType.CREATE).one().node
    else:
        raise RuntimeError(f'The node {node} is not a valid calculation')
    info_file = folder / ('aiida_info')
    info_content = f'Label: {calcjob.label}\nDescription: {calcjob.description}\nUUID: {calcjob.uuid}\n'
    info_file.write_text(info_content)
    # export the retrieved outputs  and the input files
    save_all_repository_objects(retrieved, folder, decompress)
    save_all_repository_objects(calcjob, folder, decompress)
    if include_potcar:
        export_pseudos(calcjob, folder)


@ensure_node_first_arg
@ensure_node_kwargs
def export_pseudos(calc_job_node, folder):
    """Save the pseudopotential file (POTCAR)"""
    pps = calc_job_node.get_incoming(link_label_filter='potential%').nested()['potential']
    multi_potcar = MultiPotcarIo.from_structure(calc_job_node.inputs.structure, pps)
    dst = str(folder / 'POTCAR')
    multi_potcar.write(dst)


@ensure_node_first_arg
@ensure_node_kwargs
def export_relax(workchain_node, dst, include_potcar=False, decompress=False):
    """
    Export a relaxation workflow (e.g. VaspRelaxWorkChain)

    This function exports a series of relaxation calculations in sub-folders
    """
    from aiida.orm import Node, QueryBuilder, WorkChainNode

    from aiida_vasp.workchains.relax import RelaxWorkChain
    from aiida_vasp.workchains.v2.relax import VaspRelaxWorkChain

    dst = Path(dst)
    dst.mkdir(exist_ok=True)
    if workchain_node.process_class not in (VaspRelaxWorkChain, RelaxWorkChain):
        raise ValueError(
            f'Error {workchain_node} should be `VaspRelaxWorkChain` or `RelaxWorkChain`,'
            f' but it is {workchain_node.process_class}'
        )

    q = QueryBuilder()
    q.append(Node, filters={'id': workchain_node.pk})
    q.append(WorkChainNode, tag='vaspwork', project=['id', '*'])
    q.order_by({'vaspwork': {'id': 'asc'}})  # Sort by ascending PK
    for index, (pk, node) in enumerate(q.iterall()):
        relax_folder = dst / f'relax_calc_{index:03d}'
        try:
            export_vasp_calc(node, relax_folder, decompress=decompress, include_potcar=include_potcar)
        except (ValueError, AttributeError, KeyError):
            print(f'Error exporting calculation {pk}')

    # Write POSCAR file for the input
    input_structure = workchain_node.inputs.structure
    poscar_parser = PoscarParser(data=input_structure, precision=10)
    poscar_parser.write(str(dst / 'POSCAR'))

    # Write POSCAR file for the input
    try:
        out_structure = workchain_node.outputs.relax.structure
    except AttributeError:
        print(
            'Cannot find the output structure - skipping.'
            ' This usually means that the relaxation did not finish without error.'
        )
        out_structure = None
    if out_structure:
        poscar_parser = PoscarParser(data=out_structure, precision=10)
        poscar_parser.write(str(dst / 'POSCAR_RELAXED'))

    # Write the info
    info_file = dst / ('aiida_info')
    info_content = (
        f'Label: {workchain_node.label}\nDescription: {workchain_node.description}\nUUID: {workchain_node.uuid}\n'
    )
    info_file.write_text(info_content)


@ensure_node_first_arg
def export_neb(workchain, dst, decompress=True, include_potcar=True, energy_type='energy_extrapolated'):
    """Export the neb calculation"""
    from aiida.plugins import WorkflowFactory

    energies = {key: value[energy_type] for key, value in workchain.outputs.misc['total_energies'].items()}

    # Query for the energy computed for the end structures
    q = orm.QueryBuilder()
    q.append(orm.Node, filters={'id': workchain.inputs.initial_structure.id}, tag='root')
    q.append(orm.CalcFunctionNode, with_outgoing='root', project=['attributes.function_name'])
    q.append(
        orm.StructureData,
        with_outgoing=orm.CalcFunctionNode,
        tag='relaxed',
        project=['label'],
        # edge_filters={'label': 'init_structure'},
        edge_project=['label'],
    )
    q.append(
        WorkflowFactory('vasp.v2.relax'),
        with_outgoing='relaxed',
        project=['label', 'uuid'],
        tag='relaxation',
    )
    q.append(
        orm.Dict,
        with_incoming='relaxation',
        edge_filters={'label': 'misc'},
        project=['attributes.total_energies.energy_extrapolated'],
    )
    q.append(orm.CalcJobNode, with_outgoing=orm.Dict, project=['*'])
    q.distinct()

    # First export the original calculation
    export_vasp_calc(workchain, dst, decompress=decompress, include_potcar=include_potcar)
    ends = {}
    end_id = f'{len(energies) + 1:02d}'
    for _, _, _, relax_uuid, eng, calcjob, label in q.all():
        if label.startswith('init'):
            if '00' in ends:
                print(
                    'Duplicated calculation: {relax_uuid} -> {eng} vs existing {existing}'.format(
                        relax_uuid=relax_uuid, eng=eng, existing=ends['00']
                    )
                )
            else:
                ends['00'] = calcjob

        elif label.startswith('final'):
            if end_id in ends:
                print(
                    'Duplicated calculation: {relax_uuid} -> {eng} vs existing {existing}'.format(
                        relax_uuid=relax_uuid, eng=eng, existing=ends[end_id]
                    )
                )
            else:
                ends[end_id] = calcjob
    # Export the end point calculation
    for key, value in ends.items():
        export_vasp_calc(value, Path(dst) / key, decompress=decompress, include_potcar=include_potcar)


@ensure_node_kwargs
def copy_from_aiida(name: str, node, dst: Path, decompress=False, exclude=None):
    """
    Copy objects from aiida repository.

    Args:
        name (str): The full name (including the parent path) of the object.
        node (orm.Node): Node object for which the files in the repo to be copied.
        dst (Path): Path of the destination folder.

    This is a recursive function so directory copying also works.
    """

    # For check the regex the first thing because this function will be called recursively
    if exclude and re.match(exclude, name):
        return

    obj = node.base.repository.get_object(name)

    # If it is a directory, copy the contents one by one
    if obj.file_type == FileType.DIRECTORY:
        for sub_obj in node.base.repository.list_objects(name):
            copy_from_aiida(os.path.join(name, sub_obj.name), node, dst, exclude=exclude)
    else:
        # It is a file
        with node.open(name, mode='rb') as fsource:
            # Make parent directory if needed
            frepo_path = dst / name
            Path(frepo_path.parent).mkdir(exist_ok=True, parents=True)
            # Write the file
            if name.endswith('.gz') and decompress:
                out_path = str(frepo_path)[:-3]
                out_decompress = True
            else:
                out_decompress = False
                out_path = str(frepo_path)

            if not out_decompress:
                with open(out_path, 'wb') as fdst:
                    shutil.copyfileobj(fsource, fdst)
            else:
                gobj = gzip.GzipFile(fileobj=fsource, mode='rb')
                with open(out_path, 'wb') as fdst:
                    shutil.copyfileobj(gobj, fdst)


@ensure_node_first_arg
def save_all_repository_objects(node: orm.Node, target_path: Path, decompress=False, exclude=None):
    """Copy all objects of a node saved in the repository to the disc"""
    for name in node.list_object_names():
        copy_from_aiida(name, node, target_path, decompress, exclude=exclude)
