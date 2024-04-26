"""
The ``OUTCAR`` parser interface.

----------------------------
Contains the parsing interfaces to parsevasp used to parse ``OUTCAR`` content.
"""

# pylint: disable=abstract-method
from pathlib import Path

import numpy as np
from parsevasp.outcar import Outcar

from aiida_vasp.parsers.content_parsers.base import BaseFileParser


class OutcarParser(BaseFileParser):
    """The parser interface that enables parsing of ``OUTCAR`` content.

    The parser is triggered by using the ``elastic_moduli``, ``magnetization`` or ``site-magnetization``
    ``run_stats`` or ``run_status`` quantity keys.

    """

    DEFAULT_SETTINGS = {'quantities_to_parse': ['run_status', 'run_stats']}

    PARSABLE_QUANTITIES = {
        'elastic_moduli': {'inputs': [], 'name': 'elastic_moduli', 'prerequisites': []},
        'symmetries': {'inputs': [], 'name': 'symmetries', 'prerequisites': []},
        'magnetization': {'inputs': [], 'name': 'magnetization', 'prerequisites': []},
        'site_magnetization': {'inputs': [], 'name': 'site_magnetization', 'prerequisites': []},
        'run_stats': {
            'inputs': [],
            'name': 'run_stats',
            'prerequisites': [],
        },
        'run_status': {
            'inputs': [],
            'name': 'run_status',
            'prerequisites': [],
        },
    }

    def _init_from_handler(self, handler):
        """Initialize a ``parsevasp`` object of ``Outcar`` using a file like handler.

        Parameters
        ----------
        handler : object
            A file like object that provides the necessary ``OUTCAR`` content to be parsed.

        """

        try:
            self._content_parser = Outcar(file_handler=handler, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    @property
    def run_status(self):
        """
        Fetch status of calculations.

        Returns
        -------
        status : dict
            A dictionary containing the keys ``finished``, which is True if the VASP calculation
            contain timing information in the end of the ``OUTCAR``. The key ``ionic_converged`` is
            True if the number of ionic steps detected is smaller than the supplied NSW.
            The key ``electronic_converged`` is True if the number of electronic steps is smaller than
            NELM (defaults to 60 in VASP). It is also possible to check if all the ionic steps
            did reached NELM and thus did not converged if the key ``consistent_nelm_breach`` is ``True``,
            while ``contains_nelm_breach`` is True if one or more ionic steps reached NELM and thus
            did not converge electronically.

        """
        status = self._content_parser.get_run_status()
        return status

    @property
    def run_stats(self):
        """Fetch the run statistics, which included timings and memory consumption.

        Returns
        -------
        stats : dict
            A dictionary containing timing and memory consumption information
            that are parsed from the end of the ``OUTCAR`` file. The key names are
            mostly preserved, except for the memory which is prefixed with ``mem_usage_``.
            Units are preserved from ``OUTCAR`` and there are some differences between
            VASP 5 and 6.

        """
        stats = self._content_parser.get_run_stats()
        return stats

    @property
    def symmetries(self):
        """Fetch some basic symmetry data.

        Returns
        -------
        sym : dict
            A dictionary containing the number of space group operations in the
            key ``num_space_group_operations`` and the detected supplied cell in
            ``original_cell_type``. In ``symmetrized_cell_type`` the cell on which
            VASP performs the calculation has been included. Each value in the
            dictionary is a list, where each entry represent one ionic step.

        """

        sym = self._content_parser.get_symmetry()
        return sym

    @property
    def elastic_moduli(self):
        """Fetch the elastic moduli tensor.

        Returns
        -------
        moduli : dict
            A dictionary containing ndarrays with the rigid ion elastic moduli, both symmetrized and
            non-symmetrized for the keys ``symmetrized`` and ``non_symmetrized`` respectively.
            The key ``total`` contain both the rigid ion and the ionic contributions to the
            elastic tensor for the symmetrized case.

        """

        moduli = self._content_parser.get_elastic_moduli()
        return moduli

    @property
    def site_magnetization(self):
        """Fetch the site dependent magnetization.

        Returns
        -------
        magnetization : dict
            A dictionary containing the key ``sphere`` which contains the integrated
            magnetization in units of Bohr magneton. Additional keys under ``sphere`` are
            given for each direction and for non-collinear calculations all of them are used.
            The ``site_moment`` yields the magnetization per site, with a key describing the
            site number and then the ``s``, ``p``, ``d`` etc. the projections of the site magnetization
            and ``tot`` containing the total magnetization for that site.
            The ``total_magnetization`` gives the sum of each magnetization projection and
            magnetization total for each site.
            The ``full_cell`` key yields the magnetization from the electronic part of the last
            electronic step in a list.

        """
        magnetization = self._content_parser.get_magnetization()
        return magnetization

    @property
    def magnetization(self):
        """Fetch the full cell magnetization.

        Returns
        -------
        magnetization : list
            A list containing an entry that is the total magnetization in the cell in unit of
            Bohr magneton. The magnetization returned is the one associated with the electrons for the
            last electronic step.

        """
        magnetization = self.site_magnetization
        if magnetization is not None:
            magnetization = magnetization['full_cell']
        return magnetization


class VtstNebOutcarParser(OutcarParser):
    """
    Parser for processing OUTCAR generated by VASP with VTST
    """

    DEFAULT_SETTINGS = {'quantities_to_parse': ['run_status', 'run_stats', 'neb_data']}
    PARSABLE_QUANTITIES = {
        'neb_data': {'inputs': [], 'name': 'neb_data', 'prerequisites': []},
        'outcar-forces': {'inputs': [], 'name': 'forces', 'prerequisites': []},
        'outcar-positions': {'inputs': [], 'name': 'positions', 'prerequisites': []},
        'outcar-cell': {'inputs': [], 'name': 'cell', 'prerequisites': []},
        **OutcarParser.PARSABLE_QUANTITIES,
    }

    def __init__(self, *args, **kwargs):
        """
        Instantiate the parser
        """
        self._parsed_neb_data = {}
        super().__init__(*args, **kwargs)

    def _init_from_handler(self, handler):
        """Initial fro the handler"""
        super()._init_from_handler(handler)
        # Parse the NEB results from the handle and store in a dictionary
        self._parsed_neb_data = _parse_neb_outputs(handler)

    @property
    def neb_data(self):
        """
        Parsed NEB results
        """
        return self._parsed_neb_data.get('neb_data')

    @property
    def forces(self):
        """Parsed forces"""
        return self._parsed_neb_data.get('outcar-forces')

    @property
    def positions(self):
        """Parsed forces"""
        return self._parsed_neb_data.get('outcar-positions')

    @property
    def cell(self):
        """Parsed cell vectors"""
        return self._parsed_neb_data.get('outcar-cell')


def _parse_force_block(lines):
    """
    Parse the block of total forces from the OUTCAR file

    :param lines: A list of lines containing lines including the TOTAL-FORCE block

    :returns: A tuple of position and forces
    """
    forces = []
    positions = []
    istart = len(lines)
    for idx, line in enumerate(lines):
        if 'TOTAL-FORCE (eV/Angst)' in line:
            istart = idx
        elif idx > istart + 1:
            if not line.startswith(' -----'):  # Still in the block
                values = list(map(float, line.split()))
                positions.append(values[:3])
                forces.append(values[3:])
            else:
                # Reached the end of the block
                break
    return positions, forces


def _parse_neb_outputs(path, inputs=None):  # pylint: disable=too-many-branches,too-many-statements
    """
    Scan for NEB output in the OUTCAR content

    :param path: Input path or fileobj
    :param inputs(dict): Dictionary where the parsed data should be placed

    :returns dict: A dictionary of the parsed data
    """
    inputs = {} if inputs is None else inputs

    if isinstance(path, (str, Path)):
        with open(path, 'r', encoding='utf8') as fobj:
            lines = fobj.readlines()
    # A file-like object
    elif hasattr(path, 'readlines'):
        # Reset seek
        path.seek(0)
        lines = path.readlines()
    else:
        raise ValueError(f"'path' variable is not supported: {path}")

    vtst_data = {'neb_converged': False}
    for idx, line in enumerate(lines):
        if 'NIONS' in line:
            nions = int(line.split()[-1])

        elif 'VTST: version' in line:
            vtst_data['version'] = line.split(':')[1].strip()

        elif 'NEB: Tangent' in line:
            tangents = []
            for isub in range(idx + 2, idx + 99999):
                subline = lines[isub]
                if subline.strip():
                    tangents.append([float(tmp) for tmp in subline.split()])
                else:
                    break
            vtst_data['tangents'] = tangents
        elif 'NEB: forces' in line:
            forces = [float(tmp) for tmp in line.split()[-3:]]
            vtst_data['force_par_spring'] = forces[0]
            vtst_data['force_prep_real'] = forces[1]
            vtst_data['force_dneb'] = forces[2]
        elif 'stress matrix after NEB project' in line:
            stress = []
            for isub in range(idx + 1, idx + 4):
                stress.append([float(tmp) for tmp in lines[isub].split()])
            vtst_data['stress_matrix'] = stress
        elif 'FORCES: max atom' in line:
            forces = [float(tmp) for tmp in line.split()[-2:]]
            vtst_data['force_max_atom'] = forces[0]
            vtst_data['force_rms'] = forces[1]
        elif 'FORCE total and by dimension' in line:
            forces = [float(tmp) for tmp in line.split()[-2:]]
            vtst_data['force_total'] = forces[0]
            vtst_data['force_by_dimension'] = forces[1]
        elif 'Stress total and by dimension' in line:
            forces = [float(tmp) for tmp in line.split()[-2:]]
            vtst_data['stress_total'] = forces[0]
            vtst_data['stress_by_dimension'] = forces[1]
        elif 'OPT: skip step - force has converged' in line:
            vtst_data['neb_converged'] = True
        elif 'energy(sigma->0)' in line:
            tokens = line.split()
            vtst_data['energy_extrapolated'] = float(tokens[-1])
            vtst_data['energy_without_entropy'] = float(tokens[-4])
        elif 'free  energy   TOTEN' in line:
            vtst_data['free_energy'] = float(line.split()[-2])
        elif 'TOTAL-FORCE' in line:
            positions, forces = _parse_force_block(lines[idx : idx + nions + 10])
            inputs['outcar-forces'] = np.array(forces)
            inputs['outcar-positions'] = np.array(positions)
        elif 'direct lattice vectors' in line:
            cell = []
            for subline in lines[idx + 1 : idx + 4]:
                cell.append([float(tmp) for tmp in subline.split()[:3]])
            inputs['outcar-cell'] = np.array(cell)

    inputs['neb_data'] = vtst_data
    return inputs
