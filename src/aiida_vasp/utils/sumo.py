"""
Module for plotting the AiiDA BandsData using sumo.
"""

import warnings
from importlib.util import find_spec
from typing import Optional, Union

import numpy as np
from aiida.orm import BandsData, CalcJobNode, StructureData
from pymatgen.core import Lattice
from pymatgen.electronic_structure.bandstructure import (
    BandStructureSymmLine,
    Spin,
)

if find_spec('pymatgen') is None or find_spec('sumo') is None:
    raise ImportError('This module requires the sumo and pymatgen packages to be installed.')

from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from sumo.electronic_structure.dos import load_dos
from sumo.electronic_structure.effective_mass import (
    fit_effective_mass,
    get_fitting_data,
)
from sumo.plotting import dos_plotter
from sumo.plotting.bs_plotter import SBSPlotter
from sumo.plotting.phonon_bs_plotter import SPhononBSPlotter

from .pmg import PymatgenAdapator


def get_sumo_dos_plotter(scf_node: CalcJobNode, **kwargs) -> dos_plotter.SDOSPlotter:
    """Get density of state by reading directly from the vasprun.xml file.

    :param scf_node: A node with `retrieved` output attached.
    :type scf_node: ProcessNode
    :param kwargs: additional parameters passed to `load_dos` function from sumo
    :returns: A `SDOSPlotter` object to be used for plotting the density of states.
    :rtype: SDOSPlotter
    """
    adapt = PymatgenAdapator(scf_node)
    vasprun = adapt.vasprun
    tdos, pdos = load_dos(vasprun, **kwargs)
    dp = dos_plotter.SDOSPlotter(tdos, pdos)
    return dp


def get_pmg_bandstructure(
    bands_node: BandsData, structure: StructureData = None, efermi: Optional[float] = None, **kwargs
) -> BandStructureSymmLine:
    """Return a pymatgen `BandStructureSymmLine` object from BandsData.

    :param bands_node: A BandsData object
    :type bands_node: BandsData
    :param structure: a StructureData object, required if `bands_node`
                      does not have information about the cell.
    :type structure: StructureData, optional
    :param efermi: Explicit value of the fermi energy.
    :type efermi: float, optional
    :param kwargs: additional keyword arguments
    :returns: A `BandStructureSymmLine` object
    :rtype: BandStructureSymmLine
    """
    if not isinstance(bands_node, BandsData):
        raise ValueError('The input argument must be a BandsData')
    # Load the data
    bands = bands_node.get_array('bands')  # In (num_spin, kpoints, bands) or just (kpoints, bands)
    kpoints = bands_node.get_array('kpoints')  # in (num_kpoints, 3)
    try:
        occupations = bands_node.get_array('occupations')
    except (KeyError, AttributeError):
        occupations = None

    labels = bands_node.base.attributes.get('labels')
    label_numbers = bands_node.base.attributes.get('label_numbers')

    # Construct the band_dict
    bands_shape = bands.shape
    if len(bands_shape) == 3:
        if bands_shape[0] == 2:
            bands_dict = {
                Spin.up: bands[0].T,  # Have to be (bands, kpoints)
                Spin.down: bands[1].T,  # Have to be (bands, kpoints)
            }
        else:
            bands_dict = {
                Spin.up: bands[0].T,  # Have to be (bands, kpoints)
            }
    else:
        bands_dict = {Spin.up: bands.T}

    if 'cell' in bands_node.base.attributes.keys():
        lattice = Lattice(bands_node.base.attributes.get('cell'))
    else:
        lattice = Lattice(structure.cell)

    # Constructure the label dictionary
    labels_dict = {}
    for label, number in zip(make_latex_labels(labels), label_numbers):
        labels_dict[label] = kpoints[number]

    # Try to use the fermi level from the node
    if efermi is None:
        get_efermi_from_band(bands_node)

    if efermi is None:
        if occupations is not None:
            # Use the middle of the CBM and VBM as the fermi energy....
            efermi = (find_vbm(bands, occupations) + find_cbm(bands, occupations)) / 2
        else:
            efermi = 0
            warnings.warn('Cannot find fermi energy - setting it to 0, this is probably wrong!')

    bands_structure = BandStructureSymmLine(
        kpoints,
        bands_dict,
        lattice.reciprocal_lattice,
        efermi=efermi,
        labels_dict=labels_dict,
        **kwargs,
    )
    return bands_structure


def get_sumo_bands_plotter(
    bands: BandsData, efermi: Optional[float] = None, structure: Optional[StructureData] = None, **kwargs
) -> SBSPlotter:
    """
    Return a sumo `SBSPlotter` object

    :param bands_node: A BandsData object
    :param structure (optional): a StructureData object, required if `bands_node`
      does not have information about the cell.
      efermi (float): Explicit value of the fermi energy.

    :returns: A `SBSPlotter` object
    """
    bands_structure = get_pmg_bandstructure(bands, efermi=efermi, structure=structure, **kwargs)
    return SBSPlotter(bands_structure)


def find_vbm(bands: np.ndarray, occupations: np.ndarray, tol: float = 1e-4) -> float:
    """
    Find the fermi energy, put it at the top of VBM
    NOTE: this differs from the fermi energy reported in VASP when there is any
    electronic smearing.
    """
    return bands[occupations > tol].max()


def find_cbm(bands: np.ndarray, occupations: np.ndarray, tol: float = 1e-4) -> float:
    """
    Find the fermi energy, put it at the top of VBM
    NOTE: this differs from the fermi energy reported in VASP when there is any
    electronic smearing.
    """
    return bands[occupations < tol].min()


def make_latex_labels(labels: list) -> list:
    """Convert labels to laxtex style"""
    label_mapping = {
        'GAMMA': r'\Gamma',
        'LAMBDA': r'\Lambda',
        'SIGMA': r'\Sigma',
    }
    out_labels = []
    for label in labels:
        label_ = label
        for tag, replace in label_mapping.items():
            if tag in label:
                label_ = label.replace(tag, replace)
                break
        out_labels.append(f'{label_}')
    return out_labels


def get_pymatgen_phonon_bands(
    band_structure: BandsData, input_structure: StructureData, has_nac: bool = False
) -> PhononBandStructureSymmLine:
    """
    Obtain a pymatgen phonon bandstructure plotter
    """
    qpoints = band_structure.get_kpoints()
    freq = np.transpose(band_structure.get_bands())  # Pymatgen uses (3 * natoms, number qpoints) for frequency
    structure = input_structure.get_pymatgen()
    lattice = structure.lattice.reciprocal_lattice
    idx, labels = zip(*band_structure.labels)
    labels = make_latex_labels(labels)
    labels_dict = {label: qpoints[idx] for idx, label in zip(idx, labels)}
    pbs = PhononBandStructureSymmLine(
        qpoints,
        freq,
        lattice,
        labels_dict=labels_dict,
        structure=structure,
        has_nac=has_nac,
    )
    return pbs


def get_sumo_phonon_plotter(
    band_structure: BandsData,
    input_structure: StructureData,
    has_nac: bool = False,
    imag_tol: float = -5e-2,
) -> SPhononBSPlotter:
    """
    Obtain a sumo phonon plotter object
    """
    bs = get_pymatgen_phonon_bands(band_structure, input_structure, has_nac)
    return SPhononBSPlotter(bs, imag_tol)


def bandstats(
    bs: Union[BandStructureSymmLine, BandsData],
    num_sample_points: int = 3,
    temperature: Optional[float] = None,
    degeneracy_tol: float = 1e-4,
    parabolic: bool = True,
    structure: Optional[StructureData] = None,
    efermi: Optional[float] = None,
    **kwargs,
):
    """Extract fitting data for band extrema based on spin, kpoint and band.

    NOTE: This function is modified based on sumo.cli.bandstats.band_stats

    Searches forward and backward from the extrema point, but will only sample
    there data if there are enough points in that direction.
    """
    if isinstance(bs, BandsData):
        bs = get_pmg_bandstructure(bs, structure=structure, efermi=efermi, **kwargs)

    if bs.is_metal():
        raise RuntimeError('ERROR: System is metallic!')

    vbm_data = bs.get_vbm()
    cbm_data = bs.get_cbm()

    if temperature:
        raise RuntimeError('ERROR: This feature is not yet supported!')

    else:
        # Work out where the hole and electron band edges are.
        # Fortunately, pymatgen does this for us. Points at which to calculate
        # the effective mass are identified as a tuple of:
        # (spin, band_index, kpoint_index)
        hole_extrema = []
        for spin, bands in vbm_data['band_index'].items():
            hole_extrema.extend([(spin, band, kpoint) for band in bands for kpoint in vbm_data['kpoint_index']])

        elec_extrema = []
        for spin, bands in cbm_data['band_index'].items():
            elec_extrema.extend([(spin, band, kpoint) for band in bands for kpoint in cbm_data['kpoint_index']])

        # extract the data we need for fitting from the band structure
        hole_data = []
        for extrema in hole_extrema:
            hole_data.extend(get_fitting_data(bs, *extrema, num_sample_points=num_sample_points))

        elec_data = []
        for extrema in elec_extrema:
            elec_data.extend(get_fitting_data(bs, *extrema, num_sample_points=num_sample_points))

    # calculate the effective masses and log the information
    for data in hole_data:
        eff_mass = fit_effective_mass(data['distances'], data['energies'], parabolic=parabolic)
        data['effective_mass'] = eff_mass

    for data in elec_data:
        eff_mass = fit_effective_mass(data['distances'], data['energies'], parabolic=parabolic)
        data['effective_mass'] = eff_mass

    return {'hole_data': hole_data, 'electron_data': elec_data}


def get_efermi_from_band(bands_node: BandsData) -> Optional[float]:
    """Get the fermi energy from a BandsData node"""
    efermi = bands_node.base.attributes.get('efermi', None)
    if efermi is None:
        efermi = bands_node.base.attributes.get('fermi_level', None)
    return efermi
