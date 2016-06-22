from aiida.orm.calculation.inline import make_inline
from aiida.orm import DataFactory


BandsData = DataFactory('array.bands')
ArrayData = DataFactory('array')


@make_inline
def make_refrerence_bands(wannier_bands, vasp_bands, efermi):
    '''
    Compare bandstructure results from wannier and vasp.
    Takes two input array.bands nodes, stores them if they're not already
    stored.
    Takes the relevant bands from the vasp bandstructure and stores and outputs
    them in a node with linkname 'bandcmp'.
    Also returns a parameter data node with linkname 'bandinfo' containing
    fermi energy, bandgap etc of the reference bandstructure.
    '''
    import numpy as np
    assert(isinstance(wannier_bands, BandsData))
    assert(isinstance(vasp_bands, BandsData))
    assert(hasattr(wannier_bands, 'labels'))
    assert(hasattr(vasp_bands, 'labels'))
    if vasp_bands.labels:
        assert(vasp_bands.labels == wannier_bands.labels)
    kpcomp = vasp_bands.get_kpoints() == wannier_bands.get_kpoints()
    assert(kpcomp.all())

    try:
        calc = wannier_bands.inp.bands
        wset = calc.inp.settings.get_dict()
        owindow = (
            wset['dis_win_min'],
            wset['dis_win_max']
        )
        # ~ iwindow = (
            # ~ wset['dis_froz_min'],
            # ~ wset['dis_froz_max']
        # ~ )
    except KeyError as e:
        msg = ('Missing window parameters in input to '
               'wannier_bands parent calculation:\n') + e.message
        raise KeyError(msg)
    except AttributeError as e:
        msg = ('input wannier_bands must have a WannierCalculation parent!'
               + e.message)
        raise AttributeError(msg)

    comparison_node = ArrayData()

    wbands = wannier_bands.get_bands()
    vbands, vocc = vasp_bands.get_bands(also_occupations=True)

    # throw away spin dimension if appropriate
    if vbands.ndim == 3:
        vbands = vbands[0]
        vocc = vocc[0]

    # grab the vbands within the outer_window
    # find wich bands within the window match
    # by filling up from the bottom (as w90 seems to do so)
    vbands_window = np.empty(wbands.shape)
    vocc_window = np.empty(wbands.shape)
    count = 0
    for b in range(vbands.shape[1]):
        if count < wbands.shape[1]:
            band = vbands[:, b]
            if band.min() >= owindow[0] and band.max <= owindow[1]:
                vbands_window[:, count] = band
                vocc_window[:, count] = vocc[:, b]
                count += 1

    # TODO: find each band's index (s, px, py, ...)
    # TODO: store the legend with the comparison node

    # TODO: find fermi energy from vasp_bands parent
    if not efermi:
        efermi = vasp_bands.inp.bands.out.results.get_dict()['efermi']

    comparison_node.set_array('wannier_bands', wbands)
    comparison_node.set_array('reference', vbands_window)

    return {'bandcmp': comparison_node}


def band_gap(bands, occ, efermi, get_pos=False):
    assert(bands.shape == occ.shape)
    nbands = bands.shape[1]
    occupied = [i for i in range(nbands) if occ[:, i].any()]
    unoccupied = [i for i in range(nbands) if not occ[:, i].any()]
    # highest band with any occupation
    homo = bands[:, max(occupied)]
    # lowest completely unoccupied band
    lumo = bands[:, min(unoccupied)]
    gap_lower = homo.max()
    gap_lower_k = homo.argmax()
    # if homo crosses efermi, there is no band gap
    if gap_lower > efermi:
        return 0, None, []
    gap_upper = lumo.min()
    gap_upper_k = lumo.argmin()
    # if lumo crosses efermi, something is wrong
    if gap_upper < efermi:
        raise ValueError(('The given E_fermi was higher than '
                          'the lowest point of the lowest unoccupied band'))
    gap = gap_upper - gap_lower
    direct = gap_upper_k == gap_lower_k
    # TODO: check wether the two closest points are at the same kpoint (direct)
    # TODO: or not (indirect)
    return gap, direct, [gap_lower, gap_upper]
