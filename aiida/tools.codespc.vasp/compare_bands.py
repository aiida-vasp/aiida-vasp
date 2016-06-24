from aiida.orm.calculation.inline import make_inline
from aiida.orm import DataFactory


BandsData = DataFactory('array.bands')
ArrayData = DataFactory('array')


@make_inline
def make_refrerence_bands(wannier_bands, vasp_bands, efermi=None):
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
    assert(kpcomp.all(), 'kpoints may not differ')

    owindow = get_outer_window(wannier_bands)

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
    w_nbands = wbands.shape[1]
    ref_nbands = vbands.shape[1]
    count = 0
    for b in range(ref_nbands):
        if count < w_nbands:
            band = vbands[:, b]
            occ = vocc[:, b]
            if band.min() >= owindow[0] and band.max() <= owindow[1]:
                vbands_window[:, count] = band
                vocc_window[:, count] = occ
                count += 1

    # TODO: find each band's index (s, px, py, ...)
    # TODO: store the legend with the comparison node

    # find fermi energy from vasp_bands parent or work without it
    if not efermi:
        try:
            efermi = vasp_bands.inp.bands.out.results.get_dict()['efermi']
        except:
            pass

    ref_gap_info = band_gap(vbands_window, vocc_window, efermi)
    ref_info = DataFactory('parameter')()
    ref_info.update_dict({'bandgap': ref_gap_info})
    ref_info.update_dict({'efermi': efermi})

    ref_bands = DataFactory('array.bands')()
    ref_bands.set_kpointsdata(wannier_bands)
    ref_bands.set_bands(vbands_window, occupations=vocc_window)

    return {'bands': ref_bands, 'info': ref_info}


def get_outer_window(bands_node, silent=False):
    '''
    Check if bands_node is a child of a calculation and that calculation
    has a parameter data input node with linkname settings and that
    node has the keys 'dis_win_min' and 'dis_win_max'.
    If that is the case, output outer_window = (min, max).
    '''
    owindow = None
    try:
        calc = bands_node.inp.bands
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
        if not silent:
            msg = ('Missing window parameters in input to '
                   'parent calculation:\n') + e.message
            raise KeyError(msg)
    except AttributeError as e:
        if not silent:
            msg = ('bands_node is not an output of an appropriate calc node.'
                   + e.message)
            raise AttributeError(msg)
    return owindow


def band_gap(bands, occ, efermi=None):
    '''
    find the band gap in a bandstructure
    :param numpy.array bands:
        2D bands array (as from BandsData.get_bands())
    :param numpy.array occ:
        2D occupations array matching bands
    :param float efermi:
        optional: Fermi Energy level
    :return:
        {'gap': gap energy difference,
         'direct': bool wether it's a direct band gap
         'vector': start and end points to draw the
         band gap as an arrow in the bandstructure plot.
         the points consist of (k, Dispersion) coordinates
         }
    '''
    assert(bands.shape == occ.shape)
    result = {'gap': None, 'direct': None, 'vector': []}
    nbands = bands.shape[1]
    occupied = [i for i in range(nbands) if occ[:, i].any()]
    unoccupied = [i for i in range(nbands) if not occ[:, i].any()]
    # if either homo or lumo is not included, no info can be given
    if len(occupied) == 0 or len(unoccupied) == 0:
        return result
    # highest band with any occupation
    homo = bands[:, max(occupied)]
    # lowest completely unoccupied band
    lumo = bands[:, min(unoccupied)]
    gap_lower = homo.max()
    gap_lower_k = homo.argmax()
    # if homo crosses efermi, there is no band gap
    if efermi:
        if gap_lower > efermi:
            result['gap'] = 0
            result['direct'] = False
            return result
    gap_upper = lumo.min()
    gap_upper_k = lumo.argmin()
    # if lumo crosses efermi, something is wrong
    if efermi:
        if gap_upper < efermi:
            raise ValueError(
                ('The given E_fermi was higher than '
                 'the lowest point of the lowest unoccupied band'))
    result['gap'] = gap_upper - gap_lower
    result['direct'] = gap_upper_k == gap_lower_k
    # check wether the two closest points are at the same kpoint (direct)
    # or not (indirect)
    result['vector'] = [(gap_lower_k, gap_lower), (gap_upper_k, gap_upper)]
    return result
