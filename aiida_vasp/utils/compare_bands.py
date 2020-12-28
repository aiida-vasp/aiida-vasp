"""
Utils for comparing band structures.

------------------------------------
Utilities for comparing band structures. Mostly present for legacy purposes. Will be rewritten
or moved in the future.
"""
# pylint: disable=import-outside-toplevel
from aiida.plugins import DataFactory
from aiida.engine import calcfunction

BANDS_CLS = DataFactory('array.bands')


def _firstspin(bands):
    """Get only the bands for the first spin if multiple are contained."""
    if bands.ndim not in [2, 3]:
        raise ValueError('invalid input')
    if bands.ndim == 3:
        bands = bands[0]
    return bands


@calcfunction
# pylint: disable=too-many-locals
def make_reference_bands_inline(wannier_bands, vasp_bands, efermi=None):
    """
    Compare bandstructure results from wannier and vasp.

    Takes two input array.bands nodes, stores them if they're not already
    stored. Takes the relevant bands from the vasp bandstructure and stores and outputs
    them in a node with linkname 'bandcmp'.

    Also returns a parameter data node with linkname 'bandinfo' containing
    fermi energy, bandgap etc of the reference bandstructure.
    """
    import numpy as np
    assert isinstance(wannier_bands, BANDS_CLS)
    assert isinstance(vasp_bands, BANDS_CLS)
    assert hasattr(wannier_bands, 'labels')
    assert hasattr(vasp_bands, 'labels')
    if vasp_bands.labels:
        assert vasp_bands.labels == wannier_bands.labels
    kpcomp = vasp_bands.get_kpoints() == wannier_bands.get_kpoints()
    assert kpcomp.all(), 'kpoints may not differ'

    owindow = get_outer_window(wannier_bands)

    wbands = wannier_bands.get_bands()
    vbands, vocc = vasp_bands.get_bands(also_occupations=True)

    # throw away spin dimension if appropriate
    if vbands.ndim == 3:
        vbands = vbands[0]
        vocc = vocc[0]

    # grab the vbands within the outer_window
    # find wich bands within the window match
    # by searching for the best fit using the sum of square errors
    vbands_window = np.empty(wbands.shape)
    vocc_window = np.empty(wbands.shape)
    w_nbands = wbands.shape[1]
    ref_nbands = vbands.shape[1]
    for band_idx in range(w_nbands):
        errs = [band_error(wbands[:, band_idx], vbands[:, i]) for i in range(ref_nbands)]
        minerr = np.argmin(errs)
        vbands_window[:, band_idx] = vbands[:, minerr]
        vocc_window[:, band_idx] = vocc[:, minerr]

    # For the future:
    # * find each band's index (s, px, py, ...)
    # * store the legend with the comparison node

    # find fermi energy from vasp_bands parent or work without it
    if not efermi:
        try:
            efermi = vasp_bands.inp.bands.out.results.get_dict()['efermi']
        except Exception:  # pylint: disable=broad-except
            pass

    ref_gap_info = band_gap(vbands_window, vocc_window, efermi)
    ref_info = DataFactory('parameter')()
    ref_info.update_dict({'bandgap': ref_gap_info})
    ref_info.update_dict({'efermi': efermi})
    ref_info.update_dict({'outer_window': owindow})

    ref_bands = DataFactory('array.bands')()
    ref_bands.set_kpointsdata(wannier_bands)
    ref_bands.set_bands(vbands_window, occupations=vocc_window)

    return {'bands': ref_bands, 'info': ref_info}


def get_outer_window(bands_node, silent=False):
    """
    Get the ``outer_window`` parameter as a tuple (min, max), if it was given.

    Check if bands_node

    * is a child of a calculation and
    * that calculation has a parameter data input node with linkname parameters and
    * that node has the keys 'dis_win_min' and 'dis_win_max'.

    If that is the case, output outer_window = (min, max).
    """
    owindow = None
    try:
        calc = bands_node.inp.bands
        wset = calc.inp.parameters.get_dict()
        owindow = (wset['dis_win_min'], wset['dis_win_max'])
    except KeyError as err:
        if not silent:
            raise KeyError('Missing window parameters in input to ' 'parent calculation:\n' + str(err))
    except AttributeError as err:
        if not silent:
            raise AttributeError('bands_node is not an output of an appropriate calc node.' + str(err))
    return owindow


def band_gap(bands, occ, efermi=None):
    """
    Find the band gap in a bandstructure.

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
    """
    assert bands.shape == occ.shape
    result = {'gap': None, 'direct': None, 'vector': []}
    nbands = bands.shape[1]
    occupied = [i for i in range(nbands) if occ[:, i].any()]
    unoccupied = [i for i in range(nbands) if not occ[:, i].any()]
    # if either homo or lumo is not included, no info can be given
    if not occupied or not unoccupied:  # if any of them is empty
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
            raise ValueError(('The given E_fermi was higher than ' 'the lowest point of the lowest unoccupied band'))
    result['gap'] = gap_upper - gap_lower
    result['direct'] = bool(gap_upper_k == gap_lower_k)
    # check wether the two closest points are at the same kpoint (direct)
    # or not (indirect)
    result['vector'] = [(gap_lower_k, gap_lower), (gap_upper_k, gap_upper)]
    return result


def band_error(band1, band2):
    import numpy as np
    return np.square(band1 - band2).sum()


def bands_error(bands1, bands2):
    """
    Band for band rms error sqrt((|B1_i - B2_i|^2)/n) where BX_i is the i-th band of Band Structure Node X.

    Only works for BandsData nodes with 2d band arrays.
    """
    import numpy as np
    bands_1 = bands1.get_bands()
    bands_2 = bands2.get_bands()
    assert bands_1.shape == bands_2.shape
    nbands = bands_1.shape[1]
    err = np.empty(nbands)
    for band_idx in range(nbands):
        b1_i = bands_1[:, band_idx]
        b2_i = bands_2[:, band_idx]
        sq_err = np.square(b1_i - b2_i)  # pylint: disable=assignment-from-no-return
        m_err = sq_err.sum() / len(sq_err)
        err[band_idx] = np.sqrt(m_err)
    return err


# pylint: disable=too-many-locals
def compare_bands(vasp_bands, wannier_bands_list, plot_folder=None):
    """
    Compare a band structure from vasp with different ones from wannier90 obtained for different window parameters.

    :param vasp_bands: band structure output node from vasp calculation
    :param wannier_bands_list: list of band structure output nodes from wannier90 calculations
    :param plot_folder: if given, create a plot for each comparison in that folder
    :return:
    """
    import numpy as np
    import aiida_vasp.utils.bands as btool
    owindows = {get_outer_window(b): b for b in wannier_bands_list}
    ref_bands = {k: make_reference_bands_inline(wannier_bands=b, vasp_bands=vasp_bands) for k, b in owindows.items()}
    info = {}
    for wannier_bands in wannier_bands_list:
        owindow = get_outer_window(wannier_bands)
        reference = ref_bands[owindow]['bands']
        refinfo = ref_bands[owindow]['info'].get_dict()
        wannier_calc = wannier_bands.inp.bands
        wannier_param = wannier_calc.inp.parameters.get_dict()
        iwindow = [wannier_param['dis_froz_min'], wannier_param['dis_froz_max']]
        wannier_gap = band_gap(wannier_bands.get_bands(), reference.get_array('occupations'))
        ref_gap = refinfo['bandgap']
        if wannier_gap['vector']:
            wannier_k_gap = np.array([wannier_gap['vector'][0][0], wannier_gap['vector'][1][0]])
            ref_k_gap = np.array([wannier_gap['vector'][0][0], wannier_gap['vector'][1][0]])
            error_k_gap = np.abs(ref_k_gap - wannier_k_gap)
        else:
            error_k_gap = []
        info[wannier_bands.pk] = {
            'calc': wannier_calc.pk,
            'outer_window': owindow,
            'inner_window': iwindow,
            'error_per_band': bands_error(reference, wannier_bands),
            'error_e_gap': wannier_gap['gap'] and abs(wannier_gap['gap'] - ref_gap['gap']),
            'error_direct': wannier_gap['direct'] != ref_gap['direct'],
            'error_k_gap': error_k_gap
        }
        if plot_folder:
            import os
            colors = ['r', 'b', 'g', 'm', 'c', 'y', 'k']
            title = 'Vasp-Wannier comparison for window {}'.format([owindow, iwindow])
            fig = btool.plot_bstr(reference, efermi=refinfo['efermi'], colors=colors, title=title)
            btool.plot_bands(wannier_bands, colors=colors, figure=fig, ls=':')
            xlim = btool.plt.xlim()
            btool.plt.hlines(iwindow, xlim[0], xlim[1], color='k')
            btool.plt.hlines(refinfo['efermi'], xlim[0], xlim[1], color='k', linestyles='dashed')
            btool.plt.yticks(
                list(btool.plt.yticks()[0]) + [refinfo['efermi']], [str(line) for line in btool.plt.yticks()[0]] + [r'$E_{fermi}$'])
            pdf = os.path.join(plot_folder, 'comparison_%s.pdf' % wannier_calc.pk)
            fig.savefig(pdf)
            info[wannier_bands.pk]['plot'] = pdf

    return info


def compare_from_window_wf(workflow, **kwargs):
    """Find the relevant bands in the window workflow and compare them."""
    wblist = [v for k, v in workflow.get_results().iteritems() if 'bands_' in k]
    vbands = workflow.get_result('reference_bands')
    return compare_bands(vasp_bands=vbands, wannier_bands_list=wblist, **kwargs)


def plot_errors_vs_iwsize(comparison_info):
    """Plot Band structure errors versus size of the inner window parameter for wannier90."""
    import numpy as np
    import aiida_vasp.utils.bands as btool
    ows = []
    iws = []
    data = []
    for value in comparison_info.values():
        outer_window = value['outer_window']
        ows.append(outer_window[1] - outer_window[0])
        inner_window = value['inner_window']
        iws.append(inner_window[1] - inner_window[0])
        err = value['error_per_band'].sum()
        data.append(err)
    ows_i = np.sort(list(set(ows))).tolist()
    iws_i = np.sort(list(set(iws))).tolist()
    plot_data = np.empty((len(iws_i), len(ows_i)))
    for data_i, i in enumerate(data):
        owi = ows_i.index(ows[i])
        iwi = iws_i.index(iws[i])
        plot_data[iwi, owi] = data_i
    fig = btool.plt.figure()
    lines = btool.plt.plot(iws_i, plot_data, figure=fig)
    btool.plt.legend(lines, ows_i)
    return fig, zip(ows, iws, data)
