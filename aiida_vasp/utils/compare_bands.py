from aiida.orm.calculation.inline import optional_inline
from aiida.orm import DataFactory

BandsData = DataFactory('array.bands')


def _firstspin(bands):
    if bands.ndim not in [2, 3]:
        raise ValueError('invalid input')
    if bands.ndim == 3:
        bands = bands[0]
    return bands


@optional_inline
def make_reference_bands_inline(wannier_bands, vasp_bands, efermi=None):
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
    assert (isinstance(wannier_bands, BandsData))
    assert (isinstance(vasp_bands, BandsData))
    assert (hasattr(wannier_bands, 'labels'))
    assert (hasattr(vasp_bands, 'labels'))
    if vasp_bands.labels:
        assert (vasp_bands.labels == wannier_bands.labels)
    kpcomp = vasp_bands.get_kpoints() == wannier_bands.get_kpoints()
    assert (kpcomp.all(), 'kpoints may not differ')

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
    count = 0
    for b in range(w_nbands):
        errs = [
            band_error(wbands[:, b], vbands[:, i]) for i in range(ref_nbands)
        ]
        minerr = np.argmin(errs)
        vbands_window[:, b] = vbands[:, minerr]
        vocc_window[:, b] = vocc[:, minerr]

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
    ref_info.update_dict({'outer_window': owindow})

    ref_bands = DataFactory('array.bands')()
    ref_bands.set_kpointsdata(wannier_bands)
    ref_bands.set_bands(vbands_window, occupations=vocc_window)

    return {'bands': ref_bands, 'info': ref_info}


def get_outer_window(bands_node, silent=False):
    '''
    Check if bands_node is a child of a calculation and that calculation
    has a parameter data input node with linkname parameters and that
    node has the keys 'dis_win_min' and 'dis_win_max'.
    If that is the case, output outer_window = (min, max).
    '''
    owindow = None
    try:
        calc = bands_node.inp.bands
        wset = calc.inp.parameters.get_dict()
        owindow = (wset['dis_win_min'], wset['dis_win_max'])
        # ~ iwindow = (
        # ~     wset['dis_froz_min'],
        # ~     wset['dis_froz_max']
        # ~ )
    except KeyError as e:
        if not silent:
            msg = ('Missing window parameters in input to '
                   'parent calculation:\n') + e.message
            raise KeyError(msg)
    except AttributeError as e:
        if not silent:
            msg = ('bands_node is not an output of an appropriate calc node.' +
                   e.message)
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
    assert (bands.shape == occ.shape)
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
    result['direct'] = bool(gap_upper_k == gap_lower_k)
    # check wether the two closest points are at the same kpoint (direct)
    # or not (indirect)
    result['vector'] = [(gap_lower_k, gap_lower), (gap_upper_k, gap_upper)]
    return result


def band_error(band1, band2):
    import numpy as np
    return np.square(band1 - band2).sum()


def bands_error(bands1, bands2):
    '''
    band for band rms error sqrt((|B1_i - B2_i|^2)/n) where BX_i is the i-th band of
    Band Structure Node X.
    Only works for BandsData nodes with 2d band arrays.
    '''
    import numpy as np
    b1 = bands1.get_bands()
    b2 = bands2.get_bands()
    assert (b1.shape == b2.shape)
    nbands = b1.shape[1]
    err = np.empty(nbands)
    for bi in range(nbands):
        b1_i = b1[:, bi]
        b2_i = b2[:, bi]
        sq_err = np.square(b1_i - b2_i)
        m_err = sq_err.sum() / len(sq_err)
        err[bi] = np.sqrt(m_err)
    return err


def compare_bands(vasp_bands, wannier_bands_list, plot_folder=None):
    import numpy as np
    import bands as btool
    owindows = {get_outer_window(b): b for b in wannier_bands_list}
    ref_bands = {
        k: make_reference_bands_inline(wannier_bands=b, vasp_bands=vasp_bands)
        for k, b in owindows.iteritems()
    }
    info = {}
    for wannier_bands in wannier_bands_list:
        owindow = get_outer_window(wannier_bands)
        reference = ref_bands[owindow]['bands']
        refinfo = ref_bands[owindow]['info'].get_dict()
        wannier_calc = wannier_bands.inp.bands
        wannier_param = wannier_calc.inp.parameters.get_dict()
        iwindow = [
            wannier_param['dis_froz_min'], wannier_param['dis_froz_max']
        ]
        wannier_gap = band_gap(wannier_bands.get_bands(),
                               reference.get_array('occupations'))
        ref_gap = refinfo['bandgap']
        if wannier_gap['vector']:
            wannier_k_gap = np.array(
                [wannier_gap['vector'][0][0], wannier_gap['vector'][1][0]])
            ref_k_gap = np.array(
                [wannier_gap['vector'][0][0], wannier_gap['vector'][1][0]])
            error_k_gap = np.abs(ref_k_gap - wannier_k_gap)
        else:
            error_k_gap = []
        info[wannier_bands.pk] = {
            'calc':
            wannier_calc.pk,
            'outer_window':
            owindow,
            'inner_window':
            iwindow,
            'error_per_band':
            bands_error(reference, wannier_bands),
            'error_e_gap':
            wannier_gap['gap'] and abs(wannier_gap['gap'] - ref_gap['gap']),
            'error_direct':
            wannier_gap['direct'] != ref_gap['direct'],
            'error_k_gap':
            error_k_gap
        }
        if plot_folder:
            import os
            colors = ['r', 'b', 'g', 'm', 'c', 'y', 'k']
            title = 'Vasp-Wannier comparison for window {}'.format(
                [owindow, iwindow])
            fig = btool.plot_bstr(
                reference,
                efermi=refinfo['efermi'],
                colors=colors,
                title=title)
            btool.plot_bands(wannier_bands, colors=colors, figure=fig, ls=':')
            xlim = btool.plt.xlim()
            btool.plt.hlines(iwindow, xlim[0], xlim[1], color='k')
            btool.plt.hlines(
                refinfo['efermi'],
                xlim[0],
                xlim[1],
                color='k',
                linestyles='dashed')
            btool.plt.yticks(
                list(btool.plt.yticks()[0]) + [refinfo['efermi']],
                [str(l) for l in btool.plt.yticks()[0]] + [r'$E_{fermi}$'])
            pdf = os.path.join(plot_folder,
                               'comparison_%s.pdf' % wannier_calc.pk)
            fig.savefig(pdf)
            info[wannier_bands.pk]['plot'] = pdf

    return info


def compare_from_window_wf(wf, **kwargs):
    wblist = [v for k, v in wf.get_results().iteritems() if 'bands_' in k]
    vbands = wf.get_result('reference_bands')
    return compare_bands(
        vasp_bands=vbands, wannier_bands_list=wblist, **kwargs)


def plot_errors_vs_iwsize(comparison_info):
    import numpy as np
    import bands as btool
    ows = []
    iws = []
    data = []
    for k, v in comparison_info.iteritems():
        ow = v['outer_window']
        ows.append(ow[1] - ow[0])
        iw = v['inner_window']
        iws.append(iw[1] - iw[0])
        err = v['error_per_band'].sum()
        data.append(err)
    ows_i = np.sort(list(set(ows))).tolist()
    iws_i = np.sort(list(set(iws))).tolist()
    plot_data = np.empty((len(iws_i), len(ows_i)))
    for i in range(len(data)):
        owi = ows_i.index(ows[i])
        iwi = iws_i.index(iws[i])
        plot_data[iwi, owi] = data[i]
    fig = btool.plt.figure()
    lines = btool.plt.plot(iws_i, plot_data, figure=fig)
    btool.plt.legend(lines, ows_i)
    return fig, zip(ows, iws, data)
