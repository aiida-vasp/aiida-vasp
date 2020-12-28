"""
Utils for bands structures.

---------------------------
Utilities for working with band structures. Currently this is legacy and will be
rewritten or moved.
"""
# pylint: disable=import-outside-toplevel
try:
    import matplotlib
    matplotlib.use('TKAgg')
    from matplotlib import pyplot as plt
except ImportError:
    raise ImportError('Error: matplotlib must be ' + 'installed to use this functionality')


def get_bs_dims(bands_array):
    """
    Get the dimensions from the bands array of a BandsData node.

    :param numpy.array bands_array:
        an array with bands as stored in an array.bands data node
    :return: a tuple containing num_bands, num_kp, num_spins.
        if the array is only 2d, num_spins = 0
    :rtype tuple:
    """
    bshape = bands_array.shape
    nbd = nkp = nsp = 0
    if len(bshape) == 2:
        nbd = bshape[1]
        nkp = bshape[0]
    elif len(bshape) == 3:
        nbd = bshape[2]
        nkp = bshape[1]
        nsp = bshape[0]
    return nbd, nkp, nsp


def get_kp_labels(bands_node, kpoints_node=None):
    """
    Get Kpoint labels with their x-positions in matplotlib compatible format.

    A KpointsData node can optionally be given to fall back to if no labels
    are found on the BandsData node. The caller is responsible for ensuring
    the nodes match.  This should be the case if you take the kpoints from
    the input and the bands from the
    output of a calculation node.

    :param BandsData bands_node:
        The BandsData node will be searched labels first
    :param KpointsData kpoints_node:
        The optional KpointsData node will be searched only if no labels are
        present on the BandsData node. No consistency checks are performed.
    :return:  (kpx, kpl), the x-coordinates and text labels
    :rtype: tuple(list[int], list[unicode])
    :raises AttributeError: if neither of the given nodes have a labels
        attribute
    """
    kplabs = None
    kpx = []
    kpl = []
    try:
        kplabs = bands_node.labels
    except AttributeError as err:
        if kpoints_node:
            kplabs = kpoints_node.labels
        else:
            raise err
    if kplabs:
        kpx = [i[0] for i in kplabs]
        kpl = [i[1] for i in kplabs]
        for i, kpoints in enumerate(kpl):
            if kpoints == 'G':
                kpl[i] = r'$\Gamma$'
    return kpx, kpl


def get_efermi(calc):
    """Get the fermi energy from a finished calculation."""
    efermi = None
    if calc:
        p_res = calc.get_outputs_dict().get('results')
        efermi = p_res and p_res.get_dict().get('efermi')
    return efermi


def get_kp_node(calc):
    kpoints_node = None
    if calc:
        kpoints_node = calc.get_inputs_dict().get('kpoints')
    return kpoints_node


def plot_bstr(bands_node, kpoints_node=None, title=None, efermi=None, use_parent_calc=False, **kwargs):
    """
    Use matplotlib to plot the bands stored in a BandsData node.

    A KpointsData node can optionally be given as a fallback for
    kpoint labels. The caller is responsible for giving a node
    with matching labels (as in they are in/out nodes of the same
    calculation).

    :param BandsData bands_node:
        The BandsData node will be searched labels first
    :param KpointsData kpoints_node:
        The optional KpointsData node will be searched only if no labels are
        present on the BandsData node. No consistency checks are performed.
    :return: the matplotlib figure containing the plot
    """
    fig = plt.figure()
    title = title or 'Band Structure (pk=%s)' % bands_node.pk
    bands = bands_node.get_bands()
    _, nkp, _ = get_bs_dims(bands)
    plot_bands(bands_node, **kwargs)

    parent_calc = None
    if use_parent_calc:
        inputs = bands_node.get_inputs()
        parent_calc = inputs[0] if inputs else None

    efermi = get_efermi(parent_calc)
    kpoints_node = get_kp_node(parent_calc)

    if efermi:
        plt.hlines(efermi, plt.xlim()[0], nkp - 1, linestyles='dashed')
        plt.yticks(list(plt.yticks()[0]) + [efermi], [str(line) for line in plt.yticks()[0]] + [r'$E_{fermi}$'])

    try:
        kpx, kpl = get_kp_labels(bands_node, kpoints_node)
        plt.xticks(kpx, kpl)
        plt.vlines(kpx, plt.ylim()[0], plt.ylim()[1])
    except Exception:  # pylint: disable=broad-except
        pass

    plt.ylabel('Dispersion')
    plt.suptitle(title)
    return fig


def plot_bands(bands_node, **kwargs):
    """Plot a bandstructure node using matplotlib."""
    import numpy as np

    bands = bands_node.get_bands()
    nbands, nkp, nspin = get_bs_dims(bands)
    if nspin > 0:
        allbands = np.empty((nkp, nbands * nspin))
        for i in range(nspin):
            allbands[:, i * nbands:(i + 1) * nbands] = bands[i]
        bands = allbands

    if 'colors' in kwargs:
        import itertools
        colors = itertools.cycle(kwargs.pop('colors'))
        for b_idx in range(bands.shape[1]):
            plt.plot(bands[:, b_idx], color=colors.next(), **kwargs)  # pylint: disable=no-member
    else:
        plt.plot(bands, **kwargs)
