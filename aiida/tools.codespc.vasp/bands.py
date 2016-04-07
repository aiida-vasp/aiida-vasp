def get_bs_dims(bands_array):
    '''
    py:function:: get_bs_dims(bands_array)

    get the dimensions from the bands array of a BandsData node

    :param numpy.array bands_array:
        an array with bands as stored in an array.bands data node
    :return: a tuple containing num_bands, num_kp, num_spins.
        if the array is only 2d, num_spins = 0
    :rtype tuple:
    '''
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
    '''
    py:function:: get_kp_labels(bands_node[, kpoints_node=None])

    get Kpoint labels with their x-positions in matplotlib compatible format.
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
    '''
    kplabs = None
    kpx = None
    kpl = None
    try:
        kplabs = bands_node.labels
    except AttributeError as e:
        if kpoints_node:
            kplabs = kpoints_node.labels
        else:
            raise e
    if kplabs is not None:
        kpx = [i[0] for i in kplabs]
        kpl = [i[1] for i in kplabs]
        for i, kp in enumerate(kpl):
            if kp == 'G':
                kpl[i] = r'$\Gamma$'
    return kpx, kpl


def plot_bstr(bands_node, kpoints_node=None, title=None,
              use_parent_calc=False):
    '''
    py:function:: plot_bstr(bands_node[, kpoints_node=None])

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
    '''

    from matplotlib import pyplot as plt
    fig = plt.figure()
    if not title:
        title = 'Band Structure (pk=%s)' % bands_node.pk
    bands = bands_node.get_bands()
    nbands, nkp, nspin = get_bs_dims(bands)
    for iband in range(nbands):
        if nspin > 0:
            for ispin in range(nspin):
                plt.plot(bands[ispin, :, iband])
        else:
            plt.plot(bands[:, iband])
    if use_parent_calc:
        inputs = bands_node.get_inputs()
        parent_calc = inputs and inputs[0] or None
        kpoints_node = parent_calc.get_inputs_dict().get('kpoints')
        p_res = parent_calc.get_outputs_dict().get('results')
        efermi = p_res and p_res.get_dict().get('efermi')
        if efermi:
            plt.hlines(efermi, plt.xlim()[0], nkp-1, linestyles='dashed')
            plt.yticks(list(plt.yticks()[0]) + [efermi],
                       [str(l) for l in plt.yticks()[0]] + [r'$E_{fermi}$'])
    kpx, kpl = get_kp_labels(bands_node, kpoints_node)
    if kpx and kpl:
        plt.xticks(kpx, kpl)
        plt.vlines(kpx, plt.ylim()[0], plt.ylim()[1])
    plt.ylabel('Dispersion')
    plt.suptitle(title)
    return fig
