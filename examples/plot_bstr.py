#!runaiida
from aiida.orm import load_node
from aiida_vasp.utils.bands import *
from matplotlib import pyplot as plt
import sys


# ~ def get_bs_dims(bands_array):
    # ~ '''
    # ~ py:function:: get_bs_dims(bands_array)
# ~
    # ~ get the dimensions from the bands array of a BandsData node
# ~
    # ~ :param numpy.array bands_array:
        # ~ an array with bands as stored in an array.bands data node
    # ~ :return: a tuple containing num_bands, num_kp, num_spins.
        # ~ if the array is only 2d, num_spins = 0
    # ~ :rtype tuple:
    # ~ '''
    # ~ bshape = bands_array.shape
    # ~ nbd = nkp = nsp = 0
    # ~ if len(bshape) == 2:
        # ~ nbd = bshape[1]
        # ~ nkp = bshape[0]
    # ~ elif len(bshape) == 3:
        # ~ nbd = bshape[2]
        # ~ nkp = bshape[1]
        # ~ nsp = bshape[0]
    # ~ return nbd, nkp, nsp

# ~
# ~ def get_kp_labels(bands_node, kpoints_node=None):
    # ~ '''
    # ~ py:function:: get_kp_labels(bands_node[, kpoints_node=None])
# ~
    # ~ get Kpoint labels with their x-positions in matplotlib compatible format.
    # ~ A KpointsData node can optionally be given to fall back to if no labels
    # ~ are found on the BandsData node. The caller is responsible for ensuring
    # ~ the nodes match.  This should be the case if you take the kpoints from
    # ~ the input and the bands from the
    # ~ output of a calculation node.
# ~
    # ~ :param BandsData bands_node:
        # ~ The BandsData node will be searched labels first
    # ~ :param KpointsData kpoints_node:
        # ~ The optional KpointsData node will be searched only if no labels are
        # ~ present on the BandsData node. No consistency checks are performed.
    # ~ :return:  (kpx, kpl), the x-coordinates and text labels
    # ~ :rtype: tuple(list[int], list[unicode])
    # ~ :raises AttributeError: if neither of the given nodes have a labels
        # ~ attribute
    # ~ '''
    # ~ kplabs = None
    # ~ kpx = None
    # ~ kpl = None
    # ~ try:
        # ~ kplabs = bands_node.labels
    # ~ except AttributeError as e:
        # ~ if kpoints_node:
            # ~ kplabs = kpoints_node.labels
        # ~ else:
            # ~ raise e
    # ~ if kplabs is not None:
        # ~ kpx = [i[0] for i in kplabs]
        # ~ kpl = [i[1] for i in kplabs]
    # ~ return kpx, kpl
# ~
# ~
# ~ def plot_bstr(bands_node, kpoints_node=None):
    # ~ '''
    # ~ py:function:: plot_bstr(bands_node[, kpoints_node=None])
# ~
    # ~ Use matplotlib to plot the bands stored in a BandsData node.
    # ~ A KpointsData node can optionally be given as a fallback for
    # ~ kpoint labels. The caller is responsible for giving a node
    # ~ with matching labels (as in they are in/out nodes of the same
    # ~ calculation).
# ~
    # ~ :param BandsData bands_node:
        # ~ The BandsData node will be searched labels first
    # ~ :param KpointsData kpoints_node:
        # ~ The optional KpointsData node will be searched only if no labels are
        # ~ present on the BandsData node. No consistency checks are performed.
    # ~ :return: the matplotlib figure containing the plot
    # ~ '''
# ~
    # ~ fig = plt.figure()
    # ~ bands = bands_node.get_bands()
    # ~ nbands, nkp, nspin = get_bs_dims(bands)
    # ~ for iband in range(nbands):
        # ~ if nspin > 0:
            # ~ for ispin in range(nspin):
                # ~ plt.plot(bands[ispin, :, iband])
        # ~ else:
            # ~ plt.plot(bands[:, iband])
    # ~ kpx, kpl = get_kp_labels(bands_node, kpoints_node)
    # ~ if kpx and kpl:
        # ~ plt.xticks(kpx, kpl)
        # ~ plt.vlines(kpx, plt.ylim()[0], plt.ylim()[1])
    # ~ return fig


def plot_from_calc(calc, bands_lname='bands', kp_lname='kpoints'):
    '''
    py:function:: plot_from_calc(calc[, **kwargs])

    Plots the bandstructure from a calculation node.
    If available, draws vertical lines for special kpoints
    and lables them, as well as a dashed hline for E_fermi.

    Kpoint labels are searched in the output node named by bands_lname,
    if not found there and an input node is present at kp_lname, this
    is searched.

    E_fermi is expected to be found in a 'results' output node.

    :param JobCalculation calc:
        calculation node with an output BandsData node.
    :param str [bands_lname]: default='bands', the output linkname
        for the BandsData node
    :param str [kp_lname]: default='kpoints', the input linkname for
        the KpointsData node

    Example:
        calc = load_node(<pk>)
        fig = plot_from_calc(calc)
        fig.savefig('example.pdf')
    '''
    outs = calc.get_outputs_dict()
    inps = calc.get_inputs_dict()
    bnode = outs[bands_lname]
    knode = outs.get(kp_lname)
    if not knode:
        knode = inps.get(kp_lname)
    res = outs.get('results')
    nkp = get_bs_dims(bnode.get_bands())[1]
    if res:
        resd = res.get_dict()
        efermi = resd.get('efermi')
        if efermi:
            plt.hlines(efermi, 0, nkp-1, linestyles='dashed')
    return plot_bstr(bnode, kpoints_node=knode)

if __name__ == '__main__':
    pk = int(sys.argv[1])
    pfile = sys.argv[2]

    fig = plot_from_calc(calc=load_node(pk))
    fig.savefig(pfile)
