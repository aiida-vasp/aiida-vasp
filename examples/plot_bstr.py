#!runaiida
from aiida.orm import load_node
from matplotlib import pyplot as plt
import numpy as np
import sys

def plot_bstr(pk):
    calc =  load_node(pk)
    fig = plt.figure()
    bands = calc.out.bands.get_bands()
    nbands = bands.shape[2]
    nkp = bands.shape[1]
    nspin = bands.shape[0]
    for ispin in range(nspin):
        for iband in range(nbands):
            plt.plot(bands[ispin, :, iband])
    kplabs = calc.inp.kpoints.labels
    if kplabs:
        kpx = [i[0] for i in kplabs]
        kpl = [i[1] for i in kplabs]
        plt.xticks(kpx, kpl)
        plt.vlines(kpx, plt.ylim()[0], plt.ylim()[1])
    res = calc.out.results.get_dict()
    efermi = res.get('efermi')
    if efermi:
        plt.hlines(efermi, 0, nkp-1, linestyles='dashed')
    return fig

if __name__ == '__main__':
    pk = int(sys.argv[1])
    pfile = sys.argv[2]

    fig = plot_bstr(pk)
    fig.savefig(pfile)
