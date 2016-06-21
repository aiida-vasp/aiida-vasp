from aiida.orm.calculation.inline import make_inline
from aiida.orm import DataFactory


BandsData = DataFactory('array.bands')
ArrayData = DataFactory('array')


@make_inline
def compare_bands(wannier_bands, vasp_bands):
    assert(isinstance(wannier_bands, BandsData))
    assert(isinstance(vasp_bands, BandsData))
    assert(hasattr(wannier_bands, 'labels'))
    assert(hasattr(vasp_bands, 'labels'))
    assert(vasp_bands.labels == wannier_bands.labels)

    try:
        calc = wannier_bands.inp.bands
        wset = calc.inp.settings.get_dict()
        owindow = (
            wset.get('dis_window_min'),
            wset.get('dis_window_max')
        )
        iwindow = (
            wset.get('dis_froz_min'),
            wset.get('dis_froz_max')
        )
        assert(len(filter(None, owindow)) == 2)
        assert(len(filter(None, iwindow)) == 2)
    except Exception as e:
        raise e

    comparison_node = ArrayData()

    wbands = wannier_bands.get_bands()
    vbands = vasp_bands.get_bands()

    # TODO: get v and w bands into same shape

    # TODO: grab the vbands within the outer_window
    vbands_window = '???'

    # TODO: find wich bands within the window match

    # TODO: find each band's index (s, px, py, ...)

    # TODO: store the legend with the comparison node

    comparison_node.set_array('tb', wbands)
    comparison_node.set_array('dft', vbands_window)

    return comparison_node
