"""
Utils for setting symbols on chemical elements.

-----------------------------------------------
Utilities for choosing appropriate element symbols for a chemical element.
"""
from lxml import html
import requests

VERSION = {
    'latest': {
        'version': '5.2',
        'url': 'http://cms.mpi.univie.ac.at/vasp/vasp/Recommended_PAW_potentials_DFT_calculations_using_vasp_5_2.html',
        'gw-url': 'http://cms.mpi.univie.ac.at/vasp/vasp/Recommended_GW_PAW_potentials_vasp_5_2.html'
    }
}


def get_recommendations(version_nr='latest', use_gw=False):
    """
    Get recommendations for a certain type of PAW.

    :param version_nr: VASP version number
    :param use_gw: get recommendations for GW instead LDA pseudopotentials
    :return: recommendations dict
    """
    urlkey = 'gw-url' if use_gw else 'url'
    page = requests.get(VERSION[version_nr][urlkey])
    tree = html.fromstring(page.text)
    tags = tree.xpath('//table//td[1]/b')
    rec = {}
    for tag in tags:
        item = tag.text.strip().split(' ')
        element = item[0]
        symbol = '_'.join([i for i in item if i])
        rec[element] = symbol
    return rec


# pylint: disable=too-few-public-methods
class PawInfo(object):  # pylint: disable=useless-object-inheritance
    """Simple class to bundle and pass around info about a PAW."""

    def __init__(self, symbol, default_enmax, valency):
        self.symbol = symbol
        self.default_enmax = default_enmax
        self.valency = valency

    def __str__(self):
        return self.symbol

    def __repr__(self):
        return '<paw: {} at {}>'.format(self.symbol, hex(id(self)))


def get_all(version_nr='latest', use_gw=False):
    """
    Get recommendations for all symbols.

    :param version_nr: VASP version number
    :param use_gw:  Get recommendations for GW (default: LDA)
    :return: recommendations dict
    """
    urlkey = 'gw-url' if use_gw else 'url'
    page = requests.get(VERSION[version_nr][urlkey])
    tree = html.fromstring(page.text)
    tags = tree.xpath('//table/tr')
    syms = {}
    for tag in tags:
        row = tag.text_content().strip().split('\n')
        if row[1].isdigit():
            symboll = [i for i in row[0].strip().split(' ') if i]
            element = symboll[0]
            symbol = '_'.join(symboll)
            suffix = symboll[1] if len(symboll) > 1 else '_'
            if suffix == 'GW':
                suffix = '_'
            if not syms.get(element):
                syms[element] = {}
            row[0] = symbol
            row[1] = int(row[1])
            row[2] = float(row[2])
            syms[element][suffix] = (PawInfo(*row))
    return syms


if __name__ == '__main__':
    DEF_PAW = get_recommendations()
    DEF_GW = get_recommendations(use_gw=True)
    with open('default_paws.py', 'w') as defaults:
        defaults.write('lda = {\n')
        defaults.writelines(['"{}": "{}",\n'.format(k, v) for k, v in DEF_PAW.items()])
        defaults.write('}\n\n')
        defaults.write('gw = {\n')
        defaults.writelines(['"{}": "{}",\n'.format(k, v.replace('_GW', '')) for k, v in DEF_GW.items()])
        defaults.write('}\n\n')
