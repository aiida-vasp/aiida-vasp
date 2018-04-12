# pylint: disable=no-self-use
"""Utils for VASP KPOINTS format"""
import logging

from parsevasp.kpoints import Kpoints
from aiida_vasp.io.parser import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class


class KpParser(BaseFileParser):
    """
    Parser for VASP KPOINTS format

    This is a wrapper for the parsevasp.kpoints parser. It will read a given KPOINTS file and
    parse it into a Aiida KpointsData object.
    """

    PARSABLE_ITEMS = {
        'kpoints': {
            'inputs': [],
            'parsers': ['EIGENVAL', 'IBZKPT'],
            'nodeName': 'kpoints',
            'prerequisites': []
        },
    }

    def __init__(self, *args, **kwargs):
        super(KpParser, self).__init__(*args, **kwargs)
        self._parsable_items = KpParser.PARSABLE_ITEMS
        self._parsed_data = {}

    def _parse_file(self, inputs):
        """Create a DB Node from a KPOINTS file"""

        result = inputs
        result = {}

        try:
            parsed_kpoints = Kpoints(None, None, self._filepath, None)
        except SystemExit:
            return {'kpoints': None}

        mode = parsed_kpoints.entries.get('mode')
        result['kpoints'] = getattr(self, '_get_kpointsdata_' + mode)(parsed_kpoints.entries)

        return result

    def _get_kpointsdata_explicit(self, kpoints_dict):
        """Turn an 'explicit' kpoints dictionary into Aiida KpointsData"""
        kpout = get_data_class('array.kpoints')()

        kpoints = kpoints_dict.get('points')
        cartesian = not kpoints[0].get_direct()
        kpoint_list = []
        weights = []
        for kpoint in kpoints:
            kpoint_list.append(kpoint.get_point().tolist())
            weights.append(kpoint.get_weight())

        if weights[0] is None:
            weights = None

        kpout.set_kpoints(kpoint_list, weights=weights, cartesian=cartesian)

        return kpout

    def _get_kpointsdata_automatic(self, kpoints_dict):
        """Turn an 'automatic' kpoints dictionary into Aiida KpointsData"""
        kpout = get_data_class('array.kpoints')()

        mesh = kpoints_dict.get('divisions')
        shifts = kpoints_dict.get('shifts')
        kpout.set_kpoints_mesh(mesh, offset=shifts)

        return kpout
