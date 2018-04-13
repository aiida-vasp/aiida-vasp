# pylint: disable=no-self-use
"""Utils for VASP KPOINTS format"""
import logging

from parsevasp.kpoints import Kpoints, Kpoint
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
        self.init_with_kwargs(**kwargs)

    def _init_with_path(self, filepath):
        self._filepath = filepath
        self._parsable_items = KpParser.PARSABLE_ITEMS
        self._parsed_data = {}

    def _init_with_kpointsdata(self, kpointsdata):
        """Initialise with a given kpointsData object"""
        logging.basicConfig()

        if kpointsdata.get_attrs().get('mesh'):
            mode = 'automatic'
        elif kpointsdata.get_attrs().get('array|kpoints'):
            mode = 'explicit'

        kpoints_dict = {}
        for keyword in ['comment', 'divisions', 'shifts', 'points', 'tetra', 'tetra_volume', 'mode', 'centering', 'num_kpoints']:
            kpoints_dict[keyword] = None

        kpoints_dict.update(getattr(self, '_get_kpointsdict_' + mode)(kpointsdata))

        try:
            self._data_obj = Kpoints(kpoints_dict=kpoints_dict)
        except SystemExit:
            self._data_obj = None

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
        """Turn an 'automatic' kpoints dictionary into Aiida KpointsData."""
        kpout = get_data_class('array.kpoints')()

        mesh = kpoints_dict.get('divisions')
        shifts = kpoints_dict.get('shifts')
        kpout.set_kpoints_mesh(mesh, offset=shifts)

        return kpout

    def _get_kpointsdict_explicit(self, kpointsdata):
        """Turn Aiida KpointData into an 'explicit' kpoints dictionary."""
        dictionary = {}

        kpts = []
        try:
            points, weights = kpointsdata.get_kpoints(also_weights=True)
        except AttributeError:
            points = kpointsdata.get_kpoints()
            weights = None
        for index, point in enumerate(points):
            if weights is not None:
                kpt = Kpoint(point, weight=weights[index])
            else:
                # no weights supplied, so set them to 1.0
                kpt = Kpoint(point, weight=1.0)
            kpts.append(kpt)
        dictionary["points"] = kpts
        dictionary["mode"] = "explicit"
        dictionary["num_kpoints"] = len(kpts)

        return dictionary

    def _get_kpointsdict_automatic(self, kpointsdata):
        """Turn Aiida KpointData into an 'explicit' kpoints dictionary."""
        dictionary = {}
        # automatic mode
        mesh = kpointsdata.get_kpoints_mesh()
        dictionary["divisions"] = mesh[0]
        dictionary["shifts"] = mesh[1]
        dictionary["mode"] = "automatic"
        # here we need to make a choice, so should
        # add more to Aiida to make this better
        # defined
        dictionary["centering"] = "Gamma"
        dictionary["num_kpoints"] = 0

        return dictionary
