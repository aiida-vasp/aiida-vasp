"""
KPOINTS parser.

---------------
The object parser that handles the parsing of KPOINTS.
"""
# pylint: disable=no-self-use

from parsevasp.kpoints import Kpoints, Kpoint
from aiida_vasp.parsers.object_parsers.parser import BaseFileParser
from aiida_vasp.parsers.node_composer import NodeComposer, get_node_composer_inputs_from_object_parser
from aiida_vasp.utils.aiida_utils import get_data_class


class KpointsParser(BaseFileParser):
    """
    Parser for VASP KPOINTS.

    This is a wrapper for the parsevasp.kpoints parser. It will convert
    KPOINTS representatation from parsevasp to AiiDA KpointsData objects and vice versa.

    The Parsing direction depends on whether the KpointsParser is initialised with
    'handler = ...' (use handler object) or 'data = ...' (read from AiiDA data structure).

    """

    PARSABLE_ITEMS = {
        'kpoints-kpoints': {
            'inputs': [],
            'name': 'kpoints',
            'prerequisites': [],
        },
    }

    def __init__(self, *args, **kwargs):
        """
        Initialize INCAR parser

        data : KpointsArray

        """
        super(KpointsParser, self).__init__(*args, **kwargs)
        self._kpoints = None
        if 'data' in kwargs:
            self._init_kpoints(kwargs['data'])

    def _init_kpoints(self, data):
        """Initialize with a given AiiDA KpointsData instance."""
        if isinstance(data, get_data_class('array.kpoints')):
            self._data_obj = data
        else:
            self._logger.warning('Please supply an AiiDA KpointsData datatype for `data`.')
            self._data_obj = None
        self._kpoints = data
        self._parsable_items = self.__class__.PARSABLE_ITEMS
        self._parsed_data = {}

    @property
    def _parsed_object(self):
        """
        Return an instance of parsevasp.Kpoints.

        Corresponds to the stored KpointsData, but with a different representation.
        """

        if isinstance(self._data_obj, get_data_class('array.kpoints')):
            # The KpointsData has not been successfully parsed yet. So let's parse it.
            try:
                _ = self._data_obj.get_attribute('mesh')
                mode = 'automatic'
            except AttributeError:
                pass

            try:
                _ = self._data_obj.get_attribute('array|kpoints')
                mode = 'explicit'
            except AttributeError:
                pass

            kpoints_dict = {}
            for keyword in ['comment', 'divisions', 'shifts', 'points', 'tetra', 'tetra_volume', 'mode', 'centering', 'num_kpoints']:
                kpoints_dict[keyword] = None

            kpoints_dict.update(getattr(self, '_get_kpointsdict_' + mode)(self._data_obj))

            try:
                return Kpoints(kpoints_dict=kpoints_dict, logger=self._logger)
            except SystemExit:
                return None

        # _data_obj is SingleFile:
        return self._data_obj

    def _parse_object(self, inputs):
        """Create a DB Node from KPOINTS."""

        result = inputs
        result = {}

        if isinstance(self._data_obj, get_data_class('array.kpoints')):
            return {'kpoints-kpoints': self._data_obj}

        try:
            parsed_kpoints = Kpoints(file_handler=self._data_obj.handler, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exitited abnormally. Returning None.')
            return {'kpoints-kpoints': None}

        if parsed_kpoints.entries.get('mode') == 'line':
            self._logger.warning('The read KPOINTS contained line mode which is' 'not supported. Returning None.')
            return {'kpoints-kpoints': None}
        result['kpoints-kpoints'] = parsed_kpoints.entries

        return result

    @property
    def kpoints(self):
        if self._kpoints is None:
            inputs = get_node_composer_inputs_from_object_parser(self, quantity_keys=['kpoints-kpoints'])
            self._kpoints = NodeComposer.compose('array.kpoints', inputs)
        return self._kpoints

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
                kpt = Kpoint(point, weight=weights[index], logger=self._logger)
            else:
                # no weights supplied, so set them to 1.0
                kpt = Kpoint(point, weight=1.0, logger=self._logger)
            kpts.append(kpt)
        dictionary['points'] = kpts
        dictionary['mode'] = 'explicit'
        dictionary['num_kpoints'] = len(kpts)

        return dictionary

    @staticmethod
    def _get_kpointsdata_explicit(kpoints_dict):
        """Turn an 'explicit' kpoints dictionary into Aiida KpointsData."""
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

    @staticmethod
    def _get_kpointsdata_automatic(kpoints_dict):
        """Turn an 'automatic' kpoints dictionary into Aiida KpointsData."""
        kpout = get_data_class('array.kpoints')()

        mesh = kpoints_dict.get('divisions')
        shifts = kpoints_dict.get('shifts')
        kpout.set_kpoints_mesh(mesh, offset=shifts)

        return kpout

    @staticmethod
    def _get_kpointsdict_automatic(kpointsdata):
        """Turn Aiida KpointData into an 'automatic' kpoints dictionary."""
        dictionary = {}
        # automatic mode
        mesh = kpointsdata.get_kpoints_mesh()
        dictionary['divisions'] = mesh[0]
        dictionary['shifts'] = mesh[1]
        dictionary['mode'] = 'automatic'
        # here we need to make a choice, so should
        # add more to Aiida to make this better
        # defined
        dictionary['centering'] = 'Gamma'
        dictionary['num_kpoints'] = 0

        return dictionary
