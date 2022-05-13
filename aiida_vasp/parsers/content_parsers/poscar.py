"""
The ``POSCAR``/``CONTCAR`` parser interface.

------------------------------------
Contains the parsing interfaces to ``parsevasp`` used to parse ``POSCAR``/``CONTCAR`` content.
"""
# pylint: disable=no-self-use
import numpy as np

from aiida.common.constants import elements
from parsevasp.poscar import Poscar, Site
from aiida_vasp.parsers.content_parsers.base import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class


class PoscarParser(BaseFileParser):
    """The parser interface that enables parsing of ``POSCAR``/``CONTCAR`` content.

    The parser is triggered by using the ``poscar-structure`` quantity key. The quantity key ``structure``
    will on the other hand parse the structure using the XML parser.

    Parameters
    ----------
    precision : int, optional
        An integer specifying the number of digits for floating point numbers that will be written
        to ``POSCAR``/``CONTCAR``. Defaults to 12.

    """

    DEFAULT_SETTINGS = {'quantities_to_parse': ['poscar-structure']}

    PARSABLE_QUANTITIES = {
        'poscar-structure': {
            'inputs': [],
            'name': 'structure',
            'prerequisites': [],
        },
    }

    def __init__(self, *, precision=12, **kwargs):
        """Initialize an instance of this class."""

        self._precision = precision
        super().__init__(**kwargs)

    def _init_from_handler(self, handler):
        """Initialize a ``parsevasp`` object of ``Poscar`` using a file like handler.

        Parameters
        ----------
        handler : object
            A file like object that provides the necessary ``POSCAR`/``CONTCAR``` content to be parsed.

        """

        try:
            self._content_parser = Poscar(file_handler=handler, prec=self._precision, conserve_order=True, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    def _init_from_data(self, data):
        """Initialize using an AiiDA ``StructureData`` instance.

        Parameters
        ----------
        data : object
            A valid AiiDA ``StructureData`` object.

        """

        if isinstance(data, get_data_class('structure')):
            self._content_data = data
        else:
            raise TypeError('The supplied AiiDA data structure is not a StructureData.')

    @property
    def structure(self):
        """
        Return a structure that is ready to be consumed by the the AiiDA ``StructureData``.

        Returns
        -------
        aiida_structure : dict
            A dict that contain keys ``comment``, ``unitcell`` and ``sites``, which are compatible
            with consumption of the initialization of the AiiDA ``StructureData``.

        """

        aiida_structure = parsevasp_to_aiida(self._content_parser)

        return aiida_structure

    def _content_data_to_content_parser(self):
        """
        Convert an AiiDA ``StructureData`` to a content parser instance of ``Poscar`` from ``parsevasp``.

        Returns
        -------
        content_parser : object
            An instance of ``Poscar`` from ``parsevasp``.

        """
        dictionary = {}
        dictionary['comment'] = self._content_data.label or self._content_data.get_formula()
        dictionary['unitcell'] = np.asarray(self._content_data.cell)
        # As for now all Aiida-structures are in Cartesian coordinates.
        direct = False
        sites = []
        _transform_to_bool = np.vectorize(self.transform_to_bool)
        for index, site in enumerate(self._content_data.sites):
            if self._options is None:
                _selective = [True, True, True]
            else:
                try:
                    _selective = _transform_to_bool(np.array(self._options['positions_dof'])[index, :])
                except KeyError:
                    _selective = [True, True, True]
            sites.append(Site(site.kind_name, site.position, selective=_selective, direct=direct))

        dictionary['sites'] = sites

        # We brake hard if ``parsevasp`` fail here. If we can not write we will not try another parser.
        content_parser = Poscar(poscar_dict=dictionary, prec=self._precision, conserve_order=True, logger=self._logger)

        return content_parser

    def transform_to_bool(self, value):
        """Helper function to transform the dictionary from strings or integers to bools"""
        if value in [0, 'F', 'f']:
            return False
        if value in [1, 'T', 't']:
            return True
        return True


def parsevasp_to_aiida(poscar):
    """``parsevasp`` to AiiDA conversion.

    Generate an AiiDA structure that can be consumed by ``StructureData`` on initialization
    from the ``parsevasp`` instance of the ``Poscar`` class.

    Parameters
    ----------
    poscar : object
        An instance of the ``Poscar`` class in ``parsevasp``.

    Returns
    -------
    poscar_dict : dict
        A dictionary representation which are ready to be used when creating an
        AiiDA ``StructureData`` instance.

    """

    # Fetch a dictionary containing the entries, make sure all coordinates are
    # cartesian
    poscar_dict = poscar.get_dict(direct=False)

    for site in poscar_dict['sites']:
        specie = site['specie']
        # User can specify whatever they want for the elements, but
        # the symbols entries in AiiDA only support the entries defined
        # in aiida.common.constants.elements{}

        # Strip trailing _ in case user specifies potential
        symbol = specie.split('_')[0].capitalize()
        # Check if leading entry is part of
        # aiida.common.constants.elements{}, otherwise set to X, but first
        # invert
        symbols = fetch_symbols_from_elements(elements)
        try:
            symbols[symbol]
        except KeyError:
            symbol = 'X'

        site['symbol'] = symbol
        site['kind_name'] = specie

    return poscar_dict


def fetch_symbols_from_elements(elmnts):
    """Fetch the symbol entry in the elements dictionary in Aiida."""

    new_dict = {}
    for key, value in elmnts.items():
        new_dict[value['symbol']] = key
    return new_dict
