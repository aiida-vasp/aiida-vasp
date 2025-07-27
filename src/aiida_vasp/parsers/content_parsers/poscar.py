"""
The ``POSCAR``/``CONTCAR`` parser interface.

Contains the parsing interfaces to ``parsevasp`` used to parse ``POSCAR``/``CONTCAR`` content.
"""

from typing import Any, TextIO

import numpy as np
from aiida import orm
from aiida.common.constants import elements
from parsevasp.poscar import Poscar, Site

from aiida_vasp.parsers.content_parsers.base import BaseFileParser


class PoscarParser(BaseFileParser):
    """The parser interface that enables parsing of ``POSCAR``/``CONTCAR`` content.

    The parser is triggered by using the ``poscar-structure`` quantity key. The quantity key ``structure``
    will on the other hand parse the structure using the XML parser.

    :param precision: An integer specifying the number of digits for floating point numbers that will be written
                      to ``POSCAR``/``CONTCAR``. Defaults to 12.
    :type precision: int, optional
    """

    DEFAULT_SETTINGS = {'quantities_to_parse': ['structure']}

    PARSABLE_QUANTITIES = {
        'structure': {
            'inputs': [],
            'name': 'structure',
            'prerequisites': [],
        },
    }

    def __init__(self, *, precision=12, **kwargs):
        """Initialize an instance of this class."""

        self._precision = precision
        super().__init__(**kwargs)

    def _init_from_handler(self, handler: TextIO) -> None:
        """Initialize using a file like handler.

        :param handler: A file like object that provides the necessary ``POSCAR``/``CONTCAR`` content to be parsed.
        :type handler: file-like object
        """

        try:
            self._content_parser = Poscar(
                file_handler=handler, prec=self._precision, conserve_order=True, logger=self._logger
            )
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    def _init_from_data(self, data: orm.StructureData) -> None:
        """Initialize using AiiDA ``StructureData`` instance.

        :param data: A valid AiiDA ``StructureData`` object.
        :type data: object
        """

        if isinstance(data, orm.StructureData):
            self._content_data = data
        else:
            raise TypeError('The supplied AiiDA data structure is not a StructureData.')

    @property
    def structure(self) -> dict[str, Any]:
        """Return structure from POSCAR.

        :returns: A dict that contain keys ``comment``, ``unitcell`` and ``sites``, which are compatible
                  with consumption of the initialization of the AiiDA ``StructureData``.
        :rtype: dict
        """

        aiida_structure = parsevasp_to_aiida(self._content_parser)

        return aiida_structure

    def _content_data_to_content_parser(self) -> 'PoscarParser':
        """Convert an AiiDA ``StructureData`` to a content parser instance of ``Poscar`` from ``parsevasp``.

        :returns: An instance of ``Poscar`` from ``parsevasp``.
        :rtype: object
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

    def transform_to_bool(self, value: str | int) -> bool:
        """Helper function to transform the dictionary from strings or integers to bools"""
        if value in [0, 'F', 'f']:
            return False
        if value in [1, 'T', 't']:
            return True
        return True


def parsevasp_to_aiida(poscar: Poscar) -> dict[str, Any]:
    """``parsevasp`` to AiiDA conversion.

    Generate an AiiDA structure that can be consumed by ``StructureData`` on initialization
    from the ``parsevasp`` instance of the ``Poscar`` class.

    :param poscar: An instance of the ``Poscar`` class in ``parsevasp``.
    :type poscar: object
    :returns: A dictionary representation which are ready to be used when creating an
              AiiDA ``StructureData`` instance.
    :rtype: dict
    """

    # Fetch a dictionary containing the entries, make sure all coordinates are
    # cartesian
    poscar_dict = poscar.get_dict(direct=False)
    # Inverted dictionary with element names being the keys and numbers being the values
    symbols = {value['symbol']: key for key, value in elements.items()}

    for site in poscar_dict['sites']:
        specie = site['specie']
        # User can specify whatever they want for the elements, but
        # the symbols entries in AiiDA only support the entries defined
        # in aiida.common.constants.elements{}

        # Strip trailing _ in case user specifies potential
        symbol = specie.split('_')[0].capitalize()
        # Check if leading entry is part of
        # aiida.common.constants.elements{}, otherwise set to X
        if symbol not in symbols:
            symbol = 'X'

        site['symbol'] = symbol
        site['kind_name'] = specie

    return poscar_dict
