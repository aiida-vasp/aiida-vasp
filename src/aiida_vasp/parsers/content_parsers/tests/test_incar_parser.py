"""Test the INCAR parser."""
# pylint: disable=redefined-outer-name, unused-wildcard-import, unused-argument, wildcard-import

import pytest
from aiida_vasp.parsers.content_parsers.incar import IncarParser
from aiida_vasp.utils.aiida_utils import get_data_node
from aiida_vasp.utils.fixtures import *

compare_incar = {'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'magmom': '30 * 2*0.', 'sigma': 0.5}


@pytest.mark.parametrize(['incar_parser'], [('incar',)], indirect=True)
def test_parse_incar(incar_parser):
    """Load a reference INCAR parser.

    We check that it parses and provides the correct content for the default INCAR.

    """

    # The structure for the INCAR parser should have the key `incar-structure`
    result = incar_parser.get_quantity('incar')
    assert result == compare_incar


@pytest.mark.parametrize(['incar_parser'], [('phonondb',)], indirect=True)
def test_parse_incar_phonon(incar_parser):
    """Load a reference INCAR parser.

    We check that it parses and provides the correct content for an INCAR used for
    phonon calculations.

    """

    incar = incar_parser.incar
    assert incar['prec'] == 'Accurate'
    assert incar['ibrion'] == -1
    assert incar['encut'] == pytest.approx(359.7399)
    assert incar['lreal'] is False


@pytest.mark.parametrize(['incar_parser'], [('incar',)], indirect=True)
def test_parse_incar_write(incar_parser, tmpdir):
    """Load a reference INCAR parser and check that the write functionality works.

    Here we make sure the write function of the content parser works.

    """
    # Write the loaded structure to file
    temp_path = str(tmpdir.join('INCAR'))
    incar_parser.write(temp_path)

    # Load the written structure using a new content parser instance
    content = None
    with open(temp_path, 'r', encoding='utf8') as handler:
        content = handler.readlines()
    ref_content = ['GGA = PE\n', 'GGA_COMPAT = .FALSE.\n', 'LORBIT = 11\n', 'MAGMOM = 30 * 2*0.\n', 'SIGMA = 0.5\n']
    assert content == ref_content


def test_parse_incar_data(vasp_params, tmpdir):
    """Load a reference AiiDA Dict and check that the parser can
    initialize using the data.

    Using the Dict sitting in the initialized parser it should
    write that content to an INCAR file when issuing write which is also tested,
    file is reloaded and content checked.

    """

    # Initialize parser with an existing reference Dict
    incar_parser = IncarParser(data=vasp_params)

    # Check that get_quantity return the same Dict instance
    assert vasp_params == incar_parser.get_quantity('key_does_not_matter')

    # Write the loaded Dict to file, which behind the scenes convert it
    # to a INCAR format
    temp_path = str(tmpdir.join('INCAR'))
    incar_parser.write(temp_path)

    # Load the written INCAR using a new content parser instance
    parser = None
    with open(temp_path, 'r', encoding='utf8') as handler:
        parser = IncarParser(handler=handler)
    result = parser.get_quantity('incar')
    assert vasp_params.get_dict() == result


@pytest.mark.parametrize(['incar_parser'], [(['incar', 'INCAR.nohash'],)], indirect=True)
def test_parse_incar_nohash(incar_parser):
    """Load a reference INCAR parser.

    Using parsevasp. Returned content should be None
    since parsevasp refuse to parse an INCAR where the
    comments does not start with hashtags.

    """

    result = incar_parser.incar
    assert result is None


def test_parse_incar_invalid_tag(vasp_params, tmpdir):
    """Test parsing an INCAR with an invalid tag."""
    params = vasp_params.get_dict()
    params.update(foo='bar')
    vasp_params_modified = get_data_node('core.dict', dict=params)
    parser = IncarParser(data=vasp_params_modified)
    temp_path = str(tmpdir.join('INCAR'))
    with pytest.raises(SystemExit):
        parser.write(temp_path)
