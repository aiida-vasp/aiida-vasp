import pytest
from aiida import orm

from aiida_vasp.utils.opthold import BandOptions, ConvOptions, OptionContainer, RelaxOptions


def test_option_container():
    """Test the base option container"""

    class TestHolder(OptionContainer):
        x: int = 1
        y: str = 'test'
        z: float

    holder = TestHolder(z=1.0)
    d = holder.aiida_dict()
    holder2 = TestHolder(**d.get_dict())
    assert holder == holder2
    assert holder.model_dump() == holder2.model_dump()

    assert holder.aiida_validate({'x': 2, 'y': 'test2', 'z': 2.0}) is None
    assert holder.aiida_validate({'x': 2, 'y': 2.0, 'z': '2'}) is not None

    node = holder.aiida_serialize({'x': 2, 'z': 3})
    assert node['x'] == 2
    assert node['y'] == 'test'
    assert node['z'] == 3


@pytest.mark.parametrize('cls', [RelaxOptions, BandOptions, ConvOptions])
def test_workchain_options_roundtrip(cls):
    """Test the relax option container"""
    opt = cls()
    out = opt.aiida_dict()
    assert cls(**out.get_dict()) == opt


@pytest.mark.parametrize('cls', [RelaxOptions, BandOptions, ConvOptions])
def test_workchain_options_validate(cls):
    """Test the relax option container"""
    opt = cls()
    out = opt.model_dump()
    assert cls.aiida_validate(out) is None


@pytest.mark.parametrize('cls', [RelaxOptions, BandOptions, ConvOptions])
def test_workchain_options_serialize(cls):
    """Test the relax option container"""
    opt = cls()
    out = opt.model_dump()
    assert isinstance(cls.aiida_serialize(out), orm.Dict)


@pytest.mark.parametrize('cls', [RelaxOptions, BandOptions, ConvOptions])
def test_workchain_options_description(cls):
    """Test the relax option container"""
    opt = cls()
    desc = opt.aiida_description()
    assert desc
