"""Utilities to modify wannier parameters in a proveniency-trackable way"""
from aiida.orm.calculation.inline import make_inline
from aiida.orm import DataFactory


@make_inline
def modify_wannier_parameters_inline(original, modifications):  # pylint: disable=invalid-name
    """
    An InlineCalculation for modifying wannier parameters ('.win' file).

    :key ParameterData original: base parameters, can be overridden
    :key ParameterData modifications: additional parameters and overrides

    if original comes from a VASP2WANNIER setup and num_wann is overriden,
    num_bands will automatically be set accordingly.
    No consistency checks are performed.
    """
    result = DataFactory('parameter')()
    orig_dict = original.get_dict()
    mod_dict = modifications.get_dict()

    set_num_bands = not bool(orig_dict.get('num_bands'))
    set_num_bands &= bool(mod_dict.get('num_wann'))
    set_num_bands &= bool(orig_dict.get('num_wann'))
    if set_num_bands:
        orig_dict['num_bands'] = orig_dict['num_wann']

    result.set_dict(orig_dict)
    result.update_dict(mod_dict)
    return {'wannier_parameters': result}
