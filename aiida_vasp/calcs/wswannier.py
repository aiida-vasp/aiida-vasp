# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""Wannier - Calculation: uses wswannier parser by default"""
from .wannier import WannierCalculation


class WswannierCalculation(WannierCalculation):
    """Same as WannierCalculation but uses vasp.wswannier as default parser"""
    default_parser = 'vasp.wswannier'
