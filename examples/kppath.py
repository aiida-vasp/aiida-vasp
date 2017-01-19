from aiida_vasp.calcs.base import BasicCalculation, Input


class KppathCalculation(BasicCalculation):
    kpoints = Input(types='parameter')
