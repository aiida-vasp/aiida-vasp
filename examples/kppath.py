from aiida.orm.calculation.job.vasp.base import BasicCalculation, Input


class KppathCalculation(BasicCalculation):
    kpoints = Input(types='parameter')
