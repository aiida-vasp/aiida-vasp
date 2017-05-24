from aiida.orm.calculation.job.vasp.wannier import WannierCalculation

class WswannierCalculation(WannierCalculation):
    default_parser = 'vasp.wswannier'
