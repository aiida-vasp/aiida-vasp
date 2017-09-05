from .wannier import WannierCalculation


class WswannierCalculation(WannierCalculation):
    default_parser = 'vasp.wswannier'
