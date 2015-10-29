#encoding: utf-8
from aiida.orm.calculation.job.vasp import VaspCalculation
from aiida.orm import DataFactory
from aiida.parsers.parser import Parser
from aiida.parsers.exceptions import OutputParsingError, InvalidOperation
from aiida.common.datastructures import calc_states
from pymatgen.io import vasp

__copyright__ = u'Copyright (c), 2015, Rico Häuselmann'

class VaspParser(Parser):
    '''
    Parses the output files of a Vasp run.
    '''
    def __init__(self, calculation):
        '''
        initializes with a Calculation object
        '''
        super(VaspParser, self).__init__(calculation)
        if not isinstance(calculation, VaspCalculation):
            raise OutputParsingError('Input calculation must be a VaspCalculation')
        self._calc = calculation

    def parse_with_retrieved(self, retrieved):
        '''
        parses the datafolder into results that then get stored in AiiDA.
        '''
        # check status
        state = self._calc.get_state()
        if state != calc_states.PARSING:
            raise InvalidOperation("Calculation not in {} state".format(calc_states.PARSING))
        # get the retrieved folder if exists
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
        except KeyError:
            self.logger.error('No retrieved folder found')
            return False, ()
        # ls retrieved
        list_of_files = out_folder.get_folder_list()
        # check for OUTCAR presence
        if not self._calc._OUTPUT_FILE_NAME in list_of_files:
            self.logger.error('OUTCAR not found - probably something went wrong while reading the input files!')
            return False, ()
        new_nodes_list = []
        makepar = lambda f, p : self._make_param_data(out_folder, f, p)
        new_nodes_list.append(makepar('CONTCAR', vasp.Poscar.from_file))
        new_nodes_list.append(makepar('OSZICAR', vasp.Oszicar))
        new_nodes_list.append(makepar('OUTCAR', vasp.Outcar))
        new_nodes_list.append(makepar('vasprun.xml', vasp.Vasprun, lname='vasprun'))
        # TODO as single files: chg, chgcar (additional), doscar, eigenval, pcdat,
        # procar, wavecar, xdatcar

    def _make_param_data(self, out_folder, filename, pmgparser, lname=None):
        ParameterData = DataFactory('parameter')
        list_of_files = out_folder.get_folder_list()
        getf = lambda f : out_folder.get_abs_path(f)
        if not self._calc.filename in list_of_files:
            self.logger.error('{} not found')
        outf = pmgparser(getf(filename))
        outf_ln = lname or filename.lower()
        outf_dt = ParameterData(dict=outf.as_dict())
        return (outf_ln, outf_dt)
