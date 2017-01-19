#encoding: utf-8
from aiida_vasp.calcs.vasp import VaspCalculation
from aiida.orm import DataFactory
from aiida.parsers.parser import Parser
from aiida.parsers.exceptions import OutputParsingError
from aiida.common.exceptions import InvalidOperation
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
        def makepar(f, p, lname=None):
            return self._make_param_data(out_folder, f, p, lname)

        new_nodes_list.append(makepar('CONTCAR', vasp.Poscar.from_file))
        new_nodes_list.append(makepar('OSZICAR', vasp.Oszicar))
        new_nodes_list.append(makepar('OUTCAR', vasp.Outcar))
        new_nodes_list.append(makepar('vasprun.xml', vasp.Vasprun, lname='vasprun'))
        # TODO as single files: chg, chgcar (additional), doscar, eigenval, pcdat,
        # procar, wavecar, xdatcar
        def makesf(f, lname=None):
            return self._make_single_file(out_folder, f, lname)

        new_nodes_list.append(makesf('CHG'))
        new_nodes_list.append(makesf('CHGCAR'))
        new_nodes_list.append(makesf('DOSCAR'))
        new_nodes_list.append(makesf('EIGENVAL'))
        new_nodes_list.append(makesf('PCDAT'))
        new_nodes_list.append(makesf('PROCAR'))
        new_nodes_list.append(makesf('WAVECAR'))
        new_nodes_list.append(makesf('XDATCAR'))

        logger.error(new_nodes_list)

        return True, new_nodes_list

    def _make_param_data(self, out_folder, filename, pmgparser, lname=None):
        ParameterData = DataFactory('parameter')
        list_of_files = out_folder.get_folder_list()
        if not filename in list_of_files:
            self.logger.error('{} not found')
        outf = pmgparser(out_folder.get_abs_path(filename))
        outf_ln = lname or filename.lower()
        outf_dt = ParameterData(dict=outf.as_dict())
        return (outf_ln, outf_dt)

    def _make_single_file(self, out_folder, filename, lname=None):
        SinglefileData = DataFactory('singlefile')
        outf_n = out_folder.get_abs_path(filename)
        outf_ln = lname or filename.lower()
        outf_dt = SinglefileData(file=outf_n)
        return (outf_ln, outf_dt)
