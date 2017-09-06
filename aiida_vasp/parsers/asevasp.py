#encoding: utf-8
"""AiiDA Parser plugin for aiida_vasp.AsevaspCalculation"""
from pymatgen.io import vasp
from aiida.orm import DataFactory
from aiida.parsers.parser import Parser
from aiida.parsers.exceptions import OutputParsingError
from aiida.common.exceptions import InvalidOperation
from aiida.common.datastructures import calc_states

from aiida_vasp.calcs.asevasp import AsevaspCalculation

__copyright__ = u'Copyright (c), 2015, Rico Häuselmann'


class AsevaspParser(Parser):
    """
    Parses the output files of a Vasp run.
    """

    def __init__(self, calculation):
        """
        initializes with a Calculation object
        """
        super(AsevaspParser, self).__init__(calculation)
        if not isinstance(calculation, AsevaspCalculation):
            raise OutputParsingError(
                'Input calculation must be a VaspCalculation')
        self._calc = calculation

    def parse_with_retrieved(self, retrieved):
        """
        parses the datafolder into results that then get stored in AiiDA.
        """
        # check status
        state = self._calc.get_state()
        if state != calc_states.PARSING:
            raise InvalidOperation(
                "Calculation not in {} state".format(calc_states.PARSING))
        # get the retrieved folder if exists
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]  # pylint: disable=protected-access
        except KeyError:
            self.logger.error('No retrieved folder found')
            return False, ()
        # ls retrieved
        list_of_files = out_folder.get_folder_list()
        # check for OUTCAR presence
        if self._calc._OUTPUT_FILE_NAME not in list_of_files:  # pylint: disable=protected-access
            self.logger.error(
                'OUTCAR not found - probably something went wrong while reading the input files!'
            )
            return False, ()
        new_nodes_list = []

        def makepar(filename, pmgparser, lname=None):
            return self._make_param_data(out_folder, filename, pmgparser,
                                         lname)

        new_nodes_list.append(makepar('CONTCAR', vasp.Poscar.from_file))
        new_nodes_list.append(makepar('OSZICAR', vasp.Oszicar))
        new_nodes_list.append(makepar('OUTCAR', vasp.Outcar))
        new_nodes_list.append(
            makepar('vasprun.xml', vasp.Vasprun, lname='vasprun'))

        def makesf(filename, lname=None):
            return self._make_single_file(out_folder, filename, lname)

        new_nodes_list.append(makesf('CHG'))
        new_nodes_list.append(makesf('CHGCAR'))
        new_nodes_list.append(makesf('DOSCAR'))
        new_nodes_list.append(makesf('EIGENVAL'))
        new_nodes_list.append(makesf('PCDAT'))
        new_nodes_list.append(makesf('PROCAR'))
        new_nodes_list.append(makesf('WAVECAR'))
        new_nodes_list.append(makesf('XDATCAR'))

        self.logger.debug(new_nodes_list)

        return True, new_nodes_list

    def _make_param_data(self, out_folder, filename, pmgparser, lname=None):
        """Create a parameter data node"""
        parameter_cls = DataFactory('parameter')
        list_of_files = out_folder.get_folder_list()
        if filename not in list_of_files:
            self.logger.error('{} not found')
        outf = pmgparser(out_folder.get_abs_path(filename))
        outf_ln = lname or filename.lower()
        outf_dt = parameter_cls(dict=outf.as_dict())
        return (outf_ln, outf_dt)

    @staticmethod
    def _make_single_file(out_folder, filename, lname=None):
        """Create a singlefile node"""
        singlefile_cls = DataFactory('singlefile')
        outf_n = out_folder.get_abs_path(filename)
        outf_ln = lname or filename.lower()
        outf_dt = singlefile_cls(file=outf_n)
        return (outf_ln, outf_dt)
