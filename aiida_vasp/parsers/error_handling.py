"""
Module for gather errors from VASP calculation.

Here we define the datastructures for storing and definiting
errors.
"""
import os
import re
from collections import namedtuple
from pathlib import Path

import yaml

from aiida.common.links import LinkType
from aiida.orm import CalcJobNode, WorkChainNode


class VaspError:
    """Class representing the errors given by VASP"""

    def __init__(self, regex, shortname, message, location='STDOUT', critical=True, suggestion=None):  # pylint: disable=too-many-arguments
        """
        Initialise a VaspError object.

        :param regex: Regex used for scanning
        :param message: Message to the user
        :param critical: Whether the error is critical
        """
        if isinstance(regex, str):
            self.regex = re.compile(regex)
        else:
            self.regex = regex
        self.message = message
        self.location = location
        self.shortname = shortname
        self.critical = critical
        self.suggestion = suggestion

    def __repr__(self):
        return f'VaspError(re={self.regex}, message={self.message}, critical={self.critical})'

    def error_report(self):
        """Generate the report for the error"""
        return f'{self.shortname}: {self.message}'

    def check_line(self, line):
        """Check the error in the line, return True the error is found"""
        mch = self.regex.search(line)
        if mch:
            return ErrorRecord(self, self.shortname, self.message, self.critical, self.suggestion)
        return None


class VaspGeneralError(VaspError):
    """Class representing an general VASP error"""

    def __init__(self):
        regex = re.compile(r'(FAIL.*$)|(WARN.*$)|(BAD.*)', re.IGNORECASE)
        super(VaspGeneralError, self).__init__(regex=regex, shortname='general', message=None, critical=False)

    def check_line(self, line):
        """Check error in a single line"""
        mch = self.regex.search(line)
        if mch:
            # Record the group as the error message
            message = mch.group(0)
            return ErrorRecord(self, self.shortname, message, self.critical, self.suggestion)
        return None


# Record instance for an error
ErrorRecord = namedtuple('ErrorRecord', ['error', 'shortname', 'message', 'critical', 'suggestion'])


class ErrorScanner:
    """Scanner for VASP erros"""

    ERRORS = [
        VaspError('LAPACK: Routine ZPOTRF', shortname='zpotrf', message='Error in Lapack ZPOTRF', location='STDOUT', critical=True),
        VaspError('WARNING in EDDRMM', shortname='rmm-diis', message='RMM-DIIS error', location='STDOUT', critical=True),
        VaspError('SBESSELITER: nicht', shortname='nicht', message='Error, try set LREAL = False', location='STDOUT', critical=True),
        VaspError('internal error in subroutine SGRCON:',
                  shortname='sgrcon',
                  message='Symmetry related Error',
                  location='STDOUT',
                  critical=True),
        VaspGeneralError(),
    ]
    DEFAULT_ERROR_FILE = 'vasp_errors.yaml'

    @classmethod
    def load_from_yml(cls, fname=None):
        """Read the errors from yml and save it as the class method"""
        if fname is None:
            fname = Path(__file__).parent / cls.DEFAULT_ERROR_FILE

        with open(fname) as fhd:
            data = yaml.safe_load(fhd)
        for entry in data:
            cls.ERRORS.append(VaspError(**entry))

        # Remove any general errors, and append one at the end
        generals = []
        for error in cls.ERRORS:
            if isinstance(error, VaspGeneralError):
                generals.append(error)
        for error in generals:
            cls.ERRORS.remove(error)
        cls.ERRORS.append(VaspGeneralError())

    def __init__(self, fhandle, ftype='STDOUT'):
        """
        Initialise the scanner with the file handle

        :param fhandle: A iteratble object for getting the lines
        :param ftype: Type of the file handle - it is STDOUT or OUTCAR?
        """
        self.fhandle = fhandle
        self.ftype = ftype
        self.errors_found = set()
        self.has_scanned = False
        self.scan()

    def scan(self):
        """Scan for errors in the stream"""
        for line in self.fhandle:
            for error in self.ERRORS:
                # Select only the relavent errors
                if error.location != self.ftype:
                    continue
                # Check line method returns a ErrorRecord object
                error_rec = error.check_line(line)
                if error_rec:
                    self.errors_found.add(error_rec)
                    # Once the error is found, no point to continue
                    break

        self.has_scanned = True

    def get_errors(self):
        """Get the errors"""
        return self.errors_found

    def get_lists(self):
        return [list(named) for named in self.errors_found]

    @property
    def has_critical(self):
        """Whether there are critical errors"""
        return any([error.critical for error in self.errors_found])

    @property
    def has_error(self):
        """Whether there are critical errors"""
        return bool(self.errors_found)

    @property
    def n_erros(self):
        """The number of errors found"""
        return len(self.errors_found)

    def __repr__(self):
        return f'<ErrorScanner with {len(self.errors_found)} found errors>'


# This loads the error defined in the vasp_errors.yaml placed in the same directory
ErrorScanner.load_from_yml(os.path.join(os.path.split(__file__)[0], 'vasp_errors.yaml'))


class NodeErrorScanner(ErrorScanner):
    """Scann Error for a Node"""

    def __init__(self, node):
        """Initialise an NodeErrorScanner"""

        if isinstance(node, CalcJobNode):
            calcjob = node
        elif isinstance(node, WorkChainNode):
            misc = node.outputs.misc
            calcjob = misc.get_incoming(link_type=LinkType.CREATE).one().node

        stdout_name = calcjob.get_option('scheduler_stdout')
        retrieved = calcjob.outputs.retrieved
        with retrieved.open(stdout_name) as fhandle:
            lines = fhandle.readlines()

        super(NodeErrorScanner, self).__init__(lines)
