"""Common code for parsers"""
from aiida.parsers.parser import Parser
from aiida.common.datastructures import calc_states as cstat


class BaseParser(Parser):
    """
    Does common tasks all parsers carry out and provides
    convenience methods.
    """

    def __init__(self, calc):
        self._new_nodes = {}
        super(BaseParser, self).__init__(calc)
        self.out_folder = None

    def parse_with_retrieved(self, retrieved):
        """
        sets the out_folder attribute for use in extending parsers.
        """
        # ~ super(BaseParser, self).parse_with_retrieved(retrieved)
        self.out_folder = self.get_folder(retrieved)
        return self.result(success=bool(self.out_folder is not None))

    def check_state(self):
        """
        log an error if the calculation being parsed is not in PARSING state
        """
        if self._calc.get_state() != cstat.PARSING:
            self.logger.error('Calculation not in parsing state')
            # ~ raise InvalidOperation('Calculation not in parsing state')

    def get_folder(self, retrieved):
        """convenient access to the retrieved folder"""
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
            return out_folder
        except KeyError:
            self.logger.error('No retrieved folder found')
            return None

    def result(self, success):
        """
        returns a success flag as well as the new output nodes added to
        the parser's internal list of output nodes during parsing.

        :param bool success: wether the parsing was successful
        :return: (sucess, new_nodes), the expected return values of a parser
        """
        return bool(success), self._new_nodes.items()

    def get_file(self, fname):
        """
        conveniend access to retrieved files
        :param fname: name of the file
        :return: absolute path to the retrieved file
        """
        try:
            ofname = self.out_folder.get_abs_path(fname)
            return ofname
        except OSError:
            self.logger.warning(fname + ' not found in retrieved')
            return None

    def add_node(self, linkname, node):
        """add a node to the internal list of output nodes"""
        self._new_nodes[linkname] = node

    @property
    def new_nodes(self):
        """holds a list of parsed output nodes"""
        return self._new_nodes.items()
