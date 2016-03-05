from aiida.parsers.parser import Parser
from aiida.common.datastructures import calc_states as cstat
from aiida.common.exceptions import InvalidOperation


class BaseParser(Parser):
    '''
    Does the boring stuff you need anyway.
    '''
    def __init__(self, calc):
        self._new_nodes = {}
        super(BaseParser, self).__init__(calc)

    def check_state(self):
        if self._calc.get_state() != cstat.PARSING:
            self.logger.error('Calculation not in parsing state')
            # ~ raise InvalidOperation('Calculation not in parsing state')

    def get_folder(self, retrieved):
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
            return out_folder
        except KeyError:
            self.logger.error('No retrieved folder found')
            return None

    def result(self, success):
        if success:
            return True, self._new_nodes.items()
        else:
            return False, self._new_nodes.items()

    def get_file(self, fname):
        try:
            ofname = self.out_folder.get_abs_path(fname)
            return ofname
        except OSError:
            self.logger.warning(fname+' not found in retrieved')
            return None

    def add_node(self, linkname, node):
        self._new_nodes[linkname] = node

    @property
    def new_nodes(self):
        return self._new_nodes.items()
