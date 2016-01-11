from aiida.parsers.parser import Parser
from aiida.common.datastructures import calc_states as cstat


class BaseParser(Parser):
    '''
    Does the boring stuff you need anyway.
    '''
    def check_state(self):
        if self._calc.get_state() != cstat.PARSING:
            raise InvalidOperation('Calculation not in parsing state')

    def get_folder(self, retrieved):
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
            return out_folder
        except KeyError as e:
            self.logger.error('No retrieved folder found')
            return None

    def result(self, success, new_nodes=()):
        if success:
            return True, new_nodes
        else:
            return False, new_nodes

    def get_file(self, fname):
        try:
            ofname = self.out_folder.get_abs_path(fname)
        except OSError:
            self.logger.warning(fname+' not found in retrieved')
            return None
