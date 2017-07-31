from vasp import VaspParser


class BasicParser(VaspParser):
    def parse_with_retrieved(self, retrieved):
        self.check_state()
        self.out_folder = self.get_folder(retrieved)
        if not self.out_folder:
            return self.result(success=False)
        outcar = self.get_file('OUTCAR')
        if not outcar:
            self.logger.error(
                'OUTCAR not found, ' +
                'look at the scheduler output for troubleshooting')
            return self.result(success=False)
        return self.result(success=True)
