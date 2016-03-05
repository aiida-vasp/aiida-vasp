from vasp2w90 import Vasp2w90Parser


class NscfParser(Vasp2w90Parser):
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

        self.vrp = self.read_run()
        self.dcp = self.read_dos()
        self.set_bands(self.read_eigenval()[0])
        self.set_dos(self.get_dos_node(self.vrp, self.dcp))
        self.set_win(self.get_win_node())
        self.set_wdat(self.get_wdat_node())
        self.add_node('results', self.get_output)
        return self.result(success=True)
