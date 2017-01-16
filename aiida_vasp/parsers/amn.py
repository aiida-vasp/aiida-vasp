from vasp2w90 import Vasp2w90Parser


class AmnParser(Vasp2w90Parser):
    def parse_with_retrieved(self, retrieved):
        super(AmnParser, self).parse_with_retrieved(retrieved)

        self.set_wdat(self.get_wdat_node())
        return self.result(success=True)

    def get_wdat_node(self):
        wdatnode = self._calc.new_wannier_data()
        for ext in ['mmn', 'amn', 'eig']:
            wfile = self.get_file('wannier90.'+ext)
            if wfile:
                wdatnode.add_file(wfile)
        return wdatnode
