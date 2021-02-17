import cqlat as c

class PointSelection:

    def __init__(self, coordinate_list = None):
        if None == coordinate_list:
            self.cdata = c.mk_psel()
        else:
            self.cdata = c.mk_psel(coordinate_list)

    def __del__(self):
        c.free_psel(self)

    def save(self, path):
        c.save_psel(self, path)

    def load(self, path):
        c.load_psel(self, path)

class FieldSelection:

    def __init__(self, total_site = None, n_per_tslice = -1, rng = None, psel = None):
        if None == total_site:
            self.cdata = c.mk_fsel()
        else:
            assert type(n_per_tslice) = int
            assert isinstance(rng, RngState)
            if None == psel:
                self.cdata = c.mk_fsel(total_site, n_per_tslice, rng)
            else:
                assert isinstance(psel, PointSelection)
                self.cdata = c.mk_fsel(total_site, n_per_tslice, rng, psel)
            self.update()
            self.update(n_per_tslice)

    def __del__(self):
        c.free_psel(self)

    def add_psel(self, psel):
        c.add_psel_fsel(self, psel)
        self.update()

    def update(self, n_per_tslice = -1):
        c.upate_fsel(self, n_per_tslice)

    def save(self, path):
        c.save_fsel(self, path)

    def load(self, path, n_per_tslice):
        c.load_lsel(self, path, n_per_tslice)
