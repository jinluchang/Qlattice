import cqlat as c

from qlat.rng import *

from qlat.geo import *

class PointSelection:

    def __init__(self, coordinate_list = None):
        if None == coordinate_list:
            self.cdata = c.mk_psel()
        else:
            self.cdata = c.mk_psel(coordinate_list)

    def __del__(self):
        c.free_psel(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, PointSelection)
        c.set_psel(self, v1)
        return self

    def save(self, path):
        c.save_psel(self, path)

    def load(self, path):
        c.load_psel(self, path)

    def list(self):
        return c.mk_list_psel(self)

class FieldSelection:

    def __init__(self, total_site = None, n_per_tslice = -1, rng = None, psel = None):
        if None == total_site:
            self.cdata = c.mk_fsel()
        else:
            assert isinstance(n_per_tslice, int)
            assert isinstance(rng, RngState)
            if None == psel:
                self.cdata = c.mk_fsel(total_site, n_per_tslice, rng)
            else:
                assert isinstance(psel, PointSelection)
                self.cdata = c.mk_fsel(total_site, n_per_tslice, rng, psel)
            self.update()
            self.update(n_per_tslice)

    def __del__(self):
        c.free_fsel(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, FieldSelection)
        c.set_fsel(self, v1)
        return self

    def add_psel(self, psel):
        c.add_psel_fsel(self, psel)
        self.update()

    def update(self, n_per_tslice = -1):
        c.update_fsel(self, n_per_tslice)

    def save(self, path):
        return c.save_fsel(self, path)

    def load(self, path, n_per_tslice):
        return c.load_fsel(self, path, n_per_tslice)

    def geo(self):
        geo = Geometry((0, 0, 0, 0))
        c.set_geo_fsel(geo, self)
        return geo

    def total_site(self):
        return c.get_total_site_fsel(self)

    def n_elems(self):
        return c.get_n_elems_fsel(self)

    def n_per_tslice(self):
        return c.get_n_per_tslice_fsel(self)

    def prob(self):
        return c.get_prob_fsel(self)
