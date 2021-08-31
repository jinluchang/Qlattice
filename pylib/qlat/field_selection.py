import cqlat as c

from qlat.rng_state import *

from qlat.geometry import *

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

    def copy(self):
        x = PointSelection()
        x @= self
        return x

    def set_rand(self, rs, total_site, n_points):
        c.set_rand_psel(self, rs, total_site, n_points)

    def save(self, path):
        c.save_psel(self, path)

    def load(self, path):
        c.load_psel(self, path)

    def list(self):
        return c.mk_list_psel(self)

class FieldSelection:

    def __init__(self, total_site = None, n_per_tslice = -1, rs = None, psel = None):
        self.cdata = c.mk_fsel()
        if total_site is not None:
            assert isinstance(rs, RngState)
            assert isinstance(n_per_tslice, int)
            c.set_rand_fsel(self, rs, total_site, n_per_tslice)
            if psel is not None:
                c.add_psel_fsel(self, psel)
            self.update()
            self.update(n_per_tslice)

    def __del__(self):
        c.free_fsel(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, FieldSelection)
        c.set_fsel(self, v1)
        return self

    def copy(self):
        x = FieldSelection()
        x @= self
        return x

    def set_uniform(self, total_site, val = 0):
        # default (val = 0) select every sites
        # val = -1 deselection everything
        c.set_uniform_fsel(self, total_site, val)

    def set_rand(self, rs, total_site, n_per_tslice):
        assert isinstance(rs, RngState)
        assert isinstance(n_per_tslice, int)
        c.set_rand_fsel(self, rs, total_site, n_per_tslice)
        self.update()
        self.update(n_per_tslice)

    def add_psel(self, psel):
        c.add_psel_fsel(self, psel)
        self.update()

    def update(self, n_per_tslice = -1):
        # if n_per_tslice < 0: only update various indices
        # if n_per_tslice >= 0: only update parameters
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
