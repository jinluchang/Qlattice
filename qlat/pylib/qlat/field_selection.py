import qlat.cqlat as c

from qlat_utils import *

from qlat.geometry import *
from qlat.utils_io import *

class PointSelection:

    # self.geo
    # self.cdata

    def __init__(self, coordinate_list = None, geo = None):
        if None == coordinate_list:
            self.cdata = c.mk_psel()
        else:
            self.cdata = c.mk_psel(coordinate_list)
        self.geo = geo

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_psel(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, PointSelection)
        c.set_psel(self, v1)
        self.geo = v1.geo
        return self

    def copy(self):
        x = PointSelection()
        x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def set_rand(self, rs, total_site, n_points):
        c.set_rand_psel(self, rs, total_site, n_points)
        self.geo = Geometry(total_site)

    def save(self, path):
        mk_file_dirs_info(path)
        c.save_psel(self, path)

    def load(self, path, geo = None):
        c.load_psel(self, path)
        self.geo = geo

    def n_points(self):
        return c.get_n_points_psel(self)

    def to_list(self):
        return c.mk_list_psel(self)

    def from_list(self, coordinate_list, geo = None):
        c.set_list_psel(self, coordinate_list)
        self.geo = geo
        return self

    def coordinate_from_idx(self, idx):
        return c.get_coordinate_from_idx_psel(self, idx)

cache_point_selection = mk_cache("point_selection")

def get_psel_tslice(total_site, *, t_dir = 3):
    # if t_dir = 3, then [ [0,0,0,0,], [0,0,0,1,], ..., [0,0,0,total_site[3]-1],]
    # if t_dir = 2, then [ [0,0,0,0,], [0,0,1,0,], ..., [0,0,total_site[2]-1],0,]
    # need total_site to set the psel.geo property
    assert 0 <= t_dir and t_dir < 4
    assert isinstance(total_site, list)
    param_tuple = (tuple(total_site), t_dir,)
    if param_tuple not in cache_point_selection:
        psel = PointSelection()
        c.set_tslice_psel(psel, total_site[t_dir], t_dir)
        psel.geo = Geometry(total_site)
        cache_point_selection[param_tuple] = psel
    return cache_point_selection[param_tuple]

class FieldSelection:

    # self.cdata

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
        assert isinstance(self.cdata, int)
        c.free_fsel(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, FieldSelection)
        c.set_fsel(self, v1)
        return self

    def copy(self):
        x = FieldSelection()
        x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

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

    def add_psel(self, psel, rank_psel = 1024 * 1024 * 1024 * 1024 * 1024):
        # Add psel points to the selection, with the rank specified as rank_psel.
        # If the point is already selected with lower rank, the rank is unchanged.
        c.add_psel_fsel(self, psel, rank_psel)
        self.update()

    def update(self, n_per_tslice = -1):
        # if n_per_tslice < 0: only update various indices
        # if n_per_tslice >= 0: only update parameters (n_per_tslice and prob)
        c.update_fsel(self, n_per_tslice)

    def select_rank_range(self, rank_start = 0, rank_stop = -1):
        # return new fsel with selected points that
        # rank_start <= rank and (rank < rank_stop or rank_stop == -1)
        # Does NOT change the n_per_tslice parameter for the new fsel
        fsel = FieldSelection()
        c.select_rank_range_fsel(fsel, self, rank_start, rank_stop)
        fsel.update()
        fsel.update(self.n_per_tslice())
        return fsel

    def select_t_range(self, rank_start = 0, rank_stop = -1):
        # return new fsel with selected points that
        # t_start <= t and (t < t_stop or t_stop == -1)
        # rank_start <= rank < rank_stop (rank_stop = -1 implies unlimited)
        # Does NOT change the n_per_tslice parameter for the new fsel
        fsel = FieldSelection()
        c.select_rank_range_fsel(fsel, self, rank_start, rank_stop)
        fsel.update()
        fsel.update(self.n_per_tslice())
        return fsel

    def to_psel(self):
        psel = PointSelection(None, self.geo())
        c.set_psel_fsel(psel, self)
        return psel

    def to_psel_local(self):
        psel = PointSelection(None, self.geo())
        c.set_psel_fsel_local(psel, self)
        return psel

    def save(self, path):
        mk_file_dirs_info(path)
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
        # return fsel.prob
        # n_per_tslice / spatial_volume
        return c.get_prob_fsel(self)

    def idx_from_coordinate(xg):
        return c.get_idx_from_coordinate_fsel(self, xg)

    def coordinate_from_idx(idx):
        return c.get_coordinate_from_idx_fsel(self, idx)

###

def is_matching_fsel(fsel1, fsel2):
    return c.is_matching_fsel(fsel1, fsel2)
