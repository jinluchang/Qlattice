# Qlat python bindings

[TOC]

## Geometry

```python
class Geometry:
    # self.cdata
    def __init__(self, total_site = None, multiplicity = None):
    def __imatmul__(self, v1):
    def copy(self, is_copying_data = True):
    def total_site(self):
    def total_volume(self):
    def local_volume(self):
    def multiplicity(self):
    def node_site(self):
    def eo(self):
    def expansion_left(self):
    def expansion_right(self):
    def id_node(self):
    def num_node(self):
    def coor_node(self):
    def size_node(self):
    def show(self):
    def show_all(self):
    def coordinate_g_from_l(self, xl):
    def coordinate_l_from_g(self, xg):
    def is_local(self, xl):
    def is_local_xg(self, xg):
    def xg_list(self):

def geo_reform(geo, multiplicity = 1, expansion_left = None, expansion_right = None):
def geo_eo(geo, eo = 0):
```

## LatData

```python
class LatData:
    # self.cdata
    def __init__(self):
    def __imatmul__(self, v1):
    def copy(self, is_copying_data = True):
    def __copy__(self):
    def __deepcopy__(self, memo):
    def save_node(self, path):
    def load_node(self, path):
    def bcast(self):
    def glb_sum(self):
    def save(self, path):
    def load(self, path):
    def show(self):
    def __iadd__(self, ld1):
    def __isub__(self, ld1):
    def __imul__(self, factor):
    def __add__(self, ld1):
    def __sub__(self, ld1):
    def __mul__(self, factor):
    __rmul__ = __mul__
    def __neg__(self):
    def __pos__(self):
    def set_zero(self):
    def qnorm(self):
    def is_match(self, ld1):
        # ld.info needs to be exactly equal
    def is_complex(self):
    def ndim(self, *, is_complex = True):
    def dim_sizes(self, *, is_complex = True):
    def dim_name(self, dim):
    def dim_size(self, dim):
    def dim_indices(self, dim):
    def set_dim_sizes(self, dim_sizes, *, is_complex = True):
    def set_dim_name(self, dim, name, indices = None):
    def dim_names(self, *, is_complex = True):
        # by default, return list can be used as the input argument for ld.from_numpy
    def to_list(self, *, is_complex = True):
    def from_list(self, val, *, is_complex = True):
    def to_numpy(self, *, is_complex = True):
    def from_numpy(self, val, dim_names = None, *, is_complex = True):
        # only set LatData shape if it is initially empty
        # otherwise only set data and ignore shape completely
        # dim_names should be a list of names for each dimension
    def __setitem__(self, idx, val):
        # use list with correct length as val
        # idx should be tuple or list of int
    def __getitem__(self, idx):
        # return a new list every call
        # idx should be tuple or list of int
    def __getstate__(self):
    def __setstate__(self, state):
    def info(self, dim = None, *, is_complex = True):
        # by default, return list can be used as the input argument for ld.set_info or mk_lat_data
    def set_info(self, info_list, *, is_complex = True):
        # info_list format:
        # [ [ dim_name, dim_size, dim_indices, ], ... ]
        # dim_indices can be optional

def mk_lat_data(info_list, *, is_complex = True):
def load_lat_data(path):
```

## RngState

```python
class RngState:
	# self.cdata
    def __init__(self, x = None, y = None):
    def __imatmul__(self, v1):
    def copy(self, is_copying_data = True):
    def split(self, seed):
    def rand_gen(self):
    def u_rand_gen(self, upper = 1.0, lower = 0.0):
    def g_rand_gen(self, center = 0.0, sigma = 1.0):
    def c_rand_gen(self, size):
        # size can be total_site of the lattice
    def select(self, l):

rng_state_root = RngState()
```

## PointSelection

```python
class PointSelection:
    # self.geo
    # self.cdata
    def __init__(self, coordinate_list = None, geo = None):
    def __imatmul__(self, v1):
    def copy(self):
    def set_rand(self, rs, total_site, n_points):
    def save(self, path):
    def load(self, path, geo = None):
    def to_list(self):
    def from_list(self, coordinate_list, geo = None):
    def coordinate_from_idx(self, idx):
        return c.get_coordinate_from_idx_psel(self, idx)

def get_psel_tslice(total_site):
    # [ [0,0,0,0], [0,0,0,1], ..., [0,0,0,total_site[3]-1], ]
    # need total_site to set the psel.geo property
```

## FieldSelection

```python
class FieldSelection:
    # self.cdata
    def __init__(self, total_site = None, n_per_tslice = -1, rs = None, psel = None):
    def __imatmul__(self, v1):
    def copy(self):
    def set_uniform(self, total_site, val = 0):
        # default (val = 0) select every sites
        # val = -1 deselection everything
    def set_rand(self, rs, total_site, n_per_tslice):
    def add_psel(self, psel, rank_psel = 1024 * 1024 * 1024 * 1024 * 1024):
        # Add psel points to the selection, with the rank specified as rank_psel.
        # If the point is already selected with lower rank, the rank is unchanged.
    def update(self, n_per_tslice = -1):
        # if n_per_tslice < 0: only update various indices
        # if n_per_tslice >= 0: only update parameters (n_per_tslice and prob)
    def select_rank_range(self, rank_start = 0, rank_stop = -1):
        # return new fsel with selected points that
        # rank_start <= rank and (rank < rank_stop or rank_stop == -1)
        # Does NOT change the n_per_tslice parameter for the new fsel
    def select_t_range(self, rank_start = 0, rank_stop = -1):
        # return new fsel with selected points that
        # t_start <= t and (t < t_stop or t_stop == -1)
        # rank_start <= rank < rank_stop (rank_stop = -1 implies unlimited)
        # Does NOT change the n_per_tslice parameter for the new fsel
    def to_psel(self):
    def to_psel_local(self):
    def save(self, path):
    def load(self, path, n_per_tslice):
    def geo(self):
    def total_site(self):
    def n_elems(self):
    def n_per_tslice(self):
    def prob(self):
        # return fsel.prob
        # n_per_tslice / spatial_volume
    def idx_from_coordinate(xg):
    def coordinate_from_idx(idx):

def is_matching_fsel(fsel1, fsel2):
```

## Field

```python
class Field:
	# self.ctype
    # self.cdata
    def __init__(self, ctype, geo = None, multiplicity = None):
    def __imatmul__(self, f1):
    def copy(self, is_copying_data = True):
    def swap(self, x):
    def total_site(self):
    def multiplicity(self):
    def sizeof_m(self):
    def geo(self):
    def __iadd__(self, f1):
    def __isub__(self, f1):
    def __imul__(self, factor):
		# factor can be float, complex, FieldM<Complex,1>
    def set_zero(self):
    def set_unit(self, coef = 1.0):
    def set_rand(self, rng, upper = 1.0, lower = 0.0):
    def set_rand_g(self, rng, center = 0.0, sigma = 1.0):
    def qnorm(self):
    def crc32(self):
    def save(self, path, *args):
    def load(self, path, *args):
    def save_64(self, path, *args):
    def save_double(self, path, *args):
    def save_float_from_double(self, path, *args):
    def load_64(self, path, *args):
    def load_double(self, path, *args):
    def load_double_from_float(self, path, *args):
    def as_complex_field(self):
		# return new Field("Complex") with the same content
    def from_complex_field(self, f):
    def __getitem__(self, i):
        # i can be (xg, m,) or xg
    def __setitem__(self, i, val):
        # i can be (xg, m,) or xg
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
    def xg_list(self):
    def glb_sum(self):
    def glb_sum_tslice(self):

def split_fields(fs, f):
def merge_fields(f, fs):
def merge_fields_ms(f, fms):
    # fms = [ (f0, m0,), (f1, m1,), ... ]
    # f.get_elem(x, m) = fms[m][0].get_elem(x, fms[m][1])
def mk_merged_fields_ms(fms):
    # fms = [ (f0, m0,), (f1, m1,), ... ]
    # f.get_elem(x, m) = fms[m][0].get_elem(x, fms[m][1])
    # return f

def field_expanded(f, expansion_left, expansion_right):
def refresh_expanded(field, comm_plan = None):
def refresh_expanded_1(field):
def make_field_expand_comm_plan(comm_marks):
    # comm_marks is of type Field("int8_t")
def mk_phase_field(geo: Geometry, lmom):
    # lmom is in lattice momentum unit
    # exp(i * 2*pi/L * lmom \cdot xg )
def mk_fft(is_forward, *, is_only_spatial = False, is_normalizing = False):
```

## SelectedField

```python
class SelectedField:
    # self.ctype
    # self.fsel
    # self.cdata
    def __init__(self, ctype, fsel, multiplicity = None):
    def __imatmul__(self, f1):
        # won't change self.fsel
    def copy(self, is_copying_data = True):
    def swap(self, x):
    def n_elems(self):
    def total_site(self):
    def multiplicity(self):
    def geo(self):
    def __iadd__(self, f1):
    def __isub__(self, f1):
    def __imul__(self, factor):
    def set_zero(self):
    def qnorm(self):
    def get_elems(self, idx):
    def get_elem(self, idx, m = None):
    def set_elems(self, idx, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
    def set_elem(self, idx, m, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
    def __getitem__(self, i):
        # i can be (idx, m,) or idx
    def __setitem__(self, i, val):
        # i can be (idx, m,) or idx
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
    def save(self, path, *args):
        # possible way to call:
        # f.save(path)
        # f.save(path, new_size_node)
        # f.save(sfw, fn)
    def load(self, path, *args):
        # possible way to call:
        # f.load(path)
        # f.load(sfr, fn)
        # if self.fsel is None, self.fsel will be set during f.load(sfr, fn)
    def save_64(self, path, *args):
    def save_double(self, path, *args):
    def save_float_from_double(self, path, *args):
    def load_64(self, path, *args):
    def load_double(self, path, *args):
    def load_double_from_float(self, path, *args):
    def glb_sum_tslice(self):
```

## SelectedPoints

```python
class SelectedPoints:
    # self.ctype
    # self.psel
    # self.cdata
    def __init__(self, ctype, psel, multiplicity = None):
    def __imatmul__(self, f1):
        # won't change self.psel
    def copy(self, is_copying_data = True):
    def swap(self, x):
    def n_points(self):
    def multiplicity(self):
    def __iadd__(self, f1):
    def __isub__(self, f1):
    def __imul__(self, factor):
    def set_zero(self):
    def qnorm(self):
    def get_elems(self, idx):
    def get_elem(self, idx, m = None):
    def set_elems(self, idx, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
    def set_elem(self, idx, m, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
    def __getitem__(self, i):
        # i can be (idx, m,) or idx
    def __setitem__(self, i, val):
        # i can be (idx, m,) or idx
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
    def save(self, path):
    def load(self, path):
    def save_complex(self, path):
    def load_complex(self, path):
    def to_lat_data(self):
    def from_lat_data(self, ld):
    def to_numpy(self):
    def from_numpy(self, arr):
        # need to be already initialized with ctype and psel
```

## Prop SelProp PselProp

```python
class Prop(Field):
    def __init__(self, geo = None, *, ctype = None, multiplicity = None):
    def copy(self, is_copying_data = True):

class SelProp(SelectedField):
    def __init__(self, fsel, *, ctype = None, multiplicity = None):
    def copy(self, is_copying_data = True):

class PselProp(SelectedPoints):
    def __init__(self, psel, *, ctype = None, multiplicity = None):
    def copy(self, is_copying_data = True):

def mk_point_src(geo, xg, value = 1.0):
def mk_wall_src(geo, tslice, lmom = None):
def mk_rand_u1_src(sel, rs):
    # return (prop_src, fu1,) where prop_src = Prop() and fu1 = Field("Complex")
    # fu1 stores the random u1 numbers (fu1.multiplicity() == 1)
    # sel can be psel or fsel
def get_rand_u1_sol(prop_sol, fu1, sel):
def mk_rand_u1_prop(inv, sel, rs):
    # interface function
    # return s_prop
    # sel can be psel or fsel
```

## GaugeField

```python
class GaugeField(Field):
    def __init__(self, geo = None, *, ctype = None, multiplicity = None):
    def copy(self, is_copying_data = True):
    def save(self, path):
    def load(self, path):
    def set_rand(self, rng, sigma = 0.5, n_step = 1):
    def unitarize(self):
    def plaq(self):
    def link_trace(self):
    def twist_boundary_at_boundary(self, lmom : float = -0.5, mu : int = 3):
        # modify in place
    def show_info(self):

def gf_show_info(gf):
def gf_avg_plaq(gf):
def gf_avg_link_trace(gf):
def gf_wilson_line_no_comm(wlf, m, gf_ext, path, path_n = None):
    # wlf = Field("ColorMatrix", geo)
    # will only modify the m'th component
    # e.g. path = [ mu, mu, nu, -mu-1, -mu-1, ]
    # e.g. path = [ mu, nu, -mu-1, ], path_n = [ 2, 1, 2, ]
def gf_wilson_lines_no_comm(gf_ext, path_list):
    # path_list = [ path_spec, ... ]
    # e.g. path_spec = [ mu, mu, nu, -mu-1, -mu-1, ]
    # e.g. path_spec = ([ mu, nu, -mu-1, ], [ 2, 1, 2, ],)
    # return wlf
def gf_avg_wilson_loop_normalized_tr(gf, l, t):
def set_g_rand_color_matrix_field(fc, rng, sigma, n_steps = 1):
def gf_twist_boundary_at_boundary(gf : GaugeField, lmom : float = -0.5, mu : int = 3):
    # modify gf in place
def mk_left_expanded_gauge_field(gf):
```

## GaugeTransform

```python
class GaugeTransform(Field):
    def __init__(self, geo = None, *, ctype = None, multiplicity = None):
    def copy(self, is_copying_data = True):
    def set_rand(self, rng, sigma = 0.5, n_step = 1):
    def unitarize(self):
    def __mul__(self, other):
        # other can be GaugeTransform, GaugeField, Prop, SelProp, PselProp, list
    def inv(self):
```

## Inverter InverterDwfFreeField InverterDomainWall InverterGaugeTransform

```python
class Inverter:

class InverterDwfFreeField(Inverter):
    # self.mass
    # self.m5
    # self.momtwist
    # self.timer
    def __init__(self, *, mass, m5 = 1.0, momtwist = None, qtimer = TimerNone()):
    def __mul__(self, prop_src):
        # prop_src: prop or [ prop, ... ]

class InverterDomainWall(Inverter):
	# self.cdata
	# self.timer
    def __init__(self, *, gf, fa, qtimer = TimerNone()):
    def __mul__(self, prop_src):
        # prop_src: prop or [ prop, ... ]
    def stop_rsd(self):
    def set_stop_rsd(self, stop_rsd):
    def max_num_iter(self):
    def set_max_num_iter(self, max_num_iter):
    def max_mixed_precision_cycle(self):
    def set_max_mixed_precision_cycle(self, max_mixed_precision_cycle):

class InverterGaugeTransform(Inverter):
    # self.inverter
    # self.gt
    # self.gt_inv
    # self.timer
    def __init__(self, *, inverter, gt, qtimer = TimerNone()):
    def __mul__(self, prop_src):
```

