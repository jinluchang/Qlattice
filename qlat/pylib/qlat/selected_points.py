import qlat.cqlat as c

from qlat_utils import *

from qlat.geometry import *
from qlat.field_selection import *
from qlat.field import *

class SelectedPoints:

    # self.ctype
    # self.psel
    # self.cdata

    def __init__(self, ctype, psel, multiplicity = None):
        assert isinstance(ctype, str)
        assert isinstance(psel, PointSelection)
        self.ctype = ctype
        self.psel = psel
        if multiplicity is None:
            self.cdata = c.mk_spfield(ctype)
        else:
            assert isinstance(multiplicity, int)
            self.cdata = c.mk_spfield_psel(ctype, psel, multiplicity)

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_spfield(self)

    def __imatmul__(self, f1):
        # won't change self.psel
        from qlat.selected_field import SelectedField
        assert f1.ctype == self.ctype
        if isinstance(f1, SelectedPoints):
            # two psel must be the same object
            if self.psel is f1.psel:
                c.set_spfield(self, f1)
            else:
                raise Exception("SelectedPoints @= psel not match")
        elif isinstance(f1, SelectedField):
            # psel must be subset of fsel
            c.set_spfield_sfield(self, f1)
        elif isinstance(f1, Field):
            c.set_spfield_field(self, f1)
        else:
            raise Exception("SelectedPoints @= type mismatch")
        return self

    def copy(self, is_copying_data = True):
        f = SelectedPoints(self.ctype, self.psel)
        if is_copying_data:
            f @= self
        return f

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def swap(self, x):
        assert isinstance(x, SelectedPoints)
        assert x.ctype == self.ctype
        assert x.psel is self.psel
        cdata = x.cdata
        x.cdata = self.cdata
        self.cdata = cdata

    def n_points(self):
        return c.get_n_points_spfield(self)

    def multiplicity(self):
        return c.get_multiplicity_spfield(self)

    def __iadd__(self, f1):
        assert isinstance(f1, SelectedPoints) and f1.ctype == self.ctype
        c.set_add_spfield(self, f1)
        return self

    def __isub__(self, f1):
        assert isinstance(f1, SelectedPoints) and f1.ctype == self.ctype
        c.set_sub_spfield(self, f1)
        return self

    def __imul__(self, factor):
        assert isinstance(factor, float)
        c.set_mul_double_spfield(self, factor)
        return self

    def set_zero(self):
        c.set_zero_spfield(self)

    def qnorm(self):
        return c.qnorm_spfield(self)

    def get_elems(self, idx):
        return np.array(c.get_elems_spfield(self, idx))

    def get_elem(self, idx, m = None):
        if m is None:
            return np.array(c.get_elem_spfield(self, idx))
        else:
            return np.array(c.get_elem_spfield(self, idx, m))

    def set_elems(self, idx, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(val, bytes):
            return c.set_elems_spfield(self, idx, val)
        elif isinstance(val, np.ndarray):
            return self.set_elems(idx, val.tobytes())
        else:
            assert False

    def set_elem(self, idx, m, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(val, bytes):
            return c.set_elem_spfield(self, idx, m, val)
        elif isinstance(val, np.ndarray):
            return self.set_elem(idx, m, val.tobytes())
        else:
            assert False

    def __getitem__(self, i):
        # i can be (idx, m,) or idx
        if isinstance(i, tuple) and len(i) == 2 and isinstance(i[0], int):
            idx, m = i
            return self.get_elem(idx, m)
        elif isinstance(i, int):
            idx = i
            return self.get_elems(idx)
        else:
            assert False
            return None

    def __setitem__(self, i, val):
        # i can be (idx, m,) or idx
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(i, tuple) and len(i) == 2 and isinstance(i[0], int):
            idx, m = i
            return self.set_elem(idx, m, val)
        elif isinstance(i, int):
            idx = i
            return self.set_elems(idx, val)
        else:
            assert False
            return None

    def save(self, path):
        assert isinstance(path, str)
        return self.save_complex(path)

    def load(self, path):
        assert isinstance(path, str)
        return self.load_complex(path)

    def save_complex(self, path):
        assert isinstance(path, str)
        mk_file_dirs_info(path)
        return c.save_complex_spfield(self, path)

    def load_complex(self, path):
        assert isinstance(path, str)
        return c.load_complex_spfield(self, path)

    def to_lat_data(self):
        assert self.ctype in field_ctypes_complex
        ld = LatData()
        c.lat_data_from_complex_spfield(ld, self)
        return ld

    def from_lat_data(self, ld):
        assert self.ctype in field_ctypes_complex
        assert isinstance(ld, LatData)
        c.complex_spfield_from_lat_data(self, ld)

    def to_numpy(self):
        n_points = self.n_points()
        return np.array([ self.get_elems(idx) for idx in range(n_points) ])

    def from_numpy(self, arr):
        # need to be already initialized with ctype and psel
        # arr.shape[0] == n_points
        n_points = self.n_points()
        assert arr.shape[0] == n_points
        for idx in range(n_points):
            self.set_elems(idx, arr[idx])

###

@timer
def set_selected_points(sp, f):
    # deprecated use @=
    from qlat.selected_field import SelectedField
    assert isinstance(sp, SelectedPoints)
    if isinstance(f, Field):
        c.set_spfield_field(sp, f)
    elif isinstance(f, SelectedField):
        c.set_spfield_sfield(sp, f)
    else:
        raise Exception("set_selected_points")
