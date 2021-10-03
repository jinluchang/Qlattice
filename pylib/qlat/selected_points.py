import cqlat as c

from qlat.timer import *
from qlat.geometry import *
from qlat.field_selection import *
from qlat.field import *
from qlat.lat_io import *

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
        c.free_spfield(self)

    def __imatmul__(self, f1):
        # won't change self.psel
        from qlat.selected_field import SelectedField
        assert f1.ctype == self.ctype
        if isinstance(f1, SelectedPoints):
            if self.psel is f1.psel:
                c.set_spfield(self, f1)
            else:
                raise Exception("SelectedPoints @= psel not match")
        elif isinstance(f1, SelectedField):
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
        return c.get_elems_spfield(self, idx)

    def get_elem(self, idx, m = None):
        if m is None:
            return c.get_elem_spfield(self, idx)
        else:
            return c.get_elem_spfield(self, idx, m)

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
        ld = LatData()
        c.lat_data_from_complex_spfield(ld, self)
        return ld

    def from_lat_data(self, ld):
        assert isinstance(ld, LatData)
        c.complex_spfield_from_lat_data(self, ld)

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
