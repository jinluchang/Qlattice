import cqlat as c

from qlat.timer import *
from qlat.geometry import *
from qlat.field_selection import *
from qlat.field import *

class SelectedPoints:

    def __init__(self, ctype, psel, multiplicity = None):
        self.ctype = ctype
        self.psel = psel
        if psel == None or multiplicity == None:
            self.cdata = c.mk_spfield(ctype)
        else:
            assert isinstance(psel, PointSelection)
            assert isinstance(multiplicity, int)
            self.cdata = c.mk_spfield_psel(ctype, psel, multiplicity)

    def __del__(self):
        c.free_spfield(self)

    def __imatmul__(self, f1):
        assert isinstance(f1, SelectedPoints) and f1.ctype == self.ctype
        self.psel = f1.psel
        c.set_spfield(self, f1)
        return self

    def copy(self):
        f = SelectedPoints(self.ctype, self.psel)
        f @= self
        return f

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

    def save(self, path):
        assert isinstance(path, str)
        return self.save_complex(path)

    def load(self, path):
        assert isinstance(path, str)
        return self.load_complex(path)

    def save_complex(self, path):
        assert isinstance(path, str)
        return c.save_complex_spfield(self, path)

    def load_complex(self, path):
        assert isinstance(path, str)
        return c.load_complex_spfield(self, path)

@timer
def set_selected_points(sp, f):
    assert isinstance(sp, SelectedPoints)
    if isinstance(f, Field):
        c.set_spfield_field(sp, f)
    elif isinstance(f, SelectedField):
        c.set_spfield_sfield(sp, f)
    else:
        raise Exception("set_selected_points")
