import cqlat as c

from qlat.timer import *
from qlat.geo import *
from qlat.field_selection import *
from qlat.field import *

class SelectedField:

    def __init__(self, ctype, v1 = None, v2 = None, v3 = None):
        assert isinstance(ctype, str)
        self.ctype = ctype
        if v1 == None:
            self.cdata = c.mk_sfield(ctype)
            self.fsel = None
        elif isinstance(v1, Geometry):
            geo = v1
            n_elems = v2
            multiplicity = v3
            assert isinstance(n_elems, int)
            if multiplicity == None:
                self.cdata = c.mk_sfield(ctype, geo, n_elems)
            else:
                assert isinstance(multiplicity, int)
                self.cdata = c.mk_sfield(ctype, geo, n_elems, multiplicity)
            self.fsel = None
        elif isinstance(v1, FieldSelection):
            fsel = v1
            multiplicity = v2
            assert isinstance(multiplicity, int)
            assert v3 == None
            self.cdata = c.mk_sfield_fsel(ctype, fsel, multiplicity)
            self.fsel = fsel
        else:
            raise Exception("SelectedField.__init__")

    def __del__(self):
        c.free_sfield(self)

    def __imatmul__(self, f1):
        assert isinstance(f1, SelectedField) and f1.ctype == self.ctype
        c.set_sfield(self, f1)
        return self

    def n_elems(self):
        return c.get_n_elems_sfield(self)

    def total_site(self):
        return c.get_total_site_field(self)

    def multiplicity(self):
        return c.get_multiplicity_sfield(self)

    def geo(self):
        geo = Geometry((0, 0, 0, 0))
        c.set_geo_sfield(geo, self)
        return geo

    def __iadd__(self, f1):
        assert isinstance(f1, Field) and f1.ctype == self.ctype
        c.set_add_sfield(self, f1)
        return self

    def __isub__(self, f1):
        assert isinstance(f1, Field) and f1.ctype == self.ctype
        c.set_sub_sfield(self, f1)
        return self

    def __imul__(self, factor):
        assert isinstance(factor, float)
        c.set_mul_double_sfield(self, factor)
        return self

    def set_zero(self):
        c.set_zero_sfield(self)

    def qnorm(self):
        return c.qnorm_sfield(self)

    def sparse(self, sel):
        if isinstance(sel, PointSelection):
            from qlat.selected_points import SelectedPoints, set_selected_points
            psel = sel
            sp = SelectedPoints(self.ctype)
            set_selected_points(sp, self, psel, self.fsel)
            return sp
        elif isinstance(sel, FieldSelection):
            fsel = sel
            sf = SelectedField(self.ctype)
            set_selected_field(sf, self, fsel, self.fsel)
        else:
            raise Exception("Field.sparse")

@timer
def set_selected_field(sf, v1, v2, v3 = None):
    if isinstance(v1, Field):
        f = v1
        fsel = v2
        assert isinstance(fsel, FieldSelection)
        assert v3 == None
        c.set_sfield_field(sf, f, fsel)
        sf.fsel = fsel
    elif isinstance(v1, SelectedField):
        sf0 = v1
        fsel = v2
        fsel0 = v3
        assert isinstance(fsel, FieldSelection)
        assert isinstance(fsel0, FieldSelection)
        c.set_sfield_sfield(sf, sf0, fsel, fsel0)
        sf.fsel = fsel
    else:
        raise Exception("set_selected_field")
