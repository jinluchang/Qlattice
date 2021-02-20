import cqlat as c

from qlat.geo import *

class Field:

    def __init__(self, ctype, geo = None, multiplicity = None):
        assert isinstance(ctype, str)
        self.ctype = ctype
        if geo == None:
            self.cdata = c.mk_field(ctype)
        elif multiplicity == None:
            assert isinstance(geo, Geometry)
            self.cdata = c.mk_field(ctype, geo)
        else:
            assert isinstance(geo, Geometry)
            assert isinstance(multiplicity, int)
            self.cdata = c.mk_field(ctype, geo, multiplicity)

    def __del__(self):
        c.free_field(self)

    def __imatmul__(self, f1):
        assert isinstance(f1, Field) and f1.ctype == self.ctype
        c.set_field(self, f1)
        return self

    def total_site(self):
        return c.get_total_site_field(self)

    def multiplicity(self):
        return c.get_multiplicity_field(self)

    def geo(self):
        geo = Geometry((0, 0, 0, 0))
        c.set_geo_field(geo, self)
        return geo

    def mview(self):
        return c.get_mview_field(self)

    def __iadd__(self, f1):
        assert isinstance(f1, Field) and f1.ctype == self.ctype
        c.set_add_field(self, f1)
        return self

    def __isub__(self, f1):
        assert isinstance(f1, Field) and f1.ctype == self.ctype
        c.set_sub_field(self, f1)
        return self

    def __imul__(self, factor):
        assert isinstance(factor, float)
        c.set_mul_double_field(self, factor)
        return self

    def set_zero(self):
        c.set_zero_field(self)

    def set_unit(self, coef = 1.0):
        c.set_unit_field(self, coef)

    def qnorm(self):
        return c.qnorm_field(self)

def split_fields(fs, f):
    nf = len(fs)
    assert nf >= 1
    ctype = f.ctype
    for i in range(nf):
        if not isinstance(fs[i], Field):
            fs[i] = Field(ctype)
        else:
            assert fs[i].ctype == ctype
    c.split_fields_field(fs, f)

def merge_fields(f, fs):
    nf = len(fs)
    assert nf >= 1
    assert isinstance(f, Field)
    assert f.ctype == fs[0].ctype
    c.merge_fields_field(f, fs)
