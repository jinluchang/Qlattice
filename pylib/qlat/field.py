import cqlat as c

from qlat.geo import *
from qlat.rng import *

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

    def copy(self):
        f = Field(self.ctype)
        f @= self
        return f

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

    def set_rand(self, rng, upper = 1.0, lower = 1.0):
        assert isinstance(rng, RngState)
        c.set_u_rand_double_field(self, rng, upper, lower)

    def qnorm(self):
        return c.qnorm_field(self)

    def crc32(self):
        return c.crc32_field(self)

    def sparse(self, sel):
        from qlat.field_selection import PointSelection, FieldSelection
        if isinstance(sel, PointSelection):
            from qlat.selected_points import SelectedPoints, set_selected_points
            psel = sel
            sp = SelectedPoints(self.ctype)
            set_selected_points(sp, self, psel)
            return sp
        elif isinstance(sel, FieldSelection):
            from qlat.selected_field import SelectedField, set_selected_field
            fsel = sel
            sf = SelectedField(self.ctype)
            set_selected_field(sf, self, fsel)
            return sf
        else:
            raise Exception("Field.sparse")

    def save(self, path, new_size_node = None):
        assert isinstance(path, str)
        if new_size_node == None:
            return c.save_field(self, path)
        else:
            return c.save_field(self, path, new_size_node)

    def load(self, path):
        assert isinstance(path, str)
        return c.load_field(self, path)

    def save_64(self, path, new_size_node = None):
        assert isinstance(path, str)
        f = self.copy()
        f.to_from_endianness("big_64")
        return f.save(path, new_size_node)

    def save_double(self, path, new_size_node = None):
        return self.save_64(path, new_size_node)

    def save_float_from_double(self, path, new_size_node = None):
        ff = Field("float")
        ff.float_from_double(self)
        ff.to_from_endianness("big_32")
        return ff.save(path, new_size_node)

    def load_64(self, path):
        assert isinstance(path, str)
        ret = self.load(path)
        self.to_from_endianness("big_64")
        return ret

    def load_double(self, path):
        return self.load_64(path)

    def load_double_fromt_float(self, path):
        ff = Field("float")
        ret = ff.load(path)
        ff.to_from_endianness("big_32")
        self.double_from_float(ff)
        return ret

    def float_from_double(self, f):
        assert isinstance(f, Field)
        assert self.ctype == "float"
        c.convert_float_from_double_field(self, f)

    def double_from_float(self, ff):
        assert isinstance(ff, Field)
        assert ff.ctype == "float"
        c.convert_double_from_float_field(self, ff)

    def to_from_endianness(self, tag):
        assert isinstance(tag, str)
        c.to_from_endianness(self, tag)

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
