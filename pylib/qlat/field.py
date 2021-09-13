import cqlat as c

from qlat.geometry import *
from qlat.rng_state import *

field_ctypes_double = [
        "ColorMatrix",
        "WilsonMatrix",
        "NonRelWilsonMatrix",
        "SpinMatrix",
        "WilsonVector",
        "Complex",
        "double",
        ]

field_ctypes_long = [
        "long",
        "int64_t",
        ]

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

    def copy(self, is_copying_data = True):
        f = Field(self.ctype)
        if is_copying_data:
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
        if isinstance(factor, float):
            c.set_mul_double_field(self, factor)
            return self
        elif isinstance(factor, complex):
            c.set_mul_complex_field(self, factor)
            return self
        elif isinstance(factor, Field):
            c.set_mul_cfield_field(self, factor)
            return self
        else:
            assert False

    def set_zero(self):
        c.set_zero_field(self)

    def set_unit(self, coef = 1.0):
        c.set_unit_field(self, coef)

    def set_rand(self, rng, upper = 1.0, lower = 0.0):
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
            sp = SelectedPoints(self.ctype, psel)
            set_selected_points(sp, self)
            return sp
        elif isinstance(sel, FieldSelection):
            from qlat.selected_field import SelectedField, set_selected_field
            fsel = sel
            sf = SelectedField(self.ctype, fsel)
            set_selected_field(sf, self)
            return sf
        else:
            raise Exception("Field.sparse")

    def save(self, path, *args):
        # possible way to call:
        # f.save(path)
        # f.save(path, new_size_node)
        # f.save(sfw, fn)
        from qlat.fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            if len(args) == 0:
                return c.save_field(self, path)
            else:
                [new_size_node] = args
                return c.save_field(self, path, new_size_node)
        elif isinstance(path, ShuffledFieldsWriter):
            sfw = path
            [fn] = args
            return sfw.write(fn, self)
        else:
            raise Exception("Field.save")

    def load(self, path, *args):
        # possible way to call:
        # f.load(path)
        # f.load(sfr, fn)
        from qlat.fields_io import ShuffledFieldsReader
        if isinstance(path, str):
            return c.load_field(self, path)
        elif isinstance(path, ShuffledFieldsReader):
            sfr = path
            [fn] = args
            return sfr.read(fn, self)
        else:
            raise Exception("Field.load")

    def save_64(self, path, *args):
        f = self.copy()
        from qlat.fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            f.to_from_endianness("big_64")
        elif isinstance(path, ShuffledFieldsWriter):
            f.to_from_endianness("little_64")
        return f.save(path, *args)

    def save_double(self, path, *args):
        return self.save_64(path, *args)

    def save_float_from_double(self, path, *args):
        ff = Field("float")
        ff.float_from_double(self)
        from qlat.fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            ff.to_from_endianness("big_32")
        elif isinstance(path, ShuffledFieldsWriter):
            ff.to_from_endianness("little_32")
        return ff.save(path, *args)

    def load_64(self, path, *args):
        ret = self.load(path, *args)
        if ret > 0:
            from qlat.fields_io import ShuffledFieldsReader
            if isinstance(path, str):
                self.to_from_endianness("big_64")
            elif isinstance(path, ShuffledFieldsReader):
                self.to_from_endianness("little_64")
        return ret

    def load_double(self, path, *args):
        return self.load_64(path, *args)

    def load_double_from_float(self, path, *args):
        ff = Field("float")
        ret = ff.load(path, *args)
        if ret > 0:
            from qlat.fields_io import ShuffledFieldsReader
            if isinstance(path, str):
                ff.to_from_endianness("big_32")
            elif isinstance(path, ShuffledFieldsReader):
                ff.to_from_endianness("little_32")
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
        c.to_from_endianness_field(self, tag)

    def as_complex_field(self):
        f = Field("Complex")
        c.assign_as_complex_field(f, self)
        return f

    def from_complex_field(self, f):
        c.assign_from_complex_field(self, f)
        return f

    def get_elems(self, xg):
        return c.get_elems_field(self, xg)

    def get_elem(self, xg, m = None):
        if m is None:
            return c.get_elem_field(self, xg)
        else:
            return c.get_elem_field(self, xg, m)

    def glb_sum(self):
        if self.ctype in field_ctypes_double:
            return c.glb_sum_double_field(self)
        elif self.ctype in field_ctypes_long:
            return c.glb_sum_long_field(self)
        else:
            assert False

    def glb_sum_tslice(self):
        from qlat.field_selection import PointSelection
        from qlat.selected_points import SelectedPoints, set_selected_points
        psel = PointSelection()
        psel.set_tslice(self.total_site())
        sp = SelectedPoints(self.ctype, psel)
        if self.ctype in field_ctypes_double:
            c.glb_sum_tslice_double_field(sp, self)
            return sp
        elif self.ctype in field_ctypes_long:
            c.glb_sum_tslice_long_field(sp, self)
            return sp
        else:
            assert False

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
