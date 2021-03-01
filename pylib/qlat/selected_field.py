import cqlat as c

from qlat.timer import *
from qlat.geometry import *
from qlat.field import *
from qlat.field_selection import *

class SelectedField:

    def __init__(self, ctype, fsel, multiplicity = None):
        assert isinstance(ctype, str)
        assert isinstance(fsel, FieldSelection)
        self.ctype = ctype
        self.fsel = fsel
        if multiplicity is None:
            self.cdata = c.mk_sfield(ctype)
        else:
            assert isinstance(multiplicity, int)
            self.cdata = c.mk_sfield_fsel(ctype, fsel, multiplicity)

    def __del__(self):
        c.free_sfield(self)

    def __imatmul__(self, f1):
        assert isinstance(f1, SelectedField) and f1.ctype == self.ctype
        self.fsel = f1.fsel
        c.set_sfield(self, f1)
        return self

    def copy(self):
        f = SelectedField(self.ctype, self.fsel)
        f @= self
        return f

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
        assert isinstance(f1, SelectedField) and f1.ctype == self.ctype
        c.set_add_sfield(self, f1)
        return self

    def __isub__(self, f1):
        assert isinstance(f1, SelectedField) and f1.ctype == self.ctype
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
            sp = SelectedPoints(self.ctype, psel)
            set_selected_points(sp, self)
            return sp
        elif isinstance(sel, FieldSelection):
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
            assert len(args) == 0
            return c.save_sfield(self, path)
        elif isinstance(path, ShuffledFieldsWriter):
            sfw = path
            [fn] = args
            return sfw.write(fn, self)
        else:
            raise Exception("SelectedField.save")

    def load(self, path, *args):
        # possible way to call:
        # f.load(path)
        # f.load(sfr, fn)
        from qlat.fields_io import ShuffledFieldsReader
        if isinstance(path, str):
            return c.load_sfield(self, path)
        elif isinstance(path, ShuffledFieldsReader):
            sfr = path
            [fn] = args
            return sfr.read(fn, self)
        else:
            raise Exception("SelectedField.load")

    def save_64(self, *path):
        f = self.copy()
        f.to_from_endianness("big_64")
        return f.save(*path)

    def save_double(self, *path):
        return self.save_64(*path)

    def save_float_from_double(self, *path):
        ff = SelectedField("float", self.fsel)
        ff.float_from_double(self)
        ff.to_from_endianness("big_32")
        return ff.save(*path)

    def load_64(self, *path):
        ret = self.load(*path)
        if ret > 0:
            self.to_from_endianness("big_64")
        return ret

    def load_double(self, *path):
        return self.load_64(*path)

    def load_double_from_float(self, *path):
        ff = SelectedField("float", self.fsel)
        ret = ff.load(*path)
        if ret > 0:
            ff.to_from_endianness("big_32")
            self.double_from_float(ff)
        return ret

    def float_from_double(self, f):
        assert isinstance(f, SelectedField)
        assert self.ctype == "float"
        c.convert_float_from_double_sfield(self, f)

    def double_from_float(self, ff):
        assert isinstance(ff, SelectedField)
        assert ff.ctype == "float"
        c.convert_double_from_float_sfield(self, ff)

    def to_from_endianness(self, tag):
        assert isinstance(tag, str)
        c.to_from_endianness_sfield(self, tag)

@timer
def set_selected_field(sf, f):
    assert isinstance(sf, SelectedField)
    if isinstance(f, Field):
        c.set_sfield_field(sf, f)
    elif isinstance(f, SelectedField):
        c.set_sfield_sfield(sf, f)
    else:
        raise Exception("set_selected_field")
