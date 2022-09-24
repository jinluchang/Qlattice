import qlat.cqlat as c

from qlat_utils import *

from qlat.geometry import *
from qlat.utils_io import *
import numpy as np

field_ctypes_complex = [
        "ColorMatrix",
        "WilsonMatrix",
        "NonRelWilsonMatrix",
        "SpinMatrix",
        "WilsonVector",
        "Complex",
        ]

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

    # self.ctype
    # self.cdata

    def __init__(self, ctype, geo = None, multiplicity = None):
        assert isinstance(ctype, str)
        self.ctype = ctype
        if geo is None:
            self.cdata = c.mk_field(ctype)
        elif multiplicity is None:
            assert isinstance(geo, Geometry)
            self.cdata = c.mk_field(ctype, geo)
        else:
            assert isinstance(geo, Geometry)
            assert isinstance(multiplicity, int)
            self.cdata = c.mk_field(ctype, geo, multiplicity)

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_field(self)

    def __imatmul__(self, f1):
        # f1 can be Field, SelectedField, SelectedPoints
        # field geo does not change if already initialized
        assert f1.ctype == self.ctype
        if isinstance(f1, Field):
            c.set_field(self, f1)
        else:
            from qlat.selected_field import SelectedField
            from qlat.selected_points import SelectedPoints
            if isinstance(f1, SelectedField):
                c.set_field_sfield(self, f1)
            elif isinstance(f1, SelectedPoints):
                c.set_field_spfield(self, f1)
            else:
                raise Exception(f"Field @= type mismatch {type(self)} {type(f1)}")
        return self

    def copy(self, is_copying_data = True):
        f = Field(self.ctype)
        if is_copying_data:
            f @= self
        return f

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def swap(self, x):
        assert isinstance(x, Field)
        assert x.ctype == self.ctype
        cdata = x.cdata
        x.cdata = self.cdata
        self.cdata = cdata

    def total_site(self):
        return c.get_total_site_field(self)

    def multiplicity(self):
        return c.get_multiplicity_field(self)

    def sizeof_m(self):
        return c.get_sizeof_m_field(self)

    def geo(self):
        geo = Geometry((0, 0, 0, 0))
        c.set_geo_field(geo, self)
        return geo

    def mview(self):
        return c.get_mview_field(self)

    def __iadd__(self, f1):
        # f1 can be Field, SelectedField, SelectedPoints
        if isinstance(f1, Field):
            assert f1.ctype == self.ctype
            c.set_add_field(self, f1)
        else:
            from qlat.selected_field import SelectedField
            from qlat.selected_points import SelectedPoints
            if isinstance(f1, SelectedField):
                c.acc_field_sfield(self, f1)
            elif isinstance(f1, SelectedPoints):
                assert f1.ctype == self.ctype
                c.acc_field_spfield(self, f1)
            else:
                raise Exception(f"Field += type mismatch {type(self)} {type(f1)}")
        return self

    def __isub__(self, f1):
        # f1 can be Field, SelectedField, SelectedPoints
        if isinstance(f1, Field):
            assert f1.ctype == self.ctype
            c.set_sub_field(self, f1)
        else:
            from qlat.selected_field import SelectedField
            from qlat.selected_points import SelectedPoints
            if isinstance(f1, SelectedField):
                assert f1.ctype == self.ctype
                f1n = f1.copy()
                f1n *= -1
                c.acc_field_sfield(self, f1n)
            elif isinstance(f1, SelectedPoints):
                assert f1.ctype == self.ctype
                f1n = f1.copy()
                f1n *= -1
                c.acc_field_spfield(self, f1n)
            else:
                raise Exception(f"Field += type mismatch {type(self)} {type(f1)}")
            assert False
        return self

    def __imul__(self, factor):
		# factor can be float, complex, FieldM<Complex,1>
        if isinstance(factor, float):
            c.set_mul_double_field(self, factor)
        elif isinstance(factor, complex):
            c.set_mul_complex_field(self, factor)
        elif isinstance(factor, Field):
            assert factor.ctype == "Complex"
            assert factor.multiplicity() == 1
            c.set_mul_cfield_field(self, factor)
        else:
            assert False
        return self

    def set_zero(self):
        c.set_zero_field(self)

    def set_unit(self, coef = 1.0):
        c.set_unit_field(self, coef)

    def set_rand(self, rng, upper = 1.0, lower = 0.0):
        assert isinstance(rng, RngState)
        c.set_u_rand_double_field(self, rng, upper, lower)

    def set_rand_g(self, rng, center = 0.0, sigma = 1.0):
        assert isinstance(rng, RngState)
        c.set_g_rand_double_field(self, rng, center, sigma)

    def set_checkers(self):
        # no longer needed?
        assert self.ctype=="double"
        c.set_checkers_double_field(self)

    def set_complex_from_double(self, sf):
        assert isinstance(sf, Field)
        assert self.ctype=="Complex"
        assert sf.ctype=="double"
        c.set_complex_from_double_field(self,sf)

    def set_double_from_complex(self, cf):
        assert isinstance(cf, Field)
        assert self.ctype=="double"
        assert cf.ctype=="Complex"
        c.set_double_from_complex_field(self,cf)

    def set_abs_from_complex(self, cf):
        assert isinstance(cf, Field)
        assert self.ctype=="double"
        assert cf.ctype=="Complex"
        c.set_abs_from_complex_field(self,cf)

    def set_ratio_double(self, sf1, sf2):
        assert isinstance(sf1, Field)
        assert isinstance(sf2, Field)
        assert self.ctype=="double"
        assert sf1.ctype=="double"
        assert sf2.ctype=="double"
        c.set_ratio_double_field(self,sf1,sf2)

    def less_than_double(self, sf2, mask):
        assert isinstance(sf2, Field)
        assert isinstance(mask, Field)
        assert self.ctype=="double"
        assert sf2.ctype=="double"
        assert mask.ctype=="double"
        c.less_than_double_field(self,sf2,mask)

    def invert_double(self):
        assert self.ctype=="double"
        c.invert_double_field(self)

    def multiply_double(self, factor):
        assert isinstance(factor, Field)
        assert factor.ctype=="double"
        c.multiply_double_field(self,factor)

    def qnorm(self):
        return c.qnorm_field(self)

    def crc32(self):
        return c.crc32_field(self)

    def sparse(self, sel):
        # deprecated
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
            mk_file_dirs_info(path)
            if len(args) == 0:
                return c.save_field(self, path)
            else:
                [ new_size_node, ] = args
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
            [ fn, ] = args
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
		# return new Field("Complex") with the same content
        f = Field("Complex")
        c.assign_as_complex_field(f, self)
        return f

    def from_complex_field(self, f):
        c.assign_from_complex_field(self, f)
        return f

    def get_elems(self, xg):
        return np.array(c.get_elems_field(self, xg))

    def get_elem(self, xg, m = None):
        if m is None:
            return np.array(c.get_elem_field(self, xg))
        else:
            return np.array(c.get_elem_field(self, xg, m))

    def set_elems(self, xg, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(val, bytes):
            return c.set_elems_field(self, xg, val)
        elif isinstance(val, np.ndarray):
            return self.set_elems(xg, val.tobytes())
        else:
            assert False

    def set_elem(self, xg, m, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(val, bytes):
            return c.set_elem_field(self, xg, m, val)
        elif isinstance(val, np.ndarray):
            return self.set_elem(xg, m, val.tobytes())
        else:
            assert False

    def __getitem__(self, i):
        # i can be (xg, m,) or xg
        if isinstance(i, tuple) and len(i) == 2 and isinstance(i[0], (list, tuple)):
            xg, m = i
            return self.get_elem(xg, m)
        elif isinstance(i, (list, tuple)):
            xg = i
            return self.get_elems(xg)
        else:
            assert False
            return None

    def __setitem__(self, i, val):
        # i can be (xg, m,) or xg
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(i, tuple) and len(i) == 2 and isinstance(i[0], (list, tuple)):
            xg, m = i
            return self.set_elem(xg, m, val)
        elif isinstance(i, (list, tuple)):
            xg = i
            return self.set_elems(xg, val)
        else:
            assert False
            return None

    def xg_list(self):
        # return xg for all local sites
        return self.geo().xg_list()

    def field_shift(self, shift):
        # return new shifted Field
        # shift is the coordinate to shift the field
        f1 = self.copy(is_copying_data = False)
        c.field_shift_field(f1, self, shift)
        return f1

    def reflect(self):
        # reflect the field, return None
        return c.reflect_field(self)

    def glb_sum(self):
        if self.ctype in field_ctypes_double:
            return np.array(c.glb_sum_double_field(self))
        elif self.ctype in field_ctypes_long:
            return np.array(c.glb_sum_long_field(self))
        else:
            assert False
            return None

    def glb_sum_tslice(self, *, t_dir = 3):
        # return SelectedPoints(self.ctype, get_psel_tslice(self.total_site(), t_dir = t_dir))
        from qlat.field_selection import get_psel_tslice
        from qlat.selected_points import SelectedPoints
        psel = get_psel_tslice(self.total_site(), t_dir = t_dir)
        sp = SelectedPoints(self.ctype, psel)
        if self.ctype in field_ctypes_double:
            c.glb_sum_tslice_double_field(sp, self, t_dir)
            return sp
        elif self.ctype in field_ctypes_long:
            c.glb_sum_tslice_long_field(sp, self, t_dir)
            return sp
        else:
            assert False
            return None

###

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

def merge_fields_ms(f, fms):
    # fms = [ (f0, m0,), (f1, m1,), ... ]
    # f.get_elem(x, m) = fms[m][0].get_elem(x, fms[m][1])
    multiplicity = len(fms)
    assert multiplicity >= 1
    assert isinstance(f, Field)
    assert f.ctype == fms[0][0].ctype
    fs, ms = zip(*fms)
    c.merge_fields_ms_field(f, fs, ms)

def mk_merged_fields_ms(fms):
    # fms = [ (f0, m0,), (f1, m1,), ... ]
    # f.get_elem(x, m) = fms[m][0].get_elem(x, fms[m][1])
    # return f
    multiplicity = len(fms)
    assert multiplicity >= 1
    for m in range(multiplicity):
        assert isinstance(fms[m][0], Field)
        assert isinstance(fms[m][1], int)
    ctype = fms[0][0].ctype
    for m in range(multiplicity):
        assert ctype == fms[m][0].ctype
    f = Field(ctype)
    merge_fields_ms(f, fms)
    return f
