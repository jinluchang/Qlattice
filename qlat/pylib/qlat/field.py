import qlat.c as c

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

field_ctypes_float = [
        "ComplexF",
        "float",
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
        if self.ctype in field_ctypes_double:
            c.set_u_rand_double_field(self, rng, upper, lower)
        elif self.ctype in field_ctypes_float:
            c.set_u_rand_float_field(self, rng, upper, lower)
        else:
            assert False

    def set_rand_g(self, rng, center = 0.0, sigma = 1.0):
        assert isinstance(rng, RngState)
        if self.ctype in field_ctypes_double:
            c.set_g_rand_double_field(self, rng, center, sigma)
        else:
            assert False

    def multiply_double(self, factor):
        assert isinstance(factor, Field)
        assert factor.ctype == "double"
        c.multiply_double_field(self,factor)

    def qnorm(self):
        return c.qnorm_field(self)

    def crc32(self):
        return c.crc32_field(self)

    def sparse(self, sel):
        # deprecated
        displayln_info("Field.sparse: deprecated")
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

    def save_direct(self, path, *args):
        # save Field directly (without any conversion of endianness or precision)
        # possible way to call:
        # f.save_direct(path)
        # f.save_direct(path, new_size_node)
        # f.save_direct(sfw, fn)
        from qlat.fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            if len(args) == 0:
                n_bytes = c.save_field(self, path)
            else:
                [ new_size_node, ] = args
                n_bytes = c.save_field(self, path, new_size_node)
        elif isinstance(path, ShuffledFieldsWriter):
            sfw = path
            [fn] = args
            n_bytes = sfw.write(fn, self)
        else:
            raise Exception("Field.save_direct")
        assert n_bytes != 0
        return n_bytes

    def load_direct(self, path, *args):
        # load Field directly (without any conversion of endianness or precision)
        # possible way to call:
        # f.load(path)
        # f.load(sfr, fn)
        from qlat.fields_io import ShuffledFieldsReader
        if isinstance(path, str):
            n_bytes = c.load_field(self, path)
        elif isinstance(path, ShuffledFieldsReader):
            sfr = path
            [ fn, ] = args
            n_bytes = sfr.read(fn, self)
        else:
            raise Exception("Field.load_direct")
        assert n_bytes != 0
        return n_bytes

    def save_64(self, path, *args):
        f = self.copy()
        from qlat.fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            f.to_from_endianness("big_64")
        elif isinstance(path, ShuffledFieldsWriter):
            f.to_from_endianness("little_64")
        return f.save_direct(path, *args)

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
        return ff.save_direct(path, *args)

    def load_64(self, path, *args):
        ret = self.load_direct(path, *args)
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
        ret = ff.load_direct(path, *args)
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

    def as_field(self, ctype = "Complex"):
		# return new Field(ctype) with the same content
        f = Field(ctype)
        c.assign_as_field(f, self)
        return f

    def from_field(self, f):
		# assign from f with the same content but possibly different type
        c.assign_from_field(self, f)
        return f

    def get_elems(self, index):
        # index can also be xg
        # get_elems is collective operation when xg is coordinate
        # get_elems will be specific to a single process if xg is index
        return np.array(c.get_elems_field(self, index))

    def get_elem(self, index, m = None):
        # index can also be xg
        if m is None:
            return np.array(c.get_elem_field(self, index))
        else:
            return np.array(c.get_elem_field(self, index, m))

    def set_elems(self, index, val):
        # index can also be xg
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(val, bytes):
            return c.set_elems_field(self, index, val)
        elif isinstance(val, np.ndarray):
            return self.set_elems(index, val.tobytes())
        else:
            assert False

    def set_elem(self, index, m, val):
        # index can also be xg
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(val, bytes):
            return c.set_elem_field(self, index, m, val)
        elif isinstance(val, np.ndarray):
            return self.set_elem(index, m, val.tobytes())
        else:
            assert False

    def __getitem__(self, i):
        # i can be (index, m,) or index
        # index can also be xg
        if isinstance(i, tuple) and len(i) == 2 and isinstance(i[0], (int, list, tuple,)):
            index, m = i
            return self.get_elem(index, m)
        elif isinstance(i, (int, list, tuple)):
            index = i
            return self.get_elems(index)
        else:
            assert False
            return None

    def __setitem__(self, i, val):
        # i can be (index, m,) or index
        # index can also be xg
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(i, tuple) and len(i) == 2 and isinstance(i[0], (int, list, tuple,)):
            index, m = i
            return self.set_elem(index, m, val)
        elif isinstance(i, (int, list, tuple)):
            index = i
            return self.set_elems(index, val)
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
