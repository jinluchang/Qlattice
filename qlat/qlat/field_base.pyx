# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .geometry cimport Geometry
from .field_types cimport FieldFloat
from .selected_field_types cimport SelectedFieldFloat

from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import cqlat as c
import qlat_utils as q
import numpy as np

from .field_type_dict import *

### -------------------------------------------------------------------

def Field(type ctype, Geometry geo=None, int multiplicity=0):
    assert ctype in field_type_dict
    FieldType = field_type_dict[ctype]
    field = FieldType(geo, multiplicity)
    return field

def SelectedField(type ctype, FieldSelection fsel, int multiplicity=0):
    assert ctype in field_type_dict
    FieldType = selected_field_type_dict[ctype]
    field = FieldType(fsel, multiplicity)
    return field

def SelectedPoints(type ctype, PointsSelection psel, int multiplicity=0):
    assert ctype in field_type_dict
    FieldType = selected_points_type_dict[ctype]
    field = FieldType(psel, multiplicity)
    return field

### -------------------------------------------------------------------

cdef class FieldBase:

    ctype = ElemType

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

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
        """
        f1 can be Field, SelectedField, SelectedPoints
        """
        if isinstance(f1, FieldBase):
            assert f1.ctype is self.ctype
            c.set_add_field(self, f1)
        else:
            if isinstance(f1, SelectedFieldBase):
                c.acc_field_sfield(self, f1)
            elif isinstance(f1, SelectedPointsBase):
                assert f1.ctype is self.ctype
                c.acc_field_spfield(self, f1)
            else:
                raise Exception(f"Field += type mismatch {type(self)} {type(f1)}")
        return self

    def __isub__(self, f1):
        """
        f1 can be Field, SelectedField, SelectedPoints
        """
        if isinstance(f1, FieldBase):
            assert f1.ctype is self.ctype
            c.set_sub_field(self, f1)
        else:
            if isinstance(f1, SelectedFieldBase):
                assert f1.ctype is self.ctype
                f1n = f1.copy()
                f1n *= -1
                c.acc_field_sfield(self, f1n)
            elif isinstance(f1, SelectedPointsBase):
                assert f1.ctype is self.ctype
                f1n = f1.copy()
                f1n *= -1
                c.acc_field_spfield(self, f1n)
            else:
                raise Exception(f"Field += type mismatch {type(self)} {type(f1)}")
            assert False
        return self

    def __imul__(self, factor):
        """
		factor can be float, complex, FieldM<Complex,1>
        """
        if isinstance(factor, float):
            c.set_mul_double_field(self, factor)
        elif isinstance(factor, complex):
            c.set_mul_complex_field(self, factor)
        elif isinstance(factor, FieldBase):
            assert factor.ctype in [ ElemTypeComplex, ElemTypeDouble, ]
            c.set_mul_cfield_field(self, factor)
        else:
            assert False
        return self

    @q.timer
    def set_unit(self, coef=1.0):
        c.set_unit_field(self, coef)

    @q.timer
    def set_rand(self, rng, upper=1.0, lower=0.0):
        assert isinstance(rng, RngState)
        if self.ctype in field_ctypes_double:
            c.set_u_rand_double_field(self, rng, upper, lower)
        elif self.ctype in field_ctypes_float:
            c.set_u_rand_float_field(self, rng, upper, lower)
        else:
            assert False

    @q.timer
    def set_rand_g(self, rng, center=0.0, sigma=1.0):
        assert isinstance(rng, RngState)
        if self.ctype in field_ctypes_double:
            c.set_g_rand_double_field(self, rng, center, sigma)
        else:
            assert False

    def qnorm(self):
        return c.qnorm_field(self)

    def crc32(self):
        return c.crc32_field(self)

    def save_direct(self, path, *args):
        """
        Generic save for Field object
        save Field directly (without any conversion of endianness or precision)
        possible way to call:
        f.save_direct(path)
        f.save_direct(path, new_size_node)
        f.save_direct(sfw, fn)
        """
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
        """
        Generic load for Field object
        load Field directly (without any conversion of endianness or precision)
        Field geo and multiplicity will be determined during loading
        possible way to call:
        f.load_direct(path)
        f.load_direct(sfr, fn)
        """
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
        """
        Generic save for 64-bit size element Field object
        save 64-bit Field (do conversion of endianness)
        """
        f = self.copy()
        from .fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            f.to_from_endianness("big_64")
        elif isinstance(path, ShuffledFieldsWriter):
            f.to_from_endianness("little_64")
        return f.save_direct(path, *args)

    def save_double(self, path, *args):
        """
        Generic save for double element Field object
        save double Field as double (do conversion of endianness)
        """
        return self.save_64(path, *args)

    def save_float_from_double(self, path, *args):
        """
        Generic save for double element Field object
        save double Field as float (do conversion of endianness and precision)
        """
        ff = FieldFloat()
        ff.float_from_double(self)
        from .fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            ff.to_from_endianness("big_32")
        elif isinstance(path, ShuffledFieldsWriter):
            ff.to_from_endianness("little_32")
        return ff.save_direct(path, *args)

    def load_64(self, path, *args):
        """
        Generic load for 64-bit size element Field object
        load 64-bit Field (do conversion of endianness)
        """
        ret = self.load_direct(path, *args)
        if ret > 0:
            from .fields_io import ShuffledFieldsReader
            if isinstance(path, str):
                self.to_from_endianness("big_64")
            elif isinstance(path, ShuffledFieldsReader):
                self.to_from_endianness("little_64")
        return ret

    def load_double(self, path, *args):
        """
        Generic load for double Field object
        load double Field (do conversion of endianness)
        """
        return self.load_64(path, *args)

    def load_double_from_float(self, path, *args):
        """
        Generic load for double Field object
        load double Field from float(do conversion of endianness or precision)
        """
        ff = FieldFloat()
        ret = ff.load_direct(path, *args)
        if ret > 0:
            from .fields_io import ShuffledFieldsReader
            if isinstance(path, str):
                ff.to_from_endianness("big_32")
            elif isinstance(path, ShuffledFieldsReader):
                ff.to_from_endianness("little_32")
            self.double_from_float(ff)
        return ret

    def float_from_double(self, FieldBase f):
        """
        self needs to be FieldFloat
        """
        assert isinstance(self, FieldFloat)
        c.convert_float_from_double_field(self, f)

    def double_from_float(self, FieldFloat ff):
        """
        self can be any FieldBase subtype but need to be actually contains double precision numbers
        """
        c.convert_double_from_float_field(self, ff)

    def to_from_endianness(self, tag):
        """
        Convert between the native endianness and the endianness specified by ``tag``
        tag can be ``"big_32", "big_64", "little_32", "little_64"``
        """
        assert isinstance(tag, str)
        c.to_from_endianness_field(self, tag)

    def as_field(self, ctype=ElemTypeComplex):
        """
		return new Field(ctype) with the same content
        """
        f = Field(ctype)
        c.assign_as_field(f, self)
        return f

    def from_field(self, f):
        """
		assign from f with the same content but possibly different type
        """
        c.assign_from_field(self, f)
        return f

    def get_elems(self, index):
        """
        index can also be xg
        get_elems is collective operation when xg is coordinate
        get_elems will be specific to a single process if xg is index
        """
        return np.array(c.get_elems_field(self, index))

    def get_elem(self, index, m=None):
        """
        index can also be xg
        """
        if m is None:
            return np.array(c.get_elem_field(self, index))
        else:
            return np.array(c.get_elem_field(self, index, m))

    def set_elems(self, index, val):
        """
        index can also be xg
        val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype=complex).tobytes()
        """
        if isinstance(val, bytes):
            return c.set_elems_field(self, index, val)
        elif isinstance(val, np.ndarray):
            return self.set_elems(index, val.tobytes())
        else:
            assert False

    def set_elem(self, index, m, val):
        """
        index can also be xg
        val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype=complex).tobytes()
        """
        if isinstance(val, bytes):
            return c.set_elem_field(self, index, m, val)
        elif isinstance(val, np.ndarray):
            return self.set_elem(index, m, val.tobytes())
        else:
            assert False

    def __getitem__(self, i):
        """
        i can be (index, m,) or index
        index can also be xg
        """
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
        """
        i can be (index, m,) or index
        index can also be xg
        val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype=complex).tobytes()
        """
        if isinstance(i, tuple) and len(i) == 2 and isinstance(i[0], (int, list, tuple,)):
            index, m = i
            self.set_elem(index, m, val)
        elif isinstance(i, (int, list, tuple)):
            index = i
            self.set_elems(index, val)
        else:
            assert False

    def xg_list(self):
        """
        return xg for all local sites
        """
        return self.geo().xg_list()

    def field_shift(self, shift):
        """
        return new shifted Field
        shift is the coordinate to shift the field
        """
        f1 = self.copy(is_copying_data=False)
        c.field_shift_field(f1, self, shift)
        return f1

    def reflect(self):
        """
        reflect the field, return None
        """
        return c.reflect_field(self)

    def glb_sum(self):
        if self.ctype in field_ctypes_double:
            return np.array(c.glb_sum_double_field(self))
        elif self.ctype in field_ctypes_long:
            return np.array(c.glb_sum_long_field(self))
        else:
            assert False
            return None

    def glb_sum_tslice(self, *, t_dir=3):
        """
        return SelectedPoints(self.ctype, get_psel_tslice(self.total_site(), t_dir=t_dir))
        """
        from .field_selection_utils import get_psel_tslice
        psel = get_psel_tslice(self.total_site(), t_dir=t_dir)
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

### -------------------------------------------------------------------

def split_fields(fs, f):
    nf = len(fs)
    assert nf >= 1
    ctype = f.ctype
    for i in range(nf):
        if not isinstance(fs[i], FieldBase):
            fs[i] = Field(ctype)
        else:
            assert fs[i].ctype is ctype
    c.split_fields_field(fs, f)

def merge_fields(f, fs):
    nf = len(fs)
    assert nf >= 1
    assert isinstance(f, FieldBase)
    assert f.ctype is fs[0].ctype
    c.merge_fields_field(f, fs)

def merge_fields_ms(f, fms):
    """
    fms = [ (f0, m0,), (f1, m1,), ... ]
    f.get_elem(x, m) = fms[m][0].get_elem(x, fms[m][1])
    """
    multiplicity = len(fms)
    assert multiplicity >= 1
    assert isinstance(f, FieldBase)
    assert f.ctype is fms[0][0].ctype
    fs, ms = zip(*fms)
    c.merge_fields_ms_field(f, fs, ms)

def mk_merged_fields_ms(fms):
    """
    fms = [ (f0, m0,), (f1, m1,), ... ]
    f.get_elem(x, m) = fms[m][0].get_elem(x, fms[m][1])
    return f
    """
    multiplicity = len(fms)
    assert multiplicity >= 1
    for m in range(multiplicity):
        assert isinstance(fms[m][0], FieldBase)
        assert isinstance(fms[m][1], int)
    ctype = fms[0][0].ctype
    for m in range(multiplicity):
        assert ctype is fms[m][0].ctype
    f = Field(ctype)
    merge_fields_ms(f, fms)
    return f

### -------------------------------------------------------------------

cdef class SelectedFieldBase:

    ctype = ElemType

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def n_elems(self):
        return c.get_n_elems_sfield(self)

    def total_site(self):
        return c.get_total_site_sfield(self)

    def multiplicity(self):
        return c.get_multiplicity_sfield(self)

    def geo(self):
        geo = Geometry((0, 0, 0, 0))
        c.set_geo_sfield(geo, self)
        return geo

    def __iadd__(self, f1):
        assert isinstance(f1, SelectedFieldBase)
        assert f1.ctype is self.ctype
        c.set_add_sfield(self, f1)
        return self

    def __isub__(self, f1):
        assert isinstance(f1, SelectedFieldBase)
        assert f1.ctype is self.ctype
        c.set_sub_sfield(self, f1)
        return self

    def __imul__(self, factor):
        assert isinstance(factor, float)
        c.set_mul_double_sfield(self, factor)
        return self

    def set_rand(self, rng, upper=1.0, lower=0.0):
        assert isinstance(rng, RngState)
        if self.ctype in field_ctypes_double:
            c.set_u_rand_double_sfield(self, rng, upper, lower)
        else:
            assert False

    def qnorm(self):
        return c.qnorm_sfield(self)

    def get_elems(self, idx):
        return np.array(c.get_elems_sfield(self, idx))

    def get_elem(self, idx, m=None):
        if m is None:
            return np.array(c.get_elem_sfield(self, idx))
        else:
            return np.array(c.get_elem_sfield(self, idx, m))

    def set_elems(self, idx, val):
        """
        val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype=complex).tobytes()
        """
        if isinstance(val, bytes):
            return c.set_elems_sfield(self, idx, val)
        elif isinstance(val, np.ndarray):
            return self.set_elems(idx, val.tobytes())
        else:
            assert False

    def set_elem(self, idx, m, val):
        """
        val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype=complex).tobytes()
        """
        if isinstance(val, bytes):
            return c.set_elem_sfield(self, idx, m, val)
        elif isinstance(val, np.ndarray):
            return self.set_elem(idx, m, val.tobytes())
        else:
            assert False

    def __getitem__(self, i):
        """
        i can be (idx, m,) or idx
        """
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
        """
        i can be (idx, m,) or idx
        val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype=complex).tobytes()
        """
        if isinstance(i, tuple) and len(i) == 2 and isinstance(i[0], int):
            idx, m = i
            self.set_elem(idx, m, val)
        elif isinstance(i, int):
            idx = i
            self.set_elems(idx, val)
        else:
            assert False

    def save_direct(self, path, *args):
        """
        Generic save for SelectedField object
        possible way to call:
        f.save_direct(path) # has some limitations
        f.save_direct(sfw, fn)
        """
        from .fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            assert len(args) == 0
            n_bytes = c.save_sfield(self, path)
        elif isinstance(path, ShuffledFieldsWriter):
            sfw = path
            [fn] = args
            n_bytes = sfw.write(fn, self)
        else:
            raise Exception("SelectedField.save_direct")
        assert n_bytes != 0
        return n_bytes

    def load_direct(self, path, *args):
        """
        Generic load for SelectedField object
        possible way to call:
        f.load_direct(path) # has some limitations
        f.load_direct(sfr, fn)
        if self.fsel is None, self.fsel will be set during f.load_direct(sfr, fn)
        """
        from .fields_io import ShuffledFieldsReader
        if isinstance(path, str):
            n_bytes = c.load_sfield(self, path)
        elif isinstance(path, ShuffledFieldsReader):
            sfr = path
            [fn] = args
            n_bytes = sfr.read(fn, self)
        else:
            raise Exception("SelectedField.load")
        assert n_bytes != 0
        return n_bytes

    def save_64(self, path, *args):
        """
        Generic save for SelectedField object with conversion
        """
        f = self.copy()
        from .fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            f.to_from_endianness("big_64")
        elif isinstance(path, ShuffledFieldsWriter):
            f.to_from_endianness("little_64")
        return f.save_direct(path, *args)

    def save_double(self, path, *args):
        """
        Generic save for SelectedField object with conversion
        """
        return self.save_64(path, *args)

    def save_float_from_double(self, path, *args):
        """
        Generic save for SelectedField object with conversion
        """
        ff = SelectedFieldFloat(self.fsel)
        ff.float_from_double(self)
        from .fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            ff.to_from_endianness("big_32")
        elif isinstance(path, ShuffledFieldsWriter):
            ff.to_from_endianness("little_32")
        return ff.save_direct(path, *args)

    def load_64(self, path, *args):
        """
        Generic load for SelectedField object with conversion
        """
        ret = self.load_direct(path, *args)
        if ret > 0:
            from .fields_io import ShuffledFieldsReader
            if isinstance(path, str):
                self.to_from_endianness("big_64")
            elif isinstance(path, ShuffledFieldsReader):
                self.to_from_endianness("little_64")
        return ret

    def load_double(self, path, *args):
        """
        Generic load for SelectedField object with conversion
        """
        return self.load_64(path, *args)

    def load_double_from_float(self, path, *args):
        """
        Generic load for SelectedField object with conversion
        """
        ff = SelectedField(ElemTypeFloat, self.fsel)
        ret = ff.load_direct(path, *args)
        if ret > 0:
            from .fields_io import ShuffledFieldsReader
            if isinstance(path, str):
                ff.to_from_endianness("big_32")
            elif isinstance(path, ShuffledFieldsReader):
                ff.to_from_endianness("little_32")
            self.double_from_float(ff)
        return ret

    def float_from_double(self, SelectedFieldBase f):
        assert isinstance(self, SelectedFieldFloat)
        self.fsel = f.fsel
        c.convert_float_from_double_sfield(self, f)

    def double_from_float(self, SelectedFieldFloat ff):
        self.fsel = ff.fsel
        c.convert_double_from_float_sfield(self, ff)

    def to_from_endianness(self, tag):
        assert isinstance(tag, str)
        c.to_from_endianness_sfield(self, tag)

    def field_shift(self, shift, is_reflect=False):
        """
        return new shifted SelectedField
        shift is the coordinate to shift the field
        is_reflect determine whether to negate coordinate after shift
        """
        f1 = self.copy(is_copying_data=False)
        f1.fsel = FieldSelection()
        c.field_shift_sfield(f1, self, shift, is_reflect)
        return f1

    def glb_sum_tslice(self, *, t_dir=3):
        """
        return SelectedPoints(self.ctype, get_psel_tslice(self.total_site(), t_dir=t_dir))
        """
        from qlat.field_selection_utils import get_psel_tslice
        psel = get_psel_tslice(self.total_site(), t_dir=t_dir)
        sp = SelectedPoints(self.ctype, psel)
        if self.ctype in field_ctypes_double:
            c.glb_sum_tslice_double_sfield(sp, self, t_dir)
        elif self.ctype in field_ctypes_long:
            c.glb_sum_tslice_long_sfield(sp, self, t_dir)
        else:
            assert False
        return sp

### -------------------------------------------------------------------

cdef class SelectedPointsBase:

    ctype = ElemType

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def n_points(self):
        return c.get_n_points_spfield(self)

    def multiplicity(self):
        return c.get_multiplicity_spfield(self)

    def __iadd__(self, f1):
        assert isinstance(f1, SelectedPointsBase) and f1.ctype is self.ctype
        c.set_add_spfield(self, f1)
        return self

    def __isub__(self, f1):
        assert isinstance(f1, SelectedPointsBase) and f1.ctype is self.ctype
        c.set_sub_spfield(self, f1)
        return self

    def __imul__(self, factor):
        assert isinstance(factor, float)
        c.set_mul_double_spfield(self, factor)
        return self

    def qnorm(self):
        return c.qnorm_spfield(self)

    def get_elems(self, idx):
        return np.array(c.get_elems_spfield(self, idx))

    def get_elem(self, idx, m=None):
        if m is None:
            return np.array(c.get_elem_spfield(self, idx))
        else:
            return np.array(c.get_elem_spfield(self, idx, m))

    def set_elems(self, idx, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype=complex).tobytes()
        if isinstance(val, bytes):
            return c.set_elems_spfield(self, idx, val)
        elif isinstance(val, np.ndarray):
            return self.set_elems(idx, val.tobytes())
        else:
            assert False

    def set_elem(self, idx, m, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype=complex).tobytes()
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
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype=complex).tobytes()
        if isinstance(i, tuple) and len(i) == 2 and isinstance(i[0], int):
            idx, m = i
            self.set_elem(idx, m, val)
        elif isinstance(i, int):
            idx = i
            self.set_elems(idx, val)
        else:
            assert False

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
        return np.asarray(self).copy()

    def from_numpy(self, arr):
        """
        need to be already initialized with ctype and psel
        arr.shape[0] == n_points
        """
        n_points = self.n_points()
        assert arr.shape[0] == n_points
        np.asarray(self).ravel()[:] = arr.ravel()

### -------------------------------------------------------------------
