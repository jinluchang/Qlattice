# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .geometry cimport Geometry
from .fields_io cimport (
        ShuffledFieldsReader,
        ShuffledFieldsWriter,
        )
from .field_selection cimport (
        PointsSelection,
        FieldSelection,
        )
from .field_types cimport (
        FieldChar,
        FieldRealD,
        FieldRealF,
        FieldComplexD,
        FieldComplexF,
        )
from .selected_field_types cimport (
        SelectedFieldChar,
        SelectedFieldRealD,
        SelectedFieldRealF,
        SelectedFieldComplexD,
        SelectedFieldComplexF,
        )
from .selected_points_types cimport (
        SelectedPointsChar,
        SelectedPointsRealD,
        SelectedPointsRealF,
        SelectedPointsComplexD,
        SelectedPointsComplexF,
        SelectedPointsLong,
        SelectedPointsChar,
        )

from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import cqlat as c
import qlat_utils as q
import numpy as np

from .mpi import glb_sum
from .field_type_dict import (
        field_type_dict,
        selected_field_type_dict,
        selected_points_type_dict,
        field_ctypes_complex,
        field_ctypes_complex_f,
        field_ctypes_double,
        field_ctypes_float,
        field_ctypes_long,
        field_ctypes_char,
        )

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

    @q.timer
    def cast_from(self, FieldBase other):
        """
        other can be Field but of different type
        """
        cdef cc.Long size_per_site = other.multiplicity * other.sizeof_m
        cdef cc.Long mult = size_per_site // self.sizeof_m
        assert mult * self.sizeof_m == size_per_site
        self.__init__(other.geo, mult)
        self[:].ravel().view(dtype=np.int8)[:] = other[:].ravel().view(dtype=np.int8)

    @q.timer
    def get_data_sig(self, RngState rng):
        """
        get a signature of the real_d or complex_d field
        """
        cdef FieldComplexD fc
        cdef FieldComplexF fcf
        cdef FieldRealD fr
        cdef FieldRealF frf
        cdef FieldRealD fu
        if self.ctype in field_ctypes_complex:
            fc = FieldComplexD()
            fc.cast_from(self)
            fu = FieldRealD(fc.geo, fc.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (fc[:] * fu[:]).sum()
        elif self.ctype in field_ctypes_complex_f:
            fcf = FieldComplexF()
            fcf.cast_from(self)
            fu = FieldRealD(fcf.geo, fcf.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (fcf[:] * fu[:]).sum()
        elif self.ctype in field_ctypes_double:
            fr = FieldRealD()
            fr.cast_from(self)
            fu = FieldRealD(fr.geo, fr.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (fr[:] * fu[:]).sum()
        elif self.ctype in field_ctypes_float:
            frf = FieldRealF()
            frf.cast_from(self)
            fu = FieldRealD(frf.geo, frf.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (frf[:] * fu[:]).sum()
        else:
            raise Exception(f"get_data_sig: {self.ctype}")
        return glb_sum(sig)

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
                f1n *= -1.0
                c.acc_field_sfield(self, f1n)
            elif isinstance(f1, SelectedPointsBase):
                assert f1.ctype is self.ctype
                f1n = f1.copy()
                f1n *= -1.0
                c.acc_field_spfield(self, f1n)
            else:
                raise Exception(f"Field += type mismatch {type(self)} {type(f1)}")
            assert False
        return self

    def __imul__(self, factor):
        """
		factor can be int, float, complex, FieldM<ComplexD,1>
        """
        if isinstance(factor, (int, float,)):
            c.set_mul_double_field(self, float(factor))
        elif isinstance(factor, complex):
            c.set_mul_complex_field(self, factor)
        elif isinstance(factor, FieldBase):
            assert factor.ctype in [ ElemTypeComplexD, ElemTypeRealD, ]
            c.set_mul_cfield_field(self, factor)
        else:
            assert False
        return self

    def crc32(self):
        return c.crc32_field(self)

    def save_direct(self, path, *args, **kwargs):
        """
        Generic save for Field object
        save Field directly (without any conversion of endianness or precision)
        possible way to call:
        f.save_direct(path)
        f.save_direct(sfw, fn)
        """
        cdef cc.Long n_bytes
        if isinstance(path, str):
            assert len(args) == 0
            n_bytes = self.write_direct(path, **kwargs)
        elif isinstance(path, ShuffledFieldsWriter):
            sfw = path
            fn, = args
            n_bytes = self.write_sfw_direct(sfw, fn, **kwargs)
        else:
            raise Exception("Field.save_direct")
        if n_bytes == 0:
            q.displayln_info(f"WARNING: Field.save_direct({path},*{args},**{kwargs}) n_bytes=0")
        return n_bytes

    def load_direct(self, path, *args, **kwargs):
        """
        Generic load for Field object
        load Field directly (without any conversion of endianness or precision)
        Field geo and multiplicity will be determined during loading
        possible way to call:
        f.load_direct(path)
        f.load_direct(sfr, fn)
        """
        cdef cc.Long n_bytes
        if isinstance(path, str):
            assert len(args) == 0
            n_bytes = self.read_direct(path, **kwargs)
        elif isinstance(path, ShuffledFieldsReader):
            sfr = path
            fn, = args
            n_bytes = self.read_sfr_direct(sfr, fn, **kwargs)
        else:
            raise Exception("Field.load_direct")
        if n_bytes == 0:
            q.displayln_info(f"WARNING: Field.load_direct({path},*{args},**{kwargs}) n_bytes=0")
        return n_bytes

    def save_64(self, path, *args, **kwargs):
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
        return f.save_direct(path, *args, **kwargs)

    def save_double(self, path, *args, **kwargs):
        """
        Generic save for double element Field object
        save double Field as double (do conversion of endianness)
        """
        return self.save_64(path, *args, **kwargs)

    def save_float_from_double(self, path, *args, **kwargs):
        """
        Generic save for double element Field object
        save double Field as float (do conversion of endianness and precision)
        """
        ff = FieldRealF()
        ff.float_from_double(self)
        from .fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            ff.to_from_endianness("big_32")
        elif isinstance(path, ShuffledFieldsWriter):
            ff.to_from_endianness("little_32")
        return ff.save_direct(path, *args, **kwargs)

    def load_64(self, path, *args, **kwargs):
        """
        Generic load for 64-bit size element Field object
        load 64-bit Field (do conversion of endianness)
        """
        ret = self.load_direct(path, *args, **kwargs)
        if ret > 0:
            from .fields_io import ShuffledFieldsReader
            if isinstance(path, str):
                self.to_from_endianness("big_64")
            elif isinstance(path, ShuffledFieldsReader):
                self.to_from_endianness("little_64")
        return ret

    def load_double(self, path, *args, **kwargs):
        """
        Generic load for double Field object
        load double Field (do conversion of endianness)
        """
        return self.load_64(path, *args, **kwargs)

    def load_double_from_float(self, path, *args, **kwargs):
        """
        Generic load for double Field object
        load double Field from float(do conversion of endianness or precision)
        """
        ff = FieldRealF()
        ret = ff.load_direct(path, *args, **kwargs)
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
        self needs to be FieldRealF
        """
        assert isinstance(self, FieldRealF)
        c.convert_float_from_double_field(self, f)

    def double_from_float(self, FieldRealF ff):
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

    def as_field(self, ctype=ElemTypeComplexD):
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

    def __setitem__(self, idx, val):
        """
        Implemented in terms of ``np.asarray``
        """
        np.asarray(self)[idx] = val

    def __getitem__(self, idx):
        """
        Implemented in terms of ``np.asarray``
        """
        return np.asarray(self)[idx]

    def get_elems(self, idx):
        return self[idx]

    def get_elem(self, idx, m=0):
        return self[idx, m]

    def set_elems(self, idx, val):
        """
        val should be np.ndarray. e.g. np.array([1, 2, 3], dtype=complex)
        """
        self[idx] = val

    def set_elem(self, idx, m, val):
        """
        val should be np.ndarray. e.g. np.array([1, 2, 3], dtype=complex)
        """
        self[idx, m] = val

    @q.timer
    def set_m(self, FieldBase f1, cc.Int m, cc.Int m1):
        """
        Set components `m` from `f1`'s components `m1`.
        """
        cdef FieldChar fc = FieldChar()
        cdef FieldChar f1c = FieldChar()
        self.swap_cast(fc)
        f1.swap_cast(f1c)
        cdef cc.Int sizeof_m = self.ctype.sizeof_m
        cc.set_field_m(fc.xx, f1c.xx, m, m1, sizeof_m)
        self.swap_cast(fc)
        f1.swap_cast(f1c)

    def __getnewargs__(self):
        return ()

    def __getstate__(self):
        """
        Only work when single node (or if all nodes has the same data).
        """
        geo = self.geo
        multiplicity = self.multiplicity
        data_arr = self[:]
        return [ data_arr, geo, multiplicity, ]

    def __setstate__(self, state):
        """
        Only work when single node (or if all nodes has the same data).
        """
        if self.view_count > 0:
            raise ValueError("can't load while being viewed")
        self.__init__()
        [ data_arr, geo, multiplicity, ] = state
        self.init_from_geo(geo, multiplicity)
        self[:] = data_arr

    def __len__(self):
        return self.n_sites

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

    @q.timer
    def cast_from(self, SelectedFieldBase other):
        """
        other can be SelectedFieldBase but of different type
        """
        cdef cc.Long size_per_site = other.multiplicity * other.sizeof_m
        cdef cc.Long mult = size_per_site // self.sizeof_m
        assert mult * self.sizeof_m == size_per_site
        self.__init__(other.fsel, mult)
        self[:].ravel().view(dtype=np.int8)[:] = other[:].ravel().view(dtype=np.int8)

    @q.timer
    def get_data_sig(self, RngState rng):
        """
        get a signature of the real_d or complex_d field
        """
        cdef SelectedFieldComplexD fc
        cdef SelectedFieldComplexF fcf
        cdef SelectedFieldRealD fr
        cdef SelectedFieldRealF frf
        cdef SelectedFieldRealD fu
        if self.ctype in field_ctypes_complex:
            fc = SelectedFieldComplexD()
            fc.cast_from(self)
            fu = SelectedFieldRealD(fc.fsel, fc.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (fc[:] * fu[:]).sum()
        elif self.ctype in field_ctypes_complex_f:
            fcf = SelectedFieldComplexF()
            fcf.cast_from(self)
            fu = SelectedFieldRealD(fcf.fsel, fcf.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (fcf[:] * fu[:]).sum()
        elif self.ctype in field_ctypes_double:
            fr = SelectedFieldRealD()
            fr.cast_from(self)
            fu = SelectedFieldRealD(fr.fsel, fr.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (fr[:] * fu[:]).sum()
        elif self.ctype in field_ctypes_float:
            frf = SelectedFieldRealF()
            frf.cast_from(self)
            fu = SelectedFieldRealD(frf.fsel, frf.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (frf[:] * fu[:]).sum()
        else:
            raise Exception(f"get_data_sig: {self.ctype}")
        return glb_sum(sig)

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
        assert isinstance(factor, (int, float))
        c.set_mul_double_sfield(self, float(factor))
        return self

    def __setitem__(self, idx, val):
        """
        Implemented in terms of ``np.asarray``
        """
        np.asarray(self)[idx] = val

    def __getitem__(self, idx):
        """
        Implemented in terms of ``np.asarray``
        """
        return np.asarray(self)[idx]

    def get_elems(self, idx):
        return self[idx]

    def get_elem(self, idx, m = 0):
        return self[idx, m]

    def set_elems(self, idx, val):
        """
        val should be np.ndarray. e.g. np.array([1, 2, 3], dtype=complex)
        """
        self[idx] = val

    def set_elem(self, idx, m, val):
        """
        val should be np.ndarray. e.g. np.array([1, 2, 3], dtype=complex)
        """
        self[idx, m] = val

    def save_direct(self, path, *args, **kwargs):
        """
        Generic save for SelectedField object
        possible way to call:
        f.save_direct(path) # has some limitations
        f.save_direct(sfw, fn)
        """
        cdef cc.Long n_bytes
        if isinstance(path, str):
            assert len(args) == 0
            n_bytes = self.write_direct(path, **kwargs)
        elif isinstance(path, ShuffledFieldsWriter):
            sfw = path
            fn, = args
            n_bytes = self.write_sfw_direct(sfw, fn, **kwargs)
        else:
            raise Exception("SelectedField.save_direct")
        if n_bytes == 0:
            q.displayln_info(f"WARNING: SelectedField.load_direct({path},*{args},**{kwargs}) n_bytes=0")
        return n_bytes

    def load_direct(self, path, *args, **kwargs):
        """
        Generic load for SelectedField object
        possible way to call:
        f.load_direct(path) # has some limitations
        f.load_direct(sfr, fn)
        if self.fsel is None, self.fsel will be set during f.load_direct(sfr, fn)
        """
        cdef cc.Long n_bytes
        if isinstance(path, str):
            assert len(args) == 0
            n_bytes = self.read_direct(path, **kwargs)
        elif isinstance(path, ShuffledFieldsReader):
            sfr = path
            fn, = args
            n_bytes = self.read_sfr_direct(sfr, fn, **kwargs)
        else:
            raise Exception("SelectedField.load_direct")
        if n_bytes == 0:
            q.displayln_info(f"WARNING: SelectedField.load_direct({path},*{args},**{kwargs}) n_bytes=0")
        return n_bytes

    def save_64(self, path, *args, **kwargs):
        """
        Generic save for SelectedField object with conversion
        """
        f = self.copy()
        from .fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            f.to_from_endianness("big_64")
        elif isinstance(path, ShuffledFieldsWriter):
            f.to_from_endianness("little_64")
        return f.save_direct(path, *args, **kwargs)

    def save_double(self, path, *args, **kwargs):
        """
        Generic save for SelectedField object with conversion
        """
        return self.save_64(path, *args, **kwargs)

    def save_float_from_double(self, path, *args, **kwargs):
        """
        Generic save for SelectedField object with conversion
        """
        ff = SelectedFieldRealF(self.fsel)
        ff.float_from_double(self)
        from .fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            ff.to_from_endianness("big_32")
        elif isinstance(path, ShuffledFieldsWriter):
            ff.to_from_endianness("little_32")
        return ff.save_direct(path, *args, **kwargs)

    def load_64(self, path, *args, **kwargs):
        """
        Generic load for SelectedField object with conversion
        """
        ret = self.load_direct(path, *args, **kwargs)
        if ret > 0:
            from .fields_io import ShuffledFieldsReader
            if isinstance(path, str):
                self.to_from_endianness("big_64")
            elif isinstance(path, ShuffledFieldsReader):
                self.to_from_endianness("little_64")
        return ret

    def load_double(self, path, *args, **kwargs):
        """
        Generic load for SelectedField object with conversion
        """
        return self.load_64(path, *args, **kwargs)

    def load_double_from_float(self, path, *args, **kwargs):
        """
        Generic load for SelectedField object with conversion
        """
        ff = SelectedField(ElemTypeRealF, self.fsel)
        ret = ff.load_direct(path, *args, **kwargs)
        if ret > 0:
            from .fields_io import ShuffledFieldsReader
            if isinstance(path, str):
                ff.to_from_endianness("big_32")
            elif isinstance(path, ShuffledFieldsReader):
                ff.to_from_endianness("little_32")
            self.double_from_float(ff)
        return ret

    def float_from_double(self, SelectedFieldBase f):
        assert isinstance(self, SelectedFieldRealF)
        self.fsel = f.fsel
        c.convert_float_from_double_sfield(self, f)

    def double_from_float(self, SelectedFieldRealF ff):
        self.fsel = ff.fsel
        c.convert_double_from_float_sfield(self, ff)

    def to_from_endianness(self, tag):
        assert isinstance(tag, str)
        c.to_from_endianness_sfield(self, tag)

    def glb_sum_tslice(self, *, t_dir=3):
        """
        return SelectedPoints(self.ctype, get_psel_tslice(self.total_site, t_dir=t_dir))
        """
        from .c import get_psel_tslice
        cdef PointsSelection psel = get_psel_tslice(self.total_site, t_dir=t_dir)
        sp = SelectedPoints(self.ctype, psel)
        if self.ctype in field_ctypes_double:
            c.glb_sum_tslice_double_sfield(sp, self, t_dir)
        elif self.ctype in field_ctypes_long:
            c.glb_sum_tslice_long_sfield(sp, self, t_dir)
        else:
            assert False
        return sp

    def __getnewargs__(self):
        return ()

    def __getstate__(self):
        """
        Only work when single node (or if all nodes has the same data).
        """
        fsel = self.fsel
        multiplicity = self.multiplicity
        data_arr = self[:]
        return [ data_arr, multiplicity, fsel, ]

    def __setstate__(self, state):
        """
        Only work when single node (or if all nodes has the same data).
        """
        if self.view_count > 0:
            raise ValueError("can't load while being viewed")
        self.__init__()
        [ data_arr, multiplicity, fsel, ] = state
        self.init_from_fsel(fsel, multiplicity)
        self[:] = data_arr

    def __len__(self):
        return self.n_elems

### -------------------------------------------------------------------

cdef class SelectedPointsBase:

    ctype = ElemType

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    @q.timer
    def cast_from(self, SelectedPointsBase other):
        """
        other can be SelectedPointsBase but of different type
        """
        cdef cc.Long size_per_site = other.multiplicity * other.sizeof_m
        cdef cc.Long mult = size_per_site // self.sizeof_m
        assert mult * self.sizeof_m == size_per_site
        self.__init__(other.psel, mult)
        self[:].ravel().view(dtype=np.int8)[:] = other[:].ravel().view(dtype=np.int8)

    @q.timer
    def get_data_sig(self, RngState rng):
        """
        get a signature of the real_d or complex_d field
        """
        cdef SelectedPointsComplexD fc
        cdef SelectedPointsComplexF fcf
        cdef SelectedPointsRealD fr
        cdef SelectedPointsRealF frf
        cdef SelectedPointsLong fl
        cdef SelectedPointsChar fch
        cdef SelectedPointsRealD fu
        if self.ctype in field_ctypes_complex:
            fc = SelectedPointsComplexD()
            fc.cast_from(self)
            fu = SelectedPointsRealD(fc.psel, fc.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (fc[:] * fu[:]).sum()
        elif self.ctype in field_ctypes_complex_f:
            fcf = SelectedPointsComplexF()
            fcf.cast_from(self)
            fu = SelectedPointsRealD(fcf.psel, fcf.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (fcf[:] * fu[:]).sum()
        elif self.ctype in field_ctypes_double:
            fr = SelectedPointsRealD()
            fr.cast_from(self)
            fu = SelectedPointsRealD(fr.psel, fr.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (fr[:] * fu[:]).sum()
        elif self.ctype in field_ctypes_float:
            frf = SelectedPointsRealF()
            frf.cast_from(self)
            fu = SelectedPointsRealD(frf.psel, frf.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (frf[:] * fu[:]).sum()
        elif self.ctype in field_ctypes_long:
            fl = SelectedPointsLong()
            fl.cast_from(self)
            fu = SelectedPointsRealD(fl.psel, fl.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (fl[:] * fu[:]).sum()
        elif self.ctype in field_ctypes_char:
            fch = SelectedPointsChar()
            fch.cast_from(self)
            fu = SelectedPointsRealD(fch.psel, fch.multiplicity)
            fu.set_rand(rng, 1.0, -1.0)
            sig = (fch[:] * fu[:]).sum()
        else:
            raise Exception(f"get_data_sig: {self.ctype}")
        sig = glb_sum(sig)
        if self.points_dist_type == "g":
            return sig / self.geo.num_node
        else:
            return sig

    def __setitem__(self, idx, val):
        """
        Implemented in terms of ``np.asarray``
        """
        np.asarray(self)[idx] = val

    def __getitem__(self, idx):
        """
        Implemented in terms of ``np.asarray``
        """
        return np.asarray(self)[idx]

    def get_elems(self, idx):
        return self[idx]

    def get_elem(self, idx, m = 0):
        return self[idx, m]

    def set_elems(self, idx, val):
        """
        val should be np.ndarray. e.g. np.array([1, 2, 3], dtype=complex)
        """
        self[idx] = val

    def set_elem(self, idx, m, val):
        """
        val should be np.ndarray. e.g. np.array([1, 2, 3], dtype=complex)
        """
        self[idx, m] = val

    def save_str(self):
        return self.to_lat_data().save_str()

    def load_str(self, bytes content):
        cdef LatData ld = LatData()
        ld.load_str(content)
        self.from_lat_data(ld)

    def to_numpy(self):
        return np.asarray(self).copy()

    def from_numpy(self, arr):
        """
        need to be already initialized with ctype and psel
        arr.shape[0] == n_points
        """
        v_arr = np.asarray(self)
        assert arr.shape[0] == v_arr.shape[0]
        v_arr.ravel()[:] = arr.ravel()

    def __getnewargs__(self):
        return ()

    def __getstate__(self):
        """
        Only work when single node (or if all nodes has the same data).
        """
        psel = self.psel
        n_points = self.n_points
        multiplicity = self.multiplicity
        points_dist_type = self.points_dist_type
        data_arr = self[:]
        return [ data_arr, n_points, multiplicity, points_dist_type, psel, ]

    def __setstate__(self, state):
        """
        Only work when single node (or if all nodes has the same data).
        """
        if self.view_count > 0:
            raise ValueError("can't load while being viewed")
        self.__init__()
        cdef cc.Long n_points
        cdef cc.Int multiplicity
        [ data_arr, n_points, multiplicity, points_dist_type, psel, ] = state
        self.init_from_n_points(n_points, multiplicity, points_dist_type)
        self.psel = psel
        self[:] = data_arr

    def __len__(self):
        return self.n_points

### -------------------------------------------------------------------
