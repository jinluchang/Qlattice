# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.field_base``
==========================\n
Base classes and factory functions for lattice field types, selected fields,
and selected points — the core data containers in qlat.\n
Documentation: ``docs/qlat/qlat_field_base.md``\n
.. note:: Update the documentation when updating this source file.
"""

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
    """
    SelectedField(ctype, fsel) with the default multiplicity == 0 creates an
    *empty*, uninitialized field that keeps fsel: it is meant to be filled
    later, e.g. with load_double / float_from_double (see examples-py
    selected-convert-io.py).  Pass a positive multiplicity to allocate now.
    """
    assert ctype in field_type_dict
    FieldType = selected_field_type_dict[ctype]
    field = FieldType(fsel, multiplicity)
    return field

def SelectedPoints(type ctype, PointsSelection psel, int multiplicity=0):
    """
    SelectedPoints(ctype, psel) with the default multiplicity == 0 creates an
    *empty*, uninitialized field that keeps psel; pass a positive multiplicity
    to allocate now.
    """
    assert ctype in field_type_dict
    FieldType = selected_points_type_dict[ctype]
    field = FieldType(psel, multiplicity)
    return field

def field_check_key(idx):
    """
    Validate a field index and return the NumPy index to use.

    Field buffers are ``(local_volume, multiplicity, *elem_shape)``,
    C-contiguous, indexed by the **flat local site index** with the first
    coordinate varying fastest, matching ``geo.coordinate_from_index``.

    A ``Coordinate`` (or a tuple/list containing one) is a common mistake
    because it looks like the C++ ``get_elem`` API; reject it with an
    actionable message instead of silently doing the wrong thing.
    """
    if isinstance(idx, (q.Coordinate, q.CoordinateD)):
        raise TypeError(
            f"field indices are flat local site indices, not {type(idx).__name__}"
            "; use get_elem_xg(xg, m) for global coordinates, or "
            "geo.index_from_coordinate(xl) to convert a local coordinate to a "
            "flat local index")
    if isinstance(idx, (tuple, list)):
        for key in idx:
            if isinstance(key, (q.Coordinate, q.CoordinateD)):
                raise TypeError(
                    "field indices are flat local site indices; a tuple "
                    "containing a Coordinate would be interpreted by NumPy as a "
                    "fancy index over the site axis. Use "
                    "get_elem_xg(xg, m) for global coordinates, or "
                    "geo.index_from_coordinate(xl)")
    return idx

### -------------------------------------------------------------------

cdef class FieldBase:
    """
    Base class of lattice fields defined on every site of a ``Geometry``.
    #
    The NumPy buffer interface (``np.asarray(f)``, equivalently ``f[:]``)
    exposes a writeable, zero-copy view with shape
    ``(local_volume, multiplicity, *elem_shape)``, C-contiguous. Axis 0 is the
    **flat local site index**, with the first coordinate varying fastest
    (matching ``geo.coordinate_from_index``); a ``Coordinate`` is not an
    accepted index. See ``docs/qlat/qlat_field_indexing.md``.
    """

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
        return memoryview(np.asarray(self).reshape(-1))

    def __iadd__(self, f1):
        """
        f1 can be Field, SelectedField, SelectedPoints
        """
        if isinstance(f1, FieldBase):
            assert f1.ctype is self.ctype
            self._cc_iadd(f1)
        else:
            if isinstance(f1, SelectedFieldBase):
                self._cc_acc_field_sfield(f1, f1.fsel)
            elif isinstance(f1, SelectedPointsBase):
                assert f1.ctype is self.ctype
                self._cc_acc_field_spfield(f1, f1.psel.geo, f1.psel)
            else:
                raise Exception(f"Field += type mismatch {type(self)} {type(f1)}")
        return self

    def __isub__(self, f1):
        """
        f1 can be Field, SelectedField, SelectedPoints
        """
        if isinstance(f1, FieldBase):
            assert f1.ctype is self.ctype
            self._cc_isub(f1)
        else:
            if isinstance(f1, SelectedFieldBase):
                assert f1.ctype is self.ctype
                f1n = f1.copy()
                f1n *= -1.0
                self._cc_acc_field_sfield(f1n, f1n.fsel)
            elif isinstance(f1, SelectedPointsBase):
                assert f1.ctype is self.ctype
                f1n = f1.copy()
                f1n *= -1.0
                self._cc_acc_field_spfield(f1n, f1n.psel.geo, f1n.psel)
            else:
                raise Exception(f"Field += type mismatch {type(self)} {type(f1)}")
        return self

    def __imul__(self, factor):
        """
		factor can be int, float, complex, FieldM<ComplexD,1>
        """
        if isinstance(factor, (int, float,)):
            self._cc_imul_double(float(factor))
        elif isinstance(factor, complex):
            if self.ctype in field_ctypes_complex:
                self._cc_imul_complex(factor)
            elif factor.imag == 0.0:
                # a real valued factor is well defined for any ctype
                self._cc_imul_double(float(factor.real))
            else:
                raise ValueError(
                        f"Field.__imul__: cannot multiply {self.ctype}"
                        f" by the complex factor {factor}"
                        )
        elif isinstance(factor, FieldBase):
            assert factor.ctype in [ ElemTypeComplexD, ElemTypeRealD, ]
            if factor.ctype is ElemTypeComplexD:
                self._cc_imul_complex_field(factor)
            else:
                self._cc_imul_real_field(factor)
        else:
            assert False
        return self

    def crc32(self):
        return self._cc_crc32()

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
        f._cc_to_realf(self)

    def double_from_float(self, FieldRealF ff):
        """
        self can be any FieldBase subtype but need to be actually contains double precision numbers
        """
        self._cc_convert_double_from_float(ff)

    def to_from_endianness(self, tag):
        """
        Convert between the native endianness and the endianness specified by ``tag``
        tag can be ``"big_32", "big_64", "little_32", "little_64"``
        """
        assert isinstance(tag, str)
        self._cc_to_from_endianness(tag)

    def as_field(self, ctype=ElemTypeComplexD):
        """
		return new Field(ctype) with the same content
        """
        f = Field(ctype)
        f.cast_from(self)
        return f

    def from_field(self, f):
        """
		assign from f with the same content but possibly different type
        """
        self.cast_from(f)
        return f

    def __setitem__(self, idx, val):
        """
        Implemented in terms of ``np.asarray``.
        #
        The index is a NumPy index into the buffer
        ``(local_volume, multiplicity, *elem_shape)``, which is C-contiguous
        with the flat local site index as axis 0; see the class docstring.
        #
        .. note:: For repeated calls, use ``np.asarray(f)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        cdef object arr = np.asarray(self)
        arr[field_check_key(idx)] = val

    def __getitem__(self, idx):
        """
        Implemented in terms of ``np.asarray``.
        #
        The index is a NumPy index into the buffer
        ``(local_volume, multiplicity, *elem_shape)``, which is C-contiguous
        with the flat local site index as axis 0; see the class docstring.
        #
        .. note:: For repeated calls, use ``np.asarray(f)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        cdef object arr = np.asarray(self)
        return arr[field_check_key(idx)]

    def get_elems(self, idx):
        """
        .. note:: For repeated calls, use ``np.asarray(f)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        return self[idx]

    def get_elem(self, idx, m=0):
        """
        .. note:: For repeated calls, use ``np.asarray(f)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        return self[idx, m]

    def set_elems(self, idx, val):
        """
        val should be np.ndarray. e.g. np.array([1, 2, 3], dtype=complex)
        #
        .. note:: For repeated calls, use ``np.asarray(f)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        self[idx] = val

    def set_elem(self, idx, m, val):
        """
        val should be np.ndarray. e.g. np.array([1, 2, 3], dtype=complex)
        #
        .. note:: For repeated calls, use ``np.asarray(f)`` once and index the
           resulting array directly, to avoid creating a new view each time.
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
    f._cc_split_fields(fs)

def merge_fields(f, fs):
    nf = len(fs)
    assert nf >= 1
    assert isinstance(f, FieldBase)
    assert f.ctype is fs[0].ctype
    f._cc_merge_fields(fs)

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
    f._cc_merge_fields_ms(fs, ms)

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
        self._cc_iadd(f1)
        return self

    def __isub__(self, f1):
        assert isinstance(f1, SelectedFieldBase)
        assert f1.ctype is self.ctype
        self._cc_isub(f1)
        return self

    def __imul__(self, factor):
        assert isinstance(factor, (int, float))
        self._cc_imul_double(float(factor))
        return self

    def __setitem__(self, idx, val):
        """
        Implemented in terms of ``np.asarray``.
        #
        .. note:: For repeated calls, use ``np.asarray(sf)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        np.asarray(self)[idx] = val

    def __getitem__(self, idx):
        """
        Implemented in terms of ``np.asarray``.
        #
        .. note:: For repeated calls, use ``np.asarray(sf)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        return np.asarray(self)[idx]

    def get_elems(self, idx):
        """
        .. note:: For repeated calls, use ``np.asarray(sf)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        return self[idx]

    def get_elem(self, idx, m = 0):
        """
        .. note:: For repeated calls, use ``np.asarray(sf)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        return self[idx, m]

    def set_elems(self, idx, val):
        """
        val should be np.ndarray. e.g. np.array([1, 2, 3], dtype=complex)
        #
        .. note:: For repeated calls, use ``np.asarray(sf)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        self[idx] = val

    def set_elem(self, idx, m, val):
        """
        val should be np.ndarray. e.g. np.array([1, 2, 3], dtype=complex)
        #
        .. note:: For repeated calls, use ``np.asarray(sf)`` once and index the
           resulting array directly, to avoid creating a new view each time.
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
        f._cc_to_realf(self)

    def double_from_float(self, SelectedFieldRealF ff):
        self.fsel = ff.fsel
        self._cc_convert_double_from_float(ff)

    def to_from_endianness(self, tag):
        assert isinstance(tag, str)
        self._cc_to_from_endianness(tag)

    def glb_sum_tslice(self, *, t_dir=3):
        """
        return SelectedPoints(self.ctype, get_psel_tslice(self.total_site, t_dir=t_dir))
        """
        from .c import get_psel_tslice
        cdef PointsSelection psel = get_psel_tslice(self.total_site, t_dir=t_dir)
        sp = SelectedPoints(self.ctype, psel)
        if self.ctype in field_ctypes_double:
            self._cc_glb_sum_tslice(sp, self.fsel, t_dir)
        elif self.ctype in field_ctypes_long:
            self._cc_glb_sum_tslice(sp, self.fsel, t_dir)
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
        Implemented in terms of ``np.asarray``.
        #
        .. note:: For repeated calls, use ``np.asarray(sp)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        np.asarray(self)[idx] = val

    def __getitem__(self, idx):
        """
        Implemented in terms of ``np.asarray``.
        #
        .. note:: For repeated calls, use ``np.asarray(sp)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        return np.asarray(self)[idx]

    def get_elems(self, idx):
        """
        .. note:: For repeated calls, use ``np.asarray(sp)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        return self[idx]

    def get_elem(self, idx, m = 0):
        """
        .. note:: For repeated calls, use ``np.asarray(sp)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        return self[idx, m]

    def set_elems(self, idx, val):
        """
        val should be np.ndarray. e.g. np.array([1, 2, 3], dtype=complex)
        #
        .. note:: For repeated calls, use ``np.asarray(sp)`` once and index the
           resulting array directly, to avoid creating a new view each time.
        """
        self[idx] = val

    def set_elem(self, idx, m, val):
        """
        val should be np.ndarray. e.g. np.array([1, 2, 3], dtype=complex)
        #
        .. note:: For repeated calls, use ``np.asarray(sp)`` once and index the
           resulting array directly, to avoid creating a new view each time.
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
