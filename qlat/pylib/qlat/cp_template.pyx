# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import cqlat as c
import qlat_utils as q
import numpy as np

### -------------------------------------------------------------------

cdef class Geometry:

    def __cinit__(self):
        self.xx = cc.Geometry()
        self.cdata = <long>&(self.xx)

    def __init__(self, total_site = None, int multiplicity = 1):
        # if total_site is None: create geo uninitialized
        # elif multiplicity is None: create geo with multiplicity = 1
        self.xx.init(Coordinate(total_site).xx, multiplicity)

    def __imatmul__(self, Geometry v1):
        self.xx = v1.xx
        return self

    def copy(self, is_copying_data = True):
        x = Geometry()
        if is_copying_data:
            x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def total_site(self):
        return c.get_total_site_geo(self)

    def total_volume(self):
        return c.get_total_volume_geo(self)

    def local_site(self):
        return c.get_node_site_geo(self)

    def local_volume(self):
        return c.get_local_volume_geo(self)

    def multiplicity(self):
        return c.get_multiplicity_geo(self)

    def eo(self):
        return c.get_eo_geo(self)

    def expansion_left(self):
        return c.get_expansion_left_geo(self)

    def expansion_right(self):
        return c.get_expansion_right_geo(self)

    def id_node(self):
        return c.get_id_node_geo(self)

    def num_node(self):
        return c.get_num_node_geo(self)

    def coor_node(self):
        return c.get_coor_node_geo(self)

    def size_node(self):
        return c.get_size_node_geo(self)

    def __str__(self):
        return self.show()

    def __repr__(self):
        return self.show()

    def show(self):
        total_site = self.total_site()
        multiplicity = self.multiplicity()
        expan_left = self.expansion_left()
        expan_right = self.expansion_right()
        eo = self.eo()
        zero = [ 0, 0, 0, 0, ]
        if expan_left == zero and expan_right == zero and eo == 0:
            return f"Geometry({total_site}, {multiplicity})"
        else:
            return f"Geometry({total_site}, {multiplicity}, expan_left={expan_left}, expan_right={expan_right}, eo={eo})"

    def coordinate_g_from_l(self, xl):
        return c.coordinate_g_from_l_geo(self, xl)

    def coordinate_l_from_g(self, xg):
        return c.coordinate_l_from_g_geo(self, xg)

    def is_local(self, xl):
        return c.is_local_geo(self, xl)

    def is_local_xg(self, xg):
        # return a global coordinate is inside the local volume or not
        return c.is_local_xg_geo(self, xg)

    def xg_list(self):
        # return xg for all local sites
        return c.get_xg_list(self)

### -------------------------------------------------------------------

def geo_reform(Geometry geo, int multiplicity = 1, expansion_left = None, expansion_right = None):
    cdef Coordinate el, er
    if expansion_left is None:
        el = Coordinate()
    elif isinstance(expansion_left, int):
        e = expansion_left
        el = Coordinate([ e, e, e, e, ])
    elif isinstance(expansion_left, list):
        el = Coordinate(expansion_left)
    else:
        assert isinstance(expansion_left, Coordinate)
        el = <Coordinate>expansion_left
    if expansion_right is None:
        er = Coordinate()
    elif isinstance(expansion_right, int):
        e = expansion_right
        er = Coordinate([ e, e, e, e, ])
    elif isinstance(expansion_right, list):
        er = Coordinate(expansion_right)
    else:
        assert isinstance(expansion_left, Coordinate)
        er = <Coordinate>expansion_left
    cdef Geometry geo_new = Geometry()
    geo_new.xx = cc.geo_reform(geo.xx, multiplicity, el.xx, er.xx)
    return geo_new

def geo_eo(Geometry geo, int eo = 0):
    cdef Geometry geo_new = Geometry()
    geo_new.xx = cc.geo_eo(geo.xx, eo)
    return geo_new

### -------------------------------------------------------------------

field_ctypes_complex = [
        ElemTypeColorMatrix,
        ElemTypeWilsonMatrix,
        ElemTypeNonRelWilsonMatrix,
        ElemTypeIsospinMatrix,
        ElemTypeSpinMatrix,
        ElemTypeWilsonVector,
        ElemTypeComplex,
        ]

field_ctypes_double = [
        ElemTypeColorMatrix,
        ElemTypeWilsonMatrix,
        ElemTypeNonRelWilsonMatrix,
        ElemTypeIsospinMatrix,
        ElemTypeSpinMatrix,
        ElemTypeWilsonVector,
        ElemTypeComplex,
        ElemTypeDouble,
        ]

field_ctypes_float = [
        ElemTypeComplexF,
        ElemTypeFloat,
        ]

field_ctypes_long = [
        ElemTypeLong,
        ElemTypeInt64t,
        ]

### -------------------------------------------------------------------

field_type_dict = {}

selected_field_type_dict = {}

selected_points_type_dict = {}

def Field(type ctype, Geometry geo = None, int multiplicity = 0):
    assert ctype in field_type_dict
    FieldType = field_type_dict[ctype]
    field = FieldType(geo, multiplicity)
    return field

def SelectedField(type ctype, FieldSelection fsel, int multiplicity = 0):
    assert ctype in field_type_dict
    FieldType = selected_field_type_dict[ctype]
    field = FieldType(fsel, multiplicity)
    return field

def SelectedPoints(type ctype, PointSelection psel, int multiplicity = 0):
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
        # f1 can be Field, SelectedField, SelectedPoints
        if isinstance(f1, FieldBase):
            assert f1.ctype is self.ctype
            c.set_add_field(self, f1)
        else:
            from qlat.selected_field import SelectedField
            from qlat.selected_points import SelectedPoints
            if isinstance(f1, SelectedFieldBase):
                c.acc_field_sfield(self, f1)
            elif isinstance(f1, SelectedPointsBase):
                assert f1.ctype is self.ctype
                c.acc_field_spfield(self, f1)
            else:
                raise Exception(f"Field += type mismatch {type(self)} {type(f1)}")
        return self

    def __isub__(self, f1):
        # f1 can be Field, SelectedField, SelectedPoints
        if isinstance(f1, FieldBase):
            assert f1.ctype is self.ctype
            c.set_sub_field(self, f1)
        else:
            from qlat.selected_field import SelectedField
            from qlat.selected_points import SelectedPoints
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
		# factor can be float, complex, FieldM<Complex,1>
        if isinstance(factor, float):
            c.set_mul_double_field(self, factor)
        elif isinstance(factor, complex):
            c.set_mul_complex_field(self, factor)
        elif isinstance(factor, FieldBase):
            assert factor.ctype is ElemTypeComplex
            assert factor.multiplicity() == 1
            c.set_mul_cfield_field(self, factor)
        else:
            assert False
        return self

    @q.timer
    def set_unit(self, coef = 1.0):
        c.set_unit_field(self, coef)

    @q.timer
    def set_rand(self, rng, upper = 1.0, lower = 0.0):
        assert isinstance(rng, RngState)
        if self.ctype in field_ctypes_double:
            c.set_u_rand_double_field(self, rng, upper, lower)
        elif self.ctype in field_ctypes_float:
            c.set_u_rand_float_field(self, rng, upper, lower)
        else:
            assert False

    @q.timer
    def set_rand_g(self, rng, center = 0.0, sigma = 1.0):
        assert isinstance(rng, RngState)
        if self.ctype in field_ctypes_double:
            c.set_g_rand_double_field(self, rng, center, sigma)
        else:
            assert False

    def multiply_double(self, factor):
        assert isinstance(factor, FieldBase)
        assert factor.ctype is ElemTypeDouble
        c.multiply_double_field(self,factor)

    def qnorm(self):
        return c.qnorm_field(self)

    def crc32(self):
        return c.crc32_field(self)

    def sparse(self, sel):
        # deprecated
        q.displayln_info("Field.sparse: deprecated")
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
        ff = FieldFloat()
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
        ff = FieldFloat()
        ret = ff.load_direct(path, *args)
        if ret > 0:
            from qlat.fields_io import ShuffledFieldsReader
            if isinstance(path, str):
                ff.to_from_endianness("big_32")
            elif isinstance(path, ShuffledFieldsReader):
                ff.to_from_endianness("little_32")
            self.double_from_float(ff)
        return ret

    def float_from_double(self, FieldBase f):
        assert isinstance(self, FieldFloat)
        c.convert_float_from_double_field(self, f)

    def double_from_float(self, FieldFloat ff):
        c.convert_double_from_float_field(self, ff)

    def to_from_endianness(self, tag):
        assert isinstance(tag, str)
        c.to_from_endianness_field(self, tag)

    def as_field(self, ctype = ElemTypeComplex):
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
            self.set_elem(index, m, val)
        elif isinstance(i, (int, list, tuple)):
            index = i
            self.set_elems(index, val)
        else:
            assert False

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
    # fms = [ (f0, m0,), (f1, m1,), ... ]
    # f.get_elem(x, m) = fms[m][0].get_elem(x, fms[m][1])
    multiplicity = len(fms)
    assert multiplicity >= 1
    assert isinstance(f, FieldBase)
    assert f.ctype is fms[0][0].ctype
    fs, ms = zip(*fms)
    c.merge_fields_ms_field(f, fs, ms)

def mk_merged_fields_ms(fms):
    # fms = [ (f0, m0,), (f1, m1,), ... ]
    # f.get_elem(x, m) = fms[m][0].get_elem(x, fms[m][1])
    # return f
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

cdef class PointSelection:

    def __cinit__(self):
        self.cdata = <long>&(self.xx)
        self.geo = None

    def __init__(self, coordinate_list = None, Geometry geo = None):
        cdef long n_points
        cdef long i
        self.geo = geo
        if coordinate_list is None:
            return
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        n_points = len(coordinate_list)
        self.xx = cc.PointSelection(n_points)
        for i in range(n_points):
            c = coordinate_list[i]
            if not isinstance(c, Coordinate):
                c = Coordinate(c)
            self.xx[i] = (<Coordinate>c).xx

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef long n_points = self.xx.size()
        cdef int ndim = 2
        cdef char* fmt = "i"
        cdef Buffer buf = Buffer(self, ndim, sizeof(cc.Int32t))
        cdef Py_ssize_t* shape = &buf.shape_strides[0]
        cdef Py_ssize_t* strides = &buf.shape_strides[buf.ndim]
        shape[0] = n_points
        shape[1] = 4
        buf.set_strides()
        buffer.buf = <char*>(self.xx.data())
        if flags & PyBUF_FORMAT:
            buffer.format = fmt
        else:
            buffer.format = NULL
        buffer.internal = NULL
        buffer.itemsize = buf.itemsize
        buffer.len = buf.get_len()
        buffer.ndim = buf.ndim
        buffer.obj = buf
        buffer.readonly = 0
        buffer.shape = shape
        buffer.strides = strides
        buffer.suboffsets = NULL
        self.view_count += 1

    def __releasebuffer__(self, Py_buffer *buffer):
        self.view_count -= 1

    def __imatmul__(self, PointSelection v1):
        self.geo = v1.geo
        self.xx = v1.xx
        return self

    def copy(self):
        x = PointSelection()
        x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def set_rand(self, rs, total_site, n_points):
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        c.set_rand_psel(self, rs, total_site, n_points)
        self.geo = Geometry(total_site)

    def save(self, path):
        c.save_psel(self, path)

    def load(self, path, geo):
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        c.load_psel(self, path)
        self.geo = geo

    def n_points(self):
        return c.get_n_points_psel(self)

    def to_list(self):
        return c.mk_list_psel(self)

    def from_list(self, coordinate_list, geo = None):
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        c.set_list_psel(self, coordinate_list)
        self.geo = geo
        return self

    def coordinate_from_idx(self, idx):
        return c.get_coordinate_from_idx_psel(self, idx)

### -------------------------------------------------------------------

cdef class FieldSelection:

    def __cinit__(self):
        self.cdata = <long>&(self.xx)

    def __init__(self, total_site = None, n_per_tslice = -1, rs = None, psel = None):
        if total_site is not None:
            assert isinstance(rs, RngState)
            assert isinstance(n_per_tslice, int)
            c.set_rand_fsel(self, rs, total_site, n_per_tslice)
            if psel is not None:
                c.add_psel_fsel(self, psel)
            self.update()
            self.update(n_per_tslice)

    def __imatmul__(self, FieldSelection v1):
        self.xx = v1.xx
        return self

    def copy(self):
        x = FieldSelection()
        x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def set_uniform(self, total_site, val = 0):
        # default (val = 0) select every sites
        # val = -1 deselection everything
        c.set_uniform_fsel(self, total_site, val)

    def set_rand(self, rs, total_site, n_per_tslice):
        assert isinstance(rs, RngState)
        assert isinstance(n_per_tslice, int)
        c.set_rand_fsel(self, rs, total_site, n_per_tslice)
        self.update()
        self.update(n_per_tslice)

    def add_psel(self, psel, rank_psel = 1024 * 1024 * 1024 * 1024 * 1024):
        # Add psel points to the selection, with the rank specified as rank_psel.
        # If the point is already selected with lower rank, the rank is unchanged.
        c.add_psel_fsel(self, psel, rank_psel)
        self.update()

    def update(self, n_per_tslice = -1):
        # if n_per_tslice < 0: only update various indices
        # if n_per_tslice >= 0: only update parameters (n_per_tslice and prob)
        c.update_fsel(self, n_per_tslice)

    def select_rank_range(self, rank_start = 0, rank_stop = -1):
        # return new fsel with selected points that
        # rank_start <= rank and (rank < rank_stop or rank_stop == -1)
        # Does NOT change the n_per_tslice parameter for the new fsel
        fsel = FieldSelection()
        c.select_rank_range_fsel(fsel, self, rank_start, rank_stop)
        fsel.update()
        fsel.update(self.n_per_tslice())
        return fsel

    def select_t_range(self, rank_start = 0, rank_stop = -1):
        # return new fsel with selected points that
        # t_start <= t and (t < t_stop or t_stop == -1)
        # rank_start <= rank < rank_stop (rank_stop = -1 implies unlimited)
        # Does NOT change the n_per_tslice parameter for the new fsel
        fsel = FieldSelection()
        c.select_rank_range_fsel(fsel, self, rank_start, rank_stop)
        fsel.update()
        fsel.update(self.n_per_tslice())
        return fsel

    def to_psel(self):
        psel = PointSelection(None, self.geo())
        c.set_psel_fsel(psel, self)
        return psel

    def to_psel_local(self):
        psel = PointSelection(None, self.geo())
        c.set_psel_fsel_local(psel, self)
        return psel

    def save(self, path):
        return c.save_fsel(self, path)

    def load(self, path, n_per_tslice):
        return c.load_fsel(self, path, n_per_tslice)

    def geo(self):
        geo = Geometry((0, 0, 0, 0))
        c.set_geo_fsel(geo, self)
        return geo

    def total_site(self):
        return c.get_total_site_fsel(self)

    def n_elems(self):
        return c.get_n_elems_fsel(self)

    def n_per_tslice(self):
        return c.get_n_per_tslice_fsel(self)

    def prob(self):
        # return fsel.prob
        # n_per_tslice / spatial_volume
        return c.get_prob_fsel(self)

    def idx_from_coordinate(self, xg):
        return c.get_idx_from_coordinate_fsel(self, xg)

    def coordinate_from_idx(self, idx):
        return c.get_coordinate_from_idx_fsel(self, idx)

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

    def set_rand(self, rng, upper = 1.0, lower = 0.0):
        assert isinstance(rng, RngState)
        if self.ctype in field_ctypes_double:
            c.set_u_rand_double_sfield(self, rng, upper, lower)
        else:
            assert False

    def qnorm(self):
        return c.qnorm_sfield(self)

    def get_elems(self, idx):
        return np.array(c.get_elems_sfield(self, idx))

    def get_elem(self, idx, m = None):
        if m is None:
            return np.array(c.get_elem_sfield(self, idx))
        else:
            return np.array(c.get_elem_sfield(self, idx, m))

    def set_elems(self, idx, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(val, bytes):
            return c.set_elems_sfield(self, idx, val)
        elif isinstance(val, np.ndarray):
            return self.set_elems(idx, val.tobytes())
        else:
            assert False

    def set_elem(self, idx, m, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(val, bytes):
            return c.set_elem_sfield(self, idx, m, val)
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
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(i, tuple) and len(i) == 2 and isinstance(i[0], int):
            idx, m = i
            self.set_elem(idx, m, val)
        elif isinstance(i, int):
            idx = i
            self.set_elems(idx, val)
        else:
            assert False

    def sparse(self, sel):
        # deprecated
        q.displayln_info("SelectedField.sparse: deprecated")
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
        # possible way to call:
        # f.save_direct(path)
        # f.save_direct(path, new_size_node)
        # f.save_direct(sfw, fn)
        from qlat.fields_io import ShuffledFieldsWriter
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
        # possible way to call:
        # f.load_direct(path)
        # f.load_direct(sfr, fn)
        # if self.fsel is None, self.fsel will be set during f.load_direct(sfr, fn)
        from qlat.fields_io import ShuffledFieldsReader
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
        ff = SelectedFieldFloat(self.fsel)
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
        ff = SelectedField(ElemTypeFloat, self.fsel)
        ret = ff.load_direct(path, *args)
        if ret > 0:
            from qlat.fields_io import ShuffledFieldsReader
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

    def field_shift(self, shift, is_reflect = False):
        # return new shifted SelectedField
        # shift is the coordinate to shift the field
        # is_reflect determine whether to negate coordinate after shift
        f1 = self.copy(is_copying_data = False)
        f1.fsel = FieldSelection()
        c.field_shift_sfield(f1, self, shift, is_reflect)
        return f1

    def glb_sum_tslice(self, *, t_dir = 3):
        # return SelectedPoints(self.ctype, get_psel_tslice(self.total_site(), t_dir = t_dir))
        from qlat.selected_points import SelectedPoints, get_psel_tslice
        psel = get_psel_tslice(self.total_site(), t_dir = t_dir)
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

    def get_elem(self, idx, m = None):
        if m is None:
            return np.array(c.get_elem_spfield(self, idx))
        else:
            return np.array(c.get_elem_spfield(self, idx, m))

    def set_elems(self, idx, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(val, bytes):
            return c.set_elems_spfield(self, idx, val)
        elif isinstance(val, np.ndarray):
            return self.set_elems(idx, val.tobytes())
        else:
            assert False

    def set_elem(self, idx, m, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
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
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
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
        ld = q.LatData()
        c.lat_data_from_complex_spfield(ld, self)
        return ld

    def from_lat_data(self, ld):
        assert self.ctype in field_ctypes_complex
        assert isinstance(ld, q.LatData)
        c.complex_spfield_from_lat_data(self, ld)

    def to_numpy(self):
        n_points = self.n_points()
        return np.array([ self.get_elems(idx) for idx in range(n_points) ])

    def from_numpy(self, arr):
        # need to be already initialized with ctype and psel
        # arr.shape[0] == n_points
        n_points = self.n_points()
        assert arr.shape[0] == n_points
        for idx in range(n_points):
            self.set_elems(idx, arr[idx])

### -------------------------------------------------------------------

def qremove_info(path):
    return cc.qremove_info(path)

def qremove_all_info(path):
    return cc.qremove_all_info(path)

### -------------------------------------------------------------------
