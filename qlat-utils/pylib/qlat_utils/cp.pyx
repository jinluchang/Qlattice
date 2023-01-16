# cython: c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT
cimport numpy

import cqlat_utils as c
import numpy as np
import functools

### -------------------------------------------------------------------

cdef class ElemType:

    name = ""

### -------------------------------------------------------------------

# ColorMatrix WilsonMatrix NonRelWilsonMatrix SpinMatrix WilsonVector Complex ComplexF Double Float Long Int64t Int8t Char

cdef class ElemTypeColorMatrix(ElemType):
    name = "ColorMatrix"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
    @staticmethod
    cdef int ndim():
        return 2
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t](2, 3)
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.ColorMatrix)
cdef class ElemTypeWilsonMatrix(ElemType):
    name = "WilsonMatrix"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
    @staticmethod
    cdef int ndim():
        return 2
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t](2, 12)
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.WilsonMatrix)
cdef class ElemTypeNonRelWilsonMatrix(ElemType):
    name = "NonRelWilsonMatrix"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
    @staticmethod
    cdef int ndim():
        return 2
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t](2, 6)
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.NonRelWilsonMatrix)
cdef class ElemTypeIsospinMatrix(ElemType):
    name = "IsospinMatrix"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
    @staticmethod
    cdef int ndim():
        return 2
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t](2, 2)
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.IsospinMatrix)
cdef class ElemTypeSpinMatrix(ElemType):
    name = "SpinMatrix"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
    @staticmethod
    cdef int ndim():
        return 2
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t](2, 4)
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.SpinMatrix)
cdef class ElemTypeWilsonVector(ElemType):
    name = "WilsonVector"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
    @staticmethod
    cdef int ndim():
        return 1
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t](1, 12)
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.WilsonVector)
cdef class ElemTypeComplex(ElemType):
    name = "Complex"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Complex)
cdef class ElemTypeComplexF(ElemType):
    name = "ComplexF"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zf'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.ComplexF)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.ComplexF)
cdef class ElemTypeDouble(ElemType):
    name = "Double"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Double)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Double)
cdef class ElemTypeFloat(ElemType):
    name = "Float"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'f'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Float)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Float)
cdef class ElemTypeLong(ElemType):
    name = "Long"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'q'
        if not sizeof(cc.Long) == 8:
            assert sizeof(cc.Long) == 4
            fmt = 'l'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Long)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Long)
cdef class ElemTypeInt64t(ElemType):
    name = "Int64t"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'q'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Int64t)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Int64t)
cdef class ElemTypeInt8t(ElemType):
    name = "Int8t"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'b'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Int8t)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Int8t)
cdef class ElemTypeChar(ElemType):
    name = "Char"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'c'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Char)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Char)

### -------------------------------------------------------------------

def displayln(level, *args):
    # interpret the first argument as verbose level if it is int
    if isinstance(level, int):
        if level <= verbose_level():
            print(*args, flush=True)
    else:
        print(level, *args, flush=True)

def displayln_info(*args):
    if cc.get_id_node() == 0:
        displayln(*args)

### -------------------------------------------------------------------

cdef class Timer:

    def __cinit__(self, const cc.std_string& fname, cc.bool is_verbose = False):
        self.xx = cc.Timer(fname)

    def start(self):
        self.xx.start(self.is_verbose)

    def stop(self):
        self.xx.stop(self.is_verbose)

### -------------------------------------------------------------------

cdef class TimerNone:

    def __cinit__(self):
        pass

    def start(self):
        pass

    def stop(self):
        pass

### -------------------------------------------------------------------

def get_id_node():
    return cc.get_id_node()

def get_num_node():
    return cc.get_num_node()

def verbose_level(level = None):
    if level is None:
        return cc.verbose_level()
    cdef long* p_ret = &cc.verbose_level()
    p_ret[0] = level
    assert cc.verbose_level() == level
    return level

def get_time():
    return cc.get_time()

def get_start_time(time = None):
    if time is None:
        return cc.get_start_time()
    cdef double* p_ret = &cc.get_start_time()
    p_ret[0] = time
    assert cc.get_start_time() == time
    return time

def get_actual_start_time(time = None):
    if time is None:
        return cc.get_actual_start_time()
    cdef double* p_ret = &cc.get_actual_start_time()
    p_ret[0] = time
    assert cc.get_actual_start_time() == time
    return time

def get_total_time():
    return cc.get_total_time()

def get_actual_total_time():
    return cc.get_actual_total_time()

### -------------------------------------------------------------------

def timer_display(str tag = ""):
    cc.Timer.display(tag)
    cc.flush()

def timer_autodisplay():
    cc.Timer.autodisplay()
    cc.flush()

def timer_display_stack():
    cc.Timer.display_stack()
    cc.flush()

def timer_display_stack_always():
    cc.Timer.display_stack_always()
    cc.flush()

def timer_reset(long max_call_times_for_always_show_info = -1):
    cc.Timer.reset(max_call_times_for_always_show_info)

def timer_fork(long max_call_times_for_always_show_info = -1):
    cc.Timer.fork(max_call_times_for_always_show_info)

def timer_merge():
    cc.Timer.merge()

def flush():
    cc.flush()

### -------------------------------------------------------------------

def timer(func):
    cdef cc.std_string fname = "py:" + func.__name__
    cdef cc.Timer qtimer = cc.Timer(fname)
    @functools.wraps(func)
    def qtimer_func(*args, **kwargs):
        qtimer.start()
        ret = func(*args, **kwargs)
        qtimer.stop()
        return ret
    return qtimer_func

def timer_verbose(func):
    cdef cc.std_string fname = "py:" + func.__name__
    cdef cc.Timer qtimer = cc.Timer(fname)
    @functools.wraps(func)
    def qtimer_func(*args, **kwargs):
        cdef cc.bool is_verbose = True
        qtimer.start(is_verbose)
        ret = func(*args, **kwargs)
        qtimer.stop(is_verbose)
        return ret
    return qtimer_func

def timer_flops(func):
    cdef cc.std_string fname = "py:" + func.__name__
    cdef cc.Timer qtimer = cc.Timer(fname)
    @functools.wraps(func)
    def qtimer_func(*args, **kwargs):
        qtimer.start()
        flops, ret = func(*args, **kwargs)
        qtimer.flops += flops
        qtimer.stop()
        return ret
    return qtimer_func

def timer_verbose_flops(func):
    cdef cc.std_string fname = "py:" + func.__name__
    cdef cc.Timer qtimer = cc.Timer(fname)
    @functools.wraps(func)
    def qtimer_func(*args, **kwargs):
        cdef cc.bool is_verbose = True
        qtimer.start(is_verbose)
        flops, ret = func(*args, **kwargs)
        qtimer.flops += flops
        qtimer.stop(is_verbose)
        return ret
    return qtimer_func

### -------------------------------------------------------------------

cdef class Buffer:

    def __cinit__(self, object obj = None, int ndim = 1, int itemsize = 1):
        self.obj = obj
        self.ndim = ndim
        self.itemsize = itemsize
        self.shape_strides.resize(ndim * 2)

    cdef Py_ssize_t get_len(self):
        cdef int i
        cdef Py_ssize_t ret = 1
        for i in range(self.ndim):
            ret *= self.shape_strides[i]
        return ret

    cdef void set_strides(self):
        cdef Py_ssize_t* shapes = &self.shape_strides[0]
        cdef Py_ssize_t* strides = &self.shape_strides[self.ndim]
        cdef int i
        cdef Py_ssize_t stride = self.itemsize
        for i in range(self.ndim - 1, -1, -1):
            strides[i] = stride
            stride *= shapes[i]

### -------------------------------------------------------------------

cdef class Coordinate:

    def __cinit__(self):
        pass

    def __init__(self, x = None):
        if isinstance(x, Coordinate):
            self.xx = (<Coordinate>x).xx
        elif isinstance(x, (list, tuple,)):
            assert len(x) == 4
            self.xx = cc.Coordinate(x[0], x[1], x[2], x[3])
        else:
            assert x is None

    def __imatmul__(self, Coordinate v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.bool is_copying_data = True):
        cdef Coordinate x = Coordinate()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def __repr__(self):
        return f"Coordinate({self.list()})"

    def list(self):
        return [ self.xx[0], self.xx[1], self.xx[2], self.xx[3], ]

### -------------------------------------------------------------------

cdef class RngState:

    def __cinit__(self):
        self.cdata = <long>&(self.xx)

    def __init__(self, x = None, y = None):
        cdef cc.std_string seed
        if x is None:
            assert y is None
            # make a new rng
            self.xx = cc.RngState()
        elif isinstance(x, RngState):
            if y is None:
                # make a copy of x
                self.xx = (<RngState>x).xx
            else:
                # split into a new rng
                seed = str(y)
                self.xx = cc.RngState((<RngState>x).xx, seed)
        else:
            assert y is None
            # seed a new rng
            seed = str(x)
            self.xx = cc.RngState(seed)

    def __imatmul__(self, RngState v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.bool is_copying_data = True):
        cdef RngState x = RngState()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def split(self, const cc.std_string& seed):
        cdef RngState x = RngState()
        x.xx = self.xx.split(seed)
        return x

    def rand_gen(self):
        return cc.rand_gen(self.xx)

    def u_rand_gen(self, double upper = 1.0, double lower = 0.0):
        return cc.u_rand_gen(self.xx, upper, lower)

    def g_rand_gen(self, double center = 0.0, double sigma = 1.0):
        return cc.g_rand_gen(self.xx, center, sigma)

    def c_rand_gen(self, Coordinate size):
        # size can be total_site of the lattice
        cdef Coordinate x = Coordinate()
        x.xx = cc.c_rand_gen(self.xx, size.xx)
        return x

    def select(self, list l):
        ri = self.rand_gen() % len(l)
        return l[ri]

    def u_rand_fill(self, arr, double upper = 1.0, double lower = 0.0):
        arr = arr.ravel()
        assert arr.base is not None
        arr = arr.view(np.float64)
        assert arr.base is not None
        return self.u_rand_fill_double(arr, upper, lower)

    def u_rand_fill_double(self, double[:] arr, double upper = 1.0, double lower = 0.0):
        cdef long size = arr.size
        cdef long i
        for i in range(size):
            arr[i] = cc.u_rand_gen(self.xx, upper, lower)

    def g_rand_fill(self, arr, double center = 0.0, double sigma = 1.0):
        arr = arr.ravel()
        assert arr.base is not None
        arr = arr.view(np.float64)
        assert arr.base is not None
        return self.g_rand_fill_double(arr, center, sigma)

    def g_rand_fill_double(self, double[:] arr, double center = 0.0, double sigma = 1.0):
        cdef long size = arr.size
        cdef long i
        for i in range(size):
            arr[i] = cc.g_rand_gen(self.xx, center, sigma)

### -------------------------------------------------------------------

@timer
def get_double_sig(x, RngState rs):
    # get a signature of data viewed as a collection of double numbers
    if isinstance(x, LatData):
        arr = np.asarray(x).ravel()
        arr_rand = np.zeros(arr.shape, arr.dtype)
        rs.u_rand_fill(arr_rand, 1.0, -1.0)
        return np.sum(arr * arr_rand)
    else:
        return None

### -------------------------------------------------------------------

cdef class LatData:

    def __cinit__(self):
        self.cdata = <long>&(self.xx)

    def __imatmul__(self, LatData v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.bool is_copying_data = True):
        x = LatData()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef cc.bool is_complex = self.xx.is_complex()
        cdef char* fmt
        cdef int itemsize
        if is_complex:
            fmt = 'Zd'
            itemsize = sizeof(cc.Complex)
        else:
            fmt = 'd'
            itemsize = sizeof(cc.Double)
        cdef Buffer buf = Buffer(self, self.xx.ndim(), itemsize)
        cdef Py_ssize_t* shape = &buf.shape_strides[0]
        cdef Py_ssize_t* strides = &buf.shape_strides[buf.ndim]
        cdef int i
        for i in range(buf.ndim):
            shape[i] = self.xx.info[i].size
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

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

    def save_node(self, const cc.std_string& path):
        self.xx.save(path)

    def load_node(self, const cc.std_string& path):
        self.xx.load(path)

    def save(self, const cc.std_string& path):
        if cc.get_id_node() == 0:
            self.save_node(path)

    def load(self, const cc.std_string& path):
        if cc.get_id_node() == 0:
            self.load_node(path)
        self.bcast()

    def bcast(self):
        if cc.get_num_node() != 1:
            import cqlat as c
            c.bcast_lat_data(self)
        return self

    def glb_sum_in_place(self):
        if cc.get_num_node() != 1:
            import cqlat as c
            c.glb_sum_lat_data(self)
        return self

    def glb_sum(self):
        ld = self.copy()
        return ld.glb_sum_in_place()

    def __str__(self):
        return cc.show(self.xx)

    def show(self):
        return str(self)

    def set_zero(self):
        cc.set_zero(self.xx)

    def qnorm(self):
        return cc.qnorm(self.xx)

    def is_match(self, LatData ld1):
        # ld.info needs to be exactly equal
        return cc.is_matching(self.xx, ld1.xx)

    def is_complex(self):
        return self.xx.is_complex()

    def ndim(self):
        return self.xx.ndim()

    def dim_name(self, int dim):
        assert 0 <= dim
        assert dim < self.xx.info.size()
        return self.xx.info[dim].name

    def dim_size(self, int dim):
        assert 0 <= dim
        assert dim < self.xx.info.size()
        return self.xx.info[dim].size

    def dim_indices(self, int dim):
        assert 0 <= dim
        assert dim < self.xx.info.size()
        return self.xx.info[dim].indices

    def dim_sizes(self):
        cdef int ndim = self.xx.ndim()
        cdef int i
        return [ self.xx.info[i].size for i in range(ndim) ]

    def set_dim_sizes(self, list dim_sizes, *, cc.bool is_complex = True):
        cdef int ndim = len(dim_sizes)
        cdef int ndim_real = ndim
        cdef int i
        if is_complex:
            ndim_real += 1
        self.xx.info.resize(ndim_real)
        for i in range(ndim):
            self.xx.info[i].size = dim_sizes[i]
        if is_complex:
            self.xx.info[ndim] = cc.lat_dim_re_im()
        cc.lat_data_alloc(self.xx)

    def set_dim_name(self, int dim, const cc.std_string& name, list indices = None):
        assert 0 <= dim
        assert dim < self.xx.info.size()
        cdef int size
        cdef int i
        self.xx.info[dim].name = name
        if indices is None:
            self.xx.info[dim].indices.resize(0)
        else:
            indices = [ str(idx).replace("\n", "  ") for idx in indices ]
            size = len(indices)
            self.xx.info[dim].indices.resize(size)
            for i in range(size):
                self.xx.info[dim].indices[i] = indices[i]

    def dim_names(self):
        # by default, return list can be used as the input argument for ld.from_numpy
        cdef int ndim = self.xx.ndim()
        cdef int i
        return [ self.xx.info[i].name for i in range(ndim) ]

    def to_numpy(self):
        return np.asarray(self).copy()

    def from_numpy(self, numpy.ndarray val, list dim_names = None, *, cc.bool is_complex = True):
        # only set LatData shape if it is initially empty
        # otherwise only set data and ignore shape completely
        # dim_names should be a list of names for each dimension
        cdef int ndim = val.ndim
        cdef int dim
        cdef list shape = [ val.shape[dim] for dim in range(ndim) ]
        if self.ndim() == 0:
            self.set_dim_sizes(shape, is_complex = is_complex)
            if dim_names is None:
                dim_names = [ n for n in "ijklmnopqrstuvwxyz" ]
            assert ndim <= len(dim_names)
            for dim in range(ndim):
                name = dim_names[dim]
                self.set_dim_name(dim, name)
        np.asarray(self).ravel()[:] = val.ravel()[:]
        return self

    def to_list(self):
        return np.asarray(self).ravel().tolist()

    def from_list(self, list val, *, cc.bool is_complex = True):
        if self.ndim() == 0:
            self.set_dim_sizes([ len(val), ], is_complex = is_complex)
            self.set_dim_name(0, "i")
        np.asarray(self).ravel()[:] = np.array(val)
        return self

    def __getstate__(self):
        is_complex = self.is_complex()
        ndim = self.ndim()
        dim_sizes = self.dim_sizes()
        assert len(dim_sizes) == ndim
        dim_names = [ self.dim_name(dim) for dim in range(ndim) ]
        dim_indices = [ self.dim_indices(dim) for dim in range(ndim) ]
        data_list = self.to_list()
        return [ is_complex, dim_sizes, dim_names, dim_indices, data_list ]

    def __setstate__(self, state):
        [ is_complex, dim_sizes, dim_names, dim_indices, data_list ] = state
        self.__init__()
        self.set_dim_sizes(dim_sizes, is_complex = is_complex)
        ndim = len(dim_sizes)
        for dim in range(ndim):
            self.set_dim_name(dim, dim_names[dim], dim_indices[dim])
        self.from_list(data_list)

    def info(self, dim = None):
        # by default, return list can be used as the input argument for ld.set_info or mk_lat_data
        if dim is None:
            ndim = self.ndim()
            return [ self.info(i) for i in range(ndim) ]
        else:
            dim_name = self.dim_name(dim)
            dim_size = self.dim_size(dim)
            dim_indices = self.dim_indices(dim)
            return [ dim_name, dim_size, dim_indices, ]

    def set_info(self, list info_list, *, cc.bool is_complex = True):
        # info_list format:
        # [ [ dim_name, dim_size, dim_indices, ], ... ]
        # dim_indices can be optional
        for info in info_list:
            assert len(info) >= 2
        dim_sizes = [ info[1] for info in info_list ]
        self.set_dim_sizes(dim_sizes, is_complex = is_complex)
        ndim = len(dim_sizes)
        for dim in range(ndim):
            info = info_list[dim]
            if len(info) == 2:
                dim_name, dummy_dim_size = info
                self.set_dim_name(dim, dim_name)
            elif len(info) == 3:
                dim_name, dummy_dim_size, dim_indices = info
                self.set_dim_name(dim, dim_name, dim_indices)
            else:
                raise Exception(f"LatData setinfo info_list={info_list}")
        self.set_zero()

    def __iadd__(self, ld1):
        assert isinstance(ld1, LatData)
        c.set_add_lat_data(self, ld1)
        return self

    def __isub__(self, ld1):
        assert isinstance(ld1, LatData)
        c.set_sub_lat_data(self, ld1)
        return self

    def __imul__(self, factor):
        assert isinstance(factor, float)
        c.set_mul_double_lat_data(self, factor)
        return self

    def __add__(self, ld1):
        ld = self.copy()
        ld += ld1
        return ld

    def __sub__(self, ld1):
        ld = self.copy()
        ld -= ld1
        return ld

    def __mul__(self, factor):
        ld = self.copy()
        ld *= factor
        return ld

    __rmul__ = __mul__

    def __neg__(self):
        ld = self.copy()
        ld *= -1
        return ld

    def __pos__(self):
        return self

    def __setitem__(self, idx, val):
        # use list with correct length as val
        # idx should be tuple or list of int
        if isinstance(idx, int):
            idx = [ idx, ]
        if isinstance(val, np.ndarray):
            val = val.flatten()
        c.poke_lat_data(self, idx, list(val))

    def __getitem__(self, idx):
        # return a new list every call
        # idx should be tuple or list of int
        if isinstance(idx, int):
            idx = [ idx, ]
        shape = self.dim_sizes()[len(idx):]
        return np.array(c.peek_lat_data(self, idx)).reshape(shape)

### -------------------------------------------------------------------

def mk_lat_data(list info_list, *, cc.bool is_complex = True):
    # info_list format:
    # [ [ dim_name, dim_size, dim_indices, ], ... ]
    # dim_indices can be optional
    ld = LatData()
    ld.set_info(info_list, is_complex = is_complex)
    return ld

def load_lat_data(const cc.std_string& path):
    ld = LatData()
    ld.load(path)
    return ld

### -------------------------------------------------------------------

cdef class WilsonMatrix:

    def __cinit__(self):
        self.cdata = <long>&(self.xx)

    def __imatmul__(self, WilsonMatrix v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.bool is_copying_data = True):
        cdef WilsonMatrix x = WilsonMatrix()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def set_zero(self):
        cc.set_zero(self.xx)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Buffer buf = Buffer(self, ElemTypeWilsonMatrix.ndim(), ElemTypeWilsonMatrix.itemsize())
        cdef char* fmt = ElemTypeWilsonMatrix.format()
        cdef cc.std_vector[Py_ssize_t] vec = ElemTypeWilsonMatrix.shape()
        assert vec.size() == buf.ndim
        cdef Py_ssize_t* shape = &buf.shape_strides[0]
        cdef Py_ssize_t* strides = &buf.shape_strides[buf.ndim]
        cdef int i
        for i in range(buf.ndim):
            shape[i] = vec[i]
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

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

    def g5_herm(self):
        self.xx = cc.g5_herm(self.xx)

### -------------------------------------------------------------------

cdef class SpinMatrix:

    def __cinit__(self):
        self.cdata = <long>&(self.xx)

    def __imatmul__(self, SpinMatrix v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.bool is_copying_data = True):
        cdef SpinMatrix x = SpinMatrix()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def set_zero(self):
        cc.set_zero(self.xx)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Buffer buf = Buffer(self, ElemTypeSpinMatrix.ndim(), ElemTypeSpinMatrix.itemsize())
        cdef char* fmt = ElemTypeWilsonMatrix.format()
        cdef cc.std_vector[Py_ssize_t] vec = ElemTypeSpinMatrix.shape()
        assert vec.size() == buf.ndim
        cdef Py_ssize_t* shape = &buf.shape_strides[0]
        cdef Py_ssize_t* strides = &buf.shape_strides[buf.ndim]
        cdef int i
        for i in range(buf.ndim):
            shape[i] = vec[i]
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

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

### -------------------------------------------------------------------

gamma_matrix_0 = cc.get_gamma_matrix(0)
gamma_matrix_1 = cc.get_gamma_matrix(1)
gamma_matrix_2 = cc.get_gamma_matrix(2)
gamma_matrix_3 = cc.get_gamma_matrix(3)
gamma_matrix_5 = cc.get_gamma_matrix(5)

### -------------------------------------------------------------------

def as_wilson_matrix(x):
    cdef WilsonMatrix wm
    if isinstance(x, WilsonMatrix):
        return x
    elif x == 0:
        wm = WilsonMatrix()
        cc.set_zero(wm.xx)
        return wm

def as_wilson_matrix_g5_herm(x):
    cdef WilsonMatrix wm = WilsonMatrix()
    if isinstance(x, WilsonMatrix):
        wm.xx = cc.g5_herm((<WilsonMatrix>x).xx)
    elif x == 0:
        cc.set_zero(wm.xx)
    return wm

### -------------------------------------------------------------------
