# cython: c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc
import functools

### -------------------------------------------------------------------

def flush():
    cc.flush()

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

### -------------------------------------------------------------------

def timer(func):
    cdef cc.string fname = "py:" + func.__name__
    cdef cc.Timer qtimer = cc.Timer(fname)
    @functools.wraps(func)
    def qtimer_func(*args, **kwargs):
        qtimer.start()
        ret = func(*args, **kwargs)
        qtimer.stop()
        return ret
    return qtimer_func

def timer_verbose(func):
    cdef cc.string fname = "py:" + func.__name__
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
    cdef cc.string fname = "py:" + func.__name__
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
    cdef cc.string fname = "py:" + func.__name__
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

cdef class Timer:

    def __cinit__(self, const cc.string& fname, cc.bool is_verbose = False):
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

cdef class Buffer2D:

    def __cinit__(self, object obj = None):
        self.obj = obj

### -------------------------------------------------------------------

cdef class Coordinate:

    def __cinit__(self):
        self.xx = cc.Coordinate()

    def __init__(self, list x = None):
        if x is not None:
            assert len(x) == 4
            self.xx = cc.Coordinate(x[0], x[1], x[2], x[3])

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

    def list(self):
        return [ self.xx[0], self.xx[1], self.xx[2], self.xx[3], ]

### -------------------------------------------------------------------

cdef class RngState:

    def __cinit__(self):
        self.xx = cc.RngState()
        self.cdata = <long>&(self.xx)

    def __init__(self, x = None, y = None):
        cdef cc.string seed
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

    def split(self, const cc.string& seed):
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

### -------------------------------------------------------------------

cdef class WilsonMatrix:

    def __cinit__(self):
        self.xx = cc.WilsonMatrix()
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
        cdef Py_ssize_t itemsize = sizeof(cc.Complex)
        cdef Buffer2D buf = Buffer2D(self)
        buf.shape[0] = 12
        buf.shape[1] = 12
        buf.strides[1] = itemsize
        buf.strides[0] = buf.strides[1] * buf.shape[1]
        buffer.buf = <char*>&(self.xx.p)
        buffer.format = 'Zd'
        buffer.internal = NULL
        buffer.itemsize = itemsize
        buffer.len = buf.shape[0] * buf.shape[1]
        buffer.ndim = 2
        buffer.obj = buf
        buffer.readonly = 0
        buffer.shape = buf.shape
        buffer.strides = buf.strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

    def g5_herm(self):
        self.xx = cc.g5_herm(self.xx)

### -------------------------------------------------------------------

cdef class SpinMatrix:

    def __cinit__(self):
        self.xx = cc.SpinMatrix()
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
        cdef Py_ssize_t itemsize = sizeof(cc.Complex)
        cdef Buffer2D buf = Buffer2D(self)
        buf.shape[0] = 4
        buf.shape[1] = 4
        buf.strides[1] = itemsize
        buf.strides[0] = buf.strides[1] * buf.shape[1]
        buffer.buf = <char*>&(self.xx.p)
        buffer.format = 'Zd'
        buffer.internal = NULL
        buffer.itemsize = itemsize
        buffer.len = buf.shape[0] * buf.shape[1]
        buffer.ndim = 2
        buffer.obj = buf
        buffer.readonly = 0
        buffer.shape = buf.shape
        buffer.strides = buf.strides
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

def as_wilson_matrix_from_numpy(cc.Complex[:, :] x):
    cdef WilsonMatrix wm = WilsonMatrix()
    cdef cc.Complex[:, :] view = wm
    view[:, :] = x[:, :]
    return wm

def as_spin_matrix_from_numpy(cc.Complex[:, :] x):
    cdef SpinMatrix sm = SpinMatrix()
    cdef cc.Complex[:, :] view = sm
    view[:, :] = x[:, :]
    return sm

### -------------------------------------------------------------------
