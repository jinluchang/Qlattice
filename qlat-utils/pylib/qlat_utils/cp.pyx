# cython: c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cp

import functools

### -------------------------------------------------------------------

def flush():
    cp.flush()

def timer_display(str tag = ""):
    cp.Timer.display(tag)
    cp.flush()

def timer(func):
    fname = "py:" + func.__name__
    cdef cp.Timer qtimer = cp.Timer(fname)
    @functools.wraps(func)
    def qtimer_func(*args, **kwargs):
        qtimer.start()
        ret = func(*args, **kwargs)
        qtimer.stop()
        return ret
    return qtimer_func

def timer_verbose(func):
    fname = "py:" + func.__name__
    cdef cp.Timer qtimer = cp.Timer(fname)
    @functools.wraps(func)
    def qtimer_func(*args, **kwargs):
        cdef cp.bool is_verbose = True
        qtimer.start(is_verbose)
        ret = func(*args, **kwargs)
        qtimer.stop(is_verbose)
        return ret
    return qtimer_func

### -------------------------------------------------------------------

cdef class Coordinate:

    cdef cp.Coordinate c

    def __cinit__(self):
        self.c = cp.Coordinate()

    def __init__(self, list x):
        if x is not None:
            assert len(x) == 4
            self.c = cp.Coordinate(x[0], x[1], x[2], x[3])

    def __imatmul__(self, Coordinate v1):
        self.c = v1.c
        return self

    cpdef Coordinate copy(self, cp.bool is_copying_data = True):
        cdef Coordinate x = Coordinate()
        if is_copying_data:
            x.c = self.c
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def list(self):
        return [ self.c[0], self.c[1], self.c[2], self.c[3], ]

### -------------------------------------------------------------------

cdef class RngState:

    cdef cp.RngState rs
    cdef readonly long cdata

    def __cinit__(self):
        self.rs = cp.RngState()
        self.cdata = <long>&(self.rs)

    def __init__(self, x = None, y = None):
        cdef cp.string seed
        if x is None:
            assert y is None
            # make a new rng
            self.rs = cp.RngState()
        elif isinstance(x, RngState):
            if y is None:
                # make a copy of x
                self.rs = (<RngState>x).rs
            else:
                # split into a new rng
                seed = str(y)
                self.rs = cp.RngState((<RngState>x).rs, seed)
        else:
            assert y is None
            # seed a new rng
            seed = str(x)
            self.rs = cp.RngState(seed)

    def __imatmul__(self, RngState v1):
        self.rs = v1.rs
        return self

    cpdef RngState copy(self, cp.bool is_copying_data = True):
        cdef RngState x = RngState()
        if is_copying_data:
            x.rs = self.rs
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    cpdef RngState split(self, const cp.string& seed):
        cdef RngState x = RngState()
        x.rs = self.rs.split(seed)
        return x

    cpdef cp.uint64_t rand_gen(self):
        return cp.rand_gen(self.rs)

    cpdef u_rand_gen(self, double upper = 1.0, double lower = 0.0):
        return cp.u_rand_gen(self.rs, upper, lower)

    cpdef g_rand_gen(self, double center = 0.0, double sigma = 1.0):
        return cp.g_rand_gen(self.rs, center, sigma)

    cpdef c_rand_gen(self, Coordinate size):
        # size can be total_site of the lattice
        cdef Coordinate x = Coordinate()
        x.c = cp.c_rand_gen(self.rs, size.c)
        return x

    def select(self, list l):
        ri = self.rand_gen() % len(l)
        return l[ri]

### -------------------------------------------------------------------
