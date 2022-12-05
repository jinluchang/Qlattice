# cython: c_string_type=unicode, c_string_encoding=utf8

from . cimport __init__ as cp

import functools

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
