# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc

import functools
import sys

### -------------------------------------------------------------------

def get_verbose_level():
    """
    Return the current verbosity level as integer.\n
    Usage::\n
        get_verbose_level() # return the current verbosity level
    Default get_verbose_level() depends on the ``q_verbose`` environment variable. If unset, the default is ``-1``.
    """
    return cc.get_verbose_level()

def set_verbose_level(level=None):
    """
    Set the current verbosity level as integer.\n
    Usage::\n
        set_verbose_level() # set the verbosity level to be the default value.
        set_verbose_level(level) # set the verbosity level
    """
    if level is None:
        level = cc.get_verbose_level_default()
    assert isinstance(level, int)
    cdef cc.Long* p_ret = &cc.get_verbose_level()
    p_ret[0] = level
    assert cc.get_verbose_level() == level

set_verbose_level(-1)

### -------------------------------------------------------------------

cdef void display_py_stdout(const cc.std_string& msg):
    sys.stdout.write(msg)
    sys.stdout.flush()

def set_display_method(method=None):
    if method is None:
        cc.set_display_ptr()
    elif method == "py_stdout":
        cc.set_display_ptr(display_py_stdout)
    else:
        raise Exception(f"set_display_method(method='{method}')")

### -------------------------------------------------------------------

def display(level, *args):
    """
    Print all the arguments.
    Interpret the first argument as verbose level if it is ``int``.
    Only print if ``level <= get_verbose_level()``.
    If the first argument is not integer, will always print all the arguments.
    """
    if isinstance(level, int):
        if level <= cc.get_verbose_level():
            print(*args, end='', flush=True)
    else:
        print(level, *args, end='', flush=True)

def display_info(*args):
    """
    Same as ``display`` but only print if ``get_id_node() == 0``.
    """
    if cc.get_id_node() == 0:
        display(*args)

def displayln(level, *args):
    """
    Print all the arguments and then print a newline.
    Interpret the first argument as verbose level if it is ``int``.
    Only print if ``level <= get_verbose_level()``.
    If the first argument is not integer, will always print all the arguments.
    """
    if isinstance(level, int):
        if level <= cc.get_verbose_level():
            print(*args, flush=True)
    else:
        print(level, *args, flush=True)

def displayln_info(*args):
    """
    Same as ``displayln`` but only print if ``get_id_node() == 0``.
    """
    if cc.get_id_node() == 0:
        displayln(*args)

### -------------------------------------------------------------------

cdef class Timer:

    def __cinit__(self, const cc.std_string& fname, cc.Bool is_verbose = False):
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
    """
    Return the node (MPI process) id as ``int``.
    Also works without initializing MPI, in which case will always return 0.
    """
    return cc.get_id_node()

def get_num_node():
    """
    Return number of nodes (MPI processes) as ``int``.
    Also works without initializing MPI, in which case will always return 1.
    """
    return cc.get_num_node()

def sync_node():
    """
    Perform a glb_sum to make sure all nodes are synced.
    """
    cc.sync_node()

def get_time():
    """
    Return current time in seconds since epoch.
    """
    return cc.get_time()

def get_start_time():
    """
    Return start time in seconds since epoch. Does reset by ``timer_reset``
    """
    return cc.get_start_time()

def set_start_time(time):
    cdef double* p_ret = &cc.get_start_time()
    p_ret[0] = time
    assert cc.get_start_time() == time

def get_actual_start_time():
    """
    Return start time in seconds since epoch. Does not reset by ``timer_reset``
    """
    return cc.get_actual_start_time()

def set_actual_start_time(time):
    cdef double* p_ret = &cc.get_actual_start_time()
    p_ret[0] = time
    assert cc.get_actual_start_time() == time

def get_total_time():
    """
    Return total time in seconds. Does reset by ``timer_reset``
    """
    return cc.get_total_time()

def get_actual_total_time():
    """
    Return total time in seconds. Does not reset by ``timer_reset``
    """
    return cc.get_actual_total_time()

def get_time_limit():
    """
    Return time limit of the program in seconds.
    Used in check_time_limit, get_remaining_time
    """
    return cc.get_time_limit()

def set_time_limit(time_limit=None):
    """
    Set time limit of the program in seconds.
    Usage::\n
        set_time_limit() # set to be the default value.
        set_time_limit(time_limit) # set the time limit
    """
    if time_limit is None:
        time_limit = cc.get_time_limit_default()
    cdef double* p_ret = &cc.get_time_limit()
    p_ret[0] = time_limit
    assert cc.get_time_limit() == time_limit

def get_remaining_time():
    """
    Return remaining time in seconds.
    """
    return cc.get_remaining_time()

def get_time_budget():
    return cc.get_time_budget()

def set_time_budget(time_budget=None):
    if time_budget is None:
        time_budget = cc.get_time_budget_default()
    cdef double* p_ret = &cc.get_time_budget()
    p_ret[0] = time_budget
    assert cc.get_time_budget() == time_budget

### -------------------------------------------------------------------

def timer_display(str tag = ""):
    cc.Timer.display(tag)

def timer_autodisplay():
    cc.Timer.autodisplay()

def timer_display_stack():
    cc.Timer.display_stack()

def timer_display_stack_always():
    cc.Timer.display_stack_always()

def timer_reset(cc.Long max_call_times_for_always_show_info = -1):
    """
    Reset all timers, ``get_total_time``, ``get_start_time``.
    But does not reset ``get_actual_start_time`` or ``get_actual_total_time``
    """
    cc.Timer.reset(max_call_times_for_always_show_info)

def timer_fork(cc.Long max_call_times_for_always_show_info = -1):
    cc.Timer.fork(max_call_times_for_always_show_info)

def timer_merge():
    cc.Timer.merge()

### -------------------------------------------------------------------

def timer_builder(object func, cc.Bool is_verbose, cc.Bool is_flops) -> object:
    cdef cc.std_string fname = "py:" + func.__name__
    cdef cc.Timer qtimer = cc.Timer(fname)
    if not is_flops:
        @functools.wraps(func)
        def qtimer_func(*args, **kwargs):
            qtimer.start(is_verbose)
            ret = func(*args, **kwargs)
            qtimer.stop(is_verbose)
            return ret
    else:
        @functools.wraps(func)
        def qtimer_func(*args, **kwargs):
            qtimer.start(is_verbose)
            flops, ret = func(*args, **kwargs)
            qtimer.flops += flops
            qtimer.stop(is_verbose)
            return ret
    return qtimer_func

def timer(object func) -> object:
    """
    Timing functions.\n
    Usage::\n
        @q.timer
        def function(args):
            ...
    """
    return timer_builder(func, False, False)

def timer_verbose(func):
    """
    Timing functions. Always show output if ``get_verbose_level() > 0``\n
    Usage::\n
        @q.timer_verbose
        def function(args):
            ...
    """
    return timer_builder(func, True, False)

def timer_flops(func):
    """
    Timing functions with flops.\n
    Usage::\n
        @q.timer_flops
        def function(args):
            ...
            return flops, ret
    Modified function will only return ``ret`` in above example.
    """
    return timer_builder(func, False, True)

def timer_verbose_flops(func):
    """
    Timing functions with flops. Always show output if ``get_verbose_level() > 0``\n
    Usage::\n
        @q.timer_flops
        def function(args):
            ...
            return flops, ret
    Modified function will only return ``ret`` in above example.
    """
    return timer_builder(func, True, True)
