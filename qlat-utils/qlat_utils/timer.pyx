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
            print(f"Qlat {get_actual_total_time():15.3f}:", *args, end='', flush=True)
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
            print(f"Qlat {get_actual_total_time():15.3f}:", *args, flush=True)
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

def timer_reset(cc.Long max_call_times_for_always_show_info=-1):
    """
    Reset all timers, ``get_total_time``, ``get_start_time``.
    But does not reset ``get_actual_start_time`` or ``get_actual_total_time``
    """
    cc.Timer.reset(max_call_times_for_always_show_info)

def timer_fork(cc.Long max_call_times_for_always_show_info=-1):
    cc.Timer.fork(max_call_times_for_always_show_info)

def timer_merge():
    cc.Timer.merge()

### -------------------------------------------------------------------

class TimerFork:

    """
    with TimerFork():
        ...
    #
    Will always restore `verbose` level to initial value.
    """

    def __init__(self,
                 *,
                 show_display=10,
                 display_tag="TimerFork",
                 verbose=None,
                 max_call_times_for_always_show_info=-1,
                 ):
        self.max_call_times_for_always_show_info = max_call_times_for_always_show_info
        self.show_display = show_display
        self.orig_verbose_level = get_verbose_level()
        self.verbose_level = verbose
        self.display_tag = display_tag

    def __enter__(self):
        timer_fork(self.max_call_times_for_always_show_info)
        if self.verbose_level is not None:
            set_verbose_level(self.verbose_level)

    def __exit__(self, exc_type, exc_value, traceback):
        assert exc_type is None
        assert exc_value is None
        assert traceback is None
        set_verbose_level(self.orig_verbose_level)
        assert isinstance(self.display_tag, str)
        if isinstance(self.show_display, bool):
            if self.show_display:
                timer_display(self.display_tag)
        elif isinstance(self.show_display, (int, float,)):
            if get_total_time() >= self.show_display:
                timer_display(self.display_tag + "(auto)")
        else:
            raise Exception(f"TimerFork: show_display='{self.show_display}'.")
        timer_merge()

### -------------------------------------------------------------------

default_timer_kwargs = dict(
        fname=None,
        is_verbose=False,
        is_flops=False,
        is_timer_fork=False,
        timer_fork_verbose=None,
        timer_fork_max_call_times_for_always_show_info=-1,
        timer_fork_show_display=10,
        )

def timer(object func=None, **kwargs) -> object:
    """
    Timing functions.\n
    Usage::\n
        @q.timer
        def function(args):
            ...\n
        @q.timer(is_verbose=True)
        def function(args):
            ...\n
        @q.timer(is_timer_fork=True)
        def function(args):
            ...\n
    if is_flops:
        The function `func` need to return `(flops, ret,)`
    """
    if func is not None:
        assert kwargs == dict()
        return timer()(func)
    kwargs = dict(default_timer_kwargs, **kwargs)
    fname = kwargs["fname"]
    is_verbose = kwargs["is_verbose"]
    is_flops = kwargs["is_flops"]
    is_timer_fork = kwargs["is_timer_fork"]
    timer_fork_verbose = kwargs["timer_fork_verbose"]
    timer_fork_max_call_times_for_always_show_info = kwargs["timer_fork_max_call_times_for_always_show_info"]
    timer_fork_show_display = kwargs["timer_fork_show_display"]
    #
    if is_timer_fork:
        is_verbose = True
    #
    def f(object func):
        cdef cc.std_string fname_str
        if fname is None:
            fname_str = "py:" + func.__qualname__
        else:
            fname_str = fname
        cdef cc.Timer qtimer = cc.Timer(fname_str)
        if not is_timer_fork:
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
        else:
            if not is_flops:
                @functools.wraps(func)
                def qtimer_func(*args, **kwargs):
                    with TimerFork(
                            verbose=timer_fork_verbose,
                            max_call_times_for_always_show_info=timer_fork_max_call_times_for_always_show_info,
                            show_display=timer_fork_show_display,
                            ):
                        qtimer.start(is_verbose)
                        ret = func(*args, **kwargs)
                        qtimer.stop(is_verbose)
                    return ret
            else:
                @functools.wraps(func)
                def qtimer_func(*args, **kwargs):
                    with TimerFork(
                            verbose=timer_fork_verbose,
                            max_call_times_for_always_show_info=timer_fork_max_call_times_for_always_show_info,
                            show_display=timer_fork_show_display,
                            ):
                        qtimer.start(is_verbose)
                        flops, ret = func(*args, **kwargs)
                        qtimer.flops += flops
                        qtimer.stop(is_verbose)
                    return ret
        return qtimer_func
    return f

def timer_fname(str fname) -> object:
    """
    Timing functions.\n
    Usage::\n
        @q.timer_fname("fname")
        def function(args):
            ...
    """
    return timer(fname=fname)

def timer_verbose(object func):
    """
    Timing functions. Always show output if ``get_verbose_level() > 0``\n
    Usage::\n
        @q.timer_verbose
        def function(args):
            ...
    """
    return timer(is_verbose=True)(func)

def timer_verbose_fname(str fname):
    """
    Timing functions. Always show output if ``get_verbose_level() > 0``\n
    Usage::\n
        @q.timer_verbose_fname("fname")
        def function(args):
            ...
    """
    return timer(fname=fname, is_verbose=True)

def timer_flops(object func):
    """
    Timing functions with flops.\n
    Usage::\n
        @q.timer_flops
        def function(args):
            ...
            return flops, ret
    Modified function will only return ``ret`` in above example.
    """
    return timer(is_flops=True)(func)

def timer_flops_fname(str fname):
    """
    Timing functions with flops.\n
    Usage::\n
        @q.timer_flops_fname("fname")
        def function(args):
            ...
            return flops, ret
    Modified function will only return ``ret`` in above example.
    """
    return timer(fname=fname, is_flops=True)

def timer_verbose_flops(object func):
    """
    Timing functions with flops. Always show output if ``get_verbose_level() > 0``\n
    Usage::\n
        @q.timer_flops
        def function(args):
            ...
            return flops, ret
    Modified function will only return ``ret`` in above example.
    """
    return timer(is_verbose=True, is_flops=True)(func)

def timer_verbose_flops_fname(str fname):
    """
    Timing functions with flops. Always show output if ``get_verbose_level() > 0``\n
    Usage::\n
        @q.timer_verbose_flops_fname("fname")
        def function(args):
            ...
            return flops, ret
    Modified function will only return ``ret`` in above example.
    """
    return timer(fname=fname, is_verbose=True, is_flops=True)
