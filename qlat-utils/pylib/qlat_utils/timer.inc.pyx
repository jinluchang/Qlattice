def displayln(level, *args):
    """
    Print all the arguments and then print a newline.
    Interpret the first argument as verbose level if it is ``int``.
    Only print if ``level <= verbose_level()``.
    If the first argument is not integer, will always print all the arguments.
    """
    if isinstance(level, int):
        if level <= verbose_level():
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

def verbose_level(level = None):
    """
    Return or set the current verbosity level as integer.\n
    Usage::\n
        verbose_level(level) # set the verbosity level
        verbose_level() # return the current verbosity level
    Default verbose_level() depends on the ``q_verbose`` environment variable. If unset, the default is ``0``.
    """
    if level is None:
        return cc.verbose_level()
    cdef long* p_ret = &cc.verbose_level()
    p_ret[0] = level
    assert cc.verbose_level() == level
    return level

def get_time():
    """
    Return current time in seconds since epoch.
    """
    return cc.get_time()

def get_start_time(time = None):
    """
    Return start time in seconds since epoch. Does reset by ``timer_reset``
    """
    if time is None:
        return cc.get_start_time()
    cdef double* p_ret = &cc.get_start_time()
    p_ret[0] = time
    assert cc.get_start_time() == time
    return time

def get_actual_start_time(time = None):
    """
    Return start time in seconds since epoch. Does not reset by ``timer_reset``
    """
    if time is None:
        return cc.get_actual_start_time()
    cdef double* p_ret = &cc.get_actual_start_time()
    p_ret[0] = time
    assert cc.get_actual_start_time() == time
    return time

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
    """
    Reset all timers, ``get_total_time``, ``get_start_time``.
    But does not reset ``get_actual_start_time`` or ``get_actual_total_time``
    """
    cc.Timer.reset(max_call_times_for_always_show_info)

def timer_fork(long max_call_times_for_always_show_info = -1):
    cc.Timer.fork(max_call_times_for_always_show_info)

def timer_merge():
    cc.Timer.merge()

### -------------------------------------------------------------------

def timer(func):
    """
    Timing functions.\n
    Usage::\n
        @q.timer
        def function(args):
            ...
    """
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
    """
    Timing functions. Always show output if ``verbose_level() > 0``\n
    Usage::\n
        @q.timer_verbose
        def function(args):
            ...
    """
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
    """
    Timing functions with flops.\n
    Usage::\n
        @q.timer_flops
        def function(args):
            ...
            return flops, ret
    Modified function will only return ``ret`` in above example.
    """
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
    """
    Timing functions with flops. Always show output if ``verbose_level() > 0``\n
    Usage::\n
        @q.timer_flops
        def function(args):
            ...
            return flops, ret
    Modified function will only return ``ret`` in above example.
    """
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
