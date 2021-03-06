import cqlat as c

from qlat.mpi import *

from cqlat import timer_display
from cqlat import timer_display_stack, timer_display_stack_always

import functools

timer_total_flops = {}

def acc_timer_flops(fname, flops):
    if fname in timer_total_flops:
        timer_total_flops[fname] += flops
    else:
        timer_total_flops[fname] = flops

class Timer:

    def __init__(self, fname, is_verbose = False):
        self.cdata = c.mk_timer(fname)
        self.is_verbose = is_verbose
        acc_timer_flops(fname, 0)

    def __del__(self):
        c.free_timer(self)

    def start(self, is_verbose = None):
        if is_verbose is None:
            is_verbose = self.is_verbose
        c.start_timer(self, is_verbose)

    def stop(self, is_verbose = None):
        if is_verbose is None:
            is_verbose = self.is_verbose
        c.stop_timer(self, is_verbose)

    def set_flops(self, flops):
        c.set_flops_timer(self, flops)

def timer(func):
    fname = "py:" + func.__name__
    timer = Timer(fname)
    @functools.wraps(func)
    def qtimer_func(*args, **kwargs):
        timer.start()
        start_flops = timer_total_flops[fname]
        ret = func(*args, **kwargs)
        stop_flops = timer_total_flops[fname]
        timer.set_flops(stop_flops - start_flops)
        timer.stop()
        return ret
    return qtimer_func

def timer_verbose(func):
    fname = "py:" + func.__name__
    timer = Timer(fname, True)
    @functools.wraps(func)
    def qtimer_func(*args, **kwargs):
        timer.start()
        start_flops = timer_total_flops[fname]
        ret = func(*args, **kwargs)
        stop_flops = timer_total_flops[fname]
        timer.set_flops(stop_flops - start_flops)
        timer.stop()
        return ret
    return qtimer_func

def displayln(*args):
    print(*args)

def displayln_info(*args):
    if get_id_node() == 0:
        displayln(*args)

class TimerNone(Timer):

    def __init__(self):
        pass

    def __del__(self):
        pass

    def start(self):
        pass

    def stop(self):
        pass
