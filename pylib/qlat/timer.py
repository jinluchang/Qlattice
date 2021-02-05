import cqlat as c

from qlat.mpi import *

from cqlat import timer_display
from cqlat import timer_display_stack, timer_display_stack_always

import functools

class Timer:

    def __init__(self, fname):
        self.cdata = c.mk_timer(fname)

    def __del__(self):
        c.free_timer(self)

    def start(self, is_verbose = False):
        c.start_timer(self, is_verbose)

    def stop(self, is_verbose = False):
        c.stop_timer(self, is_verbose)

    def set_flops(self, flops):
        c.set_flops_timer(self, flops)

timer_total_flops = {}

def acc_timer_flops(func, flops):
    fname = "py:" + func.__name__
    if fname in timer_total_flops:
        timer_total_flops[fname] += flops
    else:
        timer_total_flops[fname] = flops

def timer(func):
    fname = "py:" + func.__name__
    timer = Timer("py:" + func.__name__)
    acc_timer_flops(func, 0)
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
    timer = Timer(fname)
    acc_timer_flops(func, 0)
    @functools.wraps(func)
    def qtimer_func(*args, **kwargs):
        timer.start(True)
        start_flops = timer_total_flops[fname]
        ret = func(*args, **kwargs)
        stop_flops = timer_total_flops[fname]
        timer.set_flops(stop_flops - start_flops)
        timer.stop(True)
        return ret
    return qtimer_func

def displayln(*args):
    print(*args)

def displayln_info(*args):
    if get_id_node() == 0:
        displayln(*args)

