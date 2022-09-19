import cqlat_utils as cu

from cqlat_utils import get_id_node, get_num_node
from cqlat_utils import timer_reset, timer_fork, timer_merge
from cqlat_utils import verbose_level
from cqlat_utils import get_actual_start_time, get_start_time, get_time
from cqlat_utils import flush

import functools

timer_total_flops = {}

def acc_timer_flops(fname, flops):
    if fname in timer_total_flops:
        timer_total_flops[fname] += flops
    else:
        timer_total_flops[fname] = flops

class Timer:

    def __init__(self, fname, is_verbose = False):
        if fname is None:
            self.cdata = None
        else:
            self.cdata = cu.mk_timer(fname)
            self.is_verbose = is_verbose
            acc_timer_flops(fname, 0)

    def __del__(self):
        if self.cdata is None:
            return
        assert isinstance(self.cdata, int)
        cu.free_timer(self)

    def start(self, is_verbose = None):
        if self.cdata is None:
            return
        if is_verbose is None:
            is_verbose = self.is_verbose
        cu.start_timer(self, is_verbose)

    def stop(self, is_verbose = None):
        if self.cdata is None:
            return
        if is_verbose is None:
            is_verbose = self.is_verbose
        cu.stop_timer(self, is_verbose)

    def set_flops(self, flops):
        if self.cdata is None:
            return
        cu.set_flops_timer(self, flops)

def timer(func):
    fname = "py:" + func.__name__
    qtimer = Timer(fname)
    @functools.wraps(func)
    def qtimer_func(*args, **kwargs):
        qtimer.start()
        start_flops = timer_total_flops[fname]
        ret = func(*args, **kwargs)
        stop_flops = timer_total_flops[fname]
        qtimer.set_flops(stop_flops - start_flops)
        qtimer.stop()
        return ret
    return qtimer_func

def timer_verbose(func):
    fname = "py:" + func.__name__
    qtimer = Timer(fname, True)
    @functools.wraps(func)
    def qtimer_func(*args, **kwargs):
        qtimer.start()
        start_flops = timer_total_flops[fname]
        ret = func(*args, **kwargs)
        stop_flops = timer_total_flops[fname]
        qtimer.set_flops(stop_flops - start_flops)
        qtimer.stop()
        return ret
    return qtimer_func

def displayln(level, *args):
    if isinstance(level, int):
        if level <= verbose_level():
            print(*args, flush=True)
    else:
        print(level, *args, flush=True)

def displayln_info(*args):
    if get_id_node() == 0:
        displayln(*args)

class TimerNone(Timer):

    def __init__(self):
        super().__init__(None)

def get_total_time():
    return get_time() - get_start_time()

def get_actual_total_time():
    return get_time() - get_actual_start_time()

def timer_display(*params):
    cu.timer_display(*params)
    flush()

def timer_autodisplay():
    cu.timer_autodisplay()
    flush()

def timer_display_stack():
    cu.timer_display_stack()
    flush()

def timer_display_stack_always():
    cu.timer_display_stack_always()
    flush()
