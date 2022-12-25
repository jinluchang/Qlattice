import qlat_utils.c as cu

from qlat_utils.c import get_id_node, get_num_node
from qlat_utils.c import timer_reset, timer_fork, timer_merge
from qlat_utils.c import verbose_level
from qlat_utils.c import get_actual_start_time, get_start_time, get_time
from qlat_utils.c import flush
from qlat_utils.c import timer, timer_verbose
from qlat_utils.c import timer_flops, timer_verbose_flops
from qlat_utils.c import timer_display

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

def timer_autodisplay():
    cu.timer_autodisplay()
    flush()

def timer_display_stack():
    cu.timer_display_stack()
    flush()

def timer_display_stack_always():
    cu.timer_display_stack_always()
    flush()
