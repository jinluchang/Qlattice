import qlat_utils.c as cu

from qlat_utils.c import flush

from qlat_utils.c import get_id_node, get_num_node
from qlat_utils.c import verbose_level
from qlat_utils.c import get_time, get_start_time, get_actual_start_time
from qlat_utils.c import get_total_time, get_actual_total_time

from qlat_utils.c import timer_display, timer_autodisplay
from qlat_utils.c import timer_display_stack, timer_display_stack_always
from qlat_utils.c import timer_reset, timer_fork, timer_merge

from qlat_utils.c import timer, timer_verbose
from qlat_utils.c import timer_flops, timer_verbose_flops

from qlat_utils.c import Timer, TimerNone

def displayln(level, *args):
    if isinstance(level, int):
        if level <= verbose_level():
            print(*args, flush=True)
    else:
        print(level, *args, flush=True)

def displayln_info(*args):
    if get_id_node() == 0:
        displayln(*args)
