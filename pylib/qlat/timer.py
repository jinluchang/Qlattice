import cqlat as c

from qlat.mpi import *

from cqlat import timer_display
from cqlat import timer_display_stack

def displayln(*args):
    print(*args)

def displayln_info(*args):
    if get_id_node() == 0:
        displayln(*args)

