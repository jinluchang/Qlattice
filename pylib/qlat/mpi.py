import cqlat as c

from cqlat import get_id_node, get_num_node, get_size_node, get_coor_node
from cqlat import sync_node

from qlat.cache import *
from qlat.timer import *

def begin(*args):
    c.begin(*args)

def end():
    q.clean_cache()
    c.end()

@timer
def glb_sum(x):
    if isinstance(x, float):
        return c.glb_sum_double(x)
    elif isinstance(x, complex):
        return c.glb_sum_complex(x)
    elif isinstance(x, int):
        return c.glb_sum_long(x)
    else:
        return x.glb_sum()

@timer_verbose
def show_machine():
    displayln(f"id_node: {get_id_node():4} / {get_num_node()}"
            f" ; coor_node: {str(get_coor_node()):9}"
            f" / {str(get_size_node())}")
