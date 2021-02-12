import cqlat as c

from cqlat import get_id_node, get_num_node, get_size_node, get_coor_node

from cqlat import sync_node

def glb_sum(x):
    if isinstance(x, float):
        return c.glb_sum_double(x)
    elif isinstance(x, complex):
        return c.glb_sum_complex(x)
    elif isinstance(x, int):
        return c.glb_sum_long(x)
    else:
        raise Exception("glb_sum")

