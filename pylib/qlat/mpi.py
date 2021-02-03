import cqlat as c

from cqlat import get_id_node, get_num_node, get_size_node, get_coor_node

def glb_sum(x):
    if type(x) == float:
        return c.glb_sum_double(x)
    elif type(x) == complex:
        return c.glb_sum_complex(x)
    elif type(x) == int:
        return c.glb_sum_long(x)
    else:
        raise Exception("glb_sum")

