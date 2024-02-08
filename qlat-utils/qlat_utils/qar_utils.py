from .c import qcat
from .c import qcat_bytes
from .timer import *

def qcat_sync_node(path):
    if get_num_node() != 1:
        import qlat.c as c
        return c.qcat_sync_node(path)
    return qcat(path)

def qcat_bytes_sync_node(path):
    if get_num_node() != 1:
        import qlat.c as c
        return c.qcat_bytes_sync_node(path)
    return qcat_bytes(path)
