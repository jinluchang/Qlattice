from qlat_utils import *
import qlat.c as c

from qlat.field import *
from qlat.coordinate import *

import numpy as np

def glb_sum_list(ret):
    displayln_info("glb_sum_list: deprecated")
    # deprecated (use glb_sum instead)
    # ret = [ va, vb, ... ]
    # return [ glb_sum(va), glb_sum(vb), ... ]
    return glb_sum(ret)
