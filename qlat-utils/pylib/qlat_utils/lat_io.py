import qlat_utils.c as c

import numpy as np

from qlat_utils.utils_io import *

from qlat_utils.c import LatData

def mk_lat_data(info_list, *, is_complex = True):
    ld = LatData()
    ld.set_info(info_list, is_complex = is_complex)
    return ld

def load_lat_data(path):
    ld = LatData()
    ld.load(path)
    return ld
