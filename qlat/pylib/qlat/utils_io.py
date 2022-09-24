import qlat.cqlat as c

import os
import pickle

from qlat.cqlat import check_time_limit, check_stop
from qlat.cqlat import get_time_limit, get_default_budget

from qlat_utils import *

from qlat.mpi import *

def obtain_lock(path):
    mk_file_dirs_info(path)
    return c.obtain_lock(path)

def release_lock():
    return c.release_lock()

def qquit(msg):
    # clean python cache and then call c.qquit(msg) (which clear all the C++ level cache and then quit)
    clean_cache()
    return c.qquit(msg)

def get_remaining_time():
    return get_time_limit() - get_actual_total_time()
