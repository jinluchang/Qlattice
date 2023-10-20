from qlat_utils import *
from . import c

import os
import pickle

from .c import check_time_limit, check_stop

from .mpi import *

def obtain_lock(path):
    mk_file_dirs_info(path)
    return c.obtain_lock(path)

def release_lock():
    return c.release_lock()

def qquit(msg):
    """
    clean python cache and then call c.qquit(msg) (which clear all the C++ level cache and then quit)
    """
    clean_cache()
    return c.qquit(msg)
