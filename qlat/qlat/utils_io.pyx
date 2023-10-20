# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc

from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import cqlat as c
import qlat_utils as q
import numpy as np

def release_lock():
    cc.release_lock()

def obtain_lock(const cc.std_string& path):
    return cc.obtain_lock(path)

def qquit(const cc.std_string& msg):
    """
    clean python cache and then call c.qquit(msg) (which clear all the C++ level cache and then quit)
    """
    q.clean_cache()
    return cc.qquit(msg)
