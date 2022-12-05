# cython: c_string_type=unicode, c_string_encoding=utf8

from . cimport __init__ as cp

def qremove_info(path):
    return cp.qremove_info(path)

def qremove_all_info(path):
    return cp.qremove_all_info(path)
