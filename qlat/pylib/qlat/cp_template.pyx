# cython: c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

### -------------------------------------------------------------------

cdef class ElemType:

    pass

cdef class Field:

    pass

cdef class SelectedField:

    pass

cdef class SelectedPoints:

    pass

### -------------------------------------------------------------------

def qremove_info(path):
    return cc.qremove_info(path)

def qremove_all_info(path):
    return cc.qremove_all_info(path)

### -------------------------------------------------------------------
