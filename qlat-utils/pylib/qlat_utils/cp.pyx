# cython: c_string_type=unicode, c_string_encoding=utf8

from . cimport __init__ as cp

def timer_display(tag = ""):
    cp.Timer.display(tag)
