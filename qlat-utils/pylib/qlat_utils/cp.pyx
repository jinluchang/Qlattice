# cython: c_string_type=unicode, c_string_encoding=utf8

from .cp cimport *

def timer_display(tag = ""):
    Timer.display(tag)
