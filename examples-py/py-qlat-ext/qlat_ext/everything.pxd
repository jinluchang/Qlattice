# all C/C++ level import

from qlat.everything cimport *

cdef extern from "qlat-ext/core.h" namespace "qlat":

    int hello()
