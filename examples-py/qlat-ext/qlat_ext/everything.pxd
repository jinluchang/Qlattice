from qlat.everything cimport *

cdef extern from "qlat-ext/hello.h" namespace "qlat":

    int hello()
