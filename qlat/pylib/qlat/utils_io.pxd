from libcpp.string cimport string

cdef extern from "qlat/utils-io.h" namespace "qlat":
    cdef int qremove_info(const string& path)
    cdef int qremove_all_info(const string& path)
