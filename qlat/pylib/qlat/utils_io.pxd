from libcpp.string cimport string

cdef extern from "qlat/utils-io.h" namespace "qlat":

    int qremove_info(const string& path)
    int qremove_all_info(const string& path)
