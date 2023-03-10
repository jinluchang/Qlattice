from libcpp.string cimport string as std_string

cdef extern from "qlat/utils-io.h" namespace "qlat":

    int qremove_info(const std_string& path)

    int qremove_all_info(const std_string& path)
