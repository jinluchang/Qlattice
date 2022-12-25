from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "qlat-utils/core.h" namespace "qlat":

    cdef cppclass Timer:
        long long flops
        Timer()
        Timer(const string& fname)
        void start()
        void start(bool is_verbose)
        void stop()
        void stop(bool is_verbose)
        @staticmethod
        void display(const string& tag)
