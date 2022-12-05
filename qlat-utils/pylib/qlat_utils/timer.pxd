from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "qlat-utils/timer.h" namespace "qlat":
    cdef cppclass Timer:
        Timer()
        Timer(string)
        void start(bool verbose = False)
        void stop(bool verbose = False)
        @staticmethod
        void display(string& tag)
