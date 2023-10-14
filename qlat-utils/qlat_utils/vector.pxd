from libcpp.vector cimport vector as std_vector

cdef extern from "qlat-utils/handle.h" namespace "qlat":

    cdef cppclass Vector[T]:
        Vector()
        Vector(const T* p, const long n)
        T* data()
        long size()
        long data_size()
