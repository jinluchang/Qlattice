from qlat_utils.everything cimport *
from .geometry cimport *

cdef extern from "qlat/field-utils.h" namespace "qlat":

    cdef cppclass Field[T]:
        Field()
        void init()
        void init(const Geometry& geo)
        void init(const Geometry& geo, int multiplicity)
        void init(const Field[T]& field)
        const Field[T]& operator=(const Field[T]& field)
        Geometry& get_geo()
        T* data()

    Vector[T] get_data[T](const Field[T]& x)
    void qswap[T](Field[T]& x, Field[T]& y)
