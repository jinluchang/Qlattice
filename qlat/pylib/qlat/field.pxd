from .geometry cimport *

cdef extern from "qlat/field-utils.h" namespace "qlat":

    cdef cppclass Field[T]:
        Field()
        void init()
        void init(const Geometry& geo)
        void init(const Geometry& geo, int multiplicity)
        void init(const Field[T]& field)
        const Field[T]& operator=(const Field[T]& field)

