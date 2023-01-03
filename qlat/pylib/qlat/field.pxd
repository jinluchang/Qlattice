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
        const Geometry& get_geo()

    Vector[T] get_data[T](const Field[T]& x)
    void set_zero[T](Field[T]& x)
    void qswap[T](Field[T]& x, Field[T]& y)

cdef extern from "qlat/selected-field-io.h" namespace "qlat":

    cdef cppclass FieldSelection:
        FieldSelection()
        const Geometry& get_geo()

    cdef cppclass PointSelection:
        PointSelection()
        PointSelection(const long n_points)
        long size()
        Coordinate* data()
        Coordinate& operator[](long i)

    cdef cppclass SelectedField[T]:
        long n_elems;
        SelectedField()
        void init()
        void init(const Geometry& geo, const long n_elems, const int multiplicity)
        void init(const FieldSelection& fsel, const int multiplicity)
        const SelectedField[T]& operator=(const SelectedField[T]& field)
        const Geometry& get_geo()

    Vector[T] get_data[T](const SelectedField[T]& x)
    void set_zero[T](SelectedField[T]& x)
    void qswap[T](SelectedField[T]& x, SelectedField[T]& y)

    cdef cppclass SelectedPoints[T]:
        int multiplicity
        long n_points
        SelectedPoints()
        void init()
        void init(const long n_points, const int multiplicity)
        void init(const PointSelection& psel, const int multiplicity)
        const SelectedPoints[T]& operator=(const SelectedPoints[T]& field)

    Vector[T] get_data[T](const SelectedPoints[T]& x)
    void set_zero[T](SelectedPoints[T]& x)
    void qswap[T](SelectedPoints[T]& x, SelectedPoints[T]& y)
