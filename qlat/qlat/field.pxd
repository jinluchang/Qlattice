from qlat_utils.everything cimport *
from .geometry cimport *

cdef extern from "qlat/field.h" namespace "qlat":

    cdef cppclass Field[T]:
        Field()
        void init()
        void init(const Geometry& geo) except +
        void init(const Geometry& geo, int multiplicity) except +
        void init(const Field[T]& field) except +
        const Field[T]& operator=(const Field[T]& field) except +
        const Geometry& get_geo()

    Vector[T] get_data[T](const Field[T]& x)
    void set_zero[T](Field[T]& x)
    void qswap[T](Field[T]& x, Field[T]& y) except +

cdef extern from "qlat/selected-field-io.h" namespace "qlat":

    cdef cppclass FieldSelection:
        FieldSelection()
        const Geometry& get_geo()

    cdef cppclass PointsSelection:
        PointsSelection()
        PointsSelection(const long n_points) except +
        long size()
        Coordinate* data()
        Coordinate& operator[](long i)

    cdef cppclass SelectedField[T]:
        long n_elems;
        SelectedField()
        void init()
        void init(const Geometry& geo, const long n_elems, const int multiplicity) except +
        void init(const FieldSelection& fsel, const int multiplicity) except +
        const SelectedField[T]& operator=(const SelectedField[T]& field) except +
        const Geometry& get_geo()

    Vector[T] get_data[T](const SelectedField[T]& x)
    void set_zero[T](SelectedField[T]& x)
    void qswap[T](SelectedField[T]& x, SelectedField[T]& y)

    cdef cppclass SelectedPoints[T]:
        int multiplicity
        long n_points
        SelectedPoints()
        void init()
        void init(const long n_points, const int multiplicity) except +
        void init(const PointsSelection& psel, const int multiplicity) except +
        const SelectedPoints[T]& operator=(const SelectedPoints[T]& field) except +

    Vector[T] get_data[T](const SelectedPoints[T]& x)
    void set_zero[T](SelectedPoints[T]& x)
    void qswap[T](SelectedPoints[T]& x, SelectedPoints[T]& y) except +
