from qlat_utils.rng_state cimport *

cdef extern from "qlat-utils/coordinate.h" namespace "qlat":

    cdef cppclass Coordinate:
        Coordinate()
        Coordinate(int x, int y, int z, int t)
        int& operator[](unsigned long i)

    cdef Coordinate coordinate_from_index(long index, const Coordinate& size)
    cdef long index_from_coordinate(const Coordinate& x, const Coordinate& size)
    cdef int eo_from_coordinate(const Coordinate& xl)

    cdef Coordinate mod(const Coordinate& x, const Coordinate& size)
    cdef Coordinate smod(const Coordinate& x, const Coordinate& size)
    cdef Coordinate middle_mod(const Coordinate& x, const Coordinate& size)

    cdef Coordinate c_rand_gen(RngState& rs, const Coordinate& size)

    cdef cppclass CoordinateD:
        CoordinateD()
        CoordinateD(const Coordinate& x)
        CoordinateD(double x, double y, double z, double t)
        double& operator[](unsigned long i)
