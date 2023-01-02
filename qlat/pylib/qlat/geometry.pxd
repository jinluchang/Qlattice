from qlat_utils.everything cimport *

cdef extern from "qlat/geometry.h" namespace "qlat":

    cdef cppclass Geometry:
        int multiplicity
        Geometry()
        const Geometry& operator=(const Geometry& geo)
        void init()
        void init(Coordinate& total_site, int multiplicity)
        Coordinate total_site()
        Coordinate local_site()
        Coordinate local_volume()
        Coordinate total_volume()

    std_string show(const Geometry& geo)

    Geometry geo_resize(const Geometry& geo, int thick)
    Geometry geo_resize(const Geometry& geo,
                        const Coordinate& expansion_left, const Coordinate& expansion_right)
    Geometry geo_reform(const Geometry& geo, int multiplicity, int thick)
    Geometry geo_reform(const Geometry& geo, int multiplicity,
                        const Coordinate& expansion_left, const Coordinate& expansion_right)
    Geometry geo_eo(const Geometry& geo, int eo)
