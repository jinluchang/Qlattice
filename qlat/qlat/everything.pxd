# all cpp level import

from qlat_utils.everything cimport *

cdef extern from "qlat/mpi.h" namespace "qlat":

    int begin(const int id_node, const Coordinate& size_node, const int color) except +
    int end(const bool is_preserving_cache) except +

cdef extern from "qlat/geometry.h" namespace "qlat":

    cdef cppclass Geometry:
        int multiplicity
        Geometry()
        const Geometry& operator=(const Geometry& geo) except +
        void init()
        void init(Coordinate& total_site, int multiplicity) except +
        Coordinate total_site()
        Coordinate local_site()
        Coordinate local_volume()
        Coordinate total_volume()
    std_string show(const Geometry& geo) except +
    Geometry geo_resize(const Geometry& geo, int thick) except +
    Geometry geo_resize(const Geometry& geo,
                        const Coordinate& expansion_left, const Coordinate& expansion_right) except +
    Geometry geo_reform(const Geometry& geo, int multiplicity, int thick) except +
    Geometry geo_reform(const Geometry& geo, int multiplicity,
                        const Coordinate& expansion_left, const Coordinate& expansion_right) except +
    Geometry geo_eo(const Geometry& geo, int eo) except +

cdef extern from "qlat/utils-io.h" namespace "qlat":

    void release_lock() except +
    bool obtain_lock(const std_string& path) except +
    void qquit(const std_string& msg) except +
    void check_time_limit(const double budget) except +
    void check_stop(const std_string& fn) except +
    bool does_file_exist_sync_node(const std_string& fn) except +
    bool does_regular_file_exist_qar_sync_node(const std_string& fn) except +
    bool does_file_exist_qar_sync_node(const std_string& fn) except +
    bool is_directory_sync_node(const std_string& fn) except +
    bool is_regular_file_sync_node(const std_string& fn) except +
    int qmkdir_sync_node(const std_string& path) except +
    int qmkdir_p_sync_node(const std_string& path) except +
    std_vector[std_string] qls_sync_node(const std_string& path) except +
    std_vector[std_string] qls_all_sync_node(const std_string& path, const bool is_folder_before_files) except +
    int qremove_info(const std_string& path) except +
    int qremove_all_info(const std_string& path) except +
    int qar_create_info(const std_string& path_qar, const std_string& path_folder_, const bool is_remove_folder_after) except +
    int qar_extract_info(const std_string& path_qar, const std_string& path_folder_, const bool is_remove_qar_after) except +
    int qcopy_file_info(const std_string& path_src, const std_string& path_dst) except +
    std_string qcat_sync_node(const std_string& path) except +
    DataTable qload_datatable_sync_node(const std_string& path, const bool is_par) except +

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
    void set_xg_field(Field[long]& f, const Geometry& geo)

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
    cdef cppclass SelectedPoints[T]:
        int multiplicity
        long n_points
        SelectedPoints()
        void init()
        void init(const long n_points, const int multiplicity) except +
        void init(const PointsSelection& psel, const int multiplicity) except +
        const SelectedPoints[T]& operator=(const SelectedPoints[T]& field) except +
    Vector[T] get_data[T](const SelectedField[T]& x)
    void set_zero[T](SelectedField[T]& x)
    void qswap[T](SelectedField[T]& x, SelectedField[T]& y)
    Vector[T] get_data[T](const SelectedPoints[T]& x)
    void set_zero[T](SelectedPoints[T]& x)
    void qswap[T](SelectedPoints[T]& x, SelectedPoints[T]& y) except +
