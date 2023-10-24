# all cpp level import

from qlat_utils.everything cimport *

cdef extern from "qlat/mpi.h" namespace "qlat":

    int begin(const int id_node, const Coordinate& size_node, const int color) except +
    int end(const bool is_preserving_cache) except +

cdef extern from "qlat/geometry.h" namespace "qlat":

    cdef cppclass GeometryNode:
        bool initialized
        int num_node
        int id_node
        Coordinate size_node
        Coordinate coor_node
        GeometryNode()
        void init()
        void init(const int id_node, const Coordinate& size_node)

    cdef cppclass Geometry:
        bool initialized
        GeometryNode geon
        int eo
        int multiplicity
        Coordinate node_site
        Coordinate expansion_left
        Coordinate expansion_right
        Coordinate node_site_expanded
        bool is_only_local
        Geometry()
        void init()
        void init(Coordinate& total_site, int multiplicity) except +
        Coordinate total_site()
        Coordinate local_site()
        long local_volume()
        long total_volume()
        bool is_local(const Coordinate& x)
        long index_from_coordinate(const Coordinate& xl)
        Coordinate coordinate_from_index(const long index)
        Coordinate coordinate_g_from_l(const Coordinate& xl)
        Coordinate coordinate_l_from_g(const Coordinate& xg)
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
        const Geometry& get_geo()
    cdef cppclass GaugeField:
        pass
    cdef cppclass GaugeTransform:
        pass
    cdef cppclass Prop:
        pass
    cdef cppclass SelProp:
        pass
    cdef cppclass PselProp:
        pass
    cdef cppclass FermionField4d:
        pass
    Vector[T] get_data[T](const Field[T]& x)
    void set_zero[T](Field[T]& x)
    void qswap[T](Field[T]& x, Field[T]& y) except +
    void set_xg_field(Field[int]& f, const Geometry& geo)

cdef extern from "qlat/field-expand.h" namespace "qlat":

    cdef cppclass CommPlan:
        CommPlan()

cdef extern from "qlat/core.h" namespace "qlat":

    cdef cppclass FieldSelection:
        FieldSelection()
        const Geometry& get_geo()
    cdef cppclass PointsSelection:
        PointsSelection()
        PointsSelection(const long n_points) except +
        long size()
        void resize(const long n_points)
        Coordinate* data()
        Coordinate& operator[](long i)
    cdef cppclass SelectedField[T]:
        long n_elems;
        SelectedField()
        void init()
        void init(const Geometry& geo, const long n_elems, const int multiplicity) except +
        void init(const FieldSelection& fsel, const int multiplicity) except +
        const Geometry& get_geo()
    cdef cppclass SelectedPoints[T]:
        int multiplicity
        long n_points
        SelectedPoints()
        void init()
        void init(const long n_points, const int multiplicity) except +
        void init(const PointsSelection& psel, const int multiplicity) except +
    Vector[T] get_data[T](const SelectedField[T]& x)
    void set_zero[T](SelectedField[T]& x)
    void qswap[T](SelectedField[T]& x, SelectedField[T]& y) except +
    Vector[T] get_data[T](const SelectedPoints[T]& x)
    void set_zero[T](SelectedPoints[T]& x)
    void qswap[T](SelectedPoints[T]& x, SelectedPoints[T]& y) except +

cdef extern from "qlat/selected-points.h" namespace "qlat":

    PointsSelection mk_random_point_selection(const Coordinate& total_site,
                                              const long num, const RngState& rs)
    void save_point_selection_info(const PointsSelection& psel,
                                   const std_string& path)
    PointsSelection load_point_selection_info(const std_string& path)

cdef extern from "qlat/selected-field-io.h" namespace "qlat":

    void set_selected_field[t](SelectedField[t]& sf, const Field[t]& f,
                               const FieldSelection& fsel) except +
    void set_selected_field[t](SelectedField[t]& sf, const SelectedField[t] sf0,
                               const FieldSelection& fsel, const FieldSelection& fsel0) except +
    void set_selected_field[t](SelectedField[t]& sf, const SelectedPoints[t] sp,
                               const FieldSelection& fsel, const PointsSelection& psel) except +
    void set_selected_points[t](SelectedPoints[t]& sp, const Field[t] sf,
                                const PointsSelection& psel) except +
    void set_selected_points[t](SelectedPoints[t]& sp, const SelectedField[t] sf,
                                const PointsSelection& psel, const FieldSelection& fsel) except +
    void set_field_selected[t](Field[t]& sf, const SelectedField[t]& f,
                               const FieldSelection& fsel) except +
    void set_field_selected[t](Field[t]& sf, const SelectedPoints[t]& f,
                               const Geometry& geo, const PointsSelection& psel) except +

cdef extern from "qlat/qcd-prop.h" namespace "qlat":

    void set_wall_src(Prop& prop, const Geometry& geo_input,
                      const int tslice, const CoordinateD& lmom) except +

    void set_point_src(Prop& prop, const Geometry& geo_input,
                       const Coordinate& xg, const Complex& value)

cdef extern from "qlat/qcd-smear.h" namespace "qlat":

    void gf_ape_smear(GaugeField& gf, const GaugeField& gf0,
                      const double alpha, const long steps) except +
    void gf_spatial_ape_smear(GaugeField& gf, const GaugeField& gf0,
                              const double alpha, const long steps) except +
    void gf_hyp_smear(GaugeField& gf, const GaugeField& gf0,
                      const double alpha1, const double alpha2, const double alpha3) except +
    void prop_smear(Prop& prop, const GaugeField& gf1,
                    const double coef, const int step,
                    const CoordinateD& mom,
                    const bool smear_in_time_dir) except +

cdef extern from "qlat/vector_utils/utils_smear_vecs.h" namespace "qlat":

    void prop_smear_qlat_convension(Prop& prop, const GaugeField& gf1,
                                    const double coef, const int step,
                                    const CoordinateD& mom,
                                    const bool smear_in_time_dir,
                                    const int mode_smear) except +


