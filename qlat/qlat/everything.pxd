# all cpp level import

from qlat_utils.everything cimport *

cdef extern from "qlat/mpi.h" namespace "qlat":

    Long& mpi_level_count() except +
    Int begin(const Int id_node, const Coordinate& size_node, const Int color) except +
    Int end(const bool is_preserving_cache) except +
    const Coordinate& get_size_node() except +
    const Coordinate& get_coor_node() except +
    Int bcast(Long& x, const Int root) except +
    Int bcast(RealD& x, const Int root) except +
    Int bcast(RealF& x, const Int root) except +
    Int bcast(ComplexD& x, const Int root) except +
    Int bcast(ComplexF& x, const Int root) except +
    Int bcast(std_string& recv, const Int root) except +
    Int bcast(Coordinate& x, const Int root) except +
    Int bcast(LatData& ld, const Int root) except +
    Int bcast_any(Long& x, const bool b) except +
    Int bcast_any(RealD& x, const bool b) except +
    Int bcast_any(RealF& x, const bool b) except +
    Int bcast_any(ComplexD& x, const bool b) except +
    Int bcast_any(ComplexF& x, const bool b) except +
    Long f_bcast_any(const Long& x, const bool b) except +
    RealD f_bcast_any(const RealD& x, const bool b) except +
    RealF f_bcast_any(const RealF& x, const bool b) except +
    ComplexD f_bcast_any(const ComplexD& x, const bool b) except +
    ComplexF f_bcast_any(const ComplexF& x, const bool b) except +
    std_vector[Int]& get_id_node_list_for_shuffle() except +
    std_vector[Int]& get_id_node_in_shuffle_list() except +
    Int glb_sum(Long& x) except +
    Int glb_sum(RealD& x) except +
    Int glb_sum(RealF& x) except +
    Int glb_sum(ComplexD& x) except +
    Int glb_sum(ComplexF& x) except +
    Int glb_sum(LatData& ld) except +
    bool glb_all(const bool b) except +
    bool glb_any(const bool b) except +
    Long f_glb_sum(const Long& ld) except +
    RealD f_glb_sum(const RealD& ld) except +
    RealF f_glb_sum(const RealF& ld) except +
    ComplexD f_glb_sum(const ComplexD& ld) except +
    ComplexF f_glb_sum(const ComplexF& ld) except +

cdef extern from "qlat/geometry.h" namespace "qlat":

    cdef cppclass GeometryNode:
        bool initialized
        Int num_node
        Int id_node
        Coordinate size_node
        Coordinate coor_node
        GeometryNode()
        void init()
        void init(const Int id_node, const Coordinate& size_node) except +
        void init(const Coordinate& coor_node, const Coordinate& size_node) except +
    cdef cppclass Geometry:
        bool initialized
        GeometryNode geon
        Int eo
        Coordinate node_site
        Coordinate expansion_left
        Coordinate expansion_right
        Coordinate node_site_expanded
        bool is_only_local
        Geometry()
        void init()
        void init(const Coordinate& total_site) except +
        void init(const Coordinate& coor_node, const Coordinate& size_node, const Coordinate& node_site) except +
        void init(const Int& id_node, const Coordinate& size_node, const Coordinate& node_site) except +
        void reset_node_site_expanded() except +
        Coordinate total_site()
        Coordinate local_site()
        Long local_volume() except +
        Long local_volume_expanded() except +
        Long total_volume() except +
        Long spatial_volume() except +
        bool is_local(const Coordinate& x)
        Coordinate coordinate_g_from_l(const Coordinate& xl)
        Coordinate coordinate_l_from_g(const Coordinate& xg)
        Long index_from_coordinate(const Coordinate& xl)
        Coordinate coordinate_from_index(const Long index)
        Long index_from_g_coordinate(const Coordinate& xg)
        Coordinate g_coordinate_from_index(const Long index)
    bool operator==(const Geometry& geo1, const Geometry& geo2) except +
    bool operator!=(const Geometry& geo1, const Geometry& geo2) except +
    std_string show(const Geometry& geo) except +
    Geometry geo_resize(const Geometry& geo, Int thick) except +
    Geometry geo_resize(const Geometry& geo, const Coordinate& expansion_left, const Coordinate& expansion_right) except +
    Geometry geo_eo(const Geometry& geo, Int eo) except +

cdef extern from "qlat/utils-io.h" namespace "qlat":

    void release_lock() except +
    bool obtain_lock(const std_string& path) except +
    void qquit(const std_string& msg) except +
    void check_time_limit(const RealD budget) except +
    void check_stop(const std_string& fn) except +

cdef extern from "qlat/core.h" namespace "qlat":

    const GeometryNode& get_geometry_node()
    cdef cppclass PointsDistType:
        pass
    std_string show(const PointsDistType points_dist_type) except +
    PointsDistType read_points_dist_type(const std_string& points_dist_type_str) except +
    cdef cppclass PointsSelection:
        bool initialized
        PointsDistType points_dist_type
        Coordinate total_site
        PointsSelection()
        PointsSelection(const Coordinate& total_site, const Long n_points) except +
        PointsSelection(const Coordinate& total_site, const Long n_points, const PointsDistType points_dist_type) except +
        void init()
        void init(const Coordinate& total_site, const Long n_points) except +
        void init(const Coordinate& total_site, const Long n_points, const PointsDistType points_dist_type) except +
        Long size()
        Coordinate* data()
        Coordinate& operator[](Long i) except +
    bool operator==(const PointsSelection& psel1, const PointsSelection& psel2) except +
    bool operator!=(const PointsSelection& psel1, const PointsSelection& psel2) except +
    void qswap(PointsSelection& x, PointsSelection& y) except +
    cdef cppclass SelectedPoints[T]:
        bool initialized;
        PointsDistType points_dist_type
        Int multiplicity
        Long n_points
        vector[T] points
        SelectedPoints()
        void init()
        void init(const Long n_points, const Int multiplicity, const PointsDistType points_dist_type) except +
        void init(const PointsSelection& psel, const Int multiplicity) except +
    Int bcast(PointsSelection& psel, const Int root) except +
    Int bcast[T](SelectedPoints[T]& sp, const Int root) except +
    cdef cppclass Field[T]:
        Int multiplicity
        vector[T] field
        Field()
        void init()
        void init(const Geometry& geo) except +
        void init(const Geometry& geo, Int multiplicity) except +
        void init(const Field[T]& field) except +
        Geometry get_geo() except +
        T& get_elem(const Coordinate& x) except +
        T& get_elem(const Coordinate& x, const Int m) except +
        T& get_elem(const Long index) except +
        T& get_elem(const Long index, const Int m) except +
        Vector[T] get_elems(const Coordinate& x) except +
        Vector[T] get_elems(const Long index) except +
    cdef cppclass GaugeField(Field[ColorMatrix]):
        pass
    cdef cppclass GaugeMomentum(Field[ColorMatrix]):
        pass
    cdef cppclass GaugeTransform(Field[ColorMatrix]):
        pass
    cdef cppclass Prop(Field[WilsonMatrix]):
        pass
    cdef cppclass FermionField4d(Field[WilsonVector]):
        pass
    cdef cppclass FieldRank(Field[Int64t]):
        pass
    cdef cppclass FieldIndex(Field[Long]):
        pass
    cdef cppclass FieldSelection:
        Long n_elems
        FieldRank f_rank
        FieldIndex f_local_idx
        vector[Int64t] ranks
        vector[Long] indices
        FieldSelection()
        void init()
        Geometry get_geo() except +
    void set_psel_from_fsel(PointsSelection& psel, const FieldSelection& fsel) except +
    void set_fsel_from_psel(FieldSelection& fsel, const PointsSelection& psel, const Geometry& geo, const Long rank_psel) except +
    void set_geo_from_psel(Geometry& geo, const PointsSelection& psel) except +
    cdef cppclass SelectedField[T]:
        Long n_elems;
        Int multiplicity
        vector[T] field
        SelectedField()
        void init()
        void init(const Geometry& geo, const Long n_elems, const Int multiplicity) except +
        void init(const FieldSelection& fsel, const Int multiplicity) except +
        Geometry get_geo() except +
    Vector[T] get_data[T](const Field[T]& x) except +
    void set_zero[T](Field[T]& x) except +
    void set_unit[T](Field[T]& f, const ComplexD& coef) except +
    void qswap[T](Field[T]& x, Field[T]& y) except +
    Vector[T] get_data[T](const SelectedField[T]& x) except +
    void set_zero[T](SelectedField[T]& x) except +
    void qswap[T](SelectedField[T]& x, SelectedField[T]& y) except +
    Vector[T] get_data[T](const SelectedPoints[T]& x) except +
    void set_zero[T](SelectedPoints[T]& x) except +
    void qswap[T](SelectedPoints[T]& x, SelectedPoints[T]& y) except +
    void qswap_cast[M](PointsSelection& f1, SelectedPoints[M]& f2, Coordinate& total_site2) except +
    void qswap_cast[M,N](Field[M]& f1, Field[N]& f2) except +
    void qswap_cast[M,N](SelectedField[M]& f1, SelectedField[N]& f2) except +
    void qswap_cast[M,N](Field[M]& f1, SelectedPoints[N]& f2, Geometry& geo2) except +
    void qswap_cast[M,N](SelectedField[M]& f1, SelectedPoints[N]& f2, Geometry& geo2) except +
    void qswap_cast[M,N](SelectedPoints[M]& f1, SelectedPoints[N]& f2) except +
    cdef cppclass SelProp(SelectedField[WilsonMatrix]):
        pass
    cdef cppclass PselProp(SelectedPoints[WilsonMatrix]):
        pass
    void set_field_m(Field[Char]& f, const Field[Char]& f1, const Int m, const Int m1, const Int sizeof_m) except +

cdef extern from "qlat/field.h" namespace "qlat":

    RealD qnorm[M](const Field[M]& f) except +
    void qnorm_field[M](Field[RealD]& f, const Field[M]& f1) except +
    void set_u_rand[M](Field[M]& sp, const RngState& rs, const RealD upper, const RealD lower) except +
    void set_g_rand[M](Field[M]& sp, const RngState& rs, const RealD center, const RealD sigma) except +
    void set_xg_field(Field[Int]& f, const Geometry& geo) except +
    void field_shift[M](Field[M]& f, const Field[M]& f1, const Coordinate& shift) except +
    void reflect_field[M](Field[M]& f) except +
    void set_sqrt_field(Field[RealD]& f, const Field[RealD]& f1) except +

cdef extern from "qlat/field-expand.h" namespace "qlat":

    cdef cppclass CommPlan:
        CommPlan()
    void refresh_expanded[M](Field[M]& f, const CommPlan& plan) except +
    void refresh_expanded_1[M](Field[M]& f) except +

cdef extern from *:
    """
    #include <qlat/field-expand.h>
    template <class M>
    inline void py_refresh_expanded_field(qlat::Field<M>& f)
    { qlat::refresh_expanded(f); }
    inline qlat::CommPlan py_make_comm_plan(const qlat::Field<int8_t>& marks)
    { return qlat::make_comm_plan(*(const qlat::CommMarks*)&marks); }
    """
    void py_refresh_expanded_field[M](Field[M]& f) except +
    CommPlan py_make_comm_plan(const Field[Int8t]& marks) except +

cdef extern from "qlat/field-io.h" namespace "qlat":

    bool is_field(const std_string& path) except +
    bool is_d_field(const std_string& path) except +
    bool dist_repartition(const Coordinate& new_size_node, const std_string& path, const std_string& new_path) except +
    crc32_t field_crc32[M](const Field[M]& f) except +
    Long write_field[M](const Field[M]& f, const std_string& path, const Coordinate& new_size_node) except +
    Long read_field[M](Field[M]& f, const std_string& path, const Coordinate& new_size_node) except +
    void convert_field_float_from_double[M, N](Field[N]& ff, const Field[M]& f) except +
    void convert_field_double_from_float[M, N](Field[N]& ff, const Field[M]& f) except +

cdef extern from "qlat/selected-field-io.h" namespace "qlat":

    Long write_field_selection(const FieldSelection& fsel, const std_string& path) except +
    Long read_field_selection(FieldSelection& fsel, const std_string& path) except +
    #
    void mk_field_selection(FieldRank& f_rank, const Geometry& geo, const Int64t val) except +
    void mk_field_selection(FieldRank& f_rank, const Coordinate& total_site, const Long n_per_tslice, const RngState& rs) except +
    void add_field_selection(FieldRank& f_rank, const PointsSelection& psel, const Long rank_psel) except +
    void add_field_selection(FieldRank& f_rank, const FieldSelection& fsel) except +
    void update_field_selection(FieldSelection& fsel) except +
    #
    PointsSelection psel_from_fsel(const FieldSelection& fsel) except +
    PointsSelection psel_from_fsel_local(const FieldSelection& fsel) except +
    bool is_matching_fsel(const FieldSelection& fsel1, const FieldSelection& fsel2) except +
    bool is_containing(const FieldSelection& fsel, const FieldSelection& fsel_small) except +
    bool is_containing(const FieldSelection& fsel, const PointsSelection& psel_small) except +
    bool is_containing(const PointsSelection& psel, const FieldSelection& fsel_small) except +
    bool is_containing(const PointsSelection& psel, const PointsSelection& psel_small) except +
    void intersect_with(FieldSelection& fsel, const FieldSelection& fsel1) except +
    PointsSelection intersect(const FieldSelection& fsel, const PointsSelection& psel) except +

cdef extern from "qlat/selected-shuffle.h" namespace "qlat":

    cdef cppclass SelectedShufflePlan:
        PointsDistType points_dist_type_send
        PointsDistType points_dist_type_recv
        Coordinate total_site
        vector[Coordinate] size_node_send
        vector[Coordinate] coor_node_send
        vector[Coordinate] size_node_recv
        vector[Coordinate] coor_node_recv
        Long num_selected_points_send
        Long num_selected_points_recv
        vector[Long] n_points_selected_points_send
        vector[Long] n_points_selected_points_recv
        SelectedPoints[Long] shuffle_idx_points_send
        SelectedPoints[Long] shuffle_idx_points_recv
        SelectedPoints[Long] shuffle_idx_points_local
        Long total_count_send
        Long total_count_recv
        Long total_count_local
        vector[Long] sendcounts
        vector[Long] recvcounts
        vector[Long] sdispls
        vector[Long] rdispls
        void init() except +
    #
    void set_selected_shuffle_plan_l_from_g(SelectedShufflePlan& ssp, const PointsSelection& psel, const Int root) except +
    void set_selected_shuffle_plan_l_from_g(SelectedShufflePlan& ssp, const std_vector[PointsSelection]& psel_vec, const std_vector[Int]& root_vec) except +
    void set_selected_shuffle_plan_g_from_l(SelectedShufflePlan& ssp, const PointsSelection& psel, const Int root, const Geometry& geo) except +
    void set_selected_shuffle_plan_g_from_l(SelectedShufflePlan& ssp, const std_vector[PointsSelection]& psel_vec, const std_vector[Int]& root_vec, const std_vector[Geometry]& geo_vec) except +
    void set_selected_shuffle_plan_r_from_l(SelectedShufflePlan& ssp, const PointsSelection& psel, const Geometry& geo, const RngState& rs) except +
    void set_selected_shuffle_plan_r_from_l(SelectedShufflePlan& ssp, const std_vector[PointsSelection]& psel_vec, const std_vector[Geometry]& geo_vec, const RngState& rs) except +
    void set_selected_shuffle_plan_dist_r_from_l(SelectedShufflePlan& ssp, const PointsSelection& psel, const Geometry& geo, const RngState& rs, const std_vector[Int]& id_node_vec) except +
    void set_selected_shuffle_plan_t_slice_from_l(SelectedShufflePlan& ssp, const std_vector[PointsSelection]& psel_vec, const std_vector[Geometry]& geo_vec) except +
    void set_selected_shuffle_plan_dist_t_slice_from_l(SelectedShufflePlan& ssp, const PointsSelection& psel, const Geometry& geo, const Int num_field) except +
    #
    void shuffle_selected_points_char(SelectedPoints[Char]& spc, const SelectedPoints[Char]& spc0, const SelectedShufflePlan& ssp) except +
    void shuffle_selected_points_back_char(SelectedPoints[Char]& spc, const SelectedPoints[Char]& spc0, const SelectedShufflePlan& ssp) except +
    void shuffle_selected_points_char(std_vector[SelectedPoints[Char]]& spc_vec, const std_vector[SelectedPoints[Char]]& spc0_vec, const SelectedShufflePlan& ssp) except +
    void shuffle_selected_points_back_char(std_vector[SelectedPoints[Char]]& spc_vec, const std_vector[SelectedPoints[Char]]& spc0_vec, const SelectedShufflePlan& ssp) except +
    #
    void shuffle_points_selection(PointsSelection& psel, const PointsSelection& psel0, const SelectedShufflePlan& ssp) except +
    void shuffle_points_selection_back(PointsSelection& psel, const PointsSelection& psel0, const SelectedShufflePlan& ssp) except +
    void shuffle_points_selection(std_vector[PointsSelection]& psel_vec, const std_vector[PointsSelection]& psel0_vec, const SelectedShufflePlan& ssp) except +
    void shuffle_points_selection_back(std_vector[PointsSelection]& psel_vec, const std_vector[PointsSelection]& psel0_vec, const SelectedShufflePlan& ssp) except +
    #
    void shuffle_selected_points[M](SelectedPoints[M]& sp, const SelectedPoints[M]& sp0, const SelectedShufflePlan& ssp) except +
    void shuffle_selected_points_back[M](SelectedPoints[M]& sp, const SelectedPoints[M]& sp0, const SelectedShufflePlan& ssp) except +
    void shuffle_selected_points[M](std_vector[SelectedPoints[M]]& sp_vec, const std_vector[SelectedPoints[M]]& sp0_vec, const SelectedShufflePlan& ssp) except +
    void shuffle_selected_points_back[M](std_vector[SelectedPoints[M]]& sp_vec, const std_vector[SelectedPoints[M]]& sp0_vec, const SelectedShufflePlan& ssp) except +
    #
    void shuffle_field_selection(PointsSelection& psel, const FieldSelection& fsel0, const SelectedShufflePlan& ssp) except +
    void shuffle_selected_field[M](SelectedPoints[M]& sp, const SelectedField[M]& sf, const SelectedShufflePlan& ssp) except +

cdef extern from "qlat/selected-points.h" namespace "qlat":

    PointsSelection mk_tslice_points_selection(const Coordinate& total_site, const Int t_dir) except +
    PointsSelection mk_random_points_selection(const Coordinate& total_site, const Long num, const RngState& rs) except +
    void set_psel_full(PointsSelection& psel, const Geometry& geo) except +
    void lat_data_from_points_selection(LatDataInt& ld, const PointsSelection& psel) except +
    void points_selection_from_lat_data(PointsSelection& psel, const LatDataInt& ld) except +
    void points_selection_from_lat_data(PointsSelection& psel, const LatDataInt& ld, const PointsDistType points_dist_type) except +
    void save_points_selection(const PointsSelection& psel, const std_string& path) except +
    void save_points_selection_info(const PointsSelection& psel, const std_string& path) except +
    PointsSelection load_points_selection(const std_string& path) except +
    PointsSelection load_points_selection_info(const std_string& path) except +
    #
    RealD qnorm[M](const SelectedPoints[M]& sp) except +
    void qnorm_field[M](SelectedPoints[RealD]& sp, const SelectedPoints[M]& sp1) except +
    void set_u_rand[M](SelectedPoints[M]& sp, const PointsSelection& psel, const RngState& rs, const RealD upper, const RealD lower) except +
    void set_g_rand[M](SelectedPoints[M]& sp, const PointsSelection& psel, const RngState& rs, const RealD center, const RealD sigma) except +
    void lat_data_from_selected_points[M](LatData& ld, const SelectedPoints[M]& sp) except +
    void selected_points_from_lat_data[M](SelectedPoints[M]& sp, const LatData& ld) except +
    void lat_data_from_selected_points[M](LatDataRealF& ld, const SelectedPoints[M]& sp) except +
    void selected_points_from_lat_data[M](SelectedPoints[M]& sp, const LatDataRealF& ld) except +
    void lat_data_from_selected_points[M](LatDataLong& ld, const SelectedPoints[M]& sp) except +
    void selected_points_from_lat_data[M](SelectedPoints[M]& sp, const LatDataLong& ld) except +
    void lat_data_from_selected_points[M](LatDataInt& ld, const SelectedPoints[M]& sp) except +
    void selected_points_from_lat_data[M](SelectedPoints[M]& sp, const LatDataInt& ld) except +
    void save_selected_points[M](const SelectedPoints[M]& sp, QFile& qfile, const bool is_sync_node) except +
    void load_selected_points[M](SelectedPoints[M]& sp, QFile& qfile, const bool is_sync_node) except +
    void save_selected_points[M](const SelectedPoints[M]& sp, const std_string& fn, const bool is_sync_node) except +
    void load_selected_points[M](SelectedPoints[M]& sp, const std_string& fn, const bool is_sync_node) except +
    std_string save_selected_points_str[M](const SelectedPoints[M]& sp, const bool is_sync_node) except +
    void load_selected_points_str[M](SelectedPoints[M]& sp, std_string& content, const bool is_sync_node) except +
    void field_glb_sum[M](SelectedPoints[M]& sp, const Field[M]& f) except +
    void field_glb_sum_tslice[M](SelectedPoints[M]& sp, const Field[M]& f, const Int t_dir) except +
    void field_glb_max[M](SelectedPoints[M]& sp, const Field[M]& f) except +
    void field_glb_min[M](SelectedPoints[M]& sp, const Field[M]& f) except +
    void set_sqrt_field(SelectedPoints[RealD]& sp, const SelectedPoints[RealD]& sp1) except +
    void acc_field[M](Field[M]& f, const SelectedPoints[M]& sp,
            const Geometry& geo, const PointsSelection& psel) except +

cdef extern from "qlat/selected-field.h" namespace "qlat":

    void field_glb_sum_tslice[M](SelectedPoints[M]& sp, const SelectedField[M]& sf, const FieldSelection& fsel, const Int t_dir) except +
    void acc_field[M](Field[M]& f, const SelectedField[M]& sf,
            const FieldSelection& fsel) except +

cdef extern from "qlat/selected-field-io.h" namespace "qlat":

    RealD qnorm[M](const SelectedField[M]& sf) except +
    void qnorm_field[M](SelectedField[RealD]& f, const SelectedField[M]& f1) except +
    void convert_field_float_from_double[M, N](SelectedField[N]& ff, const SelectedField[M]& f) except +
    void convert_field_double_from_float[M, N](SelectedField[N]& ff, const SelectedField[M]& f) except +
    void set_u_rand[M](SelectedField[M]& sp, const FieldSelection& fsel, const RngState& rs, const RealD upper, const RealD lower) except +
    void set_g_rand[M](SelectedField[M]& sp, const FieldSelection& fsel, const RngState& rs, const RealD center, const RealD sigma) except +
    void set_selected_field[t](SelectedField[t]& sf, const Field[t]& f,
                               const FieldSelection& fsel) except +
    void set_selected_field[t](SelectedField[t]& sf, const SelectedField[t] sf0,
                               const FieldSelection& fsel, const FieldSelection& fsel0) except +
    void set_selected_field[t](SelectedField[t]& sf, const SelectedPoints[t] sp,
                               const FieldSelection& fsel, const PointsSelection& psel) except +
    void set_selected_points[t](SelectedPoints[t]& sp, const Field[t] f,
                                const PointsSelection& psel) except +
    void set_selected_points[t](SelectedPoints[t]& sp, const SelectedField[t] sf,
                                const PointsSelection& psel, const FieldSelection& fsel) except +
    void set_selected_points[t](SelectedPoints[t]& sp, const SelectedPoints[t] sp0,
                                const PointsSelection& psel, const PointsSelection& psel0) except +
    void set_field_selected[t](Field[t]& f, const SelectedField[t]& sf,
                               const FieldSelection& fsel) except +
    void set_field_selected[t](Field[t]& f, const SelectedPoints[t]& sp,
                               const Geometry& geo, const PointsSelection& psel) except +
    void set_selected_points[t](SelectedPoints[t]& sp, const Field[t] f,
                                const PointsSelection& psel, const Int m) except +
    void set_field_selected[t](Field[t]& f, const SelectedPoints[t]& sp,
                               const Geometry& geo, const PointsSelection& psel,
                               const Int multiplicity, const Int m) except +
    void set_sqrt_field(SelectedField[RealD]& f, const SelectedField[RealD]& f1) except +
    #
    bool is_selected_field(const std_string& path) except +
    Long write_selected_field[M](const SelectedField[M]& sf, const std_string& path, const FieldSelection& fsel) except +
    Long read_selected_field[M](SelectedField[M]& sf, const std_string& path, const FieldSelection& fsel) except +

cdef extern from "qlat/field-shuffle.h" namespace "qlat":

    void field_shift[M](SelectedField[M]& sf, FieldSelection& fsel,
            const SelectedField[M]& sf0, const FieldSelection& fsel0,
            const Coordinate& shift, const bool is_reflect) except +
    void shuffle_field[M](std_vector[Field[M]]& fs, const Field[M]& f, const Coordinate& new_size_node) except +
    void shuffle_field_back[M](Field[M]& f, const std_vector[Field[M]]& fs, const Coordinate& new_size_node) except +

cdef extern from "qlat/gauge-action.h" namespace "qlat":

    cdef cppclass GaugeAction:
        bool initialized
        RealD beta
        RealD c1
        void init() except +
        GaugeAction() except +

cdef extern from "qlat/qcd-prop.h" namespace "qlat":

    void set_ff_vec_from_prop(std_vector[FermionField4d]& ff_vec, const Prop& prop) except +
    void set_prop_from_ff_vec(Prop& prop, const std_vector[FermionField4d]& ff_vec) except +
    void set_wall_src(Prop& prop, const Geometry& geo_input, const Int tslice, const CoordinateD& lmom) except +
    void set_point_src(Prop& prop, const Geometry& geo_input, const Coordinate& xg, const ComplexD& value) except +
    void set_rand_vol_u1(Field[ComplexD]& fu1, const Geometry& geo_input, const Int multiplicity, const RngState& rs) except +
    void set_rand_vol_u1_src(Prop& prop, const Field[ComplexD]& fu1) except +
    void free_invert(
        Prop& p_sol, const Prop& p_src,
        const RealD mass, const RealD m5,
        const CoordinateD& momtwist) except +
    void free_mom_invert(
        Prop& p_sol, const Prop& p_src,
        const RealD mass, const RealD m5,
        const CoordinateD& momtwist) except +
    void convert_wm_from_mspincolor(Prop& prop_wm,
            const Prop& prop_msc) except +
    void convert_mspincolor_from_wm(Prop& prop_msc,
            const Prop& prop_wm) except +
    void convert_wm_from_mspincolor(SelectedField[WilsonMatrix]& prop_wm,
            const SelectedField[WilsonMatrix]& prop_msc) except +
    void convert_mspincolor_from_wm(SelectedField[WilsonMatrix]& prop_msc,
            const SelectedField[WilsonMatrix]& prop_wm) except +
    void convert_wm_from_mspincolor(SelectedPoints[WilsonMatrix]& prop_wm,
            const SelectedPoints[WilsonMatrix]& prop_msc) except +
    void convert_mspincolor_from_wm(SelectedPoints[WilsonMatrix]& prop_msc,
            const SelectedPoints[WilsonMatrix]& prop_wm) except +
    void flip_tpbc_with_tslice(SelectedPoints[WilsonMatrix]& ps_prop,
            const PointsSelection& psel, const Int tslice_flip_tpbc,
            const Int t_size) except +
    void flip_tpbc_with_tslice(SelectedField[WilsonMatrix]& s_prop,
            const FieldSelection& fsel, const Int tslice_flip_tpbc) except +

cdef extern from *:
    """
    #include <qlat/qcd-prop.h>
    inline void py_set_rand_u1_src_psel(qlat::Propagator4d& prop,
        qlat::Field<qlat::ComplexD>& fu1, const qlat::PointsSelection& psel,
        const qlat::Geometry& geo, const qlat::RngState& rs)
    {
      prop.init();
      fu1.init();
      qlat::set_rand_u1_src_psel(prop,
          *(qlat::FieldM<qlat::ComplexD, 1>*)&fu1, psel, geo, rs);
    }
    inline void py_set_rand_u1_sol_psel(
        qlat::SelectedPoints<qlat::WilsonMatrix>& sp_prop,
        const qlat::Propagator4d& prop, const qlat::Field<qlat::ComplexD>& fu1,
        const qlat::PointsSelection& psel)
    {
      qlat::set_rand_u1_sol_psel(sp_prop, prop,
          *(const qlat::FieldM<qlat::ComplexD, 1>*)&fu1, psel);
    }
    inline void py_set_rand_u1_src_fsel(qlat::Propagator4d& prop,
        qlat::Field<qlat::ComplexD>& fu1, const qlat::FieldSelection& fsel,
        const qlat::RngState& rs)
    {
      prop.init();
      fu1.init();
      qlat::set_rand_u1_src_fsel(prop,
          *(qlat::FieldM<qlat::ComplexD, 1>*)&fu1, fsel, rs);
    }
    inline void py_set_rand_u1_sol_fsel(
        qlat::SelectedField<qlat::WilsonMatrix>& sf_prop,
        const qlat::Propagator4d& prop, const qlat::Field<qlat::ComplexD>& fu1,
        const qlat::FieldSelection& fsel)
    {
      qlat::set_rand_u1_sol_fsel(sf_prop, prop,
          *(const qlat::FieldM<qlat::ComplexD, 1>*)&fu1, fsel);
    }
    """
    void py_set_rand_u1_src_psel(Prop& prop, Field[ComplexD]& fu1,
            const PointsSelection& psel, const Geometry& geo,
            const RngState& rs) except +
    void py_set_rand_u1_sol_psel(PselProp& sp_prop, const Prop& prop,
            const Field[ComplexD]& fu1, const PointsSelection& psel) except +
    void py_set_rand_u1_src_fsel(Prop& prop, Field[ComplexD]& fu1,
            const FieldSelection& fsel, const RngState& rs) except +
    void py_set_rand_u1_sol_fsel(SelProp& sf_prop, const Prop& prop,
            const Field[ComplexD]& fu1, const FieldSelection& fsel) except +

cdef extern from "qlat/qcd-smear.h" namespace "qlat":

    void gf_ape_smear(GaugeField& gf, const GaugeField& gf0, const RealD alpha, const Long steps) except +
    void gf_spatial_ape_smear(GaugeField& gf, const GaugeField& gf0, const RealD alpha, const Long steps) except +
    void gf_hyp_smear(GaugeField& gf, const GaugeField& gf0, const RealD alpha1, const RealD alpha2, const RealD alpha3) except +
    void prop_smear(Prop& prop, const GaugeField& gf1, const RealD coef, const Int step, const CoordinateD& mom, const bool smear_in_time_dir) except +
    void prop_spatial_smear_no_comm(std_vector[FermionField4d]& ff_vec, const GaugeField& gf, const RealD coef, const Long step, const CoordinateD& mom) except +
    void gf_reduce_half(GaugeField& hgf, const GaugeField& gf) except +
    void prop_smear_qlat_convension(Prop& prop, const GaugeField& gf1, const RealD coef, const Int step, const CoordinateD& mom, const bool smear_in_time_dir, const Int mode_smear) except +

cdef extern from "qlat/qcd.h" namespace "qlat":

    RealD gf_avg_plaq(const GaugeField& gf) except +
    RealD gf_avg_spatial_plaq(const GaugeField& gf) except +
    RealD gf_avg_link_trace(const GaugeField& gf) except +
    void gf_plaq_field(Field[RealD]& f_plaq, const GaugeField& gf) except +
    RealD gf_local_avg_plaq(const GaugeField& gf, const Coordinate& block_site) except +
    void unitarize(Field[ColorMatrix]& gf) except +
    void make_tr_less_anti_herm_matrix(Field[ColorMatrix]& fc) except +
    Long save_gauge_field "qlat::save_gauge_field<qlat::RealD>"(
            const GaugeField& gf, const std_string& path) except +
    Long load_gauge_field "qlat::load_gauge_field<qlat::RealD>"(
            GaugeField& gf, const std_string& path) except +
    Long save_gauge_transform_cps(const GaugeTransform& gt,
            const std_string& path) except +
    Long load_gauge_transform_cps(GaugeTransform& gt,
            const std_string& path) except +
    void twist_boundary_at_boundary(GaugeField& gf, const RealD lmom,
            const Int mu) except +

cdef extern from *:
    """
    #include <qlat/qcd.h>
    inline void py_twist_boundary_at_boundary(qlat::GaugeField& gf,
        const double lmom, const qlat::Int mu)
    { qlat::twist_boundary_at_boundary(gf, lmom, mu); }
    """
    void py_twist_boundary_at_boundary(GaugeField& gf, const RealD lmom,
            const Int mu) except +

cdef extern from "qlat/qcd-utils.h" namespace "qlat":

    void set_left_expanded_gauge_field(GaugeField& gf1, const GaugeField& gf) except +
    ColorMatrix gf_avg_wilson_loop(const GaugeField& gf, const Int l, const Int t) except +
    void set_g_rand_color_matrix_field(Field[ColorMatrix]& fc, const RngState& rs, const RealD sigma, const Int n_step) except +
    void gf_wilson_line_no_comm(Field[ColorMatrix]& wilson_line_field,
            const Int wilson_line_field_m, const GaugeField& gf_ext,
            const std_vector[Int]& path) except +
    void gf_wilson_line_no_comm(Field[ColorMatrix]& wilson_line_field,
            const Int wilson_line_field_m, const GaugeField& gf_ext,
            const std_vector[Int]& path,
            const std_vector[Int]& path_n) except +

cdef extern from "qlat/qcd-gauge-transformation.h" namespace "qlat":

    void gt_apply_gauge_transformation(GaugeTransform& gt0, const GaugeTransform& gt1) except +
    void gt_apply_gauge_transformation(GaugeTransform& gt, const GaugeTransform& gt0, const GaugeTransform& gt1) except +
    void gf_apply_gauge_transformation(GaugeField& gf, const GaugeField& gf0, const GaugeTransform& gt, const bool is_dagger) except +
    void gt_invert(GaugeTransform& gt, const GaugeTransform& gt0) except +
    void ff_apply_gauge_transformation(FermionField4d& ff, const FermionField4d& ff0, const GaugeTransform& gt) except +
    void prop_apply_gauge_transformation(Prop& prop, const Prop& prop0, const GaugeTransform& gt) except +
    void prop_apply_gauge_transformation(SelectedField[WilsonMatrix]& prop, const SelectedField[WilsonMatrix]& prop0, const GaugeTransform& gt, const FieldSelection& fsel) except +
    void prop_apply_gauge_transformation(SelectedPoints[WilsonMatrix]& prop, const SelectedPoints[WilsonMatrix]& prop0, const GaugeTransform& gt, const PointsSelection& psel) except +

cdef extern from "qlat/hmc.h" namespace "qlat":

    bool metropolis_accept(RealD& accept_prob, const RealD delta_h, const int traj, const RngState& rs_) except +
    void set_rand_gauge_momentum(GaugeMomentum& gm, const RealD sigma, const RngState& rs) except +
    void set_rand_gauge_momentum(GaugeMomentum& gm, const Field[RealD]& mf, const RngState& rs) except +
    RealD gm_hamilton_node(const GaugeMomentum& gm) except +
    RealD gm_hamilton_node(const GaugeMomentum& gm, const Field[RealD]& mf) except +
    RealD gf_hamilton_node(const GaugeField& gf, const GaugeAction& ga) except +
    RealD gf_hamilton(const GaugeField& gf, const GaugeAction& ga) except +
    void gf_evolve(GaugeField& gf, const GaugeMomentum& gm, const RealD step_size) except +
    void gf_evolve_dual(GaugeField& gf, const GaugeMomentum& gm_dual, const RealD step_size) except +
    void gf_evolve(GaugeField& gf, const GaugeMomentum& gm, const Field[RealD]& mf, const RealD step_size) except +
    void gf_evolve_dual(GaugeField& gf, const GaugeMomentum& gm_dual, const Field[RealD]& mf_dual, const RealD step_size) except +
    void field_color_matrix_exp(Field[ColorMatrix]& fc, const Field[ColorMatrix]& fc1, const ComplexD coef) except +
    void field_color_matrix_mul(Field[ColorMatrix]& fc, const Field[ColorMatrix]& fc1, const Field[ColorMatrix]& fc2) except +
    void set_gm_force(GaugeMomentum& gm_force, const GaugeField& gf, const GaugeAction& ga) except +
    void set_gm_force_dual(GaugeMomentum& gm_force_dual, const GaugeField& gf, const GaugeMomentum& gm_force) except +
    RealD project_gauge_transform(GaugeMomentum& gm, GaugeMomentum& gm_dual, const Field[RealD]& mf, const Field[RealD]& mf_dual) except +
    void set_gauge_transform_momentum(GaugeMomentum& gm, GaugeMomentum& gm_dual, const Field[ColorMatrix]& gtm) except +
    void dot_gauge_momentum(Field[RealD]& f, const GaugeMomentum& gm1, const GaugeMomentum& gm2) except +
    void set_anti_hermitian_matrix_from_basis(Field[ColorMatrix]& fc, const Field[RealD]& basis) except +
    void set_basis_from_anti_hermitian_matrix(Field[RealD]& basis, const Field[ColorMatrix]& fc) except +

cdef extern from "qlat/fields-io.h" namespace "qlat":

    cdef cppclass ShuffledBitSet:
        FieldSelection fsel
        std_vector[FieldSelection] fsels
        ShuffledBitSet()
    ShuffledBitSet mk_shuffled_bitset(const FieldRank& f_rank, const Coordinate& new_size_node) except +
    ShuffledBitSet mk_shuffled_bitset(const FieldSelection& fsel, const Coordinate& new_size_node) except +
    ShuffledBitSet mk_shuffled_bitset(const Coordinate& total_site, const PointsSelection& xgs, const Coordinate& new_size_node) except +
    ShuffledBitSet mk_shuffled_bitset(const FieldRank& f_rank, const PointsSelection& xgs, const Coordinate& new_size_node) except +
    cdef cppclass ShuffledFieldsWriter:
        std_string path
        Coordinate new_size_node
        ShuffledFieldsWriter()
        void init() except +
        void init(const std_string& path_, const Coordinate& new_size_node_, const bool is_append) except +
        void close() except +
    cdef cppclass ShuffledFieldsReader:
        std_string path
        Coordinate new_size_node
        ShuffledFieldsReader()
        void init() except +
        void init(const std_string& path_, const Coordinate& new_size_node_) except +
        void close() except +
    Long write[M](ShuffledFieldsWriter& sfw, const std_string& fn, const Field[M]& field) except +
    Long write[M](ShuffledFieldsWriter& sfw, const std_string& fn, const ShuffledBitSet& sbs, const SelectedField[M]& sf) except +
    Long read[M](ShuffledFieldsReader& sfr, const std_string& fn, Field[M]& field) except +
    Long read[M](ShuffledFieldsReader& sfr, const std_string& fn, SelectedField[M]& sf, FieldSelection& fsel) except +
    Long read[M](ShuffledFieldsReader& sfr, const std_string& fn, const ShuffledBitSet& sbs, SelectedField[M]& sf) except +
    Long flush(ShuffledFieldsWriter& sfw) except +
    void read_through_sync_node(ShuffledFieldsReader& sfr) except +
    bool does_file_exist_sync_node(const ShuffledFieldsReader& sfr, const std_string& fn) except +
    bool does_file_exist_sync_node(const ShuffledFieldsWriter& sfw, const std_string& fn) except +
    bool is_sparse_field_sync_node(const ShuffledFieldsReader& sfr, const std_string& fn) except +
    std_vector[std_string] list_fields(const ShuffledFieldsReader& sfr) except +
    std_vector[std_string] list_fields(const ShuffledFieldsWriter& sfw) except +
    Int truncate_fields_sync_node(const std_string& path, const std_vector[std_string]& fns_keep, const Coordinate& new_size_node) except +
    std_vector[std_string] properly_truncate_fields_sync_node(const std_string& path, const bool is_check_all, const bool is_only_check, const Coordinate& new_size_node) except +
    bool has_duplicates(const ShuffledFieldsReader& sfr) except +
    void fields_build_index(const ShuffledFieldsReader& sfr) except +
    void fields_build_index(const std_string& path) except +
    #
    std_vector[std_string] show_all_shuffled_fields_writer() except +

cdef extern from "qlat/compressed-eigen-io.h" namespace "qlat":

    bool check_compressed_eigen_vectors(const std_string& path) except +
    bool eigen_system_repartition(const Coordinate& new_size_node, const std_string& path, const std_string& new_path) except +

cdef extern from "qlat/wilson-flow.h" namespace "qlat":

    void set_wilson_flow_z(GaugeMomentum& z, const GaugeField& gf, const RealD c1) except +
    void gf_wilson_flow_step_euler(GaugeField& gf, const RealD epsilon, const RealD c1) except +
    void gf_wilson_flow_step(GaugeField& gf, const RealD epsilon, const RealD c1) except +
    void gf_energy_density_dir_field(Field[RealD]& fd, const GaugeField& gf) except +
    void gf_energy_density_field(Field[RealD]& fd, const GaugeField& gf) except +
    RealD gf_energy_density(const GaugeField& gf) except +
    void set_plaq_flow_z(GaugeMomentum& z, const GaugeField& gf, const Field[RealD]& plaq_factor) except +
    void gf_block_stout_smear(GaugeField& gf, const GaugeField& gf0, const Coordinate& block_site, const RealD step_size) except +
    void gf_local_stout_smear(GaugeField& gf, const GaugeField& gf0, const Coordinate& block_site, const RealD step_size) except +
    void set_local_tree_gauge_f_dir(Field[Int]& f_dir, const Geometry& geo, const Coordinate& block_site, const bool is_uniform, const RngState& rs) except +
    void gt_local_tree_gauge(GaugeTransform& gt_inv, const GaugeField& gf, const Field[Int]& f_dir, const Int num_step) except +

cdef extern from "qlat/qcd-topology.h" namespace "qlat":

    RealD topology_charge_5(const GaugeField& gf) except +
    void clf_plaq_action_density_field(Field[RealD]& topf, const GaugeField& gf) except +
    void clf_spatial_plaq_action_density_field(Field[RealD]& topf, const GaugeField& gf) except +
    void clf_topology_field(Field[RealD]& topf, const GaugeField& gf) except +
    void clf_topology_field_5(Field[RealD]& topf, const GaugeField& gf) except +
    void clf_topology_field_5_terms(Field[RealD]& topf, const GaugeField& gf) except +

cdef extern from "qlat/muon-line.h" namespace "qlat":

    bool has_cuba() except +
    #
    void test_integrationMultidimensional() except +
    #
    void clear_muon_line_interpolations() except +
    Long get_number_of_muon_line_interpolations() except +
    #
    cdef cppclass IntegrationEps:
        RealD epsabs
        RealD epsrel
        Long mineval
        Long maxeval
        IntegrationEps()
    # epsabs = 1e-8
    # epsrel = 1e-3
    # mineval = 1024 * 1024
    # maxeval = 1024 * 1024 * 1024
    # dims = [ 6, 6, 6, 6, 6, ] ~ [ 16, 16, 16, 16, 16, ]
    bool compute_save_muonline_interpolation_cc(
            const std_string& path,
            const std_vector[Int]& dims,
            const IntegrationEps& eps) except +
    std_vector[RealD] muon_line_sym_py(
            const CoordinateD& x,
            const CoordinateD& y,
            const IntegrationEps& eps) except +
    #
    bool load_multiple_muonline_interpolations(const std_string& path, std_vector[Long] idx_list) except +
    #
    std_vector[std_vector[RealD]] get_muon_line_m_extra_weights_default() except +
    std_vector[std_vector[RealD]]& get_muon_line_m_extra_weights() except +
    void set_muon_line_m_extra_weights(const std_vector[std_vector[RealD]]& weight) except +
    #
    std_vector[RealD] get_muon_line_m_py(
            const CoordinateD& x, const CoordinateD& y,
            const CoordinateD& z,
            const Int idx,
            const IntegrationEps& eps) except +
    std_vector[RealD] get_muon_line_m_extra_py(const CoordinateD& x,
            const CoordinateD& y,
            const CoordinateD& z,
            const Int tag) except +
    std_vector[RealD] get_muon_line_m_extra_lat_py(
            const Coordinate& x, const Coordinate& y, const Coordinate& z,
            const Coordinate& total_site, const RealD a, const Int tag) except +

cdef extern from "qlat/hlbl-contract.h" namespace "qlat":

    cdef cppclass SlTable:
        Int s_limit
        Int l_limit
        std_vector[ComplexD] table
        void init() except +
        void init(const Int s, const Int l) except +
        void init(const Coordinate& total_site) except +
    cdef cppclass CurrentMoments[T]:
        void init() except +
        void init(const Int lsize) except +
    #
    void set_m_z_field_tag(SelectedPoints[RealD]& smf_d, const PointsSelection& psel_d, const Geometry& geo, const Coordinate& xg_x, const Coordinate& xg_y, const RealD a, const Int tag) except +
    #
    void set_local_current_from_props(SelectedPoints[WilsonMatrix]& scf, const SelectedPoints[WilsonMatrix]& sprop1, const SelectedPoints[WilsonMatrix]& sprop2, const PointsSelection& psel_d, const Geometry& geo) except +
    #
    RealD set_psel_d_prob_xy(SelectedPoints[RealD]& psel_d_prob_xy, const PointsSelection& psel, SelectedPoints[RealD]& psel_prob, const PointsSelection& psel_d, const SelectedPoints[RealD]& psel_d_prob, const Long idx_xg_x, const Long idx_xg_y) except +
    #
    void set_current_moments_from_current(CurrentMoments[WilsonMatrix]& cm, const SelectedPoints[WilsonMatrix]& current, const PointsSelection& psel_d, const SelectedPoints[RealD]& psel_d_prob_xy, const Geometry& geo) except +
    #
    void glb_sum_current_moments(CurrentMoments[WilsonMatrix]& cm) except +
    #
    std_vector[std_string] contract_four_pair_labels(const std_vector[std_string]& tags) except +
    #
    std_vector[SlTable] contract_four_pair_no_glb_sum(const ComplexD& coef, const PointsSelection& psel, const PointsSelection& psel_d, const SelectedPoints[RealD]& psel_d_prob_xy, const Geometry& geo, const Long idx_xg_x, const Long idx_xg_y, const SelectedPoints[RealD]& smf_d, const SelectedPoints[WilsonMatrix]& sc_xy, const SelectedPoints[WilsonMatrix]& sc_yx, const CurrentMoments[WilsonMatrix]& cm_xy, const CurrentMoments[WilsonMatrix]& cm_yx, const Int inv_type, const std_vector[std_string]& tags, const Long r_sq_limit, const RealD muon_mass, const RealD z_v) except +
    #
    std_vector[std_string] contract_two_plus_two_pair_labels() except +
    #
    std_vector[SlTable] contract_two_plus_two_pair_no_glb_sum(Long& n_points_selected, Long& n_points_computed, const ComplexD& coef, const Geometry& geo, const PointsSelection& psel, const SelectedPoints[RealD]& psel_prob, const PointsSelection& psel_lps, const SelectedPoints[RealD]& psel_lps_prob, const Long idx_xg_x, const SelectedPoints[ComplexD]& lps_hvp_x, const SelectedPoints[ComplexD]& edl_list_c, const Long r_sq_limit, const RealD muon_mass, const RealD z_v) except +
    #

cdef extern from "qlat/qed.h" namespace "qlat":

    cdef cppclass SpinProp(Field[SpinMatrix]):
        pass
    #
    void free_invert(
        SpinProp& sp_sol, const SpinProp& sp_src,
        const RealD mass, const RealD m5,
        const CoordinateD& momtwist) except +
    #
    void free_mom_invert(
        SpinProp& sp_sol, const SpinProp& sp_src,
        const RealD mass, const RealD m5,
        const CoordinateD& momtwist) except +
    #
    void free_scalar_mom_invert(
        Field[ComplexD]& f, const RealD mass,
        const CoordinateD& momtwist) except +
    #
    void free_scalar_deriv_mom(
        Field[ComplexD]& f,
        const vector[Int]& deriv_order,
        const CoordinateD& momtwist) except +
    #
    void invert_qed(
        SpinProp& sp_sol, const SpinProp& sp_src, const Field[ComplexD]& gf1,
        const RealD mass, const RealD m5, const Int ls,
        const vector[ComplexD]& t_wick_phase_factor_vec,
        const bool is_dagger,
        const RealD stop_rsd, const Long max_num_iter) except +
    #
    Long invert_dwf_qed(
        Field[ComplexD]& f_out4d, const Field[ComplexD]& f_in4d,
        const Field[ComplexD]& gf1, const RealD mass,
        const RealD m5, const Int ls,
        const vector[ComplexD]& t_wick_phase_factor_vec,
        const bool is_dagger,
        const RealD stop_rsd, const Long max_num_iter) except +
    #
    Long cg_with_m_dwf_qed(
        Field[ComplexD]& f_out5d, const Field[ComplexD]& f_in5d,
        const Field[ComplexD]& gf1, const RealD mass,
        const RealD m5, const Int ls,
        const vector[ComplexD]& t_wick_phase_factor_vec,
        const bool is_dagger,
        const RealD stop_rsd, const Long max_num_iter) except +
    #
    void multiply_m_dwf_qed(
        Field[ComplexD]& f_out, const Field[ComplexD]& f_in,
        const Field[ComplexD]& gf1,
        const RealD mass, const RealD m5, const Int ls,
        const vector[ComplexD]& t_wick_phase_factor_vec,
        const bool is_dagger) except +

cdef extern from "qlat/flowed-hmc.h" namespace "qlat":

    cdef cppclass FlowStepInfo:
        Int mask
        Int mu
        RealD epsilon
        Int flow_size
        FlowStepInfo()
        FlowStepInfo(const Int mask_, const Int mu_, const RealD epsilon_,
                const Int flow_size_) except +
    cdef cppclass FlowInfo:
        std_vector[FlowStepInfo] v
        FlowInfo()
        FlowInfo(const FlowInfo&) except +
    std_string show(const FlowInfo& fi) except +
    FlowInfo mk_flow_info_step(const RngState& rs, const RealD epsilon) except +
    FlowInfo mk_flow_info_step(const RngState& rs, const RealD epsilon,
            const RealD epsilon2) except +
    void gf_flow(GaugeField& gf, const GaugeField& gf0, const FlowInfo& fi) except +
    void gf_flow_inv(GaugeField& gf, const GaugeField& gf1, const FlowInfo& fi) except +
    RealD gf_hamilton_flowed_node(const GaugeField& gf0, const GaugeAction& ga,
            const FlowInfo& fi) except +
    void set_gm_force_flowed(GaugeMomentum& gm_force, const GaugeField& gf0,
            const GaugeAction& ga, const FlowInfo& fi) except +

cdef extern from "qlat/scalar-action.h" namespace "qlat":

    cdef cppclass ScalarAction:
        bool initialized
        RealD m_sq
        RealD lmbd
        RealD alpha
        ScalarAction() except +
        ScalarAction(const RealD m_sq_, const RealD lmbd_,
                const RealD alpha_) except +
        RealD action_node(const Field[RealD]& sf) except +
        RealD hmc_m_hamilton_node(const Field[ComplexD]& sm_complex,
                const Field[RealD]& masses) except +
        void hmc_set_force(Field[RealD]& sm_force,
                const Field[RealD]& sf) except +
        void hmc_field_evolve(Field[ComplexD]& sf_complex,
                const Field[ComplexD]& sm_complex,
                const Field[RealD]& masses, const RealD step_size) except +
        RealD sum_sq(const Field[RealD]& f) except +
        void axial_current_node(Field[RealD]& axial_current,
                const Field[RealD]& sf) except +
        void hmc_set_rand_momentum(Field[ComplexD]& sm_complex,
                const Field[RealD]& masses, const RngState& rs) except +
        void hmc_predict_field(Field[ComplexD]& field_ft,
                const Field[ComplexD]& momentum_ft,
                const Field[RealD]& masses, const RealD vev_sigma) except +
        void get_polar_field(Field[RealD]& polar_fields,
                const Field[RealD]& sf) except +
        void hmc_estimate_mass(Field[RealD]& masses,
                const Field[ComplexD]& field_ft,
                const Field[ComplexD]& force_ft, const RealD phi0) except +
        void to_mass_factor(Field[RealD]& sin_domega) except +

cdef extern from "qlat/qm-action.h" namespace "qlat":

    cdef cppclass QMAction:
        bool initialized
        RealD alpha
        RealD beta
        RealD FV_offset
        RealD TV_offset
        RealD center_bar
        RealD V_FV_min
        RealD barrier_strength
        RealD L
        RealD M
        RealD epsilon
        Long t_FV_out
        Long t_FV_mid
        RealD dt
        bool measure_offset_L
        bool measure_offset_M
        QMAction() except +
        QMAction(const RealD alpha_, const RealD beta_,
                const RealD V_FV_min_, const RealD FV_offset_,
                const RealD TV_offset_, const RealD barrier_strength_,
                const RealD L_, const RealD M_, const RealD epsilon_,
                const Long t_FV_out_, const Long t_FV_mid_, const RealD dt_,
                const bool measure_offset_L_,
                const bool measure_offset_M_) except +
        RealD V(const Vector[RealD]& x, const Long t) except +
        RealD dV(const Vector[RealD]& x, const Long t, const int idx) except +
        RealD action_node(Field[RealD]& f) except +
        RealD hmc_m_hamilton_node(const Field[RealD]& mf) except +
        RealD sum_sq(const Field[RealD]& f) except +
        void hmc_set_force(Field[RealD]& force, Field[RealD]& f) except +
        void hmc_field_evolve(Field[RealD]& f, const Field[RealD]& mf,
                const RealD step_size) except +
        void hmc_set_rand_momentum(Field[RealD]& mf,
                const RngState& rs) except +

cdef extern from "qlat/fermion-action.h" namespace "qlat":

    cdef cppclass FermionAction:
        bool initialized
        RealD mass
        Int ls
        RealD m5
        RealD mobius_scale
        bool is_multiplying_dminus
        bool is_using_zmobius
        Int cg_diagonal_mee
        std_vector[ComplexD] bs
        std_vector[ComplexD] cs
        FermionAction() except +
        FermionAction(const RealD mass_, const Int ls_, const RealD m5_,
                const RealD mobius_scale_,
                const bool is_multiplying_dminus_,
                const bool is_using_zmobius_) except +

cdef extern from *:
    """
    #include <qlat/fermion-action.h>
    inline std::vector<double> py_get_omega_fermion_action_re(
        const qlat::FermionAction& fa)
    {
      std::vector<double> re(fa.bs.size());
      for (std::size_t i = 0; i < re.size(); ++i) {
        re[i] = (1.0 / (fa.bs[i] + fa.cs[i])).real();
      }
      return re;
    }
    inline std::vector<double> py_get_omega_fermion_action_im(
        const qlat::FermionAction& fa)
    {
      std::vector<double> im(fa.bs.size());
      for (std::size_t i = 0; i < im.size(); ++i) {
        im[i] = (1.0 / (fa.bs[i] + fa.cs[i])).imag();
      }
      return im;
    }
    """
    std_vector[RealD] py_get_omega_fermion_action_re(
            const FermionAction& fa) except +
    std_vector[RealD] py_get_omega_fermion_action_im(
            const FermionAction& fa) except +

cdef extern from "qlat/dslash.h" namespace "qlat":

    cdef cppclass InverterDomainWall:
        InverterDomainWall() except +
        void setup() except +
        void setup(const GaugeField& gf_, const FermionAction& fa_) except +
        RealD& stop_rsd() except +
        Long& max_num_iter() except +
        Long& max_mixed_precision_cycle() except +
    void setup_inverter(InverterDomainWall& inv, const GaugeField& gf,
            const FermionAction& fa) except +
    void invert_prop_dw "qlat::invert<qlat::InverterDomainWall, qlat::RealD>"(
            Prop& sol, const Prop& src,
            const InverterDomainWall& inv) except +

cdef extern from *:
    """
    #include <qlat/dslash.h>
    inline void py_set_stop_rsd_inverter_domain_wall(
        qlat::InverterDomainWall& inv, const qlat::RealD stop_rsd)
    { inv.stop_rsd() = stop_rsd; }
    inline void py_set_max_num_iter_inverter_domain_wall(
        qlat::InverterDomainWall& inv, const qlat::Long max_num_iter)
    { inv.max_num_iter() = max_num_iter; }
    inline void py_set_max_mixed_precision_cycle_inverter_domain_wall(
        qlat::InverterDomainWall& inv,
        const qlat::Long max_mixed_precision_cycle)
    { inv.max_mixed_precision_cycle() = max_mixed_precision_cycle; }
    """
    void py_set_stop_rsd_inverter_domain_wall(InverterDomainWall& inv,
            const RealD stop_rsd) except +
    void py_set_max_num_iter_inverter_domain_wall(InverterDomainWall& inv,
            const Long max_num_iter) except +
    void py_set_max_mixed_precision_cycle_inverter_domain_wall(
            InverterDomainWall& inv,
            const Long max_mixed_precision_cycle) except +

cdef extern from "qlat/hmc-stats.h" namespace "qlat":

    void display_gm_force_magnitudes(const GaugeMomentum& gm_force,
            const Int n_elems) except +
    std_vector[RealD] get_gm_force_magnitudes(const GaugeMomentum& gm_force,
            const Int n_elems) except +
    void save_gm_force_magnitudes_list(const std_string& fn) except +
    void display_gauge_field_info_table_with_wilson_flow(
            const std_string& fn_gf_info,
            const std_string& fn_wilson_flow_energy,
            const GaugeField& gf, const RealD flow_time,
            const Int flow_steps, const Int steps, const RealD c1) except +

cdef extern from "qlat/contract-pion.h" namespace "qlat":

    LatData contract_pion(const Prop& prop, const Int tslice_src) except +
    LatData contract_pion(const SelProp& prop, const Int tslice_src,
            const FieldSelection& fsel) except +

cdef extern from "qlat/contract-hvp.h" namespace "qlat":

    LatData contract_chvp3(const SelProp& prop1, const SelProp& prop2,
            const Int tslice_src, const FieldSelection& fsel) except +

cdef extern from *:
    """
    #include <qlat/contract-field.h>
    inline void py_contract_chvp_16(qlat::Field<qlat::ComplexD>& chvp,
        const qlat::Prop& prop1, const qlat::Prop& prop2)
    {
      qlat::contract_chvp_16(*(qlat::FieldM<qlat::ComplexD, 16>*)&chvp,
          prop1, prop2);
    }
    """
    void py_contract_chvp_16(Field[ComplexD]& chvp, const Prop& prop1,
            const Prop& prop2) except +

cdef extern from "qlat/field-double.h" namespace "qlat":

    void set_complex_from_double[M](Field[M]& cf,
            const Field[RealD]& sf) except +
    void set_double_from_complex[M](Field[M]& sf,
            const Field[ComplexD]& cf) except +
    void set_abs_from_complex[M](Field[M]& sf,
            const Field[ComplexD]& cf) except +
    void set_ratio_double[M](Field[M]& sf, const Field[RealD]& sf1,
            const Field[RealD]& sf2) except +
    void less_than_double[M](Field[M]& sf1, const Field[RealD]& sf2,
            const Field[RealD]& mask) except +
    void invert_double[M](Field[M]& sf) except +
    void multiply_double[M](Field[M]& sf, const Field[RealD]& factor) except +

cdef extern from "qlat-utils/qendian.h" namespace "qlat":

    void to_from_big_endian[M](Vector[M] v, const bool is_par) except +
    void to_from_little_endian[M](Vector[M] v, const bool is_par) except +

cdef extern from "qlat-utils/handle.h" namespace "qlat":

    cdef cppclass ConstHandle[T]:
        ConstHandle()
        ConstHandle(const T& obj)
        void init()
        void init(const T& obj)
        const T& val() except +

cdef extern from "qlat/field.h" namespace "qlat":

    void split_fields[M](std_vector[Handle[Field[M]]]& vec,
            const Field[M]& f) except +
    void merge_fields[M](Field[M]& f,
            const std_vector[ConstHandle[Field[M]]]& vec) except +
    void merge_fields_ms[M](Field[M]& f,
            const std_vector[ConstHandle[Field[M]]]& vec,
            const std_vector[Int]& m_vec) except +

cdef extern from *:
    """
    #include <qlat/field-fft.h>
    template <class M>
    inline void py_fft_complex_fields(
        std::vector<qlat::Handle<qlat::Field<M>>>& vec,
        const std::vector<int>& fft_dirs,
        const std::vector<int>& fft_is_forwards, const qlat::Int mode_fft)
    {
      std::vector<bool> forwards(fft_is_forwards.size());
      for (std::size_t i = 0; i < fft_is_forwards.size(); ++i) {
        forwards[i] = (fft_is_forwards[i] != 0);
      }
      qlat::fft_complex_fields(vec, fft_dirs, forwards, mode_fft);
    }
    template <class M>
    inline void py_field_iadd(qlat::Field<M>& f0, const qlat::Field<M>& f1)
    { f0 += f1; }
    template <class M>
    inline void py_field_isub(qlat::Field<M>& f0, const qlat::Field<M>& f1)
    { f0 -= f1; }
    template <class M>
    inline void py_field_imul_double(qlat::Field<M>& f, const qlat::RealD factor)
    { f *= factor; }
    template <class M>
    inline void py_field_imul_complex(qlat::Field<M>& f, const qlat::ComplexD& factor)
    { f *= factor; }
    template <class M>
    inline void py_field_imul_real_field(qlat::Field<M>& f, const qlat::Field<qlat::RealD>& f1)
    { f *= f1; }
    template <class M>
    inline void py_field_imul_complex_field(qlat::Field<M>& f, const qlat::Field<qlat::ComplexD>& f1)
    { f *= f1; }
    template <class M>
    inline void py_sfield_isub(qlat::SelectedField<M>& f0, const qlat::SelectedField<M>& f1)
    { f0 -= f1; }
    template <class M>
    inline void py_sfield_iadd(qlat::SelectedField<M>& f0, const qlat::SelectedField<M>& f1)
    { f0 += f1; }
    template <class M>
    inline void py_sfield_imul_double(qlat::SelectedField<M>& f, const qlat::RealD factor)
    { f *= factor; }
    template <class M>
    inline void py_set_phase_field(qlat::Field<M>& f, const qlat::CoordinateD& lmom)
    { qlat::set_phase_field(*(qlat::FieldM<M, 1>*)&f, lmom); }
    """
    void py_field_iadd[M](Field[M]& f0, const Field[M]& f1) except +
    void py_field_isub[M](Field[M]& f0, const Field[M]& f1) except +
    void py_field_imul_double[M](Field[M]& f, const RealD factor) except +
    void py_field_imul_complex[M](Field[M]& f, const ComplexD& factor) except +
    void py_field_imul_real_field[M](Field[M]& f, const Field[RealD]& f1) except +
    void py_field_imul_complex_field[M](Field[M]& f, const Field[ComplexD]& f1) except +
    void py_sfield_isub[M](SelectedField[M]& f0, const SelectedField[M]& f1) except +
    void py_sfield_iadd[M](SelectedField[M]& f0, const SelectedField[M]& f1) except +
    void py_sfield_imul_double[M](SelectedField[M]& f,
            const RealD factor) except +
    void py_set_phase_field[M](Field[M]& f, const CoordinateD& lmom) except +
    void py_fft_complex_fields[M](std_vector[Handle[Field[M]]]& vec,
            std_vector[int]& fft_dirs,
            std_vector[int]& fft_is_forwards,
            const Int mode_fft) except +

