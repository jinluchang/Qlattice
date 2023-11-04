# all cpp level import

from qlat_utils.everything cimport *

cdef extern from "qlat/mpi.h" namespace "qlat":

    int begin(const int id_node, const Coordinate& size_node, const int color) except +
    int end(const bool is_preserving_cache) except +
    const Coordinate& get_size_node() except +
    const Coordinate& get_coor_node() except +
    void sync_node() except +
    void bcast(float& x, const int root) except +
    void bcast(Long& x, const int root) except +
    void bcast(ComplexD& x, const int root) except +
    void bcast(ComplexF& x, const int root) except +
    void bcast(std_string& recv, const int root) except +
    void bcast(Coordinate& x, const int root) except +
    void bcast(LatData& ld, const int root) except +
    void glb_sum(Long& ld) except +
    void glb_sum(float& ld) except +
    void glb_sum(double& ld) except +
    void glb_sum(ComplexD& ld) except +
    void glb_sum(ComplexF& ld) except +
    void glb_sum(LatData& ld) except +

cdef extern from "qlat/geometry.h" namespace "qlat":

    cdef cppclass GeometryNode:
        bool initialized
        int num_node
        int id_node
        Coordinate size_node
        Coordinate coor_node
        GeometryNode()
        void init()
        void init(const int id_node, const Coordinate& size_node) except +
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
        Long local_volume()
        Long total_volume()
        bool is_local(const Coordinate& x)
        Long index_from_coordinate(const Coordinate& xl)
        Coordinate coordinate_from_index(const Long index)
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

cdef extern from "qlat/core.h" namespace "qlat":

    cdef cppclass Field[T]:
        vector_acc[T] field
        Field()
        void init()
        void init(const Geometry& geo) except +
        void init(const Geometry& geo, int multiplicity) except +
        void init(const Field[T]& field) except +
        const Geometry& get_geo() except +
        T& get_elem(const Coordinate& x) except +
        T& get_elem(const Coordinate& x, const int m) except +
        T& get_elem(const Long index) except +
        T& get_elem(const Long index, const int m) except +
        Vector[T] get_elems(const Coordinate& x) except +
        Vector[T] get_elems(const Long index) except +
    cdef cppclass GaugeField(Field[ColorMatrix]):
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
        vector_acc[Int64t] ranks
        vector_acc[Long] indices
        FieldSelection()
        void init()
        const Geometry& get_geo()
    cdef cppclass PointsSelection:
        PointsSelection()
        PointsSelection(const Long n_points) except +
        Long size()
        void resize(const Long n_points) except +
        Coordinate* data()
        Coordinate& operator[](Long i) except +
    cdef cppclass SelectedField[T]:
        Long n_elems;
        vector_acc[T] field
        SelectedField()
        void init()
        void init(const Geometry& geo, const Long n_elems, const int multiplicity) except +
        void init(const FieldSelection& fsel, const int multiplicity) except +
        const Geometry& get_geo()
    cdef cppclass SelectedPoints[T]:
        int multiplicity
        Long n_points
        vector_acc[T] points
        SelectedPoints()
        void init()
        void init(const Long n_points, const int multiplicity) except +
        void init(const PointsSelection& psel, const int multiplicity) except +
    Vector[T] get_data[T](const Field[T]& x) except +
    void set_zero[T](Field[T]& x) except +
    void qswap[T](Field[T]& x, Field[T]& y) except +
    Vector[T] get_data[T](const SelectedField[T]& x) except +
    void set_zero[T](SelectedField[T]& x) except +
    void qswap[T](SelectedField[T]& x, SelectedField[T]& y) except +
    Vector[T] get_data[T](const SelectedPoints[T]& x) except +
    void set_zero[T](SelectedPoints[T]& x) except +
    void qswap[T](SelectedPoints[T]& x, SelectedPoints[T]& y) except +
    cdef cppclass SelProp(SelectedField[WilsonMatrix]):
        pass
    cdef cppclass PselProp(SelectedPoints[WilsonMatrix]):
        pass

cdef extern from "qlat/field.h" namespace "qlat":

    double qnorm[M](const Field[M]& f) except +
    void qnorm_field[M](Field[RealD]& f, const Field[M]& f1) except +
    void set_xg_field(Field[Int]& f, const Geometry& geo) except +
    void field_shift[M](Field[M]& f, const Field[M]& f1, const Coordinate& shift) except +
    void reflect_field[M](Field[M]& f) except +
    void set_sqrt_field(Field[RealD]& f, const Field[RealD]& f1) except +

cdef extern from "qlat/field-expand.h" namespace "qlat":

    cdef cppclass CommPlan:
        CommPlan()

cdef extern from "qlat/field-io.h" namespace "qlat":

    bool is_field(const std_string& path) except +
    bool is_d_field(const std_string& path) except +
    bool dist_repartition(const Coordinate& new_size_node, const std_string& path, const std_string& new_path) except +
    crc32_t field_crc32[M](const Field[M]& f) except +
    Long write_field[M](const Field[M]& f, const std_string& path, const Coordinate& new_size_node) except +
    Long read_field[M](Field[M]& f, const std_string& path, const Coordinate& new_size_node) except +

cdef extern from "qlat/selected-points.h" namespace "qlat":

    double qnorm[M](const SelectedPoints[M]& sp) except +
    void qnorm_field[M](SelectedPoints[RealD]& sp, const SelectedPoints[M]& sp1) except +
    PointsSelection mk_random_point_selection(const Coordinate& total_site, const Long num, const RngState& rs) except +
    void save_point_selection_info(const PointsSelection& psel, const std_string& path) except +
    PointsSelection load_point_selection_info(const std_string& path) except +
    LatData lat_data_from_selected_points[M](const SelectedPoints[M]& sp) except +
    void selected_points_from_lat_data[M](SelectedPoints[M]& sp, const LatData& ld) except +
    void save_selected_points[M](const SelectedPoints[M]& sp, const std_string& path) except +
    void load_selected_points[M](SelectedPoints[M]& sp, const std_string& path) except +
    PointsSelection mk_tslice_point_selection(const int t_size, const int t_dir) except +
    void field_glb_sum[M](SelectedPoints[M]& sp, const Field[M]& f) except +
    void field_glb_sum_tslice[M](SelectedPoints[M]& sp, const Field[M]& f, const int t_dir) except +
    void set_sqrt_field(SelectedPoints[RealD]& sp, const SelectedPoints[RealD]& sp1) except +

cdef extern from "qlat/selected-field.h" namespace "qlat":

    double qnorm[M](const SelectedField[M]& sf) except +
    void qnorm_field[M](SelectedField[RealD]& f, const SelectedField[M]& f1) except +
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
    void set_field_selected[t](Field[t]& f, const SelectedField[t]& sf,
                               const FieldSelection& fsel) except +
    void set_field_selected[t](Field[t]& f, const SelectedPoints[t]& sp,
                               const Geometry& geo, const PointsSelection& psel) except +
    void set_selected_points[t](SelectedPoints[t]& sp, const Field[t] f,
                                const PointsSelection& psel, const int m) except +
    void set_field_selected[t](Field[t]& f, const SelectedPoints[t]& sp,
                               const Geometry& geo, const PointsSelection& psel, const int m) except +
    bool is_matching_fsel(const FieldSelection& fsel1, const FieldSelection& fsel2) except +
    void set_sqrt_field(SelectedField[RealD]& f, const SelectedField[RealD]& f1) except +

cdef extern from "qlat/selected-field-io.h" namespace "qlat":

    Long write_field_selection(const FieldSelection& fsel, const std_string& path) except +
    Long read_field_selection(FieldSelection& fsel, const std_string& path) except +
    bool is_selected_field(const std_string& path) except +
    void mk_field_selection(FieldRank& f_rank, const Coordinate& total_site, const Int64t val) except +
    void mk_field_selection(FieldRank& f_rank, const Coordinate& total_site, const Long n_per_tslice, const RngState& rs) except +
    void add_field_selection(FieldRank& f_rank, const PointsSelection& psel, const Long rank_psel) except +
    void update_field_selection(FieldSelection& fsel) except +
    PointsSelection psel_from_fsel(const FieldSelection& fsel) except +
    PointsSelection psel_from_fsel_local(const FieldSelection& fsel) except +
    Long write_selected_field[M](const SelectedField[M]& sf, const std_string& path, const FieldSelection& fsel) except +
    Long read_selected_field[M](SelectedField[M]& sf, const std_string& path, const FieldSelection& fsel) except +

cdef extern from "qlat/qcd-prop.h" namespace "qlat":

    void set_wall_src(Prop& prop, const Geometry& geo_input,
                      const int tslice, const CoordinateD& lmom) except +
    void set_point_src(Prop& prop, const Geometry& geo_input,
                       const Coordinate& xg, const ComplexD& value) except +

cdef extern from "qlat/qcd-smear.h" namespace "qlat":

    void gf_ape_smear(GaugeField& gf, const GaugeField& gf0,
                      const double alpha, const Long steps) except +
    void gf_spatial_ape_smear(GaugeField& gf, const GaugeField& gf0,
                              const double alpha, const Long steps) except +
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
    bool check_file_sync_node(ShuffledFieldsReader& sfr, const std_string& fn, std_vector[Long]& final_offsets) except +
    bool does_file_exist_sync_node(ShuffledFieldsReader& sfr, const std_string& fn) except +
    std_vector[std_string] list_fields(ShuffledFieldsReader& sfr) except +
    int truncate_fields_sync_node(const std_string& path, const std_vector[std_string]& fns_keep, const Coordinate& new_size_node) except +
    std_vector[std_string] properly_truncate_fields_sync_node(const std_string& path, const bool is_check_all, const bool is_only_check, const Coordinate& new_size_node) except +

cdef extern from "qlat/compressed-eigen-io.h" namespace "qlat":

    bool check_compressed_eigen_vectors(const std_string& path) except +
    bool eigen_system_repartition(const Coordinate& new_size_node, const std_string& path, const std_string& new_path) except +


