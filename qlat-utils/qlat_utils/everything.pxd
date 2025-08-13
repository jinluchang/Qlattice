# all cpp level import

from libcpp cimport bool
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libc.string cimport memcpy

cimport libcpp.complex
from libc.stdint cimport uint64_t
from libc.stdint cimport int64_t
from libc.stdint cimport uint32_t
from libc.stdint cimport int32_t
from libc.stdint cimport uint8_t
from libc.stdint cimport int8_t

from cpython.ref cimport PyObject

### --------------------------------------------------------

cdef extern from "qlat-utils/complex.h" namespace "qlat":

    cdef cppclass ComplexD:
        ComplexD()
    cdef cppclass ComplexF:
        ComplexF()

ctypedef double complex PyComplexD
ctypedef float complex PyComplexF

# Python Complex Cast utilities

cdef inline PyComplexD pycc_d(ComplexD x):
    return (<PyComplexD*>(&x))[0]

cdef inline PyComplexF pycc_f(ComplexF x):
    return (<PyComplexF*>(&x))[0]

cdef inline ComplexD ccpy_d(PyComplexD x):
    return (<ComplexD*>(&x))[0]

cdef inline ComplexF ccpy_f(PyComplexF x):
    return (<ComplexF*>(&x))[0]

### --------------------------------------------------------

ctypedef bool Bool

ctypedef uint64_t UInt64t

ctypedef uint32_t UInt32t

ctypedef uint8_t UInt8t

ctypedef int64_t Int64t

ctypedef int32_t Int32t

ctypedef int8_t Int8t

ctypedef uint32_t crc32_t

ctypedef double RealD

ctypedef float RealF

ctypedef RealD Real

ctypedef Int64t Long

ctypedef Int32t Int

ctypedef Int8t Char

ctypedef UInt64t ULong

ctypedef UInt32t UInt

ctypedef UInt8t UChar

ctypedef std_vector[std_vector[RealD]] DataTable

### --------------------------------------------------------

ctypedef void (*DisplayPtr)(const std_string& str)

### --------------------------------------------------------

cdef extern from "qlat-utils/mat-vec.h" namespace "qlat":

    cdef cppclass ColorMatrix:
        ColorMatrix()
        ComplexD* data()
    cdef cppclass SpinMatrix:
        SpinMatrix()
        ComplexD* data()
    cdef cppclass WilsonMatrix:
        WilsonMatrix()
        ComplexD* data()
    cdef cppclass NonRelWilsonMatrix:
        NonRelWilsonMatrix()
        ComplexD* data()
    cdef cppclass IsospinMatrix:
        IsospinMatrix()
        ComplexD* data()
    cdef cppclass WilsonVector:
        WilsonVector()
        ComplexD* data()

cdef extern from "qlat-utils/handle.h" namespace "qlat":

    cdef cppclass Handle[T]:
        T* p
        Handle()
        Handle(const T& obj)
        void init()
        void init(const T& obj)
        T& val() except +
    cdef cppclass Vector[T]:
        T* p
        Long n
        Vector()
        Vector(const T* p, const Long n)
        T& operator[](const Long i) except +
        T* data()
        Long size()
        Long data_size()

cdef extern from "qlat-utils/vector.h" namespace "qlat":

    cdef cppclass vector[T]:
        bool is_copy
        bool is_acc
        Vector[T] v
        vector()
        void init()
        Long size()
        T* data()
        T& operator[](const Long i) except +
    cdef cppclass box[T]:
        bool is_copy
        bool is_acc
        Handle[T] v
        box()
        void init()
        bool null() except +

cdef extern from "qlat-utils/utils-vec.h" namespace "qlat":

    void assign_direct[M, N](M& x, const N& y) except +
    void iadd_direct[M, N](M& x, const N& y) except +
    void isub_direct[M, N](M& x, const N& y) except +
    void imul_direct[M, N](M& x, const N& y) except +

cdef extern from "qlat-utils/mat.h" namespace "qlat":

    void set_zero(RealD& x)
    void set_zero(RealF& x)
    void set_zero(ComplexD& x)
    void set_zero(ComplexF& x)
    void set_zero(ColorMatrix& x)
    void set_zero(SpinMatrix& x)
    void set_zero(WilsonMatrix& x)
    void set_zero(NonRelWilsonMatrix& x)
    void set_zero(IsospinMatrix& x)
    void set_zero(WilsonVector& x)
    Vector[ComplexD] get_data(const ColorMatrix& x)
    Vector[ComplexD] get_data(const SpinMatrix& x)
    Vector[ComplexD] get_data(const WilsonMatrix& x)
    Vector[ComplexD] get_data(const NonRelWilsonMatrix& x)
    Vector[ComplexD] get_data(const IsospinMatrix& x)
    Vector[ComplexD] get_data(const WilsonVector& x)
    ComplexD matrix_trace(const SpinMatrix& m)
    ComplexD matrix_trace(const ColorMatrix& m)
    ComplexD matrix_trace(const WilsonMatrix& m)
    ComplexD matrix_trace(const WilsonMatrix& m1, const WilsonMatrix& m2)
    ComplexD matrix_trace(const WilsonMatrix& m1, const SpinMatrix& m2)
    ComplexD matrix_trace(const SpinMatrix& m1, const WilsonMatrix& m2)
    ComplexD matrix_trace(const SpinMatrix& m1, const SpinMatrix& m2)
    ComplexD matrix_trace(const WilsonMatrix& m1, const ColorMatrix& m2)
    ComplexD matrix_trace(const ColorMatrix& m1, const WilsonMatrix& m2)
    ComplexD matrix_trace(const ColorMatrix& m1, const ColorMatrix& m2)
    const SpinMatrix& get_gamma_matrix(const int mu)
    void benchmark_matrix_functions(const Long count)
    WilsonMatrix g5_herm(const WilsonMatrix& m)
    SpinMatrix operator*(const ComplexD& a, const SpinMatrix& m)
    SpinMatrix operator*(const SpinMatrix& m, const ComplexD& a)
    SpinMatrix operator*(const SpinMatrix& m1, const SpinMatrix& m2)
    ColorMatrix operator*(const ComplexD& a, const ColorMatrix& m)
    ColorMatrix operator*(const ColorMatrix& m, const ComplexD& a)
    ColorMatrix operator*(const ColorMatrix& m1, const ColorMatrix& m2)
    WilsonMatrix operator*(const ComplexD& a, const WilsonMatrix& m)
    WilsonMatrix operator*(const WilsonMatrix& m, const ComplexD& a)
    WilsonMatrix operator*(const WilsonMatrix& m1, const WilsonMatrix& m2)
    WilsonMatrix operator*(const SpinMatrix& m1, const WilsonMatrix& m2)
    WilsonMatrix operator*(const WilsonMatrix& m1, const SpinMatrix& m2)
    WilsonMatrix operator*(const ColorMatrix& m1, const WilsonMatrix& m2)
    WilsonMatrix operator*(const WilsonMatrix& m1, const ColorMatrix& m2)
    WilsonMatrix operator+(const WilsonMatrix& m1, const WilsonMatrix& m2)
    SpinMatrix operator+(const SpinMatrix& m1, const SpinMatrix& m2)
    ColorMatrix operator+(const ColorMatrix& m1, const ColorMatrix& m2)
    WilsonMatrix operator-(const WilsonMatrix& m1, const WilsonMatrix& m2)
    SpinMatrix operator-(const SpinMatrix& m1, const SpinMatrix& m2)
    ColorMatrix operator-(const ColorMatrix& m1, const ColorMatrix& m2)
    WilsonMatrix matrix_conjugate(const WilsonMatrix& m)
    WilsonMatrix matrix_transpose(const WilsonMatrix& m)
    WilsonMatrix matrix_adjoint(const WilsonMatrix& m)
    SpinMatrix matrix_conjugate(const SpinMatrix& m)
    SpinMatrix matrix_transpose(const SpinMatrix& m)
    SpinMatrix matrix_adjoint(const SpinMatrix& m)
    ColorMatrix matrix_conjugate(const ColorMatrix& m)
    ColorMatrix matrix_transpose(const ColorMatrix& m)
    ColorMatrix matrix_adjoint(const ColorMatrix& m)
    RealD qnorm(const WilsonMatrix& x)
    RealD qnorm(const SpinMatrix& x)
    RealD qnorm(const ColorMatrix& x)
    ComplexD epsilon_contraction(const int v_s1, const int b_s1, const int v_s2,
                                 const int b_s2, const int v_s3, const int b_s3,
                                 const WilsonMatrix& wm1, const WilsonMatrix& wm2,
                                 const WilsonMatrix& wm3)

cdef extern from "qlat-utils/env.h" namespace "qlat":

    double& get_time_limit()
    Long get_time_limit_default()
    Long& get_verbose_level()
    Long get_verbose_level_default()
    double& get_time_budget()
    double get_time_budget_default()
    Long& get_qar_multi_vol_max_size()
    Long get_qar_multi_vol_max_size_default()

cdef extern from "qlat-utils/show.h" namespace "qlat":

    void display_c_stdout(const std_string& str) except +
    DisplayPtr& get_display_ptr() except +
    void set_display_ptr() except +
    void set_display_ptr(DisplayPtr f) except +
    void display(const std_string& str) except +

cdef extern from "qlat-utils/eigen.h" namespace "qlat":

    std_string get_eigen_type() except +

cdef extern from "qlat-utils/timer.h" namespace "qlat":

    int get_id_node()
    int get_num_node()
    void sync_node() except +
    void sync_node(const Long tag) except +
    double get_time()
    double& get_start_time()
    double& get_actual_start_time()
    double get_total_time()
    double get_actual_total_time()
    double get_remaining_time()
    cdef cppclass Timer:
        Long flops
        Timer()
        Timer(const std_string& fname)
        void start()
        void start(bool is_verbose)
        void stop()
        void stop(bool is_verbose)
        @staticmethod
        void display(const std_string& tag)
        @staticmethod
        void autodisplay()
        @staticmethod
        void display_stack()
        @staticmethod
        void display_stack_always()
        @staticmethod
        void reset(Long max_call_times_for_always_show_info)
        @staticmethod
        void fork(Long max_call_times_for_always_show_info)
        @staticmethod
        void merge()

cdef extern from "qlat-utils/utils-io.h" namespace "qlat":

    std_string basename(const std_string& fn) except +
    std_string dirname(const std_string& fn) except +
    std_vector[std_string] all_dirname_vec(const std_string& fn) except +
    std_string remove_trailing_slashes(const std_string& fn) except +
    std_vector[std_string] qls(const std_string& path, const bool is_sort) except +
    std_vector[std_string] qls_all(const std_string& path, const bool is_folder_before_files, const bool is_sort) except +
    bool does_file_exist(const std_string& fn) except +
    bool is_directory(const std_string& fn) except +
    bool is_regular_file(const std_string& fn) except +
    int qmkdir(const std_string& path) except +
    int qmkdir_p(const std_string& path) except +
    int qrename(const std_string& old_path, const std_string& new_path) except +
    int qremove(const std_string& path) except +
    int qremove_all(const std_string& path) except +
    int qtruncate(const std_string& path, const Long offset) except +
    #
    void clear_is_directory_cache() except +
    void remove_entry_directory_cache(const std_string& dir_) except +
    void add_entry_directory_cache(const std_string& dir_, bool is_directory) except +
    bool is_directory_cache(const std_string& dir_) except +
    bool is_regular_file_cache(const std_string& fn) except +
    bool does_file_exist_cache(const std_string& fn) except +
    #
    int qmkdir_info(const std_string& path) except +
    int qmkdir_p_info(const std_string& path) except +
    int qrename_info(const std_string& old_path, const std_string& new_path) except +
    int qremove_info(const std_string& path) except +
    int qremove_all_info(const std_string& path) except +
#
    std_vector[std_string] qls_sync_node(const std_string& path, const bool is_sort) except +
    std_vector[std_string] qls_all_sync_node(const std_string& path, const bool is_folder_before_files, const bool is_sort) except +
    bool does_file_exist_sync_node(const std_string& fn) except +
    bool is_directory_sync_node(const std_string& fn) except +
    bool is_regular_file_sync_node(const std_string& fn) except +
    int qmkdir_sync_node(const std_string& path) except +
    int qmkdir_p_sync_node(const std_string& path) except +
    int qremove_sync_node(const std_string& path) except +
    int qremove_all_sync_node(const std_string& path) except +

cdef extern from "qlat-utils/rng-state.h" namespace "qlat":

    cdef cppclass RngState:
        RngState()
        RngState(const std_string& seed)
        RngState(const RngState& rs0, const std_string& sindex)
        RngState split(const std_string& sindex)
        RngState newtype(const unsigned long type)
    uint64_t rand_gen(RngState& rs)
    double u_rand_gen(RngState& rs, const double upper, const double lower)
    double g_rand_gen(RngState& rs, const double center, const double sigma)
    void random_permute[T](std_vector[T]& vec, const RngState& rs)

cdef extern from "qlat-utils/coordinate-d.h" namespace "qlat":

    cdef cppclass Coordinate:
        Coordinate()
        Coordinate(int x, int y, int z, int t)
        int& operator[](unsigned long i)
    cdef cppclass CoordinateD:
        CoordinateD()
        CoordinateD(const Coordinate& x)
        CoordinateD(double x, double y, double z, double t)
        double& operator[](unsigned long i)
    Coordinate coordinate_from_index(Long index, const Coordinate& size)
    Long index_from_coordinate(const Coordinate& x, const Coordinate& size)
    int eo_from_coordinate(const Coordinate& xl)
    Coordinate mod(const Coordinate& x, const Coordinate& size)
    Coordinate smod(const Coordinate& x, const Coordinate& size)
    Coordinate smod_sym(const Coordinate& x, const Coordinate& size)
    Coordinate middle_mod(const Coordinate& x, const Coordinate& y, const Coordinate& size)
    Coordinate operator+(const Coordinate& x, const Coordinate& y)
    Coordinate operator-(const Coordinate& x, const Coordinate& y)
    Coordinate operator*(const Coordinate& x, const Coordinate& y)
    Coordinate operator*(const Coordinate& x, const int y)
    Coordinate operator*(const int x, const Coordinate& y)
    Coordinate operator/(const Coordinate& x, const Coordinate& y)
    Coordinate operator/(const Coordinate& x, const int y)
    bool operator==(const Coordinate& x, const Coordinate& y)
    Coordinate c_rand_gen(RngState& rs, const Coordinate& size)
    Long sqr(const Coordinate& xg)
    Long product(const Coordinate& coor)
    CoordinateD mod(const CoordinateD& x, const CoordinateD& size)
    CoordinateD smod(const CoordinateD& x, const CoordinateD& size)
    CoordinateD smod_sym(const CoordinateD& x, const CoordinateD& size)
    CoordinateD middle_mod(const CoordinateD& x, const CoordinateD& y, const CoordinateD& size)
    CoordinateD operator+(const CoordinateD& x, const CoordinateD& y)
    CoordinateD operator-(const CoordinateD& x, const CoordinateD& y)
    CoordinateD operator*(const CoordinateD& x, const CoordinateD& y)
    CoordinateD operator*(const CoordinateD& x, const double y)
    CoordinateD operator*(const double x, const CoordinateD& y)
    CoordinateD operator/(const CoordinateD& x, const CoordinateD& y)
    CoordinateD operator/(const CoordinateD& x, const double y)
    CoordinateD operator/(const double x, const CoordinateD& y)
    bool operator==(const CoordinateD& x, const CoordinateD& y)
    RealD sqr(const CoordinateD& xg)

cdef extern from "qlat-utils/lat-io.h" namespace "qlat":

    cdef cppclass LatDim:
        std_string name
        Long size
        std_vector[std_string] indices
        LatDim()
    LatDim lat_dim_re_im()
    LatDim lat_dim_number(const std_string& name, const Long start, const Long end)
    LatDim lat_dim_number(const std_string& name, const Long start, const Long end, const Long inc)
    LatDim lat_dim_string(const std_string& name, const std_vector[std_string]& indices)
    Long lat_dim_idx(const LatDim& dim, const std_string& idx) except +
    ctypedef std_vector[LatDim] LatInfo
    cdef cppclass LatDataT[T]:
        LatInfo info
        std_vector[T] res
        LatDataT()
        void load(const std_string& fn) except +
        void save(const std_string& fn) except +
        void load_str(std_string& content) except +
        std_string save_str() except +
        bool is_complex()
        Int ndim()
        T* data()
    bool is_matching[T](const LatDataT[T]& ld1, const LatDataT[T]& ld2) except +
    Long lat_data_size[T](LatDataT[T]& ld) except +
    void lat_data_alloc[T](LatDataT[T]& ld) except +
    Vector[T] get_data[T](const LatDataT[T]& x) except +
    void set_zero[T](LatDataT[T]& x) except +
    void clear[T](LatDataT[T]& x) except +
    Int glb_sum[T](LatDataT[T]& ld) except +
    Int bcast[T](LatDataT[T]& ld, const Int root) except +
    void lat_data_load_sync_node[T](LatDataT[T]& ld, const std_string& path) except +
    void lat_data_save_info[T](const std_string& path, const LatDataT[T]& ld) except +
    cdef cppclass LatDataInt(LatDataT[Int]):
        LatDataInt()
    cdef cppclass LatDataLong(LatDataT[Long]):
        LatDataLong()
    cdef cppclass LatDataRealF(LatDataT[RealF]):
        LatDataRealF()
        LatDataRealF& operator=(const LatData&)
    cdef cppclass LatData(LatDataT[RealD]):
        LatData()
        LatData& operator=(const LatDataRealF&)
    LatData operator*(const ComplexD& a, const LatData& ld) except +
    LatData operator*(const RealD a, const LatData& ld) except +
    LatData operator*(const LatData& ld, const ComplexD& a) except +
    LatData operator*(const LatData& ld, const RealD a) except +
    LatData operator+(const LatData& ld1, const LatData& ld2) except +
    LatData operator-(const LatData& ld1, const LatData& ld2) except +
    std_string show(const LatData& x)
    RealD qnorm(const LatData& x)

cdef extern from "qlat-utils/qar.h" namespace "qlat":

    cdef cppclass QFileType:
        pass
    cdef cppclass QFileMode:
        pass
    cdef cppclass QFile:
        QFile() except +
        void init() except +
        void init(const QFileType ftype, const std_string& path, const QFileMode mode) except +
        void init(const QFileType ftype, const std_string& path, const QFileMode mode, std_string& content) except +
        void init(const QFile& qfile, const Long q_offset_start, const Long q_offset_end) except +
        void close() except +
        bool null()
        const std_string& path() except +
        QFileMode mode() except +
        const std_string& content() except +
    cdef cppclass QarFile:
        std_string path
        QFileMode mode
        Long qar_index_size_saved;
        QarFile() except +
        QarFile(const std_string& path_qar, const QFileMode mode) except +
        void init() except +
        void init(const std_string& path_qar, const QFileMode mode) except +
        void close() except +
        bool null()
        bool flush() except +
        Long index_size() except +
        void save_index(const Long max_diff) except +
    #
    QFileType read_qfile_type(const std_string& ftype) except +
    std_string show(const QFileType ftype) except +
    QFileMode read_qfile_mode(const std_string& mode) except +
    std_string show(const QFileMode mode) except +
    int clean_up_qfile_map() except +
    std_vector[std_string] show_all_qfile() except +
    QFile qfopen(const std_string& path, const QFileMode mode) except +
    QFile qfopen(const QFileType ftype, const std_string& path, const QFileMode mode) except +
    QFile qfopen(const QFileType ftype, const std_string& path, const QFileMode mode, std_string& content) except +
    void qfclose(QFile& qfile) except +
    bool qfeof(const QFile& qfile) except +
    Long qftell(const QFile& qfile) except +
    int qfflush(const QFile& qfile) except +
    int qfseek_set(QFile& qfile, const Long q_offset) except +
    int qfseek_end(QFile& qfile, const Long q_offset) except +
    int qfseek_cur(QFile& qfile, const Long q_offset) except +
    Long qfile_size(QFile& qfile) except +
    Long qfile_remaining_size(QFile& qfile) except +
    #
    std_string qgetline(QFile& qfile) except +
    std_vector[std_string] qgetlines(QFile& qfile) except +
    Long qwrite_data(const std_string& line, QFile& qfile) except +
    Long qwrite_data(const std_vector[std_string]& lines, QFile& qfile) except +
    Long write_from_qfile(QFile& qfile_out, QFile& qfile_in) except +
    std_string qcat(QFile& qfile) except +
    #
    std_vector[std_string] list(const QarFile& qar) except +
    bool has_regular_file(const QarFile& qar, const std_string& fn) except +
    bool has(const QarFile& qar, const std_string& fn) except +
    QFile read(const QarFile& qar, const std_string& fn) except +
    std_string read_data(const QarFile& qar, const std_string& fn) except +
    std_string read_info(const QarFile& qar, const std_string& fn) except +
    bool verify_index(const QarFile& qar) except +
    Long write_from_qfile(QarFile& qar, const std_string& fn, const std_string& info, QFile& qfile_in) except +
    Long write_from_data(QarFile& qar, const std_string& fn, const std_string& info, const std_string& data) except +
    Long write_from_data(QarFile& qar, const std_string& fn, const std_string& info, const std_vector[std_string]& data) except +
    #
    std_vector[std_string] properly_truncate_qar_file(const std_string& path) except +
    #
    std_string show_qar_index(const QarFile& qar) except +
    int read_qar_index(const QarFile& qar, const std_string& qar_index_content) except +
    #
    bool does_regular_file_exist_qar(const std_string& path) except +
    bool does_file_exist_qar(const std_string& path) except +
    int qar_build_index(const std_string& path_qar) except +
    int qar_create(const std_string& path_qar, const std_string& path_folder,
                   const bool is_remove_folder_after) except +
    int qar_extract(const std_string& path_qar, const std_string& path_folder,
                    const bool is_remove_qar_after) except +
    int qcopy_file(const std_string& path_src, const std_string& path_dst) except +
    std_vector[std_string] list_qar(const std_string& path) except +
    #
    std_string qcat(const std_string& path) except +
    int qtouch(const std_string& path) except +
    int qtouch(const std_string& path, const std_string& content) except +
    int qtouch(const std_string& path, const std_vector[std_string]& content) except +
    int qappend(const std_string& path, const std_string& content) except +
    int qappend(const std_string& path, const std_vector[std_string]& content) except +
    #
    DataTable qload_datatable(const std_string& path, const bool is_par) except +
    #
    crc32_t compute_crc32(QFile& qfile) except +
    crc32_t compute_crc32(const std_string& path) except +
    #
    int qar_build_index_info(const std_string& path_qar) except +
    int qar_create_info(const std_string& path_qar, const std_string& path_folder_, const bool is_remove_folder_after) except +
    int qar_extract_info(const std_string& path_qar, const std_string& path_folder_, const bool is_remove_qar_after) except +
    int qcopy_file_info(const std_string& path_src, const std_string& path_dst) except +
    std_string qcat_info(const std_string& path) except +
    int qtouch_info(const std_string& path) except +
    int qtouch_info(const std_string& path, const std_string& content) except +
    int qtouch_info(const std_string& path, const std_vector[std_string]& content) except +
    int qappend_info(const std_string& path, const std_string& content) except +
    int qappend_info(const std_string& path, const std_vector[std_string]& content) except +
    void check_all_files_crc32_info(const std_string& path) except +
    #
    bool does_regular_file_exist_qar_sync_node(const std_string& fn) except +
    bool does_file_exist_qar_sync_node(const std_string& fn) except +
    int qar_create_sync_node(const std_string& path_qar, const std_string& path_folder_, const bool is_remove_folder_after) except +
    int qar_extract_sync_node(const std_string& path_qar, const std_string& path_folder_, const bool is_remove_qar_after) except +
    int qcopy_file_sync_node(const std_string& path_src, const std_string& path_dst) except +
    std_string qcat_sync_node(const std_string& path) except +
    DataTable qload_datatable_sync_node(const std_string& path, const bool is_par) except +

cdef extern from "qlat-utils/cache.h" namespace "qlat":

    std_vector[std_string] get_all_caches_info() except +
    void displayln_malloc_stats() except +

cdef extern from "qlat-utils/vector.h" namespace "qlat":

    void clear_mem_cache() except +
    void clear_all_caches() except +

### --------------------------------------------------------

cdef fused ElemType:
    ColorMatrix
    WilsonMatrix
    NonRelWilsonMatrix
    IsospinMatrix
    SpinMatrix
    WilsonVector
    ComplexD
    ComplexF
    RealD
    RealF
    Long
    Int
    Int64t
    Int32t
    Int8t
    Char

### --------------------------------------------------------
