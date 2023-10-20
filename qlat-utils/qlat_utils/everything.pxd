# all cpp level import

from libcpp cimport bool
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector

cimport libcpp.complex
from libc.stdint cimport uint64_t
from libc.stdint cimport int64_t
from libc.stdint cimport uint32_t
from libc.stdint cimport int32_t
from libc.stdint cimport uint8_t
from libc.stdint cimport int8_t

from cpython.ref cimport PyObject

### --------------------------------------------------------

ctypedef double complex PyComplex
ctypedef float complex PyComplexF

ctypedef libcpp.complex.complex[double] Complex
ctypedef libcpp.complex.complex[float] ComplexF

cdef inline PyComplex py_complex_cast(Complex x):
    return (<PyComplex*>(&x))[0]

cdef inline PyComplexF py_complex_f_cast(ComplexF x):
    return (<PyComplexF*>(&x))[0]

cdef inline Complex complex_cast(PyComplex x):
    return (<Complex*>(&x))[0]

cdef inline ComplexF complex_f_cast(PyComplexF x):
    return (<ComplexF*>(&x))[0]

ctypedef double Double

ctypedef float Float

ctypedef long Long

ctypedef int64_t Int64t

ctypedef int32_t Int32t

ctypedef int8_t Int8t

ctypedef char Char

ctypedef uint32_t crc32_t

ctypedef std_vector[std_vector[double]] DataTable

### --------------------------------------------------------

cdef extern from "qlat-utils/mat-vec.h" namespace "qlat":

    cdef cppclass ColorMatrix:
        ColorMatrix()
        const ColorMatrix& operator=(const ColorMatrix& m)
        Complex* data()
    cdef cppclass SpinMatrix:
        SpinMatrix()
        const SpinMatrix& operator=(const SpinMatrix& m)
        Complex* data()
    cdef cppclass WilsonMatrix:
        WilsonMatrix()
        const WilsonMatrix& operator=(const WilsonMatrix& m)
        Complex* data()
    cdef cppclass NonRelWilsonMatrix:
        NonRelWilsonMatrix()
        const NonRelWilsonMatrix& operator=(const NonRelWilsonMatrix& m)
        Complex* data()
    cdef cppclass IsospinMatrix:
        IsospinMatrix()
        const IsospinMatrix& operator=(const IsospinMatrix& m)
        Complex* data()
    cdef cppclass WilsonVector:
        WilsonVector()
        const WilsonVector& operator=(const WilsonVector& m)
        Complex* data()

cdef extern from "qlat-utils/handle.h" namespace "qlat":

    cdef cppclass Vector[T]:
        Vector()
        Vector(const T* p, const long n)
        T* data()
        long size()
        long data_size()

cdef extern from "qlat-utils/mat.h" namespace "qlat":

    void set_zero(ColorMatrix& x)
    void set_zero(SpinMatrix& x)
    void set_zero(WilsonMatrix& x)
    void set_zero(NonRelWilsonMatrix& x)
    void set_zero(IsospinMatrix& x)
    void set_zero(WilsonVector& x)
    Vector[Complex] get_data(const ColorMatrix& x)
    Vector[Complex] get_data(const SpinMatrix& x)
    Vector[Complex] get_data(const WilsonMatrix& x)
    Vector[Complex] get_data(const NonRelWilsonMatrix& x)
    Vector[Complex] get_data(const IsospinMatrix& x)
    Vector[Complex] get_data(const WilsonVector& x)
    Complex matrix_trace(const SpinMatrix& m)
    Complex matrix_trace(const ColorMatrix& m)
    Complex matrix_trace(const WilsonMatrix& m)
    Complex matrix_trace(const WilsonMatrix& m1, const WilsonMatrix& m2)
    Complex matrix_trace(const WilsonMatrix& m1, const SpinMatrix& m2)
    Complex matrix_trace(const SpinMatrix& m1, const WilsonMatrix& m2)
    Complex matrix_trace(const SpinMatrix& m1, const SpinMatrix& m2)
    Complex matrix_trace(const WilsonMatrix& m1, const ColorMatrix& m2)
    Complex matrix_trace(const ColorMatrix& m1, const WilsonMatrix& m2)
    Complex matrix_trace(const ColorMatrix& m1, const ColorMatrix& m2)
    const SpinMatrix& get_gamma_matrix(const int mu)
    void benchmark_matrix_functions(const long count)
    WilsonMatrix g5_herm(const WilsonMatrix& m)
    SpinMatrix operator*(const Complex& a, const SpinMatrix& m)
    SpinMatrix operator*(const SpinMatrix& m, const Complex& a)
    SpinMatrix operator*(const SpinMatrix& m1, const SpinMatrix& m2)
    ColorMatrix operator*(const Complex& a, const ColorMatrix& m)
    ColorMatrix operator*(const ColorMatrix& m, const Complex& a)
    ColorMatrix operator*(const ColorMatrix& m1, const ColorMatrix& m2)
    WilsonMatrix operator*(const Complex& a, const WilsonMatrix& m)
    WilsonMatrix operator*(const WilsonMatrix& m, const Complex& a)
    WilsonMatrix operator*(const WilsonMatrix& m1, const WilsonMatrix& m2)
    WilsonMatrix operator*(const SpinMatrix& m1, const WilsonMatrix& m2)
    WilsonMatrix operator*(const WilsonMatrix& m1, const SpinMatrix& m2)
    WilsonMatrix operator*(const ColorMatrix& m1, const WilsonMatrix& m2)
    WilsonMatrix operator*(const WilsonMatrix& m1, const ColorMatrix& m2)

cdef extern from "qlat-utils/env.h" namespace "qlat":

    double& get_time_limit()
    long get_time_limit_default()
    long& get_verbose_level()
    long get_verbose_level_default()
    double& get_time_budget()
    double get_time_budget_default()

cdef extern from "qlat-utils/timer.h" namespace "qlat":

    int get_id_node()
    int get_num_node()
    double get_time()
    double& get_start_time()
    double& get_actual_start_time()
    double get_total_time()
    double get_actual_total_time()
    double get_remaining_time()
    cdef cppclass Timer:
        long long flops
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
        void reset(long max_call_times_for_always_show_info)
        @staticmethod
        void fork(long max_call_times_for_always_show_info)
        @staticmethod
        void merge()

cdef extern from "qlat-utils/qutils-io.h" namespace "qlat":

    void flush()
    int qtouch(const std_string& path)
    int qtouch_info(const std_string& path)
    int qtouch(const std_string& path, const std_string& content)
    int qtouch_info(const std_string& path, const std_string& content)
    int qappend(const std_string& path, const std_string& content)
    int qappend_info(const std_string& path, const std_string& content)
    int qrename(const std_string& old_path, const std_string& new_path)
    int qrename_info(const std_string& old_path, const std_string& new_path)
    std_vector[std_string] qls(const std_string& path, const bool is_sort)
    std_vector[std_string] qls_all(const std_string& path,
                                      const bool is_folder_before_files,
                                      const bool is_sort)
    bool is_directory(const std_string& fn)
    bool is_regular_file(const std_string& fn)
    bool does_file_exist(const std_string& fn)
    void clear_is_directory_cache()
    bool is_directory_cache(const std_string& dir_)
    bool is_regular_file_cache(const std_string& fn)
    bool does_file_exist_cache(const std_string& fn)
    int qremove(const std_string& path)
    int qremove_all(const std_string& path)
    int qmkdir(const std_string& path)
    int qmkdir_p(const std_string& path)
    int qmkdir_info(const std_string& path)
    int qmkdir_p_info(const std_string& path)
    crc32_t compute_crc32(const std_string& path)
    void check_all_files_crc32_info(const std_string& path)
    DataTable qload_datatable(const std_string& path,
                                 const bool is_par)

cdef extern from "qlat-utils/rng-state.h" namespace "qlat":

    cdef cppclass RngState:
        RngState()
        RngState(const std_string& seed)
        RngState(const RngState& rs0, const std_string& sindex)
        const RngState& operator=(const RngState& rs)
        RngState split(const std_string& sindex)
        RngState newtype(const unsigned long type)
    uint64_t rand_gen(RngState& rs)
    double u_rand_gen(RngState& rs, const double upper, const double lower)
    double g_rand_gen(RngState& rs, const double center, const double sigma)

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
    Coordinate coordinate_from_index(long index, const Coordinate& size)
    long index_from_coordinate(const Coordinate& x, const Coordinate& size)
    int eo_from_coordinate(const Coordinate& xl)
    Coordinate mod(const Coordinate& x, const Coordinate& size)
    Coordinate smod(const Coordinate& x, const Coordinate& size)
    Coordinate middle_mod(const Coordinate& x, const Coordinate& y, const Coordinate& size)
    Coordinate operator+(const Coordinate& x, const Coordinate& y)
    Coordinate operator-(const Coordinate& x, const Coordinate& y)
    Coordinate operator*(const Coordinate& x, const Coordinate& y)
    Coordinate operator*(const Coordinate& x, const int y)
    Coordinate operator*(const int x, const Coordinate& y)
    bool operator==(const Coordinate& x, const Coordinate& y)
    Coordinate c_rand_gen(RngState& rs, const Coordinate& size)
    long sqr(const Coordinate& xg)
    CoordinateD operator+(const CoordinateD& x, const CoordinateD& y)
    CoordinateD operator-(const CoordinateD& x, const CoordinateD& y)
    CoordinateD operator*(const CoordinateD& x, const CoordinateD& y)
    CoordinateD operator*(const CoordinateD& x, const double y)
    CoordinateD operator*(const double x, const CoordinateD& y)
    bool operator==(const CoordinateD& x, const CoordinateD& y)

cdef extern from "qlat-utils/lat-io.h" namespace "qlat":

    cdef cppclass LatDim:
        std_string name
        long size
        std_vector[std_string] indices
        LatDim()
    ctypedef std_vector[LatDim] LatInfo
    cdef cppclass LatData:
        LatInfo info
        std_vector[double] res
        LatData()
        const LatData& operator=(const LatData& ld)
        void load(const std_string& fn) except +
        void save(const std_string& fn) except +
        bool is_complex()
        int ndim()
        double* data()
    LatDim lat_dim_re_im()
    LatDim lat_dim_number(const std_string& name, const long start, const long end)
    LatDim lat_dim_number(const std_string& name, const long start, const long end, const long inc)
    LatDim lat_dim_string(const std_string& name, const std_vector[std_string]& indices)
    LatData operator*(const Complex& a, const LatData& ld) except +
    LatData operator*(const double a, const LatData& ld) except +
    LatData operator*(const LatData& ld, const Complex& a) except +
    LatData operator*(const LatData& ld, const double a) except +
    LatData operator+(const LatData& ld1, const LatData& ld2) except +
    LatData operator-(const LatData& ld1, const LatData& ld2) except +
    bool is_matching(const LatData& ld1, const LatData& ld2)
    long lat_data_size(LatData& ld)
    void lat_data_alloc(LatData& ld) except +
    Vector[double] get_data(const LatData& x)
    void set_zero(LatData& x)
    std_string show(const LatData& x)
    void clear(LatData& x)
    double qnorm(const LatData& x)

cdef extern from "qlat-utils/qar-cache.h" namespace "qlat":

    bool does_regular_file_exist_qar(const std_string& path) except +
    bool does_file_exist_qar(const std_string& path) except +
    long& get_qar_multi_vol_max_size()
    void qar_build_index(const std_string& path_qar) except +
    int qar_create(const std_string& path_qar, const std_string& path_folder,
                   const bool is_remove_folder_after) except +
    int qar_extract(const std_string& path_qar, const std_string& path_folder,
                    const bool is_remove_qar_after) except +
    int qcopy_file(const std_string& path_src, const std_string& path_dst) except +
    std_vector[std_string] list_qar(const std_string& path) except +
    std_string qcat(const std_string& path) except +

cdef extern from "qlat-utils/cache.h" namespace "qlat":

    void random_permute[T](std_vector[T]& vec, const RngState& rs)
    std_vector[std_string] get_all_caches_info()
    void clear_all_caches()
    void displayln_malloc_stats()

### --------------------------------------------------------

cdef fused ElemType:
    ColorMatrix
    WilsonMatrix
    NonRelWilsonMatrix
    IsospinMatrix
    SpinMatrix
    WilsonVector
    Complex
    ComplexF
    Double
    Float
    Long
    Int64t
    Int8t
    Char

### --------------------------------------------------------
