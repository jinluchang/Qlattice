from libcpp cimport bool
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector

from .complex cimport *

ctypedef uint32_t crc32_t

ctypedef std_vector[std_vector[double]] DataTable

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

    bool is_directory(const std_string& fn);

    bool is_regular_file(const std_string& fn);

    bool does_file_exist(const std_string& fn);

    void clear_is_directory_cache();

    bool is_directory_cache(const std_string& dir_);

    bool is_regular_file_cache(const std_string& fn);

    bool does_file_exist_cache(const std_string& fn);

    int qremove(const std_string& path);

    int qremove_all(const std_string& path);

    int qmkdir(const std_string& path);

    int qmkdir_p(const std_string& path);

    int qmkdir_info(const std_string& path);

    int qmkdir_p_info(const std_string& path);

    crc32_t compute_crc32(const std_string& path);

    void check_all_files_crc32_info(const std_string& path);

    DataTable qload_datatable(const std_string& path,
                                 const bool is_par);
