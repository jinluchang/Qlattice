from libcpp cimport bool
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector

from .complex cimport *

ctypedef uint32_t crc32_t

ctypedef std_vector[std_vector[double]] DataTable

cdef extern from "qlat-utils/lib/utils-io.h" namespace "qlat":

    void flush()

    int cc_qtouch(const std_string& path)

    int cc_qtouch_info(const std_string& path)

    int cc_qtouch(const std_string& path, const std_string& content)

    int cc_qtouch_info(const std_string& path, const std_string& content)

    int cc_qappend(const std_string& path, const std_string& content)

    int cc_qappend_info(const std_string& path, const std_string& content)

    int cc_qrename(const std_string& old_path, const std_string& new_path)

    int cc_qrename_info(const std_string& old_path, const std_string& new_path)

    std_vector[std_string] cc_qls(const std_string& path, const bool is_sort)

    std_vector[std_string] cc_qls_all(const std_string& path,
                                      const bool is_folder_before_files,
                                      const bool is_sort)

    bool cc_is_directory(const std_string& fn);

    bool cc_is_regular_file(const std_string& fn);

    int cc_qremove(const std_string& path);

    int cc_qremove_all(const std_string& path);

    int cc_qmkdir(const std_string& path);

    int cc_qmkdir_p(const std_string& path);

    int cc_qmkdir_info(const std_string& path);

    int cc_qmkdir_p_info(const std_string& path);

    bool cc_does_file_exist(const std_string& fn);

    crc32_t cc_compute_crc32(const std_string& path);

    void cc_check_all_files_crc32_info(const std_string& path);

    DataTable cc_qload_datatable(const std_string& path,
                                 const bool is_par);
