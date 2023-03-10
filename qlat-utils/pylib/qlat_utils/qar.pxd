from libcpp cimport bool
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector

cdef extern from "qlat-utils/qar-cache.h" namespace "qlat":

    bool does_regular_file_exist_qar(const std_string& path)

    bool does_file_exist_qar(const std_string& path)

    long& get_qar_multi_vol_max_size()

    void qar_build_index(const std_string& path_qar)

    int qar_create(const std_string& path_qar, const std_string& path_folder,
                   const bool is_remove_folder_after)

    int qar_extract(const std_string& path_qar, const std_string& path_folder,
                    const bool is_remove_qar_after)

    int qcopy_file(const std_string& path_src, const std_string& path_dst)

    std_vector[std_string] list_qar(const std_string& path)

    std_string qcat(const std_string& path)
