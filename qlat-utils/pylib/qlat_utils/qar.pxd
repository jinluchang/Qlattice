from libcpp cimport bool
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector

cdef extern from "qlat-utils/lib/qar.h" namespace "qlat":

    bool cc_does_regular_file_exist_qar(const std_string& path)

    bool cc_does_file_exist_qar(const std_string& path)

    long& cc_get_qar_multi_vol_max_size()

    int cc_qar_create(const std_string& path_qar, const std_string& path_folder,
                      const bool is_remove_folder_after)

    int cc_qar_extract(const std_string& path_qar, const std_string& path_folder,
                       const bool is_remove_qar_after)

    int cc_qcopy_file(const std_string& path_src, const std_string& path_dst)

    std_vector[std_string] cc_list_qar(const std_string& path)

    std_string cc_qcat(const std_string& path)
