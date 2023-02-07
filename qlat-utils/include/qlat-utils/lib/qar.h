#pragma once

namespace qlat
{  //

bool cc_does_regular_file_exist_qar(const std::string& path);

bool cc_does_file_exist_qar(const std::string& path);

long& cc_get_qar_multi_vol_max_size();

int cc_qar_create(const std::string& path_qar, const std::string& path_folder,
                  const bool is_remove_folder_after = false);

int cc_qar_extract(const std::string& path_qar, const std::string& path_folder,
                   const bool is_remove_qar_after = false);

int cc_qcopy_file(const std::string& path_src, const std::string& path_dst);

std::vector<std::string> cc_list_qar(const std::string& path);

std::string cc_qcat(const std::string& path);

}  // namespace qlat
