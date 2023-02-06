#include <qlat-utils/qar-cache.h>
#include <qlat-utils/lib/qar.h>

namespace qlat
{  //

bool cc_does_regular_file_exist_qar(const std::string& path)
{
  return does_regular_file_exist_qar(path);
}

bool cc_does_file_exist_qar(const std::string& path)
{
  return does_file_exist_qar(path);
}

long& cc_get_qar_multi_vol_max_size() { return get_qar_multi_vol_max_size(); }

int cc_qar_create(const std::string& path_qar, const std::string& path_folder,
                  const bool is_remove_folder_after)
{
  return qar_create(path_qar, path_folder, is_remove_folder_after);
}

int cc_qar_extract(const std::string& path_qar, const std::string& path_folder,
                   const bool is_remove_qar_after)
{
  return qar_extract(path_qar, path_folder, is_remove_qar_after);
}

int cc_qcopy_file(const std::string& path_src, const std::string& path_dst)
{
  return qcopy_file(path_src, path_dst);
}

std::vector<std::string> cc_list_qar(const std::string& path)
{
  return list_qar(path);
}

std::string cc_qcat(const std::string& path) { return qcat(path); }

}  // namespace qlat
