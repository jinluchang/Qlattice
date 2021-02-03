#include "lib.h"

EXPORT(get_id_node, {
  using namespace qlat;
  return py_convert(get_id_node());
});

EXPORT(get_num_node, {
  using namespace qlat;
  return py_convert(get_num_node());
});

EXPORT(get_coor_node, {
  using namespace qlat;
  return py_convert(get_coor_node());
});

EXPORT(get_size_node, {
  using namespace qlat;
  return py_convert(get_size_node());
});

