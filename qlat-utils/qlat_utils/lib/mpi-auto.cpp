#include <qlat-utils/mpi-auto.h>

namespace qlat
{  //

void sync_node()
{
  TIMER("sync_node");
  RngState& rs = get_sync_node_rs();
  const Long v = rand_gen(rs) % (1024 * 1024);
  Long s = v;
  glb_sum_val(s);
  qassert(s == v * get_num_node());
}

Int bcast_val(Long& data, const int root)
{
  return bcast_byte_vec(&data, sizeof(Long), root);
}

Int bcast_val(std::string& data, const Int root)
{
  TIMER("bcast_val(str)");
  if (get_num_node() == 1) {
    return 0;
  }
  Long size = data.size();
  Int ret = bcast_val(size, root);
  if (ret != 0) {
    return ret;
  }
  if (get_id_node() == root) {
    qassert((Long)data.size() == size);
  } else {
    data.resize(size);
  }
  if (size == 0) {
    return 0;
  }
  Vector<Char> vec((Char*)&data[0], data.size());
  return bcast_vec(vec, root);
}

Int bcast_val(std::vector<std::string>& data, const Int root)
{
  TIMER("bcast_val(vec<str>)");
  if (get_num_node() == 1) {
    return 0;
  }
  Long size = data.size();
  Int ret = bcast_val(size, root);
  if (ret != 0) {
    return ret;
  }
  if (get_id_node() == root) {
    qassert((Long)data.size() == size);
  } else {
    data.resize(size);
  }
  for (Long i = 0; i < size; ++i) {
    ret = bcast_val(data[i], root);
    if (ret != 0) {
      return ret;
    }
  }
  return 0;
}

template <class T>
Int bcast_val_tt(std::vector<T>& data, const Int root)
{
  if (get_num_node() == 1) {
    return 0;
  }
  Long size = data.size();
  Int ret = bcast_val(size, root);
  if (ret != 0) {
    return ret;
  }
  if (get_id_node() == root) {
    qassert((Long)data.size() == size);
  } else {
    data.resize(size);
  }
  if (size == 0) {
    return 0;
  }
  Vector<Char> vec((Char*)&data[0], data.size() * sizeof(T));
  return bcast_vec(vec, root);
}

Int bcast_val(std::vector<Int>& data, const Int root)
{
  TIMER("bcast_val(vec<Int>)");
  return bcast_val_tt(data, root);
}

Int bcast_val(std::vector<Long>& data, const Int root)
{
  TIMER("bcast_val(vec<Long>)");
  return bcast_val_tt(data, root);
}

Int bcast_val(std::vector<RealD>& data, const Int root)
{
  TIMER("bcast_val(vec<RealD>)");
  return bcast_val_tt(data, root);
}

Int bcast_val(std::vector<RealF>& data, const Int root)
{
  TIMER("bcast_val(vec<RealF>)");
  return bcast_val_tt(data, root);
}

Int bcast_val(std::vector<std::vector<RealD>>& data, const Int root)
{
  TIMER("bcast_val(vec<vec<RealD>>)");
  if (get_num_node() == 1) {
    return 0;
  }
  Long size = data.size();
  Int ret = bcast_val(size, root);
  if (ret != 0) {
    return ret;
  }
  if (get_id_node() == root) {
    qassert((Long)data.size() == size);
  } else {
    data.resize(size);
  }
  for (Long i = 0; i < size; ++i) {
    ret = bcast_val(data[i], root);
    if (ret != 0) {
      return ret;
    }
  }
  return 0;
}

}  // namespace qlat
