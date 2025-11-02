#include <qlat-utils/mpi-auto.h>

namespace qlat
{  //

void sync_node(const Long tag)
{
  TIMER("sync_node");
  RngState& rs = get_sync_node_rs();
  const Long rand_limit = 1024L * 1024L * 1024L;
  const Long v = rand_gen(rs) % rand_limit + tag;
  Long s = v;
  glb_sum_val(s);
  qassert(s == v * get_num_node());
}

// ----------------------------

Int glb_sum_long_vec_local(void* ptr, const Long size)
{
  assert(ptr != NULL);
  assert(size >= 0);
  assert(get_num_node() == 1);
  return 0;
}

Int glb_sum_long_vec(void* ptr, const Long size)
{
  return get_glb_sum_long_vec_ptr()(ptr, size);
}

Int glb_sum_vec(Vector<Long> vec)
{
  return glb_sum_long_vec(vec.data(), vec.data_size());
}

Int glb_sum_val(Long& x) { return glb_sum_long_vec(&x, sizeof(Long)); }

// ----------------------------

Int glb_sum_int_vec_local(void* ptr, const Long size)
{
  assert(ptr != NULL);
  assert(size >= 0);
  assert(get_num_node() == 1);
  return 0;
}

Int glb_sum_int_vec(void* ptr, const Long size)
{
  return get_glb_sum_int_vec_ptr()(ptr, size);
}

Int glb_sum_vec(Vector<Int> vec)
{
  return glb_sum_int_vec(vec.data(), vec.data_size());
}

Int glb_sum_val(Int& x) { return glb_sum_int_vec(&x, sizeof(Int)); }

// ----------------------------

Int glb_sum_real_d_vec_local(void* ptr, const Long size)
{
  assert(ptr != NULL);
  assert(size >= 0);
  assert(get_num_node() == 1);
  return 0;
}

Int glb_sum_real_d_vec(void* ptr, const Long size)
{
  return get_glb_sum_real_d_vec_ptr()(ptr, size);
}

Int glb_sum_vec(Vector<RealD> vec)
{
  return glb_sum_real_d_vec(vec.data(), vec.data_size());
}

Int glb_sum_val(RealD& x) { return glb_sum_real_d_vec(&x, sizeof(RealD)); }

// -------------------

Int glb_sum_real_f_vec_local(void* ptr, const Long size)
{
  assert(ptr != NULL);
  assert(size >= 0);
  assert(get_num_node() == 1);
  return 0;
}

Int glb_sum_real_f_vec(void* ptr, const Long size)
{
  return get_glb_sum_real_f_vec_ptr()(ptr, size);
}

Int glb_sum_vec(Vector<RealF> vec)
{
  return glb_sum_real_f_vec(vec.data(), vec.data_size());
}

Int glb_sum_val(RealF& x) { return glb_sum_real_f_vec(&x, sizeof(RealF)); }

// -------------------

Int glb_sum_byte_vec_local(void* ptr, const Long size)
{
  assert(ptr != NULL);
  assert(size >= 0);
  assert(get_num_node() == 1);
  return 0;
}

Int glb_sum_byte_vec(void* ptr, const Long size)
{
  return get_glb_sum_byte_vec_ptr()(ptr, size);
}

Int glb_sum_vec(Vector<Char> vec)
{
  return glb_sum_byte_vec(vec.data(), vec.data_size());
}

Int glb_sum_val(Char& x) { return glb_sum_byte_vec(&x, sizeof(Char)); }

// ----------------------------

Int bcast_byte_vec_local(void* ptr, const Long size, const Int root)
{
  assert(ptr != NULL);
  assert(size >= 0);
  assert(root == 0);
  assert(get_num_node() == 1);
  return 0;
}

Int bcast_byte_vec(void* ptr, const Long size, const Int root)
{
  return get_bcast_byte_vec_ptr()(ptr, size, root);
}

Int bcast_vec(Vector<Char> vec, const Int root)
{
  return bcast_byte_vec(vec.data(), vec.data_size(), root);
}

Int bcast_vec(Vector<Long> vec, const Int root)
{
  return bcast_byte_vec(vec.data(), vec.data_size(), root);
}

Int bcast_vec(Vector<Int> vec, const Int root)
{
  return bcast_byte_vec(vec.data(), vec.data_size(), root);
}

Int bcast_vec(Vector<RealD> vec, const Int root)
{
  return bcast_byte_vec(vec.data(), vec.data_size(), root);
}

// ----------------------------

Int bcast_val(Long& data, const Int root)
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
