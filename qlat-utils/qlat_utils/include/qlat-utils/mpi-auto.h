// vim: set ts=2 sw=2 expandtab:

// Copyright (c) 2014 Luchang Jin
// All rights reserved.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <qlat-utils/assert.h>
#include <qlat-utils/env.h>
#include <qlat-utils/handle.h>
#include <qlat-utils/rng-state.h>
#include <qlat-utils/show.h>
#include <qlat-utils/timer.h>

namespace qlat
{  //

void sync_node();

// -------------------

API inline RngState& get_sync_node_rs_local()
{
  static RngState rs("sync_node_local");
  return rs;
}

using RngStatePtr = RngState*;

API inline RngStatePtr& get_sync_node_rs_ptr()
{
  static RngStatePtr ptr = &get_sync_node_rs_local();
  return ptr;
}

inline RngState& get_sync_node_rs() { return *get_sync_node_rs_ptr(); }

// -------------------

using MpiSumFuncPtr = Int (*)(void*, const Long);  // size is the data size.

inline Int glb_sum_long_vec_local(void* ptr, const Long size)
{
  assert(ptr != NULL);
  assert(size >= 0);
  assert(get_num_node() == 1);
  return 0;
}

API inline MpiSumFuncPtr& get_glb_sum_long_vec_ptr()
{
  static MpiSumFuncPtr ptr = glb_sum_long_vec_local;
  return ptr;
}

inline Int glb_sum_long_vec(void* ptr, const Long size)
{
  return get_glb_sum_long_vec_ptr()(ptr, size);
}

inline Int glb_sum_vec(Vector<Long> vec)
{
  return glb_sum_long_vec(vec.data(), vec.data_size());
}

inline Int glb_sum_val(Long& x) { return glb_sum_long_vec(&x, sizeof(Long)); }

// -------------------

using MpiSumFuncPtr = Int (*)(void*, const Long);  // size is the data size.

inline Int glb_sum_int_vec_local(void* ptr, const Long size)
{
  assert(ptr != NULL);
  assert(size >= 0);
  assert(get_num_node() == 1);
  return 0;
}

API inline MpiSumFuncPtr& get_glb_sum_int_vec_ptr()
{
  static MpiSumFuncPtr ptr = glb_sum_int_vec_local;
  return ptr;
}

inline Int glb_sum_int_vec(void* ptr, const Long size)
{
  return get_glb_sum_int_vec_ptr()(ptr, size);
}

inline Int glb_sum_vec(Vector<Int> vec)
{
  return glb_sum_int_vec(vec.data(), vec.data_size());
}

inline Int glb_sum_val(Int& x) { return glb_sum_int_vec(&x, sizeof(Int)); }

// -------------------

inline Int glb_sum_real_d_vec_local(void* ptr, const Long size)
{
  assert(ptr != NULL);
  assert(size >= 0);
  assert(get_num_node() == 1);
  return 0;
}

API inline MpiSumFuncPtr& get_glb_sum_real_d_vec_ptr()
{
  static MpiSumFuncPtr ptr = glb_sum_real_d_vec_local;
  return ptr;
}

inline Int glb_sum_real_d_vec(void* ptr, const Long size)
{
  return get_glb_sum_real_d_vec_ptr()(ptr, size);
}

inline Int glb_sum_vec(Vector<RealD> vec)
{
  return glb_sum_real_d_vec(vec.data(), vec.data_size());
}

inline Int glb_sum_val(RealD& x)
{
  return glb_sum_real_d_vec(&x, sizeof(RealD));
}

// -------------------

inline Int glb_sum_real_f_vec_local(void* ptr, const Long size)
{
  assert(ptr != NULL);
  assert(size >= 0);
  assert(get_num_node() == 1);
  return 0;
}

API inline MpiSumFuncPtr& get_glb_sum_real_f_vec_ptr()
{
  static MpiSumFuncPtr ptr = glb_sum_real_f_vec_local;
  return ptr;
}

inline Int glb_sum_real_f_vec(void* ptr, const Long size)
{
  return get_glb_sum_real_f_vec_ptr()(ptr, size);
}

inline Int glb_sum_vec(Vector<RealF> vec)
{
  return glb_sum_real_f_vec(vec.data(), vec.data_size());
}

inline Int glb_sum_val(RealF& x)
{
  return glb_sum_real_f_vec(&x, sizeof(RealF));
}

// -------------------

inline Int glb_sum_byte_vec_local(void* ptr, const Long size)
{
  assert(ptr != NULL);
  assert(size >= 0);
  assert(get_num_node() == 1);
  return 0;
}

API inline MpiSumFuncPtr& get_glb_sum_byte_vec_ptr()
{
  static MpiSumFuncPtr ptr = glb_sum_byte_vec_local;
  return ptr;
}

inline Int glb_sum_byte_vec(void* ptr, const Long size)
{
  return get_glb_sum_byte_vec_ptr()(ptr, size);
}

inline Int glb_sum_vec(Vector<Char> vec)
{
  return glb_sum_byte_vec(vec.data(), vec.data_size());
}

inline Int glb_sum_val(Char& x) { return glb_sum_byte_vec(&x, sizeof(Char)); }

// -------------------

using MpiBcastFuncPtr = Int (*)(void*, const Long,
                                const Int);  // size is the data size.

inline Int bcast_byte_vec_local(void* ptr, const Long size, const Int root)
{
  assert(ptr != NULL);
  assert(size >= 0);
  assert(root == 0);
  assert(get_num_node() == 1);
  return 0;
}

API inline MpiBcastFuncPtr& get_bcast_byte_vec_ptr()
{
  static MpiBcastFuncPtr ptr = bcast_byte_vec_local;
  return ptr;
}

inline Int bcast_byte_vec(void* ptr, const Long size, const Int root)
{
  return get_bcast_byte_vec_ptr()(ptr, size, root);
}

inline Int bcast_vec(Vector<Char> vec, const Int root)
{
  return bcast_byte_vec(vec.data(), vec.data_size(), root);
}

inline Int bcast_vec(Vector<Long> vec, const Int root)
{
  return bcast_byte_vec(vec.data(), vec.data_size(), root);
}

inline Int bcast_vec(Vector<Int> vec, const Int root)
{
  return bcast_byte_vec(vec.data(), vec.data_size(), root);
}

inline Int bcast_vec(Vector<RealD> vec, const Int root)
{
  return bcast_byte_vec(vec.data(), vec.data_size(), root);
}

// -------------------

Int bcast_val(Long& data, const Int root);

Int bcast_val(std::string& data, const Int root);

Int bcast_val(std::vector<std::string>& data, const Int root);

Int bcast_val(std::vector<Int>& data, const Int root);

Int bcast_val(std::vector<Long>& data, const Int root);

Int bcast_val(std::vector<RealD>& data, const Int root);

Int bcast_val(std::vector<RealF>& data, const Int root);

Int bcast_val(std::vector<std::vector<RealD>>& data, const Int root);

}  // namespace qlat
