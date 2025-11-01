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

#include <qlat-utils/qassert.h>
#include <qlat-utils/crc32.h>
#include <qlat-utils/env.h>
#include <qlat-utils/handle.h>
#include <qlat-utils/rng-state.h>
#include <qlat-utils/show.h>
#include <qlat-utils/timer.h>

namespace qlat
{  //

#define SYNC_NODE()                                                      \
  {                                                                      \
    static const qlat::Long CRC32_HASH_FOR_FILE_NAME =                   \
        (qlat::Long)qlat::crc32((void*)__FILE__, std::strlen(__FILE__)); \
    static const qlat::Long SYNC_NODE_TAG =                              \
        CRC32_HASH_FOR_FILE_NAME + __LINE__;                             \
    qlat::sync_node(SYNC_NODE_TAG);                                      \
  }

void sync_node(const Long tag = 0);

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

// -------------------

Int glb_sum_long_vec_local(void* ptr, const Long size);

API inline MpiSumFuncPtr& get_glb_sum_long_vec_ptr()
{
  static MpiSumFuncPtr ptr = glb_sum_long_vec_local;
  return ptr;
}

Int glb_sum_long_vec(void* ptr, const Long size);

Int glb_sum_vec(Vector<Long> vec);

Int glb_sum_val(Long& x);

// -------------------

Int glb_sum_int_vec_local(void* ptr, const Long size);

API inline MpiSumFuncPtr& get_glb_sum_int_vec_ptr()
{
  static MpiSumFuncPtr ptr = glb_sum_int_vec_local;
  return ptr;
}

Int glb_sum_int_vec(void* ptr, const Long size);

Int glb_sum_vec(Vector<Int> vec);

Int glb_sum_val(Int& x);

// -------------------

Int glb_sum_real_d_vec_local(void* ptr, const Long size);

API inline MpiSumFuncPtr& get_glb_sum_real_d_vec_ptr()
{
  static MpiSumFuncPtr ptr = glb_sum_real_d_vec_local;
  return ptr;
}

Int glb_sum_real_d_vec(void* ptr, const Long size);

Int glb_sum_vec(Vector<RealD> vec);

Int glb_sum_val(RealD& x);

// -------------------

Int glb_sum_real_f_vec_local(void* ptr, const Long size);

API inline MpiSumFuncPtr& get_glb_sum_real_f_vec_ptr()
{
  static MpiSumFuncPtr ptr = glb_sum_real_f_vec_local;
  return ptr;
}

Int glb_sum_real_f_vec(void* ptr, const Long size);

Int glb_sum_vec(Vector<RealF> vec);

Int glb_sum_val(RealF& x);

// -------------------

Int glb_sum_byte_vec_local(void* ptr, const Long size);

API inline MpiSumFuncPtr& get_glb_sum_byte_vec_ptr()
{
  static MpiSumFuncPtr ptr = glb_sum_byte_vec_local;
  return ptr;
}

Int glb_sum_byte_vec(void* ptr, const Long size);

Int glb_sum_vec(Vector<Char> vec);

Int glb_sum_val(Char& x);

// -------------------

using MpiBcastFuncPtr = Int (*)(void*, const Long,
                                const Int);  // size is the data size.

Int bcast_byte_vec_local(void* ptr, const Long size, const Int root);

API inline MpiBcastFuncPtr& get_bcast_byte_vec_ptr()
{
  static MpiBcastFuncPtr ptr = bcast_byte_vec_local;
  return ptr;
}

Int bcast_byte_vec(void* ptr, const Long size, const Int root);

Int bcast_vec(Vector<Char> vec, const Int root);

Int bcast_vec(Vector<Long> vec, const Int root);

Int bcast_vec(Vector<Int> vec, const Int root);

Int bcast_vec(Vector<RealD> vec, const Int root);

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
