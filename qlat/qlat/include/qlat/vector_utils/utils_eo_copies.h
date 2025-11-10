// utils_eo_copies.h
// Gen Wang
// Jun. 2023

#ifndef UTILS_EO_COPIES_H
#define UTILS_EO_COPIES_H
#pragma once

#include "utils_vector_cs.h"

namespace qlat{

inline void qlat_map_eo_site(qlat::FieldM<int8_t, 1>& eo, const Geometry& geo)
{
  TIMER("qlat_map_eo_site");
  if(eo.initialized){
    const Geometry& geo_ = eo.geo();
    if(geo_ == geo){return ;}
  }
  eo.init(geo, 1);
  int8_t* res = (int8_t*) qlat::get_data(eo).data();
  ////only bool is not write thread safe
  qacc_for(isp, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(isp);
    Int site_eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    res[isp] = site_eo;
  });
}

///////src and res can be the same pointer
template <class Ty, Int civ>
void apply_eo_signT(Ty* sP, Ty* rP, qlat::FieldM<int8_t, 1>& eo, const int8_t dir = 1)
{
  TIMER("apply_eo_sign");
  const Geometry& geo = eo.geo();
  Qassert(eo.initialized);
  Qassert(dir == 0 or dir == 1 or dir == -1);
  int8_t* eP = (int8_t*) qlat::get_data(eo).data();
  const Int Nd = civ;
  //const Int Nd = civ * multiplicity;
  ///////DATA_TYPE typenum = get_data_type<Ty >();
  qacc_for(isp, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(isp);
    const Long index    = geo.offset_from_coordinate(xl, 1);
    //qlat::ComplexD sign = qlat::ComplexD(-1.0 * dir *(eP[index]*2 - 1), 0);
    Ty sign = double(-1.0 * dir *(eP[index]*2 - 1));
    for(Int ic=0;ic<Nd;ic++){rP[index*Nd+ic] = sign * sP[index*Nd+ic];}
  });
}

template <class Ty, Int civ>
void apply_eo_signT(qlat::vector<Ty* >& sP, qlat::vector<Ty* >& rP, qlat::FieldM<int8_t, 1>& eo, const int8_t dir = 1)
{
  Qassert(sP.size() == rP.size());
  for(unsigned int si=0;si<sP.size();si++)
  {
    apply_eo_signT<Ty, civ>(sP[si], rP[si], eo, dir);
  }
}


/////src and res can be tthe same pointer
template <class Ty, Int civ>
void apply_eo_sign(qlat::FieldM<Ty , civ>& src, qlat::FieldM<Ty , civ>& res, qlat::FieldM<int8_t, 1>& eo, const int8_t dir = 1)
{
  if(!src.initialized or !res.initialized){abort_r("src should be initialized with geo!\n");}
  const Geometry& geo = src.geo();
  if(!eo.initialized){qlat_map_eo_site(eo, geo);}
  Qassert(geo == eo.geo() and src.multiplicity == res.multiplicity);
  Qassert(civ == res.multiplicity);
  //qmessage("eo sizes civ %5d, multi %5d \n", civ, res.multiplicity);
  Ty*   sP = (Ty*  ) qlat::get_data(src).data();
  Ty*   rP = (Ty*  ) qlat::get_data(res).data();
  apply_eo_signT<Ty, civ>(sP, rP, eo, dir);
}

template <class Ty, Int civ>
void apply_eo_sign(std::vector<qlat::FieldM<Ty , civ> >& src, std::vector<qlat::FieldM<Ty , civ> >& res, qlat::FieldM<int8_t, 1>& eo, const int8_t dir = 1)
{
  Qassert(src.size() == res.size());
  for(unsigned int si = 0; si<src.size(); si++){
    apply_eo_sign<Ty, civ>(src[si], res[si], eo, dir);
  }
}

///////src and res can be the same pointer
template <class Ty, Int civ>
void apply_eo_zeros(Ty* sP, Ty* rP, qlat::FieldM<int8_t, 1>& eo, const int8_t dir = 1)
{
  TIMER("apply_eo_zeros");
  const Geometry& geo = eo.geo();
  Qassert(eo.initialized);
  int8_t* eP = (int8_t*) qlat::get_data(eo).data();
  Qassert(dir == 0 or dir == 1);
  ///////DATA_TYPE typenum = get_data_type<Ty >();
  qacc_for(isp, geo.local_volume(), {
    qlat::ComplexD sign = 1.0;
    if(eP[isp] == 0 and dir == 0){sign = 0;}////even to be zero
    if(eP[isp] == 1 and dir == 1){sign = 0;}////odd  to be zero
    //qlat::ComplexD sign = qlat::ComplexD(-1.0 * dir *(eP[isp]*2 - 1), 0);
    for(Int ic=0;ic<civ;ic++){rP[isp*civ+ic] = sign * sP[isp*civ+ic];}
  });
}

/////src and res can be tthe same pointer
template <class Ty, Int civ>
void apply_eo_zeros(qlat::FieldM<Ty , civ>& src, qlat::FieldM<Ty , civ>& res, qlat::FieldM<int8_t, 1>& eo, const int8_t dir = 1)
{
  if(!src.initialized or !res.initialized){abort_r("src should be initialized with geo!\n");}
  const Geometry& geo = src.geo();
  if(!eo.initialized){qlat_map_eo_site(eo, geo);}
  Ty*   sP = (Ty*  ) qlat::get_data(src).data();
  Ty*   rP = (Ty*  ) qlat::get_data(res).data();
  apply_eo_zeros<Ty, civ>(sP, rP, eo, dir);
}

template <class Ty, Int civ>
void apply_eo_zeros(std::vector<qlat::FieldM<Ty , civ> >& src, std::vector<qlat::FieldM<Ty , civ> >& res, qlat::FieldM<int8_t, 1>& eo, const int8_t dir = 1)
{
  Qassert(src.size() == res.size());
  for(unsigned int si = 0; si<src.size(); si++){
    apply_eo_zeros<Ty, civ>(src[si], res[si], eo, dir);
  }
}

template <class Ty, Int civ, int8_t Real>
void reduce_color_Re(qlat::vector_gpu<Ty >& res, qlat::FieldM<Ty , civ>& p0, qlat::FieldM<Ty , civ>& p1, const Int sign)
{
  TIMER("reduce_color");
  Qassert(p0.initialized and p1.initialized );
  const Geometry& geo = p0.geo();
  const Long V = geo.local_volume();
  const Ty* pa = (Ty*) qlat::get_data(p0).data();
  const Ty* pb = (Ty*) qlat::get_data(p1).data();
  Qassert(Long(res.size()) == V);
  Ty* resP = res.data();
  qacc_for(isp, V , {
    if(Real ==  1)for(Int ic=0;ic<civ;ic++){resP[isp] += (sign * (qlat::qconj(pa[isp*civ + ic]) * pb[isp*civ + ic]).real());}
    if(Real ==  0)for(Int ic=0;ic<civ;ic++){resP[isp] += (sign * (qlat::qconj(pa[isp*civ + ic]) * pb[isp*civ + ic])       );}
    if(Real == -1)for(Int ic=0;ic<civ;ic++){resP[isp] += (sign * (qlat::qconj(pa[isp*civ + ic]) * pb[isp*civ + ic]).imag());}
  });
}

template <class Ty, Int civ>
void reduce_color(qlat::vector_gpu<Ty >& res, qlat::FieldM<Ty , civ>& p0, qlat::FieldM<Ty , civ>& p1, const Int sign, const int8_t Real )
{
  if(Real ==  1){reduce_color_Re<Ty, civ,  1>(res, p0, p1, sign);}
  if(Real ==  0){reduce_color_Re<Ty, civ,  0>(res, p0, p1, sign);}
  if(Real == -1){reduce_color_Re<Ty, civ, -1>(res, p0, p1, sign);}
}

template <class Ty, Int civ>
void reduce_color(std::vector<qlat::vector_gpu<Ty > >& res, std::vector<qlat::FieldM<Ty , civ> >& p0,
  std::vector<qlat::FieldM<Ty , civ> >& p1, const Int sign, const int8_t Real )
{
  const unsigned int nsrc = res.size();
  Qassert(p0.size() == nsrc and p1.size() == nsrc);
  for(unsigned int ni=0;ni<nsrc;ni++)
  {
    reduce_color(res[ni], p0[ni], p1[ni], sign, Real);
  }
}

template <class Ty>
void reduce_color(qlat::vector<Ty >& resC, qlat::vector_gpu<Ty >& p0, qlat::vector_gpu<Ty >& p1, const Int sign, const Int t0, const Int clear, const Geometry& geo)
{
  TIMER("reduce_color");
  const Long V = geo.local_volume();
  Qassert(p0.size() == p1.size() and p0.size()%V == 0);
  const Ty* pa = (Ty*) qlat::get_data(p0).data();
  const Ty* pb = (Ty*) qlat::get_data(p1).data();
  qlat::vector_gpu<Ty > res;res.resize(V);
  Ty* resP = res.data();
  const Int Nc = p0.size()/V;
  qacc_for(isp, V , {
    for(Int ic=0;ic<Nc;ic++){
      resP[isp] += ( sign * qlat::qconj(pa[isp*Nc + ic]) * pb[isp*Nc + ic] );
    }
  });

  /////reduce zero
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  vec_corrE(res.data(), resC, fd, 1, clear, Coordinate(), Ty(1.0, 0.0), t0);

  //qlat::vector_gpu<Ty > tmpT;tmpT.resize(fd.Nt);qlat::set_zero(tmpT);
  //reduce_vec(res.data(), tmpT.data(), fd.Nx*fd.Ny*fd.Nz, fd.Nt);
  //sum_all_size(tmpT.data(), tmpT.size(), true);

  //Qassert(resC.size() == nt);
  //const Int init = fd.init;
  //const Int nt = fd.nt;
  //Ty* d = resC.data();Ty* a = tmpT.data();
  //qacc_for(it, Nt, {d[(it + init - t0 + nt)%nt] += a[it]});
}

inline void get_index_mappings_reverse(qlat::vector<Long >& map, const Geometry& geo)
{
  TIMER("get_index_mappings_reverse");
  const Long V = geo.local_volume();
  const Long Vh = V / 2;

  if(map.size() == V){return ;}
  else{map.resize(V);}
  qacc_for(qlat_idx_4d, V , {
    const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
    const Int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    const Long quda_idx = eo * Vh + qlat_idx_4d / 2;
    ////map[qlat_idx_4d] = quda_idx;
    map[quda_idx] = qlat_idx_4d;
  });
}

inline void get_index_mappings(qlat::vector<Long >& map, const Geometry& geo)
{
  TIMER("get_index_mappings");
  const Long V = geo.local_volume();
  const Long Vh = V / 2;

  if(map.size() == V){return ;}
  else{map.resize(V);}
  qacc_for(qlat_idx_4d, V , {
    const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
    const Int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    const Long quda_idx = eo * Vh + qlat_idx_4d / 2;
    map[qlat_idx_4d] = quda_idx;
  });
}

/////GPU order with color to be outside even odd
template <class T1, class Ty, Int dir>
void qlat_cf_to_quda_cfT(T1*  quda_cf, Ty* src, const Int Dim, const Geometry& geo, qlat::vector<Long >& map_)
{
  TIMER("qlat_cf_to_quda_cf");
  Qassert(sizeof(T1) == sizeof(ComplexT<double>));// quda can only use double interface
  const Long V = geo.local_volume();
  Long Vh = V / 2;
  if(map_.size() != V){get_index_mappings(map_, geo);}
  qlat::vector<Long >& map = map_;
  qacc_for(qlat_idx_4d, V, {
    const Long quda_idx = map[qlat_idx_4d];
    const Long eo = quda_idx / Vh;
    const Long qi = quda_idx % Vh;
    for(Int dc = 0; dc < Dim; dc++)
    {
      //if(dir == 1){quda_cf[ quda_idx*Dim + dc] = src[qlat_idx_4d*Dim + dc];}
      //if(dir == 0){src[qlat_idx_4d*Dim + dc] = quda_cf[quda_idx*Dim + dc];}
      if(dir == 1){quda_cf[(eo*Dim + dc)*Vh + qi] = src[qlat_idx_4d*Dim + dc];}
      if(dir == 0){src[qlat_idx_4d*Dim + dc] = quda_cf[(eo*Dim + dc)*Vh + qi];}
    }
  });
}

template <class T1, class Ty>
void qlat_cf_to_quda_cf_pointer(T1*  quda_cf, Ty* src, const Int Dim, const Geometry& geo, qlat::vector<Long >& map)
{
  qlat_cf_to_quda_cfT<T1, Ty, 1>(quda_cf, src, Dim, geo, map);
}

template <class T1, class Ty>
void quda_cf_to_qlat_cf_pointer(Ty* res, T1*  quda_cf, const Int Dim, const Geometry& geo, qlat::vector<Long >& map)
{
  qlat_cf_to_quda_cfT<T1, Ty, 0>(quda_cf, res, Dim, geo, map);
}

template<typename Ty, Int dir, Int mode>
void copy_eo_cs_to_fieldMT(qlat::vector_gpu<Ty >& res, const Int civ, const Geometry& geo, vector_cs<Ty >& even, vector_cs<Ty >& odd,
  Int e0, Int e1, Int o0, Int o1, qlat::vector<Long >& map)
{
  TIMER_FLOPS("copy_eo_cs_to_fieldM");
  Qassert(even.initialized and odd.initialized);
  // //fft_desc_basic& fd = get_fft_desc_basic_plan(res.geo());
  // const bool GPU = even.GPU;
  const Long V = geo.local_volume();
  const Long Vh= V/2;
  if(dir == 1){if(Long(res.size()) != civ*V){res.resize(civ * V, QMSYNC);}}
  if(dir == 0){Qassert(Long(res.size()) == V * civ );}
  ////Qassert(check_GPU_same(get_type_mem(res.GPU), even.GPU));
  Int GPU_copy = check_GPU_multi(get_type_mem(res.GPU), even.GPU);
  Qassert(GPU_copy != -2);
  //Qassert(res.GPU == even.GPU);

  const Int DIM = 3;
  Qassert(civ % DIM == 0);
  if(map.size() == 0){get_index_mappings_reverse(map, geo);}
  Int nvec = civ / DIM;
  Qassert(e1 > e0 and o1 > o0);
  Int ne = e1 - e0;
  Int no = o1 - o0;
  Qassert(ne == no);
  Qassert(ne <= nvec and no <= nvec and e1 <= even.nvec and o1 <= odd.nvec);
  const Int b_size = even.b_size;
  const Int btotal = even.btotal;

  Qassert(btotal * b_size == DIM*V/2);
  qlat::vector<Ty** > eP;eP.resize(2*ne);
  //qlat::vector<Ty** > oP;oP.resize(no);
  ////Ty** eP = even.get_pointers(ni)
  for(Int ei=0;ei<ne;ei++){eP[ei]      = even.get_pointers(ei + e0);}
  for(Int oi=0;oi<no;oi++){eP[ne + oi] =  odd.get_pointers(oi + o0);}
  const Long* mapP = (Long*) qlat::get_data(map).data();

  ////NtNzNyNx, DIM x nvec
  Ty* r = (Ty*) qlat::get_data(res).data();
  //if(mode == 1)
  {
  qGPU_for(qi, Vh, GPU_copy, {
    //const Int ni = ci / DIM;
    //const Int c  = ci % DIM;
    for(Int ni=0;ni<ne;ni++)
    for(Int eo = 0; eo < 2;eo++)
    {
      const Long quda_idx = eo*Vh + qi;
      const Long qlat_idx_4d = mapP[quda_idx];
      Ty* rr = &r[qlat_idx_4d*civ];
      for(Int c = 0; c < DIM ; c++)
      {
        Long bv =  qi*DIM + c;////default mode with c to be inside
        ////const Long bv = qi*DIM + c; ////quda vectors this order
        if(mode == 1){bv = c*Vh + qi;} ////quda vectors this order
        const Long bi = bv / b_size;
        const Long bj = bv % b_size;
        {
          {
          if(dir == 1){rr[c*nvec + ni] = eP[eo*ne + ni][bi][bj];}
          if(dir == 0){eP[eo*ne + ni][bi][bj] = rr[c*nvec + ni];}
          }
        }
      }
    }
  });}

  timer.flops += double(V) * DIM * ne * sizeof(Ty);
}

template<typename Ty>
void copy_eo_cs_to_fieldM(qlat::vector_gpu<Ty >& res, const Int civ, const Geometry& geo, vector_cs<Ty >& even, vector_cs<Ty >& odd,
  Int e0, Int e1, Int o0, Int o1, qlat::vector<Long >& map, Int mode = 0, Int dir = 1)
{
  if(dir == 0 and mode == 0){copy_eo_cs_to_fieldMT<Ty, 0, 0>(res, civ, geo, even, odd, e0, e1, o0, o1, map);return ;}
  if(dir == 0 and mode == 1){copy_eo_cs_to_fieldMT<Ty, 0, 1>(res, civ, geo, even, odd, e0, e1, o0, o1, map);return ;}
  if(dir == 1 and mode == 0){copy_eo_cs_to_fieldMT<Ty, 1, 0>(res, civ, geo, even, odd, e0, e1, o0, o1, map);return ;}
  if(dir == 1 and mode == 1){copy_eo_cs_to_fieldMT<Ty, 1, 1>(res, civ, geo, even, odd, e0, e1, o0, o1, map);return ;}
  abort_r("Can not find cases");
}

template<typename Ty>
void copy_fieldM_to_eo_cs(vector_cs<Ty >& even, vector_cs<Ty >& odd, qlat::vector_gpu<Ty >& res, const Int civ, const Geometry& geo, 
  Int e0, Int e1, Int o0, Int o1, qlat::vector<Long >& map, Int mode = 0)
{
  copy_eo_cs_to_fieldM(res, civ, geo, even, odd, e0, e1, o0, o1, map, mode, 0);
}

}

#endif

