// utils_grid_multi.h
// Gen Wang
// Apr. 2025

#ifndef UTILS_GRID_MULTI_H
#define UTILS_GRID_MULTI_H
#pragma once

#include "general_funs.h"
#include "utils_fft_desc.h"

namespace qlat
{

/* in memeory reshape
  geoA  original geometry
  geo0 current  geometry
  geo1 new      geometry
  src to res : with civ to original geoA
*/
template <class Ty>
void grid_memory_reshape(qlat::vector<Ty* >& res, qlat::vector<Ty* >& src, const Int civ, const Geometry& geo1, const Geometry& geo0, const Geometry& geoA )
{
  TIMERA("grid_memory_reshape");
  Qassert(geoA.node_site != qlat::Coordinate());
  Qassert(geo0.node_site != qlat::Coordinate());
  Qassert(geo1.node_site != qlat::Coordinate());
  const Int biva = src.size();
  Qassert(res.size() == biva);
  Qassert(biva > 0 and civ > 0);
  const Coordinate NA = geoA.total_site();
  const Coordinate N0 = geo0.total_site();
  const Coordinate N1 = geo1.total_site();
  const Long Vol = geoA.local_volume();

  // copy if pointer the same
  if(N0 == N1){
    for(Int bi=0;bi<biva;bi++){
      if(res[bi] != src[bi]){
        cpy_GPU(res[bi], src[bi], Vol * civ, 1, 1, QFALSE);
      }
    }
    qacc_barrier(dummy);
    return ;
  }

  Coordinate nfac0;
  Coordinate nfac1;
  LInt       Vfac0 = 1;
  LInt       Vfac1 = 1;
  for(Int i=0;i<4;i++){
    Qassert(NA[i] % N0[i] == 0 and NA[i] % N1[i] == 0);
    Qassert(NA[i] >= N0[i] and NA[i] >= N1[i]);
    nfac0[i] = NA[i] / N0[i];
    nfac1[i] = NA[i] / N1[i];
    Vfac0 *= nfac0[i];
    Vfac1 *= nfac1[i];
  }

  const Int bmax = 24;
  Int bL = (bmax + civ - 1) / civ; 
  if(biva*civ < bmax){bL = biva;}

  const size_t Ndata = bL * Vol * civ*sizeof(Ty) / sizeof(char);
  VectorGPUKey gkey(0, std::string("shift_vec_buf"), 1);
  vector_gpu<char >& buf = get_vector_gpu_plan<char >(gkey);
  buf.resizeL(Ndata);
  Ty* bufP = (Ty*) buf.data();
  const double mem = 1.0 * Ndata / (1024 * 1024 * 1024.0);
  if(mem > 0.5){
    qmessage("WARNING! large memeory used for rehape! %.3f GB\n ", mem);
  }

  std::vector<Long > jobA = job_create(biva, bL);
  qlat::vector<Ty* > bP;bP.resize(bL);
  for(Int bi=0;bi<bL;bi++){
    bP[bi] = &bufP[bi * Vol * civ];
  }

  for(LInt jobi=0;jobi < jobA.size()/2; jobi++)
  {
    const Long bini = jobA[jobi*2 + 0]; const Long bcut = jobA[jobi*2+1];

    qacc_for(index, Vol, {
      const Coordinate xl = geoA.coordinate_from_index(index);
      Coordinate x0;
      Coordinate x1;
      Coordinate y0;
      Coordinate y1;
      for(Int i=0;i<4;i++)
      {
        x0[i] = xl[i] / nfac0[i];
        x1[i] = xl[i] / nfac1[i];

        y0[i] = xl[i] % nfac0[i];
        y1[i] = xl[i] % nfac1[i];
      }

      const LInt c0 = ( (y0[3] * nfac0[2] + y0[2]) * nfac0[1] +  y0[1] ) * nfac0[0] + y0[0];
      const LInt c1 = ( (y1[3] * nfac1[2] + y1[2]) * nfac1[1] +  y1[1] ) * nfac1[0] + y1[0];
      const LInt g0 = geo0.index_from_coordinate(x0);
      const LInt g1 = geo1.index_from_coordinate(x1);
      const LInt off0 = (g0*Vfac0 + c0) * civ;
      const LInt off1 = (g1*Vfac1 + c1) * civ;
      for(Int bi = 0; bi < bcut; bi++)
      {
        Ty* sP = &src[bini + bi][ off0];
        Ty* rP = &bP[ bi][        off1] ;
        for(Int ci = 0; ci < civ; ci++)
        {
          rP[ci] = sP[ci];
        }
      }
    });

    for(Int bi = 0; bi < bcut; bi++)
    {
      cpy_GPU(res[bini + bi], bP[bi], Vol * civ, 1, 1, QFALSE);
    }
    qacc_barrier(dummy);

  }
}

template <class Ty>
void grid_memory_reshape(qlat::vector<Ty* >& res, qlat::vector<Ty* >& src, const Int civ,
  const Coordinate& n1, const Coordinate& n0, const Coordinate& nA)
{
  const Geometry& geoA = get_geo_cache(nA);
  const Geometry& geo0 = get_geo_cache(n0);
  const Geometry& geo1 = get_geo_cache(n1);
  grid_memory_reshape(res, src, civ, geo1, geo0, geoA);
}

}

#endif
