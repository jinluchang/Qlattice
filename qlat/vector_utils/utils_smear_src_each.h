// utils_smear_src.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_SMEAR_SRC_Each_H
#define UTILS_SMEAR_SRC_Each_H
#pragma once

#include "general_funs.h"
#include <qlat/qcd-utils.h>
#include "utils_smear_src.h"
#include "utils_Vec_redistribute.h"
//#include <openacc.h>
//#include "mpi-ext.h" /* Needed for CUDA-aware check */
 
namespace qlat{

///////////psrc order, bfac, c, d0, t,z,y,x
template <class T, int bfac, int d0, int dirL>
__global__ void smear_global4(T* pres, const T* psrc, const T* gf, double bw, double norm,
  unsigned long Nvol, long* map_bufD)
{

  __shared__ T ls[(dirL*2)*3*3 + 1];
  __shared__ T ps[(dirL*2)*bfac*3*d0 + 1];
  ///tid should larger than 9
  unsigned int   tid   =  threadIdx.x;
  unsigned long  index =  blockIdx.y*gridDim.x + blockIdx.x;
  unsigned int   ns = blockDim.x;

  if(index < Nvol){

  ///////(2*dirL) --> c0 , c1
  ///////TODO need to check dirL is not 3
  unsigned int off = tid;
  {
    const T* gf_t = &gf[index*(2*dirL)*9];
    while(off < (2*dirL)*9){
      //int dir = off/9;int c = off%9;
      //ls[(c/3)*(2*dirL)*3 + dir*3 + c%3] = gf[index*(2*dirL)*9 + off];off += ns;
      //ls[off] = gf_t[off];
      int dir = off/9;
      int c0 = (off/3)%3;
      int c1 =  off%3;
      //ls[c0*(2*dirL)*3 + dir*3 + c1 ] = gf_t[(dir + (dirL/4))*9 + c0*3 + c1];
      ls[c0*(2*dirL)*3 + dir*3 + c1 ] = gf_t[off];
      //ls[(c0*(2*dirL)+dir)*3 + c1 ] = gf_t[off];
      //ls[c0][dir*3+c1] = gf_t[off];
      off += ns;
    }
  }

  const long res_off = index;
  T* wm = &pres[res_off*bfac*3*d0];

  ///////(2*dir) --> bi, c1 , d0
  for (int dir = -dirL; dir < dirL; ++dir){
    long src_off = map_bufD[(dir+4)*Nvol + index];
    const T* src_t = &psrc[src_off*bfac*3*d0];
    //T* res_t     = &ps[(dir+dirL)*bfac*3*d0];
    ////T* res_t     = &ps[dir+dirL];
    off = tid;
    while(off < bfac*3*d0){
      //res_t[off] = src_t[off]; off+=ns;
      unsigned int bi = off/(3*d0);
      unsigned int c  = (off%(3*d0))/d0;
      unsigned int di = off%d0;
      ps[(bi*d0 + di)*(2*dirL*3) + (dir+dirL)*3 + c] = src_t[off];
      off+=ns;
    }
  }
  __syncthreads();

  T tmp = 0.0;
  off = tid;
  while(off < bfac*3*d0){
    unsigned int c0 = off/(bfac*d0);
    unsigned int bi = (off%(bfac*d0))/d0;
    unsigned int di = off%d0;

    T* a = &ls[c0*(2*dirL*3)];
    T* b = &ps[(bi*d0 + di)*(2*dirL*3)];
    //T* b = &ps[(bi*d0 + di)];

    tmp = 0.0; 

    for(int dir=0;dir<(2*dirL)*3; dir++)
    {
      //tmp += b[dir*bfac*d0] * a[dir];
      tmp += b[dir] * a[dir];
    }

    wm[(bi*3 + c0)*d0 + di] = norm*wm[ (bi*3 + c0)*d0 + di ] + norm*bw*tmp;
    off += ns;
  }

  }

}

template <class T>
__global__ void cpy_smear_global4(T* p0, T* p1, int bfac, int Nvol, long* map_buf)
{
  unsigned int   tid   =  threadIdx.x;
  unsigned int   ns    =  blockDim.x;

  unsigned long  index =  blockIdx.y*gridDim.x + blockIdx.x;

  if(index < Nvol){
    long off0 = map_buf[index];
    long off1 = index;

    T* pres  = &p0[off0*bfac];
    T* psrc  = &p1[off1*bfac];
    unsigned int off = tid;
    while(off < bfac){pres[off] = psrc[off]; off += ns;}
  }
}

template <class T>
void cpy_vec_GPU(T* p0, T* p1, int bfac, smear_fun& smf)
{
  unsigned long Nvol = smf.Nvol;
  int  nt = 16;
  long sn = long(std::sqrt(Nvol))+1;
  dim3 dimBlock(  nt, 1, 1);
  dim3 dimGrid(   sn,sn, 1);
  cpy_smear_global4<T ><<< dimGrid, dimBlock >>>(p0, p1, bfac, smf.Nvol, &smf.map_buf[0]);
  qacc_barrier(dummy);
}

template <class T>
void smear_propagator4(T* prop, const T* gf,
                      const double width, const int step, T* prop_buf, smear_fun& smf, int bfac = 4, int d0 = 12)
{
  const double aw   = 3.0*width*width/(2*step);
  const double bw = width*width/(4.0*step - 6.0*width*width);

  unsigned long Nvol = smf.Nvol;

  //int bfac =  4;
  //int d0   = 12;
  int dirL = 3;

  int nt  = 3*3*9;
  if(bfac*d0 <= 12){ nt =        32;}
  if(bfac*d0 <=  6){ nt = 3*bfac*d0;}
  //int nt = 32;
  dim3 dimBlock(nt, 1, 1);
  long sn = Nvol;
  dim3 dimGrid( sn, 1, 1);

  for (int i = 0; i < step; ++i) {
    TIMER("Matrix multiply");
    {
    TIMER("Copy prop");
    cpy_vec_GPU(prop_buf, prop, (bfac*3)*d0, smf);
    }

    refresh_expanded_GPU(prop_buf, smf);

    bool cfind = false;
    {
    TIMER("kernal");
    #define smear_macros(ba,da) if(bfac == ba and d0 ==  da){cfind = true; \
      smear_global4<T, ba, da, 3 ><<< dimGrid, dimBlock >>>(prop, prop_buf, gf, bw, 1-aw, Nvol, &smf.map_bufD[0]);}

    ///////macros for color 3 * dirac 4 in inner prop
    smear_macros(   1,  4);
    smear_macros(   2,  4);
    smear_macros(   3,  4);
    smear_macros(   4,  4);
    smear_macros(   5,  4);
    smear_macros(   6,  4);
    smear_macros(   7,  4);
    smear_macros(   8,  4);
    smear_macros(   9,  4);
    smear_macros(  10,  4);
    smear_macros(  11,  4);
    smear_macros(  12,  4);

    ///////macros for inner color 3 and all outter prop
    smear_macros(   1,  1);
    smear_macros(   2,  1);
    smear_macros(   3,  1);
    smear_macros(   4,  1);
    smear_macros(   5,  1);
    smear_macros(   6,  1);
    smear_macros(   7,  1);
    smear_macros(   8,  1);
    smear_macros(   9,  1);
    smear_macros(  10,  1);
    smear_macros(  11,  1);
    smear_macros(  12,  1);

    smear_macros(  24,  1);

    smear_macros(   1, 48);
    smear_macros(   4, 12);
    smear_macros(   3, 16);

    qacc_barrier(dummy);
    #undef smear_macros
    }
    qassert(cfind);
  }

}

template <class T>
void rotate_Vec_prop(Propagator4dT<T>& prop, qlat::vector_acc<T > &propT, unsigned int NVmpi, unsigned int groupP, int dir = 0)
{
  TIMER("Rotate Vec prop");
  ///unsigned int NVmpi = fd.mz*fd.my*fd.mx;
  ///int groupP = (12+NVmpi-1)/NVmpi;
  qassert(groupP <= 12);qassert(groupP * NVmpi >= 12);
  qassert(groupP >   0);

  const Geometry& geo = prop.geo();
  long long Nvol =  geo.local_volume();
  if(dir == 0)propT.resize(NVmpi*Nvol*groupP*12);

  qacc_for(index, long(Nvol), {
    qlat::WilsonMatrixT<T>& v0 =  prop.get_elem(index);

    for(int c1 = 0;c1 < 3; c1++)
    for(int d1 = 0;d1 < 4; d1++)
    for(int c0 = 0;c0 < 3; c0++)
    for(int d0 = 0;d0 < 4; d0++)
    {
      int off1 = c1*4+d1;
      int n0 = off1/groupP;
      int n1 = off1%groupP;

      int off0 = c0*4+d0;
      //LInt off = (c0*3+c1)*16+d0*4 + d1;
      //long offP = n0*Nvol*groupP*12 + index*groupP*12 + n1*12 + off0;
      long offP = ((n0*Nvol + index)*groupP + n1)*12 + off0;
      if(dir == 0){propT[offP] = v0(d0*3+c0, d1*3+c1);}
      if(dir == 1){v0(d0*3+c0, d1*3+c1) = propT[offP];}
    }
  });


}

template <class T, class Tg>
void smear_propagator_gpu4(Propagator4dT<T>& prop, const GaugeFieldT<Tg >& gf, const double width, const int step, int modedis=1)
{
  long long Tfloat = 0;
  double mem       = 0.0;
  const Geometry& geo = prop.geo();

  {long long Lat = geo.local_volume();
  int nsrc = 12;
  long long vGb = Lat *nsrc*4;
  int Fcount = 3*(3*6 + 2*2); 
  int direction   = 6;
  Tfloat = step*direction*vGb*Fcount;
  mem = (Lat*nsrc*12 + Lat*4*9)*8.0;}
  ////timer.flops += Tfloat;
  print0("Memory size %.3e GB, %.3e Gflop \n", 
    mem/(1024.0*1024*1024), Tfloat/(1024.0*1024*1024));

  // set_left_expanded_gauge_field(gf1, gf)
  // prop is of qnormal size
  long Nvol = geo.local_volume();
  qlat::vector_acc<T > gfE;gfE.resize(6*Nvol*9);
  const int dir_limit = 3;
  qacc_for(index,  geo.local_volume(),{
    for (int dir = -dir_limit; dir < dir_limit; ++dir) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const ColorMatrixT<Tg > link =
          dir >= 0 ? gf.get_elem(xl, dir)
                   : (ColorMatrixT<Tg >)matrix_adjoint(
                         gf.get_elem(coordinate_shifts(xl, dir), -dir - 1));
      for(int ci=0; ci<9; ci++){
        gfE[index*(dir_limit*2)*9 + (dir+dir_limit)*9 +  ci] = link.p[ci];
      }
    }
  });

  const Geometry geo1 = geo_resize(prop.geo(), Coordinate(1, 1, 1, 0), Coordinate(1, 1, 1, 0));
  Propagator4dT<T> prop_buf;
  prop_buf.init(geo1);
  smear_fun smf;
  smf.init_mem(prop.geo(), prop_buf.geo(), sizeof(WilsonMatrixT<T>)/sizeof(T), sizeof(T));


  //EigenV pE; EigenV p_buf;
  //prop_to_EigenV(prop, pE);
  //p_buf.resize(pE.size());
  if(modedis == 0)
  {
    //fft_desc_basic fd(geo);
    //Vec_redistribute vec_large(fd, 1);
    //unsigned int NVmpi = fd.mz*fd.my*fd.mx;
    //qassert(NVmpi < 12);int groupP = (12+NVmpi-1)/NVmpi;
    //print0("====Vec setup, NVmpi %d, groupP %d \n", NVmpi, groupP);
    //qlat::vector<T > propT;qlat::vector<T > propT_buf;
    //rotate_Vec_prop(prop, propT, NVmpi, groupP, 0);
    //propT_buf.resize(propT.size());
    //vec_large.reorder((T*) &propT[0],(T*) &propT_buf[0], 1, groupP*12 ,   0);
    //vec_large.reorder((T*) &propT[0],(T*) &propT_buf[0], 1, groupP*12 , 100);
    //rotate_Vec_prop(prop, propT, NVmpi,groupP, 1);

    TIMER_FLOPS("==compute time");
    timer.flops += Tfloat;
    int bfac = 4; int d0 = 12;
    smear_propagator4((T*) prop.field[0].p, (T*) gfE.data(), width, step, (T*) prop_buf.field[0].p, smf, bfac, d0);

    return ;
  }

  ///////Redistribute the data
  fft_desc_basic fd(geo);
  Vec_redistribute vec_large(fd);

  unsigned int NVmpi = fd.mz*fd.my*fd.mx;
  qassert(NVmpi < 12);int groupP = (12+NVmpi-1)/NVmpi;
  print0("====Vec setup, NVmpi %d, groupP %d \n", NVmpi, groupP);

  size_t g_sizeT = NVmpi*size_t(gfE.size());
  qassert(g_sizeT < 2147483647);/////Check long limit reached?

  qlat::vector_acc<T > gfET;gfET.resize(NVmpi*gfE.size());
  qlat::vector_acc<T > gfET_buf;gfET_buf.resize(NVmpi*gfE.size());

  //qlat::vector<T > gfET0;gfET0.resize(NVmpi*gfE.size());
  qacc_for(index, gfE.size(),{
    for(long vi=0;vi<NVmpi;vi++)gfET[vi*gfE.size() + index] = gfE[index];
  });

  //for(long vi=0;vi<NVmpi;vi++){
  //  qthread_for(index, gfE.size(),{
  //    gfET_buf[vi*gfE.size() + index] = gfE[index];
  //    //gfET0[vi*gfE.size() + index] = gfE[index];
  //    //gfET[vi*gfE.size() + index] = gfE[index];
  //  });
  //}

  vec_large.reorder((T*) gfET.data(),(T*) gfET_buf.data(), 1, (dir_limit*2)*9 ,   0);

  //////for(long long a=0;a< gfET.size();a++){
  //////  gfET[a] = gfET_buf[a];
  //////}

  //////int groupP = 1;

  qlat::vector_acc<T > propT;
  qlat::vector_acc<T > propT_buf;

  //////rotate_Vec_prop(prop, propT, NVmpi, groupP, 0);

  //for(long long a=0;a< propT_buf.size();a++){
  //  propT[a] = propT_buf[a];
  //}

  //rotate_Vec_prop(prop, propT_buf, NVmpi,groupP, 0);
  //diff_EigenM(propT_buf, propT, "Check");

  smf.init_distribute(prop.geo());

  rotate_Vec_prop(prop, propT, NVmpi, groupP, 0);
  propT_buf.resize(propT.size());
  //rotate_prop(prop,0);
  {
    TIMER_FLOPS("==compute time");

    {TIMER("Vec prop");vec_large.reorder((T*) propT.data(),(T*) propT_buf.data(), 1, groupP*12 ,   0);}

    int bfac = groupP; int d0 = 4;
    //int bfac =1 ; int d0 = 48;
    smear_propagator4((T*) propT.data(), (T*) gfET.data(), width, step, (T*) propT_buf.data(), smf, bfac, d0);
    //smear_propagator4((T*) &propT[0], &gfET[0], width, step, (T*) prop_buf.field[0].p, smf, bfac, d0);
    //smear_propagator4((T*) &prop.field[0].p, &gfET[0], width, step, (T*) prop_buf.field[0].p, smf, bfac, d0);

    timer.flops += Tfloat;
    {TIMER("Vec prop");vec_large.reorder((T*) propT.data(),(T*) propT_buf.data(), 1, groupP*12 , 100);}
  }
  rotate_Vec_prop(prop, propT, NVmpi,groupP, 1);
  //rotate_prop(prop,1);

  //for(long long a=0;a< propT_buf.size();a++){
  //  propT_buf[a] = propT[a];
  //}


  ////rotate_Vec_prop(prop, propT, NVmpi, groupP, 1);


}






}
#endif
