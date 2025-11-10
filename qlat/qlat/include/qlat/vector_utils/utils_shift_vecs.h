// utils_FFT_GPU.h
// Gen Wang
// Sep. 2021

#ifndef UTILS_SHIFT_VECS_H
#define UTILS_SHIFT_VECS_H
#pragma once

#include "general_funs.h"
#include "utils_fft_desc.h"
#include "utils_grid_multi.h"

//////TODO
//////GPU support, Td template, shift of fieldM

namespace qlat
{

template<class Ta>
inline void free_vector_8(std::vector<Ta* >& RES)
{
  for(unsigned int i=0;i<RES.size();i++)
  {
    if(RES[i] != NULL){
      RES[i]->resize(0);
      delete RES[i];
      RES[i] = NULL;
    }
  }
  RES.resize(0);
}

template<class Ta>
inline void init_vector_8(std::vector<Ta* >& RES, const Int size = 8)
{
  if(int(RES.size()) != size){
    free_vector_8(RES);
    RES.resize(size);
    for(Int i=0;i<size;i++)
    {
      if(RES[i] == NULL){
        RES[i] = new Ta(0);
      }
    }
  }
}

struct shift_vec{
  bool initialized;
  bool flag_shift_set;
  bool GPU;

  Int rank;Int Nmpi;
  Int nx,ny,nz,nt;
  Coordinate total_site;
  LInt noden;
  LInt vol;
  LInt Nvol;
  Int Nx,Ny,Nz;

  Int N0,N1,N2,Nt;
  qlat::vector<Int> Nv,nv;

  //fft_desc_basic fd;

  //////Shift under periodic condition
  Int periodic;

  //std::vector<std::vector<Int> > sendlist;
  //std::vector<std::vector<Int> > recvlist;

  Int civ,biva;
  LInt Length;

  Int dir_cur;
  std::vector<std::vector<Int > > rank_sr;

  std::vector<qlat::vector<LInt >* > buffoffa;
  std::vector<qlat::vector<LInt >* > buffoffb;
  std::vector<qlat::vector<LInt >* > sendoffa;
  std::vector<qlat::vector<LInt >* > sendoffb;
  std::vector<qlat::vector<LInt >* > sendoffx;

  std::vector<qlat::vector_gpu<char >* > sendbufP;
  std::vector<qlat::vector_gpu<char >* > recvbufP;

  qlat::vector_gpu<char >* zeroP;
  qlat::vector_gpu<char >* bufsP;
  qlat::vector_gpu<char >* bufrP;

  std::vector<size_t > MPI_size;
  unsigned int MPI_off;
  MPI_Datatype MPI_curr;

  unsigned int bsize;
  move_index mv_civ;

  void* gauge;
  Int gauge_is_double;
  Int gbfac;int gd0;
  bool Conj;bool src_gauge;

  inline void init(fft_desc_basic &fds, bool GPU_set = true);
  inline void init(const Geometry& geo, bool GPU_set);
  inline void init(const Coordinate& site, bool GPU_set);

  shift_vec(){
    initialized = false;
    flag_shift_set = false;
    civ = 0;
    GPU = true;
    gauge = NULL;
  }
  shift_vec(fft_desc_basic &fds, bool GPU_set = true){   init(fds, GPU_set);}
  shift_vec(const Geometry& geo, bool GPU_set = true){   init(geo, GPU_set);}
  shift_vec(const Coordinate& site, bool GPU_set = true){init(site, GPU_set);}

  inline void print_info();
  ~shift_vec(){
    initialized = false;
    flag_shift_set = false;bsize = 0;
    dir_cur = 0;civ = -1;biva = -1;
    periodic = 1;

    for(Int dir=0;dir<8;dir++){clear_mem_dir(dir);}
    MPI_size.resize(0);

    zeroP->resize(0);
    bufsP->resize(0);
    bufrP->resize(0);
    delete zeroP;zeroP = NULL;
    delete bufsP;bufsP = NULL;
    delete bufrP;bufrP = NULL;

    rank_sr.resize(0);
    free_vector_8(buffoffa);
    free_vector_8(buffoffb);
    free_vector_8(sendoffa);
    free_vector_8(sendoffb);
    free_vector_8(sendoffx);

    free_vector_8(sendbufP);
    free_vector_8(recvbufP);
  }

  inline void shift_set();

  template<typename Ty>
  void set_MPI_size(Int biva_or, Int civ_or, Int dir_or = 0);
  template<typename Ty>
  void set_MPI_size(Int dir_or);

  template<typename Ty, Int flag>
  void write_send_recv(Ty* src, Ty* res);

  template<typename Ty>
  void call_MPI(Ty *src, Ty *res,Int dir_or);

  template<typename Ty>
  void shift_Evec(std::vector<qlat::vector<Ty > > &srcE,std::vector<qlat::vector<Ty > > &srcEf,std::vector<Int >& iDir,Int civ_or);

  template<typename Ty>
  void shift_vecs(std::vector<Ty* > &src,std::vector<Ty* > &res,std::vector<Int >& iDir ,Int civ_or);

  template<typename Ty>
  void shift_vecP(Ty* src,Ty* res,std::vector<Int >& iDir ,Int civ_or);

  template<typename Ty>
  void shift_vecs_dir(Ty* src, Ty* res, Int civ_, Int mu, Int sign);

  template<typename Ty>
  void shift_vecs_dir(qlat::Field<Ty>& src, qlat::Field<Ty>& res, Int mu, Int sign);

  template<typename Ty, Int civ_>
  void shift_vecs_dir(std::vector<qlat::FieldM<Ty, civ_> >& src, std::vector<qlat::FieldM<Ty, civ_> >& res, Int mu, Int sign);

  template<typename Ty>
  void shift_vecs_dirG(std::vector<qlat::FieldG<Ty> >& src, std::vector<qlat::FieldG<Ty> >& res, Int mu, Int sign);

  template<typename Ta>
  void set_gauge(Ta* gauge_, Int gbfac_, Int gd0_, bool Conj_=false, bool src_gauge_ = false)
  { 
    using D = typename IsBasicDataType<Ta>::ElementaryType;
    const Int cur = get_data_type_is_Double<D>();
    //const Int cur = Is_data_double<Ta>();
    //DATA_TYPE cur = get_data_type<Ta>();
    if(cur == 0 or cur == 1){
      gauge_is_double = cur;
      //if( get_data_type_is_double<Ta >()){gauge_is_double = 1;}
      //else{gauge_is_double = 0;}
    }else{Qassert(gauge_is_double != -1);}
    gauge = gauge_;gbfac = gbfac_;gd0 = gd0_;Conj = Conj_;src_gauge = src_gauge_;
  }

  void set_bfacs(Int gbfac_, Int gd0_, bool Conj_=false, bool src_gauge_ = false)
  {Qassert(gauge!=NULL);gbfac = gbfac_;gd0 = gd0_;Conj = Conj_;src_gauge = src_gauge_;}

  template<typename Ty, bool Conj_>
  void mult_gauge(void* pt, Int dir_or);

  inline void clear_mem_dir(Int dir){
    sendbufP[dir]->resize(0);
    recvbufP[dir]->resize(0);
  }

  inline void clear_mem()
  {
    dir_cur = 0;biva = -1;civ = -1;
    for(Int i=0;i<8;i++){MPI_size[i] = 0;}
    zeroP->resize(0);
    bufsP->resize(0);
    bufrP->resize(0);

    for(Int dir=0;dir<8;dir++)
    {
      clear_mem_dir(dir);
      MPI_size[dir] = 0;
      sendbufP[dir]->resize(0);
      recvbufP[dir]->resize(0);
    }
  }

};

inline void shift_vec::init(fft_desc_basic &fds, bool GPU_set)
{
  TIMERB("Construct shift_vec");
  (void)GPU_set;
  #ifndef QLAT_USE_ACC
  GPU = false;
  #else
  GPU = GPU_set;
  #endif

  noden = fds.noden;
  rank  = fds.rank;
  Nmpi  = fds.Nmpi;
  nx=fds.nx;ny=fds.ny;nz=fds.nz;nt=fds.nt;
  vol  = fds.vol;Nvol = fds.Nvol;
  total_site = Coordinate(nx, ny, nz, nt);

  Nx=fds.Nx;Ny=fds.Ny;Nz=fds.Nz;
  Nv = fds.Nv;nv = fds.nv;

  N0 = fds.Nv[fds.orderN[0]];N1 = fds.Nv[fds.orderN[1]];N2 = fds.Nv[fds.orderN[2]];
  Nt = fds.Nt;

  zeroP = new qlat::vector_gpu<char >(0);
  bufsP = new qlat::vector_gpu<char >(0);
  bufrP = new qlat::vector_gpu<char >(0);

  init_vector_8(sendbufP, 8);
  init_vector_8(recvbufP, 8);

  MPI_size.resize(8);
  for(Int i=0;i<8;i++){MPI_size[i] = 0;}

  flag_shift_set = false;bsize = 0;
  dir_cur = 0;biva = -1;civ = -1;
  periodic = 1;

  MPI_off = 0;MPI_curr = MPI_CHAR;

  shift_set();

  gauge = NULL;
  gauge_is_double = -1;
  gbfac = 1; gd0 = 1;Conj = false;src_gauge = false;
  initialized = true;
}

inline void shift_vec::init(const Geometry& geo, bool GPU_set)
{
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  init(fd, GPU_set);
}

inline void shift_vec::init(const Coordinate& site, bool GPU_set)
{
  Geometry geo;geo.init(site);
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  init(fd, GPU_set);
}

inline void shift_vec::shift_set()
{
  TIMERB("shift_vec::shift_set");
  if(flag_shift_set){return ;}

  Geometry geo;geo.init(total_site);
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  rank_sr.resize(8);
  init_vector_8(buffoffa, 8);
  init_vector_8(buffoffb, 8);
  init_vector_8(sendoffa, 8);
  init_vector_8(sendoffb, 8);
  init_vector_8(sendoffx, 8);

  for(Int diru=0;diru<8;diru++)
  {
    Int dir = diru%4;int sign = 1;
    if(diru >= 4){sign = -1;}

    Int s0=fd.Pos0[rank][dir];int ds = sign;int Ns = Nv[dir];

    rank_sr[diru].resize(2);
    // Int count = 0;
    for(Int ranki = 0;ranki<Nmpi;ranki++)
    {
      Int flag = 1;
      for(Int idir=0;idir<4;idir++)if(idir != dir)
      {
        if(fd.Pos0[ranki][idir] != fd.Pos0[rank][idir])flag=0;
        if(fd.Pos0[ranki][idir] != fd.Pos0[rank][idir])flag=0;
        if(fd.Pos0[ranki][idir] != fd.Pos0[rank][idir])flag=0;
      }
      if(flag == 1)
      {
        if(fd.Pos0[ranki][dir] == (s0 + ds*Ns + nv[dir])%nv[dir])
        {
          rank_sr[diru][0] = ranki;
          // count += 1;
        }
        if(fd.Pos0[ranki][dir] == (s0 - ds*Ns + nv[dir])%nv[dir])
        {
          rank_sr[diru][1] = ranki;
          // count += Nmpi;
        }
      }
    }

    Length = Nt*Nx*Ny*Nz;

    buffoffa[diru]->resize(0);
    buffoffb[diru]->resize(0);
    sendoffa[diru]->resize(0);
    sendoffb[diru]->resize(0);
    sendoffx[diru]->resize(0);

    std::vector<LInt > sendoffVa;sendoffVa.resize(0);
    std::vector<LInt > sendoffVb;sendoffVb.resize(0);
    std::vector<LInt > sendoffVx;sendoffVx.resize(0);
    std::vector<LInt > buffoffVa;buffoffVa.resize(0);
    std::vector<LInt > buffoffVb;buffoffVb.resize(0);

    for(LInt off0=0;off0<Length;off0++)
    {
      Int p[4];

      p[3] = off0/(N0*N1*N2);

      p[fd.orderN[0]] = (off0/(N1*N2))%N0;
      p[fd.orderN[1]] = (off0/(N2))%N1;
      p[fd.orderN[2]] = off0%N2;

      if(rank_sr[diru][0] == rank)p[dir] = (p[dir] + ds + Nv[dir])%Nv[dir];
      if(rank_sr[diru][0] != rank)p[dir] = p[dir] + ds;

      LInt off1 = ((p[3]*N0+p[fd.orderN[0]])*N1+p[fd.orderN[1]])*N2+p[fd.orderN[2]];

      if(p[dir] >=0 and p[dir] < Nv[dir])
      {
        buffoffVa.push_back(off0);
        buffoffVb.push_back(off1);
      }
      else{
        p[dir] = (p[dir] + Nv[dir])%Nv[dir];
        off1 = ((p[3]*N0+p[fd.orderN[0]])*N1+p[fd.orderN[1]])*N2+p[fd.orderN[2]];
        sendoffVa.push_back(off0);
        sendoffVb.push_back(off1);
      }
    }

    /////May check cotinious

    sendoffa[diru]->resize(sendoffVa.size());
    sendoffb[diru]->resize(sendoffVa.size());
    sendoffx[diru]->resize(sendoffVa.size());
    #pragma omp parallel for
    for(LInt ix=0;ix<sendoffVa.size();ix++){
      (*sendoffa[diru])[ix] = sendoffVa[ix];
      (*sendoffb[diru])[ix] = sendoffVb[ix];
      (*sendoffx[diru])[ix] = ix;
    }

    buffoffa[diru]->resize(buffoffVa.size());
    buffoffb[diru]->resize(buffoffVa.size());
    #pragma omp parallel for
    for(LInt ix=0;ix<buffoffVa.size();ix++){
      (*buffoffa[diru])[ix] = buffoffVa[ix];
      (*buffoffb[diru])[ix] = buffoffVb[ix];
    }

  }
  flag_shift_set = true;
}

template<typename Ty>
void shift_vec::set_MPI_size(Int biva_or, Int civ_or, Int dir_or )
{
  TIMER("shift_vec::set_MPI_size");
  Qassert(initialized == true and flag_shift_set == true);
  if(flag_shift_set == false){shift_set();}
  Qassert(biva_or > 0 and civ_or > 0);

  /////zeroP = NULL;bufsP = NULL;bufrP = NULL;
  /////====set up bufs for shift
  //fflush_MPI();
  Qassert(Nt > 0 and N0 > 0 and N1 > 0 and N2 > 0);

  LInt Ng = Nt*N0*N1*N2;
  Qassert(zeroP != NULL and bufsP != NULL and bufrP != NULL);

  zeroP->resize(size_t(Ng) * sizeof(Ty), GPU);
  bufsP->resize(size_t(Ng)*biva_or*civ_or * sizeof(Ty), GPU);
  bufrP->resize(size_t(Ng)*biva_or*civ_or * sizeof(Ty), GPU);
  //fflush_MPI();
  //qmessage("check point1 !\n");
  //fflush_MPI();
  /////====set up bufs for shift

  ////==assign current direction
  dir_cur = dir_or;
  ////==assign current direction
  if(biva_or == biva and civ_or == civ and bsize == sizeof(Ty)){
    if(sendoffa[dir_cur]->size() == 0){
      return ;
    }else{
      if(sendbufP[dir_cur]->size()/sizeof(Ty) == (LInt) biva_or*civ_or*(sendoffa[dir_cur]->size())){
        return  ;
      }
    }
  }

  if(sizeof(Ty) != bsize){
    bsize = sizeof(Ty);
    MPI_off = sizeof(Ty);////MPI_Datatype curr = MPI_BYTE;
    unsigned int M_size = get_mpi_type<Ty >(MPI_curr );
    Qassert(MPI_off%M_size == 0);MPI_off = MPI_off/M_size;Qassert(MPI_off != 0);
  }

  ////===assign biva and civ
  biva = biva_or;
  civ  = civ_or;
  ////===assign biva and civ

  MPI_size[dir_cur] = biva*civ*(sendoffa[dir_cur]->size());
  sendbufP[dir_cur]->resize(MPI_size[dir_cur] * sizeof(Ty), GPU);
  recvbufP[dir_cur]->resize(MPI_size[dir_cur] * sizeof(Ty), GPU);
}

template<typename Ty>
void shift_vec::set_MPI_size(Int dir_or)
{
  if(civ == -1 or biva == -1 or dir_or == -1){abort_r("Need to set up dir_cur , civ and biva first. \n");}
  if(dir_or < 0 or dir_or > 8){qmessage("dir_cur size wrong %8d. \n",dir_or);abort_r();}
  if(bsize != sizeof(Ty)){abort_r("Ty type not match with previous usage.!\n");}

  set_MPI_size<Ty >(biva, civ, dir_or);
}

inline void shift_vec::print_info()
{

  qmessage("dir_curr %d,", dir_cur);
  qmessage("biva %d, civ %d, bsize %d. \n", biva, civ, bsize);
  for(Int di=0;di<8;di++)
  {
    qmessage("dir %d, bufsize %ld, MPI_size %ld, sendsize %ld, copysize %ld \n",
           di, (long)(sendbufP[di]->size() / bsize), (long)(MPI_size[di]),
           (long)sendoffa[di]->size(), (long)buffoffa[di]->size());
  }
  fflush_MPI();

}

template<typename Ty, Int flag>
void shift_vec::write_send_recv(Ty* src, Ty* res)
{
  TIMERA("shift_vec::write_send_recv");
  if(sendoffa[dir_cur]->size() != 0 and sendbufP[dir_cur]->size() == 0){
    qmessage("Memeory not set for dir %d", dir_cur);
    abort_r();
  }
  Ty* s_tem = (Ty*) (sendbufP[dir_cur]->data());
  Ty* r_tem = (Ty*) (recvbufP[dir_cur]->data());
  const LInt writeN = sendoffa[dir_cur]->size();
  if(src == NULL){abort_r("buf not defined");}
  ////print_info();
  /////Write Send buf
  /////TODO need update for GPU
  if(flag == 0 and sendoffa[dir_cur]->size() != 0)
  {
  LInt* s0 = (LInt*) qlat::get_data(*sendoffa[dir_cur]).data();
  LInt* s1 = (LInt*) qlat::get_data(*sendoffx[dir_cur]).data();

  //for(Int bi=0;bi<biva;bi++){
  //cpy_data_from_index( &s_tem[bi*writeN*civ], &src[bi*Length*civ], 
  //     s1, s0, sendoffa[dir_cur].size(), civ, GPU, QFALSE);
  //}
  cpy_data_from_index( &s_tem[0], &src[0], 
       s1, s0, sendoffa[dir_cur]->size(), civ, GPU, QTRUE, biva, writeN*civ, Length*civ);
  }

  ////Write Result
  if(flag == 1 and sendoffb[dir_cur]->size() != 0)
  {
  LInt* s1 = (LInt*) qlat::get_data(*sendoffb[dir_cur]).data();
  LInt* s0 = (LInt*) qlat::get_data(*sendoffx[dir_cur]).data();
  //for(Int bi=0;bi<biva;bi++){
  //cpy_data_from_index( &res[bi*Length*civ], &r_tem[bi*writeN*civ], 
  //      s1, s0, sendoffb[dir_cur].size(), civ, GPU, QFALSE);
  //}
  cpy_data_from_index( &res[0], &r_tem[0], 
        s1, s0, sendoffb[dir_cur]->size(), civ, GPU, QTRUE, biva, Length*civ, writeN*civ);
  }

  if(flag == 2 and buffoffa[dir_cur]->size() != 0)
  {
  LInt* s1 = (LInt*) qlat::get_data(*buffoffb[dir_cur]).data();
  LInt* s0 = (LInt*) qlat::get_data(*buffoffa[dir_cur]).data();
  //for(Int bi=0;bi<biva;bi++){
  //cpy_data_from_index( &res[bi*Length*civ], &src[bi*Length*civ], 
  //     s1, s0, buffoffa[dir_cur].size(), civ, GPU, QFALSE);
  //}
  cpy_data_from_index( &res[0], &src[0], 
       s1, s0, buffoffa[dir_cur]->size(), civ, GPU, QTRUE, biva, Length*civ, Length*civ);
  }

  s_tem = NULL; r_tem = NULL;
}

#ifdef QLAT_USE_ACC
template <typename Ty, unsigned int gs, bool Conj>
__global__ void multiply_gauge_global(Ty* a, Ty* b, const Int dir_gauge, const Int biva, const Int gbfac)
{
  __shared__ Ty ls[9];
  __shared__ Ty ds[3*gs];
  const unsigned int nt = blockDim.x * blockDim.y;
  const unsigned int tid=  threadIdx.y*blockDim.x + threadIdx.x;
  unsigned long index =blockIdx.x;
  const Int dir_limit = 4;

  unsigned int off = tid;
  while(off < 9){
    //if(!Conj)ls[(off%3)*3 + off/3] = a[(index*dir_limit*2 + dir_gauge)*9 + off];
    //if( Conj)ls[(off%3)*3 + off/3] = qlat::qconj(a[(index*dir_limit*2 + dir_gauge)*9 + off]);
    if(!Conj)ls[off] = a[(index*dir_limit*2 + dir_gauge)*9 + off];
    if( Conj)ls[off] = qlat::qconj(a[(index*dir_limit*2 + dir_gauge)*9 + off]);
    off += nt;
  }
  __syncthreads();
  for(Int bi=0;bi<biva;bi++)
  for(Int g1=0;g1<gbfac;g1++)
  {
    Ty* offB = &b[((bi*gridDim.x + index)*gbfac + g1)*3*gs];

    unsigned int off = tid;
    while(off < 3*gs){ds[off] = offB[off]; off += nt;}
    __syncthreads();

    offB[tid] = 0;
    for(Int ic=0;ic<3;ic++){offB[tid] += ls[threadIdx.y*3 + ic] * ds[ic*gs + threadIdx.x];}
    __syncthreads();
  }
}
#endif

/////gs will be ignored if cs == -1
template<typename Cy, Int gs, Int cs, bool Conj>
void multiply_gauge(void *src, void* gauge, const Int dir_gauge,const Int biva,const Long Length, const Int gbfac, const Int gd0, const bool GPU)
{
  (void)GPU;
  const Int dir_limit = 4;
  if(cs != -1){Qassert(gd0 == gs);}
  ////convention not the same as Qlattice
  ////current x, y, z,t ,-x,-y,-z,-t; 
  //////Qlat -t,-z,-y,-x, x, y, z, t
  ///to gwu convention of shift with \psi
  ////shift vec direction opposite to gwu code

  //std::vector<Int > map_dir = {3,2,1,0,  4,5,6,7};
  //const Int dir_gauge = map_dir[dir_or];
  Int fast_eigen = 0;

  #ifdef QLAT_USE_ACC
  if(GPU and cs != -1){fast_eigen = 2;}
  #else
  if(cs != -1){        fast_eigen = 1;}
  #endif
  Qassert( fast_eigen == 0 or fast_eigen == 1 or fast_eigen == 2);
  ////fast_eigen = 0;

  ////std::vector<Int > map_dir = {4,5,6,7,  3,2,1,0};
  ////const Int dir_gauge = map_dir[dir_or];
  ///qacc cannot accept struct elements or input limit of lambda functions
  ////size_t Ndata = size_t(Length) * biva * gbfac * 3 * gs;

  ////gpu fast mode
  #ifdef QLAT_USE_ACC
  if(fast_eigen == 2){
    Qassert(gs < 512);
    dim3 dimGrid( Length, 1, 1);
    dim3 dimBlock( gs, 3, 1);
    if(!Conj)multiply_gauge_global<Cy, gs, false><<< dimGrid, dimBlock >>>((Cy*) gauge, (Cy*) src, dir_gauge, biva, gbfac);
    if( Conj)multiply_gauge_global<Cy, gs, true ><<< dimGrid, dimBlock >>>((Cy*) gauge, (Cy*) src, dir_gauge, biva, gbfac);
    qacc_barrier(dummy);
  }
  #endif

  ////cpu fast mode
  if(fast_eigen == 1){
  qthread_for(index,  Long(Length), {
    QLAT_ALIGN(QLAT_ALIGNED_BYTES) Cy buf[9];
    //if(!Conj)for(Int ci=0;ci<9;ci++){buf[(ci%3)*3 + ci/3] = ((Cy*) gauge)[(index*dir_limit*2 + dir_gauge)*9 +  ci];}
    //if( Conj)for(Int ci=0;ci<9;ci++){buf[(ci%3)*3 + ci/3] = qlat::qconj(((Cy*) gauge)[(index*dir_limit*2 + dir_gauge)*9 +  ci]);}
    if(!Conj)for(Int ci=0;ci<9;ci++){buf[ci] = ((Cy*) gauge)[(index*dir_limit*2 + dir_gauge)*9 +  ci];}
    if( Conj)for(Int ci=0;ci<9;ci++){buf[ci] = qlat::qconj(((Cy*) gauge)[(index*dir_limit*2 + dir_gauge)*9 +  ci]);}
    Eigen::Matrix<Cy, 3   , 3, Eigen::ColMajor>&     lE = *((Eigen::Matrix<Cy, 3   , 3, Eigen::ColMajor>*) buf);
    /////lE * dE0 ???
    for(Int bi=0;bi<biva;bi++)
    for(Int g1=0;g1<gbfac;g1++)
    {
      Cy* d0 = &((Cy*)src)[((bi*Length + index)*gbfac + g1)*3*gs];
      if(gs != 1){
      Eigen::Matrix<Cy, gs, 3, Eigen::ColMajor>&     dEC = *((Eigen::Matrix<Cy, gs, 3, Eigen::ColMajor>*) d0);
      dEC *= lE;}
      if(gs == 1){
      Eigen::Matrix<Cy, gs, 3, Eigen::RowMajor>&     dER = *((Eigen::Matrix<Cy, gs, 3, Eigen::RowMajor>*) d0);
      dER *= lE; }
      /////something wrong with gs x 3 Col \times 3 x 3 Col with gs = 1
      /////something wrong with mix of Row and Col with gs x 3 ColM \tims 3 x 3 RowM== 1 x 3
      /////something may be good with mix of Row and Col with gs x 3 RowM \tims 3 x 3 ColwM == 1 x 3
    }
  });}

  //////gpu and cpu sloow mode
  if(fast_eigen == 0){
  Qassert(gd0 <= 128);
  qacc_for(index,  Long(Length), {
    QLAT_ALIGN(QLAT_ALIGNED_BYTES) Cy buf[9];
    QLAT_ALIGN(QLAT_ALIGNED_BYTES) Cy res[128*3];
    if(!Conj)for(Int ci=0;ci<9;ci++){buf[ci] = ((Cy*) gauge)[(index*dir_limit*2 + dir_gauge)*9 +  ci];}
    if( Conj)for(Int ci=0;ci<9;ci++){buf[ci] = qlat::qconj(((Cy*) gauge)[(index*dir_limit*2 + dir_gauge)*9 +  ci]);}
    for(Int bi=0;bi<biva;bi++)
    for(Int g1=0;g1<gbfac;g1++)
    {
      Cy* d0 = &((Cy*)src)[((bi*Length + index)*gbfac + g1)*3*gd0];
      for(Int g0=0;g0<3*gd0;g0++){res[g0] = 0;}
      for(Int c0=0;c0<3;c0++)
      {
        for(Int c1=0;c1<3;c1++){
          /////Cy tem = buf[c1*3 + c0];
          Cy tem = buf[c0*3 + c1];
          for(Int g0=0;g0<gd0;g0++){res[c0*gd0 + g0] += d0[c1*gd0 + g0] * tem;}
        }
      }
      for(Int g0=0;g0<3*gd0;g0++){d0[g0] = res[g0];}
    }
  });}
}

template<typename Ty, bool Conj_>
void shift_vec::mult_gauge(void* pt, Int dir_gauge){
  TIMERB("Gauge multiplication");
  const Int id = Is_data_double<Ty>();
  //const bool id = get_data_type_is_double<Ty >();
  ////qmessage("civ %5d, gbfac %5d, gd0 %5d \n", int(civ), int(gbfac), int(gd0));
  Qassert(id == 0 or id == 1);
  if( id){Qassert(Long(civ*sizeof(Ty)/16) == gbfac * 3 * gd0);Qassert(gauge_is_double == 1 );}
  if(!id){Qassert(Long(civ*sizeof(Ty)/8 ) == gbfac * 3 * gd0);Qassert(gauge_is_double == 0 );}
  bool cfind = false;
  #define shift_macros(bf) if(bf == gd0){cfind = true; \
    if( id)multiply_gauge<qlat::ComplexD, bf, 1, Conj_>(pt, gauge, dir_gauge, biva,Length,gbfac,gd0,GPU); \
    if(!id)multiply_gauge<qlat::ComplexF, bf, 1, Conj_>(pt, gauge, dir_gauge, biva,Length,gbfac,gd0,GPU);}
  shift_macros(1);
  shift_macros(2);
  shift_macros(3);
  shift_macros(4);
  shift_macros(6);
  shift_macros(8);
  shift_macros(16);
  shift_macros(12);
  shift_macros(4*12);
  shift_macros(4*16);
  shift_macros(12*12);

  if(!cfind){cfind = true;
    if( id)multiply_gauge<qlat::ComplexD , 1, -1, Conj_>(pt, gauge, dir_gauge, biva,Length,gbfac,gd0,GPU);
    if(!id)multiply_gauge<qlat::ComplexF, 1, -1, Conj_>(pt,  gauge, dir_gauge, biva,Length,gbfac,gd0,GPU);}
  #undef shift_macros
  Qassert(cfind);
}

template<typename Ty>
void shift_vec::call_MPI(Ty *src, Ty *res,Int dir_or)
{
  if(bsize != sizeof(Ty)){abort_r("mem set not match!\n");}
  //TIMER("MPI shift calls ");
  if(flag_shift_set == false){qmessage("Need to set up shifts. \n");abort_r();}
  //if(dir_or != -1)
  ////if(dir_or != dir_cur)set_MPI_size<Ty >(dir_or);
  set_MPI_size<Ty >(dir_or);
  if(dir_cur < 0 or dir_cur > 8){qmessage("dir_cur size wrong %8d. \n",dir_cur);abort_r();}

  /////===set src pointer for MPI
  //resP = (void*) res_or;
  //srcP = (void*) src_or;

  ////multiply link shift
  /////gauge should be V * dir_limit(4) * 2(-+) * 9 (c3x3)
  std::vector<Int > map_dir0 = {3,2,1,0,  4,5,6,7};
  std::vector<Int > map_dir1 = {7,6,5,4,  0,1,2,3};
  Int dir_gauge = map_dir0[dir_cur];
  if(src_gauge){dir_gauge = map_dir1[dir_cur];}

  /////src_gauge need test for cases, will need to change src vectors
  if(gauge != NULL and src_gauge == true){
  if(!Conj)mult_gauge<Ty, false >((void*) src, dir_gauge);
  if( Conj)mult_gauge<Ty, true  >((void*) src, dir_gauge);}

  Ty* s_tem= (Ty*) sendbufP[dir_cur]->data();
  Ty* r_tem= (Ty*) recvbufP[dir_cur]->data();

  MPI_Request send_req;
  MPI_Request recv_req;
  const Int tags = 10240 + 777;  // AMD machine MPI have tag issues with Quda ...
  const Int tagr = 10240 + 777; 
  //const Int tags = 10;
  //const Int tagr = 10;
  //int tags = omp_get_thread_num()*Nmpi + rank;
  //int tagr = omp_get_thread_num()*Nmpi + rank_sr[dir_cur][1];
  if(MPI_size[dir_cur] == 0)write_send_recv<Ty, 2 >(src, res);//Write same node
  if(MPI_size[dir_cur] != 0)
  {
    write_send_recv<Ty, 0 >(src, res);//Write send

    ////MPI_Recv((Ftype*) &recvbuf[0] ,MPI_size,CMPI,rank_sr[dir_cur][1],tagr,comm, &status);
    MPI_Isend(s_tem, MPI_size[dir_cur]*MPI_off, MPI_curr,rank_sr[dir_cur][0],tags, get_comm(), &send_req);
    MPI_Irecv(r_tem, MPI_size[dir_cur]*MPI_off, MPI_curr,rank_sr[dir_cur][1],tagr, get_comm(), &recv_req);

    write_send_recv<Ty, 2 >(src, res);//Write same node

    MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
    MPI_Wait(&send_req, MPI_STATUS_IGNORE);
    ////qmessage("SIZE! MPI_off %d, MPI_size %d \n", int(MPI_off), int(M_size));

    //MPI_Wait(&request, &status);
    //if(omp_get_thread_num()==0)MPI_Wait(&request, MPI_STATUS_IGNORE);
    //synchronize();
    write_send_recv<Ty, 1 >(src, res);//Write from recv
  }

  if(gauge != NULL and src_gauge == false){
  if(!Conj)mult_gauge<Ty, false >((void*) res, dir_gauge);
  if( Conj)mult_gauge<Ty, true  >((void*) res, dir_gauge);}

  s_tem = NULL; r_tem = NULL;
}

inline void get_periodic(Int &dx,Int nx)
{
  dx = dx%nx;
  if(std::abs(dx) > nx/2.0)
  {
    Int sign = 1;if(dx<0)sign = -1;
    Int v = std::abs(dx);
    v = nx - v;
    dx = (-1)*sign*v;
  }
  //dx = dx;
}

template<typename Ty>
void shift_vec::shift_vecs(std::vector<Ty* > &src,std::vector<Ty* > &res,std::vector<Int >& iDir, Int civ_or)
{
  /////TODO change the use of biva, civa 
  /////dividable or change inner loop to 1
  const LInt Ng = Nt*N0*N1*N2;
  Qassert(Ng == Length);
  //#if PRINT_TIMER>4
  TIMER_FLOPS("shift_Evec");
  {
    Int count = 1; for(LInt di=0;di<iDir.size();di++){count += int(std::abs(iDir[di]));}
    timer.flops += count * src.size() * Ng*civ_or*sizeof(Ty) ;
  }
  //#endif

  Int flag_abort=0;int biva_or = src.size();
  if(iDir.size()!=4){qmessage("shift directions wrong .");flag_abort=1;}
  if(biva_or <=0 or civ_or<=0){qmessage("Cannot do it with biva_or <=0 or civ_or==0");flag_abort=1;}
  if(flag_abort==1){abort_r();}

  std::vector<Int > dir_curl,dir_numl;
  dir_curl.resize(0);
  dir_numl.resize(0);
  for(Int ir=0;ir<4;ir++)
  {
    if(iDir[ir] != 0)
    {
      Int dirc = ir;
      Int curriDir = iDir[ir];
      if(periodic == 1)get_periodic(curriDir,nv[ir]);
      //if(std::abs(iDir[ir])> nv[ir]/2.0)
      Int dirn = curriDir;
      if(curriDir < 0){dirc = dirc + 4;dirn = -1*dirn;}
      //if(periodic == 1)if(dirn >= nv[ir]){dirn = nv[ir] - dirn;}
      dir_curl.push_back(dirc);
      dir_numl.push_back(dirn);
    }
  }

  if(dir_curl.size()==0){
    LInt Nsize = Nt*N0*N1*N2*civ_or;
    VectorGPUKey gkey(size_t(Nsize)*sizeof(Ty), ssprintf("shift_vec_buf"), GPU);
    qlat::vector_gpu<char >& tem = get_vector_gpu_plan<char >(gkey);
    for(LInt vi=0;vi<(LInt) biva_or;vi++)
    {
      ///in case of a simple memory vi shift
      if(src[vi] != res[vi]){
        cpy_data_thread((Ty*) tem.data(), &src[vi][0], Nsize, GPU, QTRUE);
        cpy_data_thread(&res[vi][0], (Ty*) tem.data(), Nsize, GPU, QTRUE);
      }
    }
    return ;
  }

  Int size_vec = biva_or*civ_or;
  Ty* zero = (Ty*) zeroP->data();

  //qmessage("Flag %d, civ %d %d , biva %d \n",int(flag_shift_set),civ_or,civ, biva);
  //if(flag_shift_set == false){if(civ_or==1)set_MPI_size<Ty >(1,12);if(civ_or != 1)set_MPI_size<Ty >(1, civ_or);}
  if(civ == -1 or biva == -1){
    if(civ_or == 1){  set_MPI_size<Ty >(1,     12);}
    if(civ_or != 1){set_MPI_size<Ty >(1  , civ_or);}
  }
  if(civ_or != 1){if(civ_or != civ ){
    qmessage("civor %3d, civ %3d \n", civ_or, civ);abort_r("Configuration not equal \n");
  }}

  //Qassert(biva_or == biva); // code can be grouped with internal biva, biva_or is the outter loops
  Qassert(biva_or >= biva);

  std::vector<Ty *> ptem0;ptem0.resize(civ );

  Int count = 0;
  Int flagend = 0;

  Ty* vec_s = (Ty*) bufsP->data();
  Ty* vec_r = (Ty*) bufrP->data();

  for(Int fftn=0;fftn< (size_vec+biva*civ-1)/(biva*civ);fftn++)
  {
    if(flagend == 1){break;}
    Int start = count;

    for(Int li=0;li< biva ;li++)
    {
      if(civ_or == 1)
      {
      {
        for(Int ci=0;ci<civ ;ci++)
        {
          if(count < size_vec)
          {
            ptem0[ci] = (Ty* ) &(src[count][0]);
          }
          else{ptem0[ci] = (Ty* ) &zero[0];flagend = 1;}
          count = count + 1;
        }
        for(Int ci=0;ci<civ ;ci++)
        {
          ////memcpy(&vec_s[(li*civ +ci)*Ng+0],&ptem0[ci][0],sizeof(Ty)*Ng);
          cpy_data_thread(&vec_s[(li*civ +ci)*Ng+0],&ptem0[ci][0], Ng, GPU, QFALSE);
        }
      }
      }

      if(civ_or != 1)
      {
      if(count < size_vec){
      ////memcpy(&vec_s[li*Ng*civ_or+0],&src[count/civ_or][0],sizeof(Ty)*Ng*civ_or);
      cpy_data_thread(&vec_s[li*Ng*civ_or+0],&src[count/civ_or][0], Ng*civ_or, GPU, QFALSE);
      count = count + civ_or;
      }
      if(count == size_vec)flagend = 1;
      }
      qacc_barrier(dummy);

    }
    // move from biva * civ * Ng --> biva * Ng * civ
    if(civ_or == 1 and civ != 1)mv_civ.move_civ_in(&vec_s[0], &vec_s[0], biva, civ, Ng, 1, GPU);

    /////Shift direction kernal
    for(LInt di=0;di<dir_curl.size();di++)
    {
    for(Int si=0;si<dir_numl[di];si++)
    {
      call_MPI((Ty*) &vec_s[0], (Ty*) &vec_r[0],dir_curl[di]);
      ///memcpy(&vec_s[0],&vec_r[0],sizeof(Ty)*(biva*Ng*civ));
      cpy_data_thread(&vec_s[0],&vec_r[0], (biva*Ng*civ), GPU, QFALSE);
    }
    }
    /////Shift direction kernal

    // move from biva * Ng * civ --> biva * civ * Ng
    if(civ_or == 1 and civ != 1)mv_civ.move_civ_out(&vec_s[0], &vec_s[0], biva, Ng, civ, 1, GPU);

    //TIMER("Reorder heavy data.");
    for(Int li=0;li< biva ;li++)
    {
      //int nxi = li/biv;
      //int is  = li%biv;
      if(civ_or == 1)
      {
        for(Int ci=0;ci<civ ;ci++)
        {
          if(start < size_vec)
          {
            ptem0[ci] = (Ty*) &(res[start][0]);
          }
          else{ptem0[ci] = (Ty*) &zero[0];}
          start = start + 1;
        }
        for(Int ci=0;ci<civ ;ci++)
        {
          //memcpy(&ptem0[ci][0],&vec_s[(li*civ +ci)*Ng+0],sizeof(Ty)*Ng);
          cpy_data_thread(&ptem0[ci][0],&vec_s[(li*civ +ci)*Ng+0], Ng, GPU, QFALSE);
        }
      }
      //write_in_MPIsend(ptem0,li,1);

      if(civ_or != 1)
      {
      if(start < size_vec){
        //memcpy(&res[start/civ_or][0],&vec_s[li*Ng*civ_or+0],sizeof(Ty)*Ng*civ_or);
        cpy_data_thread(&res[start/civ_or][0],&vec_s[li*Ng*civ_or+0], Ng*civ_or, GPU, QFALSE);
        start = start + civ_or;
      }
      qacc_barrier(dummy);
      if(start == size_vec)flagend = 1;
      }
    }
  }

}

template<typename Ty>
void shift_vec::shift_vecP(Ty* src, Ty* res, std::vector<Int >& iDir, Int civ_or)
{
  std::vector<Ty* > srcL(1);srcL[0] = src;
  std::vector<Ty* > resL(1);resL[0] = res;
  shift_vecs(srcL, resL, iDir, civ_or);
}

template<typename Ty>
void shift_vec::shift_Evec(std::vector<qlat::vector<Ty > > &srcE,std::vector<qlat::vector<Ty > > &srcEf,std::vector<Int >& iDir,Int civ_or)
{
  Int flag_abort=0;
  if(srcE.size()==0){qmessage("Cannot do it with srcE.size()==0");flag_abort=1;}
  if(civ_or<=0){qmessage("Cannot do it with civ_or==0");flag_abort=1;}
  if(iDir.size()!=4){qmessage("shift directions wrong .");flag_abort=1;}
  if(srcE[0].size()!=Nt*N0*N1*N2*civ_or){
    qmessage("omp num %3d \n",omp_get_thread_num());
    qmessage("Cannot do it with srcE[0].size()!=Nt*N0*N1*N2*civ_or,srcE[0].size() %6d,Nt*N0 %6d,N1 %6d,N2 %6d \n",
    int(srcE[0].size()),Nt*N0,N1,N2*civ_or);flag_abort=1;}
  if(srcEf.size()==0){qmessage("Cannot do it with srcEf.size()==0");flag_abort=1;}
  if(srcEf[0].size()!=Nt*N0*N1*N2*civ_or){
    qmessage("omp num %3d \n",omp_get_thread_num());
    qmessage("Cannot do it with srcEf[0].size()!=N0*N1*N2*civ_or,srcEf[0].size() %6d,Nt*N0 %6d,N1 %6d,N2 %6d \n",
    int(srcEf[0].size()),Nt*N0,N1,N2*civ_or);flag_abort=1;}

  if(flag_abort==1){abort_r();}

  Int biva_or = srcE.size();
  std::vector<Ty* > src;std::vector<Ty* > res;
  src.resize(biva_or);res.resize(biva_or);
  for(Int bi=0;bi<biva_or;bi++){
    src[bi] = (Ty*) qlat::get_data(srcE[bi]).data();
    res[bi] = (Ty*) qlat::get_data(srcE[bi]).data();
  }

  shift_vecs(src, res, iDir ,civ_or);

}

template<typename Ty>
void shift_vec::shift_vecs_dir(Ty* src, Ty* res, Int civ_, Int mu, Int sign)
{
  std::vector<Int > iDir(4);for(Int i=0;i<4;i++){iDir[i] = 0;}
  iDir[mu] = sign;
  shift_vecP(src, res, iDir , civ_);
}

template<typename Ty >
void shift_vec::shift_vecs_dir(qlat::Field<Ty>& src, qlat::Field<Ty>& res, const Int mu, const Int sign){
  Qassert(IsTypeComplex<Ty>());
  Qassert(get_mem_order(src) == QLAT_DEFAULT and get_mem_order(res) == QLAT_DEFAULT);
  Qassert(src.initialized and res.initialized );
  Ty* ps0 = (Ty*) qlat::get_data(src).data();
  Ty* ps1 = (Ty*) qlat::get_data(res).data();
  const Int civ_ = src.multiplicity;
  Qassert(res.multiplicity == civ_);
  std::vector<Int > iDir(4);for(Int i=0;i<4;i++){iDir[i] = 0;}
  iDir[mu] = sign;
  shift_vecP(ps0, ps1, iDir , civ_);
}

template<typename Ty, Int civ_>
void shift_vec::shift_vecs_dir(std::vector<qlat::FieldM<Ty, civ_> >& src, std::vector<qlat::FieldM<Ty, civ_> >& res, Int mu, Int sign)
{
  Qassert(IsTypeComplex<Ty>());
  Qassert(src.size() == res.size());
  Int Nsrc = src.size();if(Nsrc == 0){return ;}
  std::vector<Ty* > Psrc; std::vector<Ty* > Pres; 
  Psrc.resize(Nsrc);Pres.resize(Nsrc);
  for(Int si = 0; si < Nsrc; si++)
  {
    Qassert(src[si].initialized and res[si].initialized);
    ////Qassert(src[si].multiplicity == civ_ and res[si].multiplicity == civ_);
    Qassert(get_mem_order(src[si]) == QLAT_DEFAULT and get_mem_order(res[si]) == QLAT_DEFAULT);
    Psrc[si] = (Ty*) qlat::get_data(src[si]).data();
    Pres[si] = (Ty*) qlat::get_data(res[si]).data();
  }
  std::vector<Int > iDir(4);for(Int i=0;i<4;i++){iDir[i] = 0;}iDir[mu] = sign;
  shift_vecs(Psrc, Pres, iDir, civ_);
}

template<typename Ty>
void shift_vec::shift_vecs_dirG(std::vector<qlat::FieldG<Ty> >& src, std::vector<qlat::FieldG<Ty> >& res, Int mu, Int sign)
{
  Qassert(IsTypeComplex<Ty>());
  Qassert(src.size() == res.size());
  Int Nsrc = src.size();if(Nsrc == 0){return ;}
  std::vector<Ty* > Psrc; std::vector<Ty* > Pres; 
  Psrc.resize(Nsrc);Pres.resize(Nsrc);
  const Int civ_ = src[0].multiplicity;
  for(Int si = 0; si < Nsrc; si++)
  {
    Qassert(src[si].initialized and res[si].initialized);
    Qassert(src[si].multiplicity == civ_ and res[si].multiplicity == civ_);
    Qassert(get_mem_order(src[si]) == QLAT_DEFAULT and get_mem_order(res[si]) == QLAT_DEFAULT);
    Psrc[si] = (Ty*) qlat::get_data(src[si]).data();
    Pres[si] = (Ty*) qlat::get_data(res[si]).data();
  }
  std::vector<Int > iDir(4);for(Int i=0;i<4;i++){iDir[i] = 0;}iDir[mu] = sign;
  shift_vecs(Psrc, Pres, iDir, civ_);
}

/////  ===shift_vec buffers related
struct ShiftVecsKey {
  Coordinate total_site;
  Int GPU;
  ShiftVecsKey(const Coordinate& total_site_, Int GPU_ = 1)
  {
    total_site = total_site_;
    GPU = GPU_;
  }
};

inline bool operator<(const ShiftVecsKey& x, const ShiftVecsKey& y)
{
  if(x.total_site < y.total_site ){  return true;}
  if(y.total_site < x.total_site ){  return false;}
  if(x.GPU < y.GPU ){  return true;}
  if(y.GPU < x.GPU ){  return false;}
  return false;
}

inline Cache<ShiftVecsKey, shift_vec >& get_shift_vec_cache()
{
  static Cache<ShiftVecsKey, shift_vec > cache("ShiftVecsKey", 16);
  return cache;
}

inline shift_vec& get_shift_vec_plan(const ShiftVecsKey& fkey)
{
  if (!get_shift_vec_cache().has(fkey)) {
    get_shift_vec_cache()[fkey].init(fkey.total_site, fkey.GPU);
  }
  shift_vec& buf = get_shift_vec_cache()[fkey];
  return buf;
}

inline shift_vec& get_shift_vec_plan(const Coordinate& total_site)
{
  ShiftVecsKey fkey(total_site);
  return get_shift_vec_plan(fkey);
}
// ========

template <class Ty, Int civ>
void shift_fieldM(shift_vec& svec, std::vector<qlat::FieldM<Ty, civ> >& src, std::vector<qlat::FieldM<Ty, civ> >& res, std::vector<Int >& iDir)
{
  if(src.size() < 1)return;
  Int biva_or = src.size();

  std::vector<Ty* > srcP;std::vector<Ty* > resP;
  srcP.resize(biva_or);resP.resize(biva_or);
  for(Int bi=0;bi<biva_or;bi++){
    srcP[bi] = (Ty*) qlat::get_data(src[bi]).data();
    resP[bi] = (Ty*) qlat::get_data(res[bi]).data();
  }

  svec.shift_vecs(srcP, resP, iDir , civ);
}


template <class Ty>
void shift_fields(shift_vec& svec, Ty* src, Ty* res, std::vector<Int >& iDir, const Int biva, const Int civ)
{
  const LInt V = svec.Nvol;
  std::vector<Ty* > srcP;std::vector<Ty* > resP;
  srcP.resize(biva);resP.resize(biva);
  for(long bi=0;bi<long(biva);bi++){
    srcP[bi] = (Ty*) &src[bi * V * civ];
    resP[bi] = (Ty*) &res[bi * V * civ];
  }
  svec.shift_vecs(srcP, resP, iDir , civ);
}

template <class Td>
void shift_fieldM(shift_vec& svec, std::vector<Propagator4dT<Td >* >& src, std::vector<Propagator4dT<Td >* >& res, std::vector<Int >& iDir)
{
  if(src.size() < 1)return;
  Int biva_or = src.size();

  std::vector<ComplexT<Td>* > srcP;std::vector<ComplexT<Td>* > resP;
  srcP.resize(biva_or);resP.resize(biva_or);
  for(Int bi=0;bi<biva_or;bi++){
    srcP[bi] = (ComplexT<Td>*) qlat::get_data(*src[bi]).data();
    resP[bi] = (ComplexT<Td>*) qlat::get_data(*res[bi]).data();
  }
  svec.shift_vecs(srcP, resP, iDir , 12*12);
}

template <class Td>
void shift_fieldM(shift_vec& svec, Propagator4dT<Td >& src, Propagator4dT<Td >& res, std::vector<Int >& iDir)
{
  std::vector<Propagator4dT<Td >* >srcP(1);srcP[0] = &src;
  std::vector<Propagator4dT<Td >* >resP(1);resP[0] = &res;
  shift_fieldM(svec, srcP, resP, iDir);
}

template <class Td>
void shift_fieldM(shift_vec& svec, std::vector<Propagator4dT<Td > > & src, std::vector< Propagator4dT<Td > >& res, std::vector<Int >& iDir)
{
  std::vector<Propagator4dT<Td >* >srcP(0);srcP.resize(src.size());
  std::vector<Propagator4dT<Td >* >resP(0);resP.resize(res.size());
  for(unsigned int i=0;i<src.size();i++){srcP[i] = &src[i];}
  for(unsigned int i=0;i<res.size();i++){resP[i] = &res[i];}
  shift_fieldM(svec, srcP, resP, iDir);
}

template <class Td>
void symmetric_shift(shift_vec& svec, std::vector<Propagator4dT<Td > >& src, std::vector< Propagator4dT<Td > >& res, std::vector< Propagator4dT<Td > >& buf, Int idir)
{
  if(src.size() == 0){res.resize(0);return ;}
  std::vector<Int > iDir(4);for(Int i=0;i<4;i++){iDir[i] = 0;}

  const Geometry& geo = src[0].geo();
  const Long Nvol = geo.local_volume();

  if(res.size() != src.size()){res.resize(src.size());}
  if(buf.size() != src.size()){buf.resize(src.size());}
  qlat::vector<qlat::ComplexT<Td>* > srcP;srcP.resize(src.size());
  qlat::vector<qlat::ComplexT<Td>* > resP;resP.resize(src.size());
  qlat::vector<qlat::ComplexT<Td>* > bufP;bufP.resize(src.size());
  for(unsigned int i=0;i<res.size();i++){
    if(!res[i].initialized){res[i].init(src[0].geo());}
    if(!buf[i].initialized){buf[i].init(src[0].geo());}
    srcP[i] = (qlat::ComplexT<Td>*) qlat::get_data(src[i]).data();
    resP[i] = (qlat::ComplexT<Td>*) qlat::get_data(res[i]).data();
    bufP[i] = (qlat::ComplexT<Td>*) qlat::get_data(buf[i]).data();
  }

  //for(unsigned int i=0;i<res.size();i++){
  //  qacc_for(isp, Nvol*12*12, {resP[i][isp]  = srcP[i][isp];});
  //}

  iDir[idir] = +1;
  shift_fieldM(svec, src, buf, iDir);
  for(unsigned int i=0;i<res.size();i++){
    qacc_for(isp, Nvol*12*12, {resP[i][isp]  = 0.5 * bufP[i][isp];});
  }

  iDir[idir] = -1;
  shift_fieldM(svec, src, buf, iDir);

  for(unsigned int i=0;i<res.size();i++){
    qacc_for(isp, Nvol*12*12, {resP[i][isp] -= 0.5 * bufP[i][isp];});
  }

}


///////src should be different than res
template<typename Ty>
void shift_vecs_dir_qpropT(std::vector<qlat::FieldM<Ty, 12*12> >& src, std::vector<qlat::FieldM<Ty, 12*12> >& res, Int mu, Int sign, shift_vec& svec)
{
  if(src.size() == 0){res.resize(0); return; }

  Qassert(src.size() == res.size());
  for(unsigned int i=0;i<src.size();i++)
  {
    qprop_move_dc_in(src[i]);
  }

  for(unsigned int iv=0;iv<src.size();iv++)
  {
    if(!res[iv].initialized){res[iv].init(src[0].geo());}
  }

  svec.shift_vecs_dir(src, res, mu, sign);

  for(unsigned int i=0;i<src.size();i++)
  {
    qprop_move_dc_out(res[i]);
    qprop_move_dc_out(src[i]);
  }
}

/*
  covariant shifts to 4 directions
*/
template<typename Ty>
void shift_vecs_cov_fieldG(std::vector< std::vector<FieldG<Ty> > >& res, std::vector< FieldG<Ty> >& s0, shift_vec& svec,
  std::vector<std::vector<FieldG<Ty > >>& buf)
{
  if(res.size() != 5){res.resize(5);}
  if(buf.size() != 2){buf.resize(2);}
  Qassert(s0.size() != 0 and s0[0].initialized);
  const Long Nsrc = s0.size();
  const bool reorder = get_mem_order(s0[0]) == QLAT_DEFAULT ? false : true;

  if(reorder)
  for(Long si=0;si<Nsrc;si++){
    switch_orders(s0[si], QLAT_DEFAULT);
  }

  init_fieldsG(buf[0], s0[0], s0.size());
  init_fieldsG(buf[1], s0[0], s0.size());
  for(unsigned int i=0;i<res.size();i++)
  {
    init_fieldsG(res[i], s0[0], s0.size());
  }

  fields_operations(res[0], res[0], s0, Ty(0.0, 0.0), Ty(0.0, 0.0), Ty(1.0, 0.0));// equal

  for(Int nu = 0; nu < 4 ; nu++)
  {
    svec.shift_vecs_dirG(s0, buf[0], nu, +1);
    svec.shift_vecs_dirG(s0, buf[1], nu, -1);
    // r = 0.5 * (b0 - b1)
    fields_operations(res[1 + nu], buf[0], buf[1], Ty(0.0, 0.0), Ty(0.5,0.0), Ty(-0.5, 0.0));
  }

  //swap src back to orginal layout and res to correct layout
  if(reorder){
    for(Long si=0;si<Nsrc;si++){
      switch_orders(s0[si], QLAT_OUTTER);
    }
    for(unsigned int i=0;i<res.size();i++)
    {
      for(Long si=0;si<Nsrc;si++){
        switch_orders(res[i][si], QLAT_OUTTER);
      }
    }
  }

}


template <class Ty, class Ta>
void shift_fields_qlat(Ty* src, Ta* res, const std::vector<Int >& iDir, const Int Nvec, const Geometry& geo, const Int move_in = 1)
{
  TIMER("shift_fields_qlat");
  //Qassert(Nvec <= Ngroup);
  const LInt V = geo.local_volume();
  Coordinate shift;
  for(Int i=0;i<4;i++){shift[i] = 1.0 * iDir[i];}

  Field<MvectorT<1, Ty>> f0;
  Field<MvectorT<1, Ty>> f1;
  f0.init(geo, Nvec);
  f1.init(geo, Nvec);
  move_index mv_civ;

  Ty* p0 = (Ty* ) qlat::get_data(f0).data();
  Ty* p1 = (Ty* ) qlat::get_data(f1).data();

  cpy_GPU(p0, src, V * Nvec);
  if(move_in == 1){mv_civ.move_civ_in(p0, p0, 1, Nvec, V, 1, true);}

  field_shift_direct(f1, f0, shift);

  if(move_in == 1){mv_civ.move_civ_out(p1, p1, 1, V, Nvec, 1, true);}
  cpy_GPU(res, p1, V * Nvec);
}

template <class Ty>
void shift_fields_gridPT(Ty** src, Ty** res, const std::vector<Int >& iDir, const Int biva, const Int civ, const Geometry& geo, const Int mode = 1)
{
  TIMER("shift_fields_grid");
  const Coordinate  total_site = geo.total_site();
  const LInt V = geo.local_volume();
  Qassert(geo.is_only_local);
  /*
    (void)mode;
    Fast mode to avoid too many shifts and copies
    Not related to field grid shifts...
  */
  if(mode == 1){
    const Int GPU = 1;
    std::vector<Field<Ty > > bs;
    std::vector<Field<Ty > > br;
    bs.resize(biva);
    br.resize(biva);
    for(unsigned int iv=0;iv<bs.size();iv++){
      set_field(bs[iv], src[iv], V * civ, geo);
      set_field(br[iv], res[iv], V * civ, geo);
    }

    const size_t Nd = size_t(V) * civ * sizeof(Ty);
    VectorGPUKey gkey(0, ssprintf("shift_vec_buf"), GPU);
    qlat::vector_gpu<char >& tem = get_vector_gpu_plan<char >(gkey);

    tem.resizeL(2 * biva * Nd);
    std::vector<vector<Ty> > to_bufL;to_bufL.resize(2 * biva);
    for(size_t i=0;i<to_bufL.size();i++){
      Vector<Ty> v( (Ty*) &tem[i * Nd], Nd / sizeof(Ty) );
      to_bufL[i].set_mem_type(MemType::Acc);
      to_bufL[i].set_view(v);
    }

    Coordinate shift;for(Int i=0;i<4;i++){shift[i] = 1.0 * iDir[i];}
    field_shift_directT(br, bs, shift, to_bufL);
    return ;
  }

  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  const Coordinate  local_site(fd.Nx, fd.Ny, fd.Nz, fd.Nt);

  for(Int bi=0;bi<biva;bi++){
    if(res[bi] != src[bi]){
      cpy_GPU(res[bi], src[bi],  V * civ);
    }
  }

  qlat::vector<Ty* > sP;sP.resize(biva);
  std::vector<Ty* > sPd;sPd.resize(biva);
  for(Int bi=0;bi<biva;bi++){
    sP[bi]  = res[bi];
    sPd[bi] = res[bi];
  }

  // three step scale
  std::vector<Coordinate > total_siteL;
  total_siteL.resize(3);

  /*
    0 to be full block small lattice
    1 3/4 Nx shifts
    2 max volume
  */
  for(Int s=0;s<3;s++){
    for(Int i=0;i<4;i++){
      if(s == 0 ){total_siteL[s][i] = total_site[i] / local_site[i];}
      if(s == 1 ){
        total_siteL[s][i] = total_site[i];
        if(local_site[i] % 2 == 0){total_siteL[s][i] = total_site[i] / 2;}
        if(local_site[i] % 3 == 0){total_siteL[s][i] = total_site[i] / 3;}
        if(local_site[i] % 4 == 0){total_siteL[s][i] = total_site[i] / 4;}
      }
      if(s == 2 ){total_siteL[s][i] = total_site[i];}
    }
  }

  //for(unsigned int s=0;s<total_siteL.size();s++){
  //  qmessage("s %d : ", s);
  //  for(Int i=0;i<4;i++){qmessage(" %d ", int(total_siteL[s][i]));}
  //  qmessage("\n");
  //}

  std::vector<Int > d0(4);
  std::vector<Int > d1(4);
  d0 = iDir;d1 = iDir;

  // from coase to finner
  Coordinate lat0 = total_site ;
  Coordinate lat1 = total_site ;

  // lat0 current geo; lat1 shift geo
  for(unsigned int s=0;s<total_siteL.size();s++)
  {
    lat1 = total_siteL[s];
    LInt small_vol = 1;
    for(Int i=0;i<4;i++){
      Int fac =  total_site[i] / lat1[i] ;
      Qassert(fac >= 1);
      d0[i] = d1[i] % fac;
      d1[i] = d1[i] / fac;
      small_vol *= fac;
    }

    bool shift=false;
    for(Int i=0;i<4;i++){
      if(d1[i] != 0){
        shift = true;
        break;
      }
    }
    // no shift, move to next without reshape
    if(shift == false){
      d1 = d0;
      continue;
    }
    //qmessage("biva %d small vol %d, civ %d, shift %3d %3d %3d %3d \n", 
    //  int(biva), int(small_vol), int(civ), d0[0], d0[1], d0[2], d0[3]);

    grid_memory_reshape(sP, sP, civ, lat1, lat0, total_site);

    shift_vec& svec_m = get_shift_vec_plan(lat1);
    svec_m.set_MPI_size<Ty >(biva, small_vol * civ);

    //void shift_fields(shift_vec& svec, Ty* src, Ty* res, std::vector<Int >& iDir, const Int biva, const Int civ)
    svec_m.shift_vecs(sPd, sPd, d1, small_vol * civ);
    //shift_fields(svec_m, res, res, d1, biva, small_vol * civ);

    d1 = d0;
    lat0 = lat1;
  }
  for(Int i=0;i<4;i++){Qassert(d1[i] == 0);}

  // reorder to ori layout if needed
  if(lat0 != total_site){
    lat1 = total_site;
    grid_memory_reshape(sP, sP, civ, lat1, lat0, total_site);
  }
}

template <class Ty>
void shift_fields_gridP(Ty** src, Ty** res, const std::vector<Int >& iDir, const Int biva, const Int civ, const Geometry& geo, const Int max_group = -1)
{
  Int max_biva = biva;
  if(max_group != -1 and max_group < biva){
    max_biva = max_group;
  }

  std::vector<Long > jobA = job_create( biva, max_biva);
  for(LInt jobi=0;jobi < jobA.size()/2; jobi++){
    const Long bini = jobA[jobi*2+0];
    const Long bcut = jobA[jobi*2+1];
    shift_fields_gridPT(&src[bini], &res[bini], iDir, bcut, civ, geo);
  }
}

template <class Ty>
void shift_fields_grid(Ty* src, Ty* res, const std::vector<Int >& iDir, const Int biva, const Int civ, const Geometry& geo, const Int max_group = -1)
{
  vector<Ty* > sP;
  vector<Ty* > rP;
  sP.resize(biva);
  rP.resize(biva);
  const LInt V = geo.local_volume();
  for(Int bi=0;bi<biva;bi++){
    sP[bi] = &src[bi * V * civ];
    rP[bi] = &res[bi * V * civ];
  }
  shift_fields_gridP(sP.data(), rP.data(), iDir, biva, civ, geo, max_group);
}

/*
  shift of FieldG 
    QLAT_OUTTER  civ = 1
    QLAT_DEFAULT civ = multiplicity
    sign = -1, shift to origin
*/
template <class Ty>
void shift_fields_grid(std::vector<FieldG<Ty > >& src, std::vector<FieldG<Ty > >& res, const Coordinate& sp, const Int sign = 1, const Int max_group = -1)
{
  if(src.size() == 0){res.resize(0);return ;}
  std::vector<Int > iDir;iDir.resize(4);
  for(Int i=0;i<4;i++){iDir[i] = sign * sp[i];}

  const Long Nsrc = src.size();
  const Geometry& geo = src[0].geo();
  const Long V = geo.local_volume();
  Int biva = 0;
  Int civ  = 0;
  vector<Ty* > sP;
  vector<Ty* > rP;

  for(Long si=0;si<Nsrc;si++){
    Qassert(src[si].initialized);
    Qassert(res[si].initialized);
    Ty* stmp = (Ty*) get_data(src[si]).data();
    Ty* rtmp = (Ty*) get_data(res[si]).data();
    if(src[si].mem_order == QLAT_OUTTER){
      if(civ != 0){Qassert(civ == 1 and biva == Nsrc * src[si].multiplicity);}
      else{
        civ = 1;
        biva = Nsrc * src[si].multiplicity;
        sP.resize(biva);
        rP.resize(biva);
      }
      for(Int ci=0;ci<src[si].multiplicity;ci++){
        const Int bi = si*src[si].multiplicity + ci;
        sP[bi] = &stmp[ci * V];
        rP[bi] = &rtmp[ci * V];
      }
    }
    if(src[si].mem_order == QLAT_DEFAULT){
      if(civ != 0){Qassert(civ == src[si].multiplicity and biva == Nsrc);}
      else{
        civ = src[si].multiplicity;
        biva = Nsrc;
        sP.resize(biva);
        rP.resize(biva);
      }
      sP[si] = &stmp[0];
      rP[si] = &rtmp[0];
    }
  }
  shift_fields_gridP(sP.data(), rP.data(), iDir, biva, civ, geo, max_group);
}

inline void clear_shift_plan_cache()
{
  get_shift_vec_cache().clear();
}



}

#endif
