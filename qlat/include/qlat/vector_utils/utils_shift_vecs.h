// utils_FFT_GPU.h
// Gen Wang
// Sep. 2021

#ifndef UTILS_SHIFT_VECS_H
#define UTILS_SHIFT_VECS_H
#pragma once

#include "general_funs.h"
#include "utils_fft_desc.h"
#include "utils_copy_data.h"

//////TODO
//////GPU support, Ty template, shift of fieldM

namespace qlat
{

struct shift_vec{
  int rank;int Nmpi;
  int nx,ny,nz,nt;
  LInt noden;
  LInt vol;
  LInt Nvol;
  int Nx,Ny,Nz;

  int N0,N1,N2,Nt;
  qlat::vector_acc<int> Nv,nv;

  fft_desc_basic fd;

  //////Shift under periodic condition
  int periodic;

  //std::vector<std::vector<int> > sendlist;
  //std::vector<std::vector<int> > recvlist;

  int civ,biva;
  LInt Length;
  bool flag_shift_set;
  bool GPU;

  int dir_cur;
  std::vector<std::vector<int > > rank_sr;

  //std::vector<std::vector<int > > sendoffa;
  //std::vector<std::vector<int > > sendoffb;
  //std::vector<std::vector<int > > sendoffx;
  //std::vector<std::vector<int > > buffoffa;
  //std::vector<std::vector<int > > buffoffb;

  std::vector<qlat::vector_acc<LInt > > sendoffa;
  std::vector<qlat::vector_acc<LInt > > sendoffb;
  std::vector<qlat::vector_acc<LInt > > sendoffx;
  std::vector<qlat::vector_acc<LInt > > buffoffa;
  std::vector<qlat::vector_acc<LInt > > buffoffb;
  std::vector<unsigned int > CON_SEND, CON_BUFF;

  std::vector<size_t > MPI_size;
  unsigned int MPI_off;
  MPI_Datatype MPI_curr;


  //void* srcP;
  //void* resP;
  unsigned int bsize;

  move_index mv_civ;

  void* zeroP;LInt zeroP_Size;
  void* bufsP;LInt bufsP_Size;
  void* bufrP;LInt bufrP_Size;

  std::vector<void* > sendbufP;
  std::vector<void* > recvbufP;

  std::vector<LInt > bufP_Size;
  ////std::vector<size_t > recvbufP_Size;

  //std::vector<Ftype > sendbuf;
  //std::vector<Ftype > recvbuf;

  //std::vector<int > currsend;
  //std::vector<int > currspls;
  //std::vector<int > currrecv;
  //std::vector<int > currrpls;

  void* gauge;
  int gbfac;int gd0;
  bool Conj;bool src_gauge;

  shift_vec(fft_desc_basic &fds, bool GPU_set = true);
  void print_info();
  ~shift_vec();

  void shift_set();

  template<typename Ty>
  void set_MPI_size(int biva_or, int civ_or, int dir_or = 0);
  template<typename Ty>
  void set_MPI_size(int dir_or);

  template<typename Ty, int flag>
  void write_send_recv(Ty* src, Ty* res);

  template<typename Ty>
  void call_MPI(Ty *src, Ty *res,int dir_or);

  template<typename Ty>
  void shift_Evec(std::vector<qlat::vector_acc<Ty > > &srcE,std::vector<qlat::vector_acc<Ty > > &srcEf,std::vector<int >& iDir,int civ_or);

  template<typename Ty>
  void shift_vecs(std::vector<Ty* > &src,std::vector<Ty* > &res,std::vector<int >& iDir ,int civ_or);

  template<typename Ty>
  void shift_vecP(Ty* src,Ty* res,std::vector<int >& iDir ,int civ_or);

  template<typename Ty, int civ_>
  void shift_vecs_dir(qlat::FieldM<Ty , civ_>& src, qlat::FieldM<Ty , civ_>& res, int mu, int sign);

  template<typename Ty, int civ_>
  void shift_vecs_dir(std::vector<qlat::FieldM<Ty , civ_> >& src, std::vector<qlat::FieldM<Ty , civ_> >& res, int mu, int sign);


  void set_gauge(void* gauge_, int gbfac_, int gd0_, bool Conj_=false, bool src_gauge_ = false)
  {gauge = gauge_;gbfac = gbfac_;gd0 = gd0_;Conj = Conj_;src_gauge = src_gauge_;}
  template<typename Ty, bool Conj_>
  void mult_gauge(void* pt, int dir_or);

  void clear_mem_dir(int dir);
  ////void clear_mem();

};

shift_vec::shift_vec(fft_desc_basic &fds, bool GPU_set)
{
  TIMERB("Construct shift_vec");
  (void)GPU_set;
  #ifndef QLAT_USE_ACC
  GPU = false;
  #else
  GPU = GPU_set;
  #endif

  copy_fft_desc(fd, fds);
  ////fd = &fds;

  noden = fds.noden;
  rank  = fds.rank;
  Nmpi  = fds.Nmpi;
  nx=fds.nx;ny=fds.ny;nz=fds.nz;nt=fds.nt;
  vol  = fds.vol;Nvol = fds.Nvol;

  Nx=fds.Nx;Ny=fds.Ny;Nz=fds.Nz;
  //mx=fds.mx;my=fds.my;mz=fds.mz;mt=fds.mt;
  Nv = fds.Nv;nv = fds.nv;
  //mv = fds.mv;
  //orderN = fds.orderN;

  N0 = fds.Nv[fds.orderN[0]];N1 = fds.Nv[fds.orderN[1]];N2 = fds.Nv[fds.orderN[2]];
  Nt = fds.Nt;

  sendbufP.resize(8);recvbufP.resize(8);bufP_Size.resize(8);MPI_size.resize(8);
  CON_SEND.resize(8);CON_BUFF.resize(8);
  for(int i=0;i<8;i++){
    sendbufP[i] = NULL;recvbufP[i] = NULL;bufP_Size[i] = 0;MPI_size[i] = 0;
    CON_SEND[i] = 1;CON_BUFF[i] = 1;
  }

  flag_shift_set = false;bsize = 0;
  dir_cur = 0;biva = -1;civ = -1;
  periodic = 1;

  MPI_off = 0;MPI_curr = MPI_CHAR;
  ///srcP = NULL;resP = NULL;

  shift_set();

  zeroP = NULL;bufsP = NULL;bufrP = NULL;
  zeroP_Size = 0;bufsP_Size = 0;bufrP_Size = 0;

  gauge = NULL;
  gbfac = 1; gd0 = 1;Conj = false;src_gauge = false;

}

void shift_vec::shift_set()
{
  TIMERB("shift_vec::shift_set");
  if(flag_shift_set){return ;}

  rank_sr.resize(8);
  buffoffa.resize(8);
  buffoffb.resize(8);
  sendoffa.resize(8);
  sendoffb.resize(8);
  sendoffx.resize(8);

  ////MPI_size.resize(8);

  for(int diru=0;diru<8;diru++)
  {
    int dir = diru%4;int sign = 1;
    if(diru >= 4){sign = -1;}

    int s0=fd.Pos0[rank][dir];int ds = sign;int Ns = Nv[dir];

    rank_sr[diru].resize(2);
    int count = 0;
    for(int ranki = 0;ranki<Nmpi;ranki++)
    {
      int flag = 1;
      for(int idir=0;idir<4;idir++)if(idir != dir)
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
          count += 1;
        }
        if(fd.Pos0[ranki][dir] == (s0 - ds*Ns + nv[dir])%nv[dir])
        {
          rank_sr[diru][1] = ranki;
          count += Nmpi;
        }
      }
    }
    //double flag = 0.0;
    ///if(count != (1+Nmpi)){print0("rank %8d. map wrong %8d. \n",rank,count);flag=1.0;}abort_sum(flag);

    Length = Nt*Nx*Ny*Nz;

    //buffoffa[diru].resize(0);
    //buffoffb[diru].resize(0);
    //sendoffa[diru].resize(0);
    //sendoffb[diru].resize(0);
    //sendoffx[diru].resize(0);

    sendoffa[diru].resize(0);
    sendoffb[diru].resize(0);
    sendoffx[diru].resize(0);
    buffoffa[diru].resize(0);
    buffoffb[diru].resize(0);

    std::vector<LInt > sendoffVa;sendoffVa.resize(0);
    std::vector<LInt > sendoffVb;sendoffVb.resize(0);
    std::vector<LInt > sendoffVx;sendoffVx.resize(0);
    std::vector<LInt > buffoffVa;buffoffVa.resize(0);
    std::vector<LInt > buffoffVb;buffoffVb.resize(0);

    for(LInt off0=0;off0<Length;off0++)
    {
      int p[4];

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
        //buffoffa[diru].push_back(off0);
        //buffoffb[diru].push_back(off1);

      }
      else{
        p[dir] = (p[dir] + Nv[dir])%Nv[dir];
        off1 = ((p[3]*N0+p[fd.orderN[0]])*N1+p[fd.orderN[1]])*N2+p[fd.orderN[2]];
        sendoffVa.push_back(off0);
        sendoffVb.push_back(off1);
        //sendoffa[diru].push_back(off0);
        //sendoffb[diru].push_back(off1);

      }
    }

    //for(LInt ix=0;ix<sendoffa[diru].size();ix++){sendoffx[diru].push_back(ix);}

    //CON_SEND = 1; CON_BUFF = 1;
    /////May check cotinious

    sendoffa[diru].resize(sendoffVa.size());
    sendoffb[diru].resize(sendoffVa.size());
    sendoffx[diru].resize(sendoffVa.size());
    #pragma omp parallel for
    for(LInt ix=0;ix<sendoffVa.size();ix++){
      sendoffa[diru][ix] = sendoffVa[ix];
      sendoffb[diru][ix] = sendoffVb[ix];
      sendoffx[diru][ix] = ix;
    }

    buffoffa[diru].resize(buffoffVa.size());
    buffoffb[diru].resize(buffoffVa.size());
    #pragma omp parallel for
    for(LInt ix=0;ix<buffoffVa.size();ix++){
      buffoffa[diru][ix] = buffoffVa[ix];
      buffoffb[diru][ix] = buffoffVb[ix];
    }


    ///reorder_civ((char*) &sendoff[diru][0], (char*) &sendoff[diru][0], 1, 2, sendoff[diru].size()/2, 1, sizeof(sendoff[diru][0]));
    /////reorder_civ((char*) &buffoff[diru][0], (char*) &buffoff[diru][0], 1, 2, buffoff[diru].size(), 1, sizeof(buffoff[diru][0]));

  }
  flag_shift_set = true;
}

template<typename Ty>
void shift_vec::set_MPI_size(int biva_or, int civ_or, int dir_or )
{
  if(flag_shift_set == false){shift_set();}

  /////zeroP = NULL;bufsP = NULL;bufrP = NULL;
  /////====set up bufs for shift
  LInt Ng = Nt*N0*N1*N2;
  if(zeroP_Size != Ng){
  free_buf(zeroP, GPU);zeroP=NULL;
  if(GPU){gpuMalloc(zeroP, Ng, Ty);}
  else{zeroP = aligned_alloc_no_acc(Ng * sizeof(Ty));}zeroP_Size = Ng;
  ////Ty* zero = (Ty*) zeroP;
  zero_Ty((Ty*) zeroP, Ng, GPU, true);
  }

  if(bufsP_Size != Ng*biva_or*civ_or){
  free_buf(bufsP, GPU);bufsP=NULL;
  if(GPU){gpuMalloc(bufsP, Ng*biva_or*civ_or, Ty);}
  else{bufsP = aligned_alloc_no_acc(Ng*biva_or*civ_or * sizeof(Ty));}bufsP_Size = Ng*biva_or*civ_or;}

  if(bufrP_Size != Ng*biva_or*civ_or){
  free_buf(bufrP, GPU);bufrP=NULL;
  if(GPU){gpuMalloc(bufrP, Ng*biva_or*civ_or, Ty);}
  else{bufrP = aligned_alloc_no_acc(Ng*biva_or*civ_or * sizeof(Ty));}bufrP_Size = Ng*biva_or*civ_or;}
  /////====set up bufs for shift

  ////==assign current direction
  dir_cur = dir_or;
  ////==assign current direction
  if(biva_or == biva and civ_or == civ and bsize == sizeof(Ty)){
    if(sendoffa[dir_cur].size() == 0){return ;}else
    {if(bufP_Size[dir_cur] == (LInt) biva_or*civ_or*(sendoffa[dir_cur].size())*CON_SEND[dir_cur]){return  ;}}
  }

  //Length = Nt*N0*N1*N2;
  //flag_set_MPI_size = 1;

  if(sizeof(Ty) != bsize){
    bsize = sizeof(Ty);
    MPI_off = sizeof(Ty);////MPI_Datatype curr = MPI_BYTE;
    unsigned int M_size = get_MPI_type<Ty >(MPI_curr );
    qassert(MPI_off%M_size == 0);MPI_off = MPI_off/M_size;qassert(MPI_off != 0);
  }

  ////===assign biva and civ
  biva = biva_or;
  civ  = civ_or;
  ////===assign biva and civ

  MPI_size[dir_cur] = biva*civ*(sendoffa[dir_cur].size())*CON_SEND[dir_cur];
  if(bufP_Size[dir_cur] != MPI_size[dir_cur])
  {
    clear_mem_dir(dir_cur);
    if(GPU){gpuMalloc(sendbufP[dir_cur], MPI_size[dir_cur], Ty);}
    else{sendbufP[dir_cur] = aligned_alloc_no_acc(MPI_size[dir_cur] * sizeof(Ty));}

    if(GPU){gpuMalloc(recvbufP[dir_cur], MPI_size[dir_cur], Ty);}
    else{recvbufP[dir_cur] = aligned_alloc_no_acc(MPI_size[dir_cur] * sizeof(Ty));}

    bufP_Size[dir_cur] = MPI_size[dir_cur];
  }
}

template<typename Ty>
void shift_vec::set_MPI_size(int dir_or)
{
  if(civ == -1 or biva == -1 or dir_or == -1){abort_r("Need to set up dir_cur , civ and biva first. \n");}
  if(dir_or < 0 or dir_or > 8){print0("dir_cur size wrong %8d. \n",dir_or);abort_r();}
  if(bsize != sizeof(Ty)){abort_r("Ty type not match with previous usage.!\n");}

  set_MPI_size<Ty >(biva, civ, dir_or);
}

void shift_vec::print_info()
{

  print0("dir_curr %d,", dir_cur);
  print0("biva %d, civ %d, bsize %d. \n", biva, civ, bsize);
  for(int di=0;di<8;di++)
  {
    print0("dir %d, bufsize %ld, MPI_size %ld, sendsize %ld, copysize %ld \n", 
            di, long(bufP_Size[di]), long(MPI_size[di]), sendoffa[di].size(), buffoffa[di].size());
  }
  fflush_MPI();

}

template<typename Ty, int flag>
void shift_vec::write_send_recv(Ty* src, Ty* res)
{
  TIMERA("shift_vec::write_send_recv");
  if(sendoffa[dir_cur].size() != 0 and bufP_Size[dir_cur] == 0){print0("Memeory not set for dir %d", dir_cur);abort_r();}
  Ty* s_tem = (Ty*) (sendbufP[dir_cur]);
  Ty* r_tem = (Ty*) (recvbufP[dir_cur]);
  LInt writeN = sendoffa[dir_cur].size();
  if(src == NULL){abort_r("buf not defined");}
  ////print_info();
  /////Write Send buf
  /////TODO need update for GPU
  if(flag == 0)
  {
  for(int bi=0;bi<biva;bi++){
  LInt* s0 = (LInt*) qlat::get_data(sendoffa[dir_cur]).data();
  LInt* s1 = (LInt*) qlat::get_data(sendoffx[dir_cur]).data();
  cpy_data_from_index( &s_tem[bi*writeN*civ], &src[bi*Length*civ], 
       s1, s0, sendoffa[dir_cur].size(), CON_SEND[dir_cur]*civ, GPU, false);
  }
  qacc_barrier(dummy);

  }

  ////Write Result
  if(flag == 1)
  {
  for(int bi=0;bi<biva;bi++){
  LInt* s1 = (LInt*) qlat::get_data(sendoffb[dir_cur]).data();
  LInt* s0 = (LInt*) qlat::get_data(sendoffx[dir_cur]).data();
  cpy_data_from_index( &res[bi*Length*civ], &r_tem[bi*writeN*civ], 
        s1, s0, sendoffb[dir_cur].size(), CON_SEND[dir_cur]*civ, GPU, false);
  }
  qacc_barrier(dummy);

  }

  if(flag == 2)
  {
  for(int bi=0;bi<biva;bi++){
  LInt* s1 = (LInt*) qlat::get_data(buffoffb[dir_cur]).data();
  LInt* s0 = (LInt*) qlat::get_data(buffoffa[dir_cur]).data();
  cpy_data_from_index( &res[bi*Length*civ], &src[bi*Length*civ], 
       s1, s0, buffoffa[dir_cur].size(), CON_BUFF[dir_cur]*civ, GPU, false);
  }
  qacc_barrier(dummy);
  }

  s_tem = NULL; r_tem = NULL;
}

#ifdef QLAT_USE_ACC
template <typename Ty, unsigned int gs, bool Conj>
__global__ void multiply_gauge_global(Ty* a, Ty* b, const int dir_gauge, const int biva, const int gbfac)
{
  __shared__ Ty ls[9];
  __shared__ Ty ds[3*gs];
  const unsigned int nt = blockDim.x * blockDim.y;
  const unsigned int tid=  threadIdx.y*blockDim.x + threadIdx.x;
  unsigned long index =blockIdx.x;
  const int dir_limit = 4;

  unsigned int off = tid;
  while(off < 9){
    //if(!Conj)ls[(off%3)*3 + off/3] = a[(index*dir_limit*2 + dir_gauge)*9 + off];
    //if( Conj)ls[(off%3)*3 + off/3] = qlat::qconj(a[(index*dir_limit*2 + dir_gauge)*9 + off]);
    if(!Conj)ls[off] = a[(index*dir_limit*2 + dir_gauge)*9 + off];
    if( Conj)ls[off] = qlat::qconj(a[(index*dir_limit*2 + dir_gauge)*9 + off]);
    off += nt;
  }
  __syncthreads();
  for(int bi=0;bi<biva;bi++)
  for(int g1=0;g1<gbfac;g1++)
  {
    Ty* offB = &b[((bi*blockDim.x + index)*gbfac + g1)*3*gs];

    unsigned int off = tid;
    while(off < 3*gs){ds[off] = offB[off]; off += nt;}
    __syncthreads();

    offB[tid] = 0;
    for(int ic=0;ic<3;ic++){offB[tid] += ls[threadIdx.y*3 + ic] * ds[ic*gs + threadIdx.x];}
  }
}
#endif

/////gs will be ignored if cs == -1
template<typename Cy, int gs, int cs, bool Conj>
void multiply_gauge(void *src, void* gauge, const int dir_gauge,const int biva,const long Length, const int gbfac, const int gd0, const bool GPU)
{
  (void)GPU;
  const int dir_limit = 4;
  if(cs != -1){qassert(gd0 == gs);}
  ////convention not the same as Qlattice
  ////current x, y, z,t ,-x,-y,-z,-t; 
  //////Qlat -t,-z,-y,-x, x, y, z, t
  ///to gwu convention of shift with \psi
  ////shift vec direction opposite to gwu code

  //std::vector<int > map_dir = {3,2,1,0,  4,5,6,7};
  //const int dir_gauge = map_dir[dir_or];
  int fast_eigen = 0;

  #ifdef QLAT_USE_ACC
  if(GPU and cs != -1){fast_eigen = 2;}
  #else
  if(cs != -1){        fast_eigen = 1;}
  #endif
  qassert( fast_eigen == 0 or fast_eigen == 1 or fast_eigen == 2);
  ////fast_eigen = 0;

  ////std::vector<int > map_dir = {4,5,6,7,  3,2,1,0};
  ////const int dir_gauge = map_dir[dir_or];
  ///qacc cannot accept struct elements or input limit of lambda functions
  ////size_t Ndata = size_t(Length) * biva * gbfac * 3 * gs;

  ////gpu fast mode
  #ifdef QLAT_USE_ACC
  if(fast_eigen == 2){
    qassert(gs < 512);
    dim3 dimGrid( Length, 1, 1);
    dim3 dimBlock( gs, 3, 1);
    if(!Conj)multiply_gauge_global<Cy, gs, false><<< dimGrid, dimBlock >>>((Cy*) gauge, (Cy*) src, dir_gauge, biva, gbfac);
    if( Conj)multiply_gauge_global<Cy, gs, true ><<< dimGrid, dimBlock >>>((Cy*) gauge, (Cy*) src, dir_gauge, biva, gbfac);
    qacc_barrier(dummy);
  }
  #endif

  ////cpu fast mode
  if(fast_eigen == 1){
  qthread_for(index,  long(Length), {
    Cy buf[9];
    //if(!Conj)for(int ci=0;ci<9;ci++){buf[(ci%3)*3 + ci/3] = ((Cy*) gauge)[(index*dir_limit*2 + dir_gauge)*9 +  ci];}
    //if( Conj)for(int ci=0;ci<9;ci++){buf[(ci%3)*3 + ci/3] = qlat::qconj(((Cy*) gauge)[(index*dir_limit*2 + dir_gauge)*9 +  ci]);}
    if(!Conj)for(int ci=0;ci<9;ci++){buf[ci] = ((Cy*) gauge)[(index*dir_limit*2 + dir_gauge)*9 +  ci];}
    if( Conj)for(int ci=0;ci<9;ci++){buf[ci] = qlat::qconj(((Cy*) gauge)[(index*dir_limit*2 + dir_gauge)*9 +  ci]);}
    Eigen::Matrix<Cy, 3   , 3, Eigen::ColMajor>&     lE = *((Eigen::Matrix<Cy, 3   , 3, Eigen::ColMajor>*) buf);
    /////lE * dE0 ???
    for(int bi=0;bi<biva;bi++)
    for(int g1=0;g1<gbfac;g1++)
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
  qassert(gd0 <= 128);
  qacc_for(index,  long(Length), {
    Cy buf[9];
    Cy res[128*3];
    if(!Conj)for(int ci=0;ci<9;ci++){buf[ci] = ((Cy*) gauge)[(index*dir_limit*2 + dir_gauge)*9 +  ci];}
    if( Conj)for(int ci=0;ci<9;ci++){buf[ci] = qlat::qconj(((Cy*) gauge)[(index*dir_limit*2 + dir_gauge)*9 +  ci]);}
    for(int bi=0;bi<biva;bi++)
    for(int g1=0;g1<gbfac;g1++)
    {
      Cy* d0 = &((Cy*)src)[((bi*Length + index)*gbfac + g1)*3*gd0];
      for(int g0=0;g0<3*gd0;g0++){res[g0] = 0;}
      for(int c0=0;c0<3;c0++)
      {
        for(int c1=0;c1<3;c1++){
          /////Cy tem = buf[c1*3 + c0];
          Cy tem = buf[c0*3 + c1];
          for(int g0=0;g0<gd0;g0++){res[c0*gd0 + g0] += d0[c1*gd0 + g0] * tem;}
        }
      }
      for(int g0=0;g0<3*gd0;g0++){d0[g0] = res[g0];}
    }
  });}
}

template<typename Ty, bool Conj_>
void shift_vec::mult_gauge(void* pt, int dir_gauge){
  TIMERB("Gauge multiplication");
  bool id = get_data_type_is_double<Ty >();
  if( id){qassert(long(civ*sizeof(Ty)/16) == gbfac * 3 * gd0);}
  if(!id){qassert(long(civ*sizeof(Ty)/8 ) == gbfac * 3 * gd0);}
  bool cfind = false;
  #define shift_macros(bf) if(bf == gd0){cfind = true; \
    if( id)multiply_gauge<qlat::Complex , bf, 1, Conj_>(pt, gauge, dir_gauge, biva,Length,gbfac,gd0,GPU); \
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
    if( id)multiply_gauge<qlat::Complex , 1, -1, Conj_>(pt, gauge, dir_gauge, biva,Length,gbfac,gd0,GPU);
    if(!id)multiply_gauge<qlat::ComplexF, 1, -1, Conj_>(pt, gauge, dir_gauge, biva,Length,gbfac,gd0,GPU);}
  #undef shift_macros
  qassert(cfind);
}

template<typename Ty>
void shift_vec::call_MPI(Ty *src, Ty *res,int dir_or)
{
  if(bsize != sizeof(Ty)){abort_r("mem set not match!\n");}
  //TIMER("MPI shift calls ");
  if(flag_shift_set == false){print0("Need to set up shifts. \n");abort_r();}
  //if(dir_or != -1)
  ////if(dir_or != dir_cur)set_MPI_size<Ty >(dir_or);
  set_MPI_size<Ty >(dir_or);
  if(dir_cur < 0 or dir_cur > 8){print0("dir_cur size wrong %8d. \n",dir_cur);abort_r();}

  /////===set src pointer for MPI
  //resP = (void*) res_or;
  //srcP = (void*) src_or;

  ////multiply link shift
  /////gauge should be V * dir_limit(4) * 2(-+) * 9 (c3x3)
  std::vector<int > map_dir0 = {3,2,1,0,  4,5,6,7};
  std::vector<int > map_dir1 = {7,6,5,4,  0,1,2,3};
  int dir_gauge = map_dir0[dir_cur];
  if(src_gauge){dir_gauge = map_dir1[dir_cur];}

  /////src_gauge need test for cases, will need to change src vectors
  if(gauge != NULL and src_gauge == true){
  if(!Conj)mult_gauge<Ty, false >((void*) src, dir_gauge);
  if( Conj)mult_gauge<Ty, true  >((void*) src, dir_gauge);}

  Ty* s_tem= (Ty*) sendbufP[dir_cur]; Ty* r_tem= (Ty*) recvbufP[dir_cur];

  MPI_Request request;
  int tags = omp_get_thread_num()*Nmpi + rank;
  int tagr = omp_get_thread_num()*Nmpi + rank_sr[dir_cur][1];
  if(MPI_size[dir_cur] == 0)write_send_recv<Ty, 2 >(src, res);//Write same node
  if(MPI_size[dir_cur] != 0)
  {
    write_send_recv<Ty, 0 >(src, res);//Write send

    ////MPI_Recv((Ftype*) &recvbuf[0] ,MPI_size,CMPI,rank_sr[dir_cur][1],tagr,comm, &status);
    MPI_Isend(s_tem, MPI_size[dir_cur]*MPI_off, MPI_curr,rank_sr[dir_cur][0],tags, get_comm(), &request);
    MPI_Recv( r_tem, MPI_size[dir_cur]*MPI_off, MPI_curr,rank_sr[dir_cur][1],tagr, get_comm(), MPI_STATUS_IGNORE);

    write_send_recv<Ty, 2 >(src, res);//Write same node

    MPI_Wait(&request, MPI_STATUS_IGNORE);
    ////print0("SIZE! MPI_off %d, MPI_size %d \n", int(MPI_off), int(M_size));

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

void get_periodic(int &dx,int nx)
{
  dx = dx%nx;
  if(std::abs(dx) > nx/2.0)
  {
    int sign = 1;if(dx<0)sign = -1;
    int v = std::abs(dx);
    v = nx - v;
    dx = (-1)*sign*v;
  }
  //dx = dx;
}

template<typename Ty>
void shift_vec::shift_vecs(std::vector<Ty* > &src,std::vector<Ty* > &res,std::vector<int >& iDir, int civ_or)
{
  /////TODO change the use of biva, civa 
  /////dividable or change inner loop to 1
  LInt Ng = Nt*N0*N1*N2;
  //#if PRINT_TIMER>4
  TIMER_FLOPS("shift_Evec");
  {
    int count = 1; for(LInt di=0;di<iDir.size();di++){count += int(std::abs(iDir[di]));}
    timer.flops += count * src.size() * Ng*civ_or*sizeof(Ty) ;
  }
  //#endif

  int flag_abort=0;int biva_or = src.size();
  if(iDir.size()!=4){print0("shift directions wrong .");flag_abort=1;}
  if(biva_or <=0 or civ_or<=0){print0("Cannot do it with biva_or <=0 or civ_or==0");flag_abort=1;}
  if(flag_abort==1){abort_r();}

  std::vector<int > dir_curl,dir_numl;
  dir_curl.resize(0);
  dir_numl.resize(0);
  for(int ir=0;ir<4;ir++)
  {
    if(iDir[ir] != 0)
    {
      int dirc = ir;
      int curriDir = iDir[ir];
      if(periodic == 1)get_periodic(curriDir,nv[ir]);
      //if(std::abs(iDir[ir])> nv[ir]/2.0)
      int dirn = curriDir;
      if(curriDir < 0){dirc = dirc + 4;dirn = -1*dirn;}
      //if(periodic == 1)if(dirn >= nv[ir]){dirn = nv[ir] - dirn;}
      dir_curl.push_back(dirc);
      dir_numl.push_back(dirn);
    }
  }

  if(dir_curl.size()==0){
    LInt Nsize = Nt*N0*N1*N2*civ_or;
    std::vector<Ty > tem;tem.resize(Nsize);
    for(LInt vi=0;vi<(LInt) biva_or;vi++)
    {
      //memcpy(&tem[0], &src[vi][0],sizeof(Ty )*Nsize);
      //memcpy(&res[vi][0] ,&tem[0],sizeof(Ty )*Nsize);
      if(src[vi] != res[vi]){
      cpy_data_thread(&tem[0], &src[vi][0], Nsize, GPU, false);
      cpy_data_thread(&res[vi][0] ,&tem[0], Nsize, GPU, true );}
    }
    return;
  }

  int size_vec = biva_or*civ_or;
  Ty* zero = (Ty*) zeroP;

  //print0("Flag %d, civ %d %d , biva %d \n",int(flag_shift_set),civ_or,civ, biva);
  //if(flag_shift_set == false){if(civ_or==1)set_MPI_size<Ty >(1,12);if(civ_or != 1)set_MPI_size<Ty >(1, civ_or);}
  if(civ == -1 or biva == -1){if(civ_or==1){set_MPI_size<Ty >(1,12);}if(civ_or != 1){set_MPI_size<Ty >(1, civ_or);}}
  if(civ_or != 1){if(civ_or != civ ){abort_r("Configuration not equal \n");}}

  std::vector<Ty *> ptem0;ptem0.resize(civ );

  int count = 0;
  int flagend = 0;

  Ty* vec_s = (Ty*) bufsP;
  Ty* vec_r = (Ty*) bufrP;

  for(int fftn=0;fftn< (size_vec+biva*civ-1)/(biva*civ);fftn++)
  {
    if(flagend == 1){break;}
    int start = count;

    for(int li=0;li< biva ;li++)
    {
      if(civ_or == 1)
      {
      {
        for(int ci=0;ci<civ ;ci++)
        {
          if(count < size_vec)
          {
            ptem0[ci] = (Ty* ) &(src[count][0]);
          }
          else{ptem0[ci] = (Ty* ) &zero[0];flagend = 1;}
          count = count + 1;
        }
        for(int ci=0;ci<civ ;ci++)
        {
          ////memcpy(&vec_s[(li*civ +ci)*Ng+0],&ptem0[ci][0],sizeof(Ty)*Ng);
          cpy_data_thread(&vec_s[(li*civ +ci)*Ng+0],&ptem0[ci][0], Ng, GPU, false);
        }
      }
      }

      if(civ_or != 1)
      {
      if(count < size_vec){
      ////memcpy(&vec_s[li*Ng*civ_or+0],&src[count/civ_or][0],sizeof(Ty)*Ng*civ_or);
      cpy_data_thread(&vec_s[li*Ng*civ_or+0],&src[count/civ_or][0], Ng*civ_or, GPU, false);
      count = count + civ_or;
      }
      if(count == size_vec)flagend = 1;
      }
      qacc_barrier(dummy);

    }
    ////if(civ_or == 1)reorder_civ((char*) &vec_s[0],(char*) &vec_s[0], civ ,biva,Ng,0, sizeof(Ty));
    if(civ_or == 1)mv_civ.dojob(&vec_s[0], &vec_s[0],civ, biva, Ng, 0, 1, GPU);

    /////Shift direction kernal
    for(LInt di=0;di<dir_curl.size();di++)
    {
    for(int si=0;si<dir_numl[di];si++)
    {
      call_MPI((Ty*) &vec_s[0], (Ty*) &vec_r[0],dir_curl[di]);
      ///memcpy(&vec_s[0],&vec_r[0],sizeof(Ty)*(biva*Ng*civ));
      cpy_data_thread(&vec_s[0],&vec_r[0], (biva*Ng*civ), GPU, true);
    }
    }
    /////Shift direction kernal

    ///if(civ_or == 1)reorder_civ((char*) &vec_s[0],(char*) &vec_s[0],civ ,biva,Ng,1, sizeof(Ty));
    if(civ_or == 1)mv_civ.dojob(&vec_s[0], &vec_s[0],civ, biva, Ng,1, 1, GPU);

    //TIMER("Reorder heavy data.");
    for(int li=0;li< biva ;li++)
    {
      //int nxi = li/biv;
      //int is  = li%biv;
      if(civ_or == 1)
      {
        for(int ci=0;ci<civ ;ci++)
        {
          if(start < size_vec)
          {
            ptem0[ci] = (Ty*) &(res[start][0]);
          }
          else{ptem0[ci] = (Ty*) &zero[0];}
          start = start + 1;
        }
        for(int ci=0;ci<civ ;ci++)
        {
          //memcpy(&ptem0[ci][0],&vec_s[(li*civ +ci)*Ng+0],sizeof(Ty)*Ng);
          cpy_data_thread(&ptem0[ci][0],&vec_s[(li*civ +ci)*Ng+0], Ng, GPU, false);
        }
      }
      //write_in_MPIsend(ptem0,li,1);

      if(civ_or != 1)
      {
      if(start < size_vec){
        //memcpy(&res[start/civ_or][0],&vec_s[li*Ng*civ_or+0],sizeof(Ty)*Ng*civ_or);
        cpy_data_thread(&res[start/civ_or][0],&vec_s[li*Ng*civ_or+0], Ng*civ_or, GPU, false);
        start = start + civ_or;
      }
      qacc_barrier(dummy);
      if(start == size_vec)flagend = 1;
      }
    }
  }

}

template<typename Ty>
void shift_vec::shift_vecP(Ty* src, Ty* res, std::vector<int >& iDir, int civ_or)
{
  std::vector<Ty* > srcL(1);srcL[0] = src;
  std::vector<Ty* > resL(1);resL[0] = res;
  shift_vecs(srcL, resL, iDir, civ_or);
}

template<typename Ty>
void shift_vec::shift_Evec(std::vector<qlat::vector_acc<Ty > > &srcE,std::vector<qlat::vector_acc<Ty > > &srcEf,std::vector<int >& iDir,int civ_or)
{
  //////#if PRINT_TIMER>4
  //TIMER_FLOPS("shift_Evec");
  //{
  //int count = 1; for(int di=0;di<iDir.size();di++){count += int(std::abs(iDir[di]));}
  //timer.flops += count * srcE.size() * get_data(srcE[0]).data_size() ;
  //}
  //////#endif

  int flag_abort=0;
  if(srcE.size()==0){print0("Cannot do it with srcE.size()==0");flag_abort=1;}
  if(civ_or<=0){print0("Cannot do it with civ_or==0");flag_abort=1;}
  if(iDir.size()!=4){print0("shift directions wrong .");flag_abort=1;}
  if(srcE[0].size()!=Nt*N0*N1*N2*civ_or){
    print0("omp num %3d \n",omp_get_thread_num());
    print0("Cannot do it with srcE[0].size()!=Nt*N0*N1*N2*civ_or,srcE[0].size() %6d,Nt*N0 %6d,N1 %6d,N2 %6d \n",
    int(srcE[0].size()),Nt*N0,N1,N2*civ_or);flag_abort=1;}
  if(srcEf.size()==0){print0("Cannot do it with srcEf.size()==0");flag_abort=1;}
  if(srcEf[0].size()!=Nt*N0*N1*N2*civ_or){
    print0("omp num %3d \n",omp_get_thread_num());
    print0("Cannot do it with srcEf[0].size()!=N0*N1*N2*civ_or,srcEf[0].size() %6d,Nt*N0 %6d,N1 %6d,N2 %6d \n",
    int(srcEf[0].size()),Nt*N0,N1,N2*civ_or);flag_abort=1;}

  if(flag_abort==1){abort_r();}

  int biva_or = srcE.size();
  std::vector<Ty* > src;std::vector<Ty* > res;
  src.resize(biva_or);res.resize(biva_or);
  for(int bi=0;bi<biva_or;bi++){
    src[bi] = (Ty*) qlat::get_data(srcE[bi]).data();
    res[bi] = (Ty*) qlat::get_data(srcE[bi]).data();
  }

  shift_vecs(src, res, iDir ,civ_or);

}

template<typename Ty, const int civ_ >
void shift_vec::shift_vecs_dir(qlat::FieldM<Ty, civ_ >& src, qlat::FieldM<Ty, civ_ >& res, const int mu, const int sign){
  qassert(src.initialized and res.initialized );
  Ty* ps0 = (Ty*) qlat::get_data(src).data();
  Ty* ps1 = (Ty*) qlat::get_data(res).data();
  std::vector<int > iDir(4);for(int i=0;i<4;i++){iDir[i] = 0;}
  iDir[mu] = sign;
  shift_vecP(ps0, ps1, iDir , civ_);
}

template<typename Ty, int civ_>
void shift_vec::shift_vecs_dir(std::vector<qlat::FieldM<Ty , civ_> >& src, std::vector<qlat::FieldM<Ty , civ_> >& res, int mu, int sign)
{
  qassert(src.size() == res.size());
  int Nsrc = src.size();if(Nsrc == 0){return ;}
  std::vector<Ty* > Psrc; std::vector<Ty* > Pres; 
  Psrc.resize(Nsrc);Pres.resize(Nsrc);
  for(int si = 0; si < Nsrc; si++)
  {
    qassert(src[si].initialized and res[si].initialized);
    Psrc[si] = (Ty*) qlat::get_data(src[si]).data();
    Pres[si] = (Ty*) qlat::get_data(res[si]).data();
  }
  std::vector<int > iDir(4);for(int i=0;i<4;i++){iDir[i] = 0;}iDir[mu] = sign;
  shift_vecs(Psrc, Pres, iDir, civ_);
}

void shift_vec::clear_mem_dir(int dir){
  free_buf(sendbufP[dir],GPU);sendbufP[dir]=NULL;
  free_buf(recvbufP[dir],GPU);recvbufP[dir]=NULL;
  bufP_Size[dir] = 0;
  //MPI_size[dir_cur] = 0;
}

shift_vec::~shift_vec(){
  flag_shift_set = false;bsize = 0;
  dir_cur = 0;civ = -1;biva = -1;
  periodic = 1;

  for(int dir=0;dir<8;dir++){clear_mem_dir(dir);}
  MPI_size.resize(0);

  free_buf(zeroP, GPU);free_buf(bufsP, GPU);free_buf(bufrP, GPU);
  zeroP_Size = 0;bufsP_Size = 0;bufrP_Size = 0;
  zeroP = NULL;bufsP = NULL;bufrP = NULL;

  rank_sr.resize(0);
  buffoffa.resize(0);
  buffoffb.resize(0);
  sendoffa.resize(0);
  sendoffb.resize(0);
  sendoffx.resize(0);

}

//template <class Ty>
//void shift_fieldM(shift_vec& svec, Ty* src, Ty* res, const int dir, const int civ)
//{
//  if(src.size() < 1)return;
//  int biva_or = src.size();
//  std::vector<int > iDir(4);
//  for(int i=0;i<4;i++){iDir[i] = 0;}
//  if(dir >= 0){iDir[  dir ] = 1;}
//  if(dir <  0){iDir[-dir-1] = -1;}
//
//  long Nvol = svec.Nvol;
//  qassert(src.size() == res.size());
//  qassert(src.size()%Nvol == 0);
//
//  std::vector<Ty* > srcP;std::vector<Ty* > resP;
//  srcP.resize(biva_or);resP.resize(biva_or);
//  for(int bi=0;bi<biva_or;bi++){
//    srcP[bi] = (Ty*) src.qlat::get_data(src[bi]).data();
//    resP[bi] = (Ty*) src.qlat::get_data(res[bi]).data();
//  }
//
//  svec.shift_vecs(srcP, resP, iDir , civ);
//}


template <class Ty, int civ>
void shift_fieldM(shift_vec& svec, std::vector<qlat::FieldM<Ty, civ> >& src, std::vector<qlat::FieldM<Ty, civ> >& res, std::vector<int >& iDir)
{
  if(src.size() < 1)return;
  int biva_or = src.size();

  std::vector<Ty* > srcP;std::vector<Ty* > resP;
  srcP.resize(biva_or);resP.resize(biva_or);
  for(int bi=0;bi<biva_or;bi++){
    srcP[bi] = (Ty*) qlat::get_data(src[bi]).data();
    resP[bi] = (Ty*) qlat::get_data(res[bi]).data();
  }

  svec.shift_vecs(srcP, resP, iDir , civ);
}

template <class Ty>
void shift_fieldM(shift_vec& svec, std::vector<Propagator4dT<Ty >* >& src, std::vector<Propagator4dT<Ty >* >& res, std::vector<int >& iDir)
{
  if(src.size() < 1)return;
  int biva_or = src.size();

  std::vector<Ty* > srcP;std::vector<Ty* > resP;
  srcP.resize(biva_or);resP.resize(biva_or);
  for(int bi=0;bi<biva_or;bi++){
    srcP[bi] = (Ty*) qlat::get_data(*src[bi]).data();
    resP[bi] = (Ty*) qlat::get_data(*res[bi]).data();
  }
  svec.shift_vecs(srcP, resP, iDir , 12*12);
}

template <class Ty>
void shift_fieldM(shift_vec& svec, Propagator4dT<Ty > & src, Propagator4dT<Ty >& res, std::vector<int >& iDir)
{
  std::vector<Propagator4dT<Ty >* >srcP(1);srcP[0] = &src;
  std::vector<Propagator4dT<Ty >* >resP(1);resP[0] = &res;
  shift_fieldM(svec, srcP, resP, iDir);
}

template <class Ty>
void shift_fieldM(shift_vec& svec, std::vector<Propagator4dT<Ty > > & src, std::vector< Propagator4dT<Ty > >& res, std::vector<int >& iDir)
{
  std::vector<Propagator4dT<Ty >* >srcP(0);srcP.resize(src.size());
  std::vector<Propagator4dT<Ty >* >resP(0);resP.resize(res.size());
  for(int i=0;i<src.size();i++){srcP[i] = &src[i];}
  for(int i=0;i<res.size();i++){resP[i] = &res[i];}
  shift_fieldM(svec, srcP, resP, iDir);
}

}

#endif
