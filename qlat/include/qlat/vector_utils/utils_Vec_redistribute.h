// utils_VEC_redistribute.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_VEC_REDISTRIBUTE
#define UTILS_VEC_REDISTRIBUTE

#pragma once
#include "utils_copy_data.h"
#include "general_funs.h"
#include "utils_fft_desc.h"

namespace qlat{


////////each MPI NT will have all nx,ny,nz
struct Vec_redistribute
{
  long noden;int rank;int Nmpi;
  int nx,ny,nz,nt,vol,Nvol;
  int Nx,Ny,Nz,Nt;
  int mx,my,mz,mt;
  bool GPU;

  qlat::vector_acc<int> Nv,nv,mv;

  qlat::vector_acc<int > orderN;
  fft_desc_basic *fd;


  int b0,civa;
  void* sendV; void* recvV;
  int tem_off;
  bool update_off;
  /////void* bufV ;

  std::vector<int > currsend;
  std::vector<int > currspls;
  std::vector<int > currrecv;
  std::vector<int > currrpls;

  std::vector<int > sendM;
  std::vector<int > splsM;
  std::vector<int > recvM;
  std::vector<int > rplsM;

  qlat::vector_acc<LInt > map_order;
  qlat::vector_acc<LInt > map_Dorder;

  /////May need to change?
  std::vector<int > secT;

  inline void set_mem(int b0_or,int civa_or);
  int flag_set_mem;

  template<typename Ty>
  void call_MPI(int flag);

  bool copy_same_node;

  template<typename Ty>
  void re_order_recv(int flag);

  template<typename Ty>
  void reorder(Ty *sendbuf,Ty *recvbuf,int b0_or,int civa_or,int flag= 0);
  ////void reorder(Ftype *sendbuf,Ftype *recvbuf,int b0_or,int civa_or,qlat::vector<int > secT_or,int flag= 0);

  std::vector<int > map_mpi_vec;

  /////std::vector<MPI_Comm > vec_comm_list;
  MPI_Comm vec_comm;
  int mode_MPI;

  //int flag_set_fft;

  ////int flag_set_fft_force;
  ////void set_fft(int civau=-1);
  ////void do_fft(bool fftdir=true);

  Vec_redistribute(fft_desc_basic &fds, bool GPU_set = true);
  ~Vec_redistribute();

};

inline Vec_redistribute::Vec_redistribute(fft_desc_basic &fds, bool GPU_set)
{
  TIMERA("Construct Vec_redistribute");
  (void)GPU_set;
  fd = &fds;
  //GPU = GPU_set;
  #ifndef QLAT_USE_ACC
  GPU = false;
  #else
  GPU = GPU_set;
  #endif

  noden = fds.noden;
  rank  = fds.rank;
  Nmpi  = fds.Nmpi;
  nx=fds.nx;ny=fds.ny;nz=fds.nz;nt=fds.nt;
  //nt=fds.nt;
  vol  = fds.vol;Nvol = fds.Nvol;

  Nx=fds.Nx;Ny=fds.Ny;Nz=fds.Nz;Nt=fds.Nt;
  mx=fds.mx;my=fds.my;mz=fds.mz;mt=fds.mt;

  Nv = fds.Nv;nv = fds.nv;
  //mv = fds.mv;
  orderN = fds.orderN;
  /////Pos0 = fds.Pos0;
  /////mi_list = fds.mi_list;
  /////map_order.set_acc(true);map_Dorder.set_acc(true);
  tem_off = -1;

  b0 =-1; civa =-1;
  copy_same_node = true;

  //////Needed for three point functions, need copy from dev
  secT.resize(fd->Nmpi);
  for(LInt i=0;i<secT.size();i++){secT[i] = fd->Nt;}

  //////if(secT_or.size()!=fd->Nmpi){print0("secT wrong %8d ! \n", int(secT_or.size()));qassert(false);}
  ////secT = secT_or;
  ////Check same number secT for MPI
  ////long mvol = mx*my*mz;
  std::vector<int > check;check.resize(0);
  for(int n=0;n<fd->Nmpi;n++)
  {
    int init = fd->Pos0[n][3];
    int ins = secT[n];
    for(int m=0;m<fd->Nmpi;m++)
    {
    if(fd->Pos0[m][3] == init)
    {
    int ns   = secT[m];
    if(ins != ns){print0("Not match, n %8d, m %8d, ins %8d, ns %8d \n",n,m,ins,ns);qassert(false);}
    }
    }
  }
  ////Check same number secT for MPI

  map_mpi_vec.resize(Nmpi);
  for(LInt mapi=0;mapi<map_mpi_vec.size();mapi++){map_mpi_vec[mapi] = 0;}

  //for(int icomm=0;icomm<vec_comm_list.size();icomm++)
  //{MPI_Comm_free(&vec_comm_list[icomm]);}
  ///vec_comm_list.resize(0);

  ////MPI_Comm vec_comm;
  int color_xyz = fd->init;
  MPI_Comm_split(get_comm() ,color_xyz,rank,&vec_comm);
  {
    int int_tem = -1;
    MPI_Comm_rank(vec_comm,&int_tem);
    map_mpi_vec[fd->rank] = int_tem;
  }
  sum_all_size((int*) (&map_mpi_vec[0]),Nmpi);

  flag_set_mem = 0;
  mode_MPI     = 1;

  //flag_set_fft = 0;
  //howmany = 0;

  /////flag_set_fft_force = 0;

}

inline void Vec_redistribute::set_mem(int b0_or,int civa_or)
{
  TIMERA("Vec redistribute set mem");

  b0 = b0_or;civa = civa_or;

  if(b0<=0 or civa<=0){
    print0("Need at least biv = 1, civ %6d !\n",civa);qassert(false);
  }

  /////set up the map on fd
  //fd->set_up_map();
  ///map i --> nz*ny*nx   to recieve position
  /////2 is the fatest
  std::vector<long > mapcur_Vtoi;
  mapcur_Vtoi.resize(nx*ny*nz/Nv[orderN[2]]);
  int Nts = secT[fd->rank];
  for(int tmi=0;tmi<mt;tmi++)
  {
    int t0 = fd->Pos0[fd->mi_list[tmi][0]][3];
    if(t0 == fd->Pos0[rank][3])
    for(LInt temmi=0;temmi<fd->mi_list[tmi].size();temmi++)
    {
      //int off = Nx*Ny*Nz*Nt*temmi;
      //for(int iv=0;iv<Nx*Ny*Nz*Nt;iv++)
      //#pragma omp parallel for
      for(LInt iv=0;iv< ((LInt) Nx)*Ny*Nz;iv++)
      {
        //LInt Vi = fd->maporder_Nitoipos[fd->mi_list[tmi][temmi]][iv];
        int rank_tem = fd->mi_list[tmi][temmi];
        LInt vi_tem = fd->index_g_from_local(iv, rank_tem);
        LInt Vi = vi_tem - fd->Pos0[rank_tem][3] * fd->vol;

        int umi = temmi;
        LInt offnode = Nx*Ny*Nz*Nts*umi;
        LInt temiv = iv%(Nv[0]*Nv[1]*Nv[2]);
        LInt inode = b0*offnode + temiv;
        if(Vi%Nv[orderN[2]] == 0)
        if(Vi/(nv[0]*nv[1]*nv[2]) == 0)
        {
          LInt temVi = Vi%(nv[0]*nv[1]*nv[2]);
          mapcur_Vtoi[temVi/Nv[orderN[2]]]    = inode/Nv[orderN[2]];
        }
      }
    }
  }


  size_t svol = fd->nx*fd->ny*fd->nz;
  size_t LoopN = svol/Nv[orderN[2]];

  map_order.resize(b0*Nts * LoopN);
  map_Dorder.resize(b0*Nts * LoopN);
  qthread_for(iv, long(b0*Nts * LoopN), {
    int bi   =  iv/(Nts*LoopN);
    int ti   = (iv%(Nts*LoopN))/LoopN;
    size_t i =  iv%(LoopN);
    size_t Aoff = (bi*Nts + ti)*svol;
    size_t Boff = (bi*Nts + ti)*Nv[0]*Nv[1]*Nv[2];
    map_order[iv]  = (Aoff +      i        *Nv[orderN[2]])/Nv[orderN[2]]; 
    map_Dorder[iv] = (Boff + mapcur_Vtoi[i]*Nv[orderN[2]])/Nv[orderN[2]]; 
  });


  currsend.resize(Nmpi   );
  currrecv.resize(Nmpi   );
  currspls.resize(Nmpi/mt);
  currrpls.resize(Nmpi/mt);

  sendM.resize(Nmpi   );
  recvM.resize(Nmpi   );
  splsM.resize(Nmpi/mt);
  rplsM.resize(Nmpi/mt);

  for(int n=0;n<Nmpi   ;n++)currsend[n] = 0;
  for(int n=0;n<Nmpi   ;n++)currrecv[n] = 0;
  for(int n=0;n<Nmpi/mt;n++)currspls[n] = 0;
  for(int n=0;n<Nmpi/mt;n++)currrpls[n] = 0;


  for(int tmi=0;tmi<mt;tmi++) //Loop of initial time on each node
  {
    int t0 = fd->Pos0[fd->mi_list[tmi][0]][3]; //Determine the initial time
    if(t0 == fd->Pos0[rank][3]) // If current node initial time is the working initial time
    {
      for(LInt temmi=0;temmi<fd->mi_list[tmi].size();temmi++) //Determine the communication nodes of this node
      {
        //int umi = mi_list[tmi][temmi];
        int umi = temmi;

        LInt rankglobal = fd->mi_list[tmi][umi];
        LInt ranklocal = map_mpi_vec[rankglobal];

        ////no 2, 2 is inside b0
        currsend[ranklocal] = b0*Nx*Ny*Nz*secT[fd->rank]*civa;
        currspls[ranklocal] = b0*Nx*Ny*Nz*secT[fd->rank]*umi*civa;
        currrecv[ranklocal] = b0*Nx*Ny*Nz*secT[fd->rank]*civa;
        currrpls[ranklocal] = b0*Nx*Ny*Nz*secT[fd->rank]*umi*civa;
      }
    }
  }

  //////update offset from data type and new ios
  update_off   = true;
  flag_set_mem = 1;
}

template<typename Ty>
void Vec_redistribute::call_MPI(int flag)
{
  #if PRINT_TIMER>4
  TIMER_FLOPS("==Vec redistribute MPI reorder");
  double Total = 0.0;
  for(long i=0;i<long(currsend.size());i++){Total += double(currsend[i]);}
  timer.flops  += Total*sizeof(Ty);
  #endif

  if(flag_set_mem==0){print0("Buf not set. \n");qassert(false);}

  Ty* src = NULL;Ty* res = NULL;
  if(flag == 0){res = (Ty*) recvV; src = (Ty*) sendV;}
  if(flag == 1){res = (Ty*) sendV; src = (Ty*) recvV;}

  unsigned int off = sizeof(Ty);MPI_Datatype curr = MPI_BYTE;
  unsigned int M_size = get_MPI_type<Ty >(curr );
  qassert(off%M_size == 0);off = off/M_size;

  if(tem_off != int(off) or update_off == true){
    ///if(findN && sizeof(Ty)== 8){curr = MPI_FLOAT ;off = off/sizeof(float)  ;findN=false;}
    ///if(findN && sizeof(Ty)==16){curr = MPI_DOUBLE;off = off/sizeof(double) ;findN=false;}
    ///////print0("Check int %d, long %d \n",sizeof(int), sizeof(long));
    #pragma omp parallel for
    for(int n=0;n<Nmpi/mt;n++)sendM[n] = off*currsend[n];
    #pragma omp parallel for
    for(int n=0;n<Nmpi/mt;n++)recvM[n] = off*currrecv[n];
    #pragma omp parallel for
    for(int n=0;n<Nmpi/mt;n++)splsM[n] = off*currspls[n];
    #pragma omp parallel for
    for(int n=0;n<Nmpi/mt;n++)rplsM[n] = off*currrpls[n];
    tem_off = off;
    update_off = false;
  }

  ////======Copy data
  if(copy_same_node){
  int ranklocal = map_mpi_vec[fd->rank];
  //int ranklocal = fd->rank;
  qassert(currsend[ranklocal] == currrecv[ranklocal]);
  if(currsend[ranklocal] != 0){
    cpy_data_thread(&res[currrpls[ranklocal]], &src[currspls[ranklocal]], currsend[ranklocal], GPU, false);
    sendM[ranklocal] = 0;
    recvM[ranklocal] = 0;
  }}

  if(mode_MPI == 0)
  {MPI_Alltoallv(src,(int*) &sendM[0],(int*) &splsM[0], curr,
                 res,(int*) &recvM[0],(int*) &rplsM[0], curr, vec_comm);}

  if(mode_MPI == 1)
  {
    std::vector<MPI_Request> send_reqs(Nmpi/mt);
    int mpi_tag = omp_get_thread_num()*Nmpi + map_mpi_vec[fd->rank];
    int c1 = 0;
    for(int n = 0; n < Nmpi/mt; n++){
      if(sendM[n]!=0){MPI_Isend(&src[currspls[n]], sendM[n], curr, n, mpi_tag + n, vec_comm, &send_reqs[c1]);c1 += 1;}
    }

    for(int n = 0; n < Nmpi/mt; n++){
      if(recvM[n]!=0){MPI_Recv( &res[currrpls[n]], recvM[n], curr, n, mpi_tag + n, vec_comm, MPI_STATUS_IGNORE);}
    }    
    MPI_Waitall(c1, send_reqs.data(), MPI_STATUS_IGNORE);
  }
  qacc_barrier(dummy);

  src = NULL;res = NULL;
}

template<typename Ty>
void Vec_redistribute::re_order_recv(int flag)
{
  TIMERB("Copy and arrange data");
  Ty* recv = (Ty*) recvV;
  Ty* send = (Ty*) sendV;

  long bfac = Nv[orderN[2]]*civa;
  LInt* m0 = (LInt*) qlat::get_data(map_order).data();
  LInt* m1 = (LInt*) qlat::get_data(map_Dorder).data();
  if(flag==0){cpy_data_from_index(&send[0],&recv[0], m0, m1, map_order.size(), bfac, GPU, true);}
  if(flag==1){cpy_data_from_index(&recv[0],&send[0], m1, m0, map_order.size(), bfac, GPU, true);}
}

//////buf size  --> b0 * Nt * (nx*ny*nz/(Nx*Ny*Nz)) * Nx*Ny*Nz * civa * sizeof(Ty)
//////Data will be modified for sendbuf and recvbuf, results on sendbuf
//////vector maps ---> b0,mz,my,mx, (Nxyz), civ --> fd.get_mi_curr()*b0 + nbi, (nxyz), civ ;
template<typename Ty>
void Vec_redistribute::reorder(Ty *sendbuf,Ty *recvbuf,int b0_or,int civa_or,int flag)
{
  TIMERB("Vec_redistribute::reorder");

  if(flag_set_mem==0){set_mem(b0_or,civa_or);}
  if(flag_set_mem==1){if(b0 != b0_or or civa != civa_or){
    set_mem(b0_or,civa_or);}
  }

  if(fd->mz*fd->my*fd->mx == 1){return ;}

  //if(flag%2 != 0)if(howmany != civa_or/2)set_fft();

  ////Set the size of civa, civa
  sendV = sendbuf;recvV = recvbuf;

  if(flag > -2 and flag < 2 )
  {
    call_MPI<Ty >(0);
    /////from recvbuf to sendbuf
    re_order_recv<Ty>(0);
  }

  //////From reorder to original
  if(flag == 100)
  {
    re_order_recv<Ty>(1);
    call_MPI<Ty >(1);
  }
}

inline Vec_redistribute::~Vec_redistribute(){
  currsend.resize(0);
  currrecv.resize(0);
  currspls.resize(0);
  currrpls.resize(0);
  //////mapcur_Vtoi.resize(0);
}

/////To T distribute or whole vectors on each node
struct Rotate_vecs{

  fft_desc_basic fd;
  fft_desc_basic fd0;
  Vec_redistribute *vec_re0;
  Vec_redistribute *vec_re1;
  bool  GPU;
  void* src;
  void* buf;
  int mode;
  //int duel_rorate;
  int b0, civa, N_extra;
  int bsize;
  size_t vol_buf,Bsize;
  ////std::vector<int > map_vecs;

  bool flag_mem_set;

  Rotate_vecs(fft_desc_basic &fd_set, bool GPU_set=true):fd(),fd0(){
    (void)GPU_set;
    #ifndef QLAT_USE_ACC
    GPU = false;
    #else
    GPU = GPU_set;
    #endif

    copy_fft_desc(fd, fd_set);
    vec_re0 = NULL;vec_re1 = NULL;buf = NULL;src = NULL;
    vol_buf = 0; Bsize = 0;
    b0 = -1;civa = -1;bsize = -1;
    flag_mem_set = false;
  }

  void check_Bsize(size_t srcB){
    if(Bsize == 0 or srcB != Bsize){
    print0("==buf size %zu, input %zu", Bsize, srcB);
    abort_r("Size not match");
    }
  }

  /////mode -1 --> no rotation, mode 0 --> rotate to T zyx, mode 1 --> every node have tzyx data
  //////buf size  --> b0 * (nx*ny*nz/(Nx*Ny*Nz)) * Nx*Ny*Nz * civa * sizeof(Ty)
  template<typename Ty>
  void set_mem(int b0_or, int civa_or, int mode_set = -2){
    TIMERA("Rotate_vecs set_mem");
    if(b0_or == -2 and civa_or == -2){qassert(b0 > 0 and civa > 0 and mode > -2 and vol_buf > 0);return;}
    if(b0_or == b0 and civa_or == civa and sizeof(Ty) == bsize){
      if(mode_set == -2 or mode_set == mode){return ;}
    }
    else{clear_mem();}

    if(mode_set != -2){mode = mode_set;}else{if(mode <= -2){abort_r("Need to set mode of rotation!\n");}}
    b0 = b0_or;civa = civa_or;bsize = sizeof(Ty);
    ///if(mode == -1){clear_mem(); return; }
    if(mode < -1 or mode > 1){abort_r("Rotate_vecs mode not supportted! \n");}

    if(mode == -1){
    vol_buf = fd.Nt * fd.Nz * fd.ny * size_t(fd.nx);
    N_extra = 1;
    }

    if(mode == 0){
    vol_buf = fd.Nt * fd.nz * fd.ny * size_t(fd.nx);
    N_extra = fd.mz * fd.my * fd.mx;
    }
    if(mode == 1){
    vol_buf = fd.nt * fd.nz * fd.ny * size_t(fd.nx);
    N_extra = fd.mt * fd.mz * fd.my * fd.mx;
    }


    size_t Bsize0 = vol_buf * b0 *  civa * sizeof(Ty);
    if(Bsize != Bsize0){
    Bsize = Bsize0;
    free_buf(buf, GPU);free_buf(src, GPU);buf=NULL;src=NULL;
    if(GPU){gpuMalloc(buf, Bsize/sizeof(Ty), Ty);gpuMalloc(src, Bsize/sizeof(Ty), Ty);}
    else{ 
      src = aligned_alloc_no_acc(Bsize);
      if(mode != -1)buf = aligned_alloc_no_acc(Bsize);
    }}
    qassert(b0 > 0 and civa > 0 and vol_buf > 0);

    //map_vecs.resize( b0);
    //for(int bi=0;bi<b0;bi++){map_vecs[bi] = fd.get_mi_curr(mode + 3)*b0 + bi;}

    /////make new rotations or not
    //if(fdp.Nv[0] != fdp.nv[0]){duel_rorate = 1;}
    //else{duel_rorate = 0;}

    if(mode == -1){vec_re0 = NULL;}

    if(mode == 0){
      vec_re0 = new Vec_redistribute(fd , GPU);
      ////vec_large.reorder((T*) &gfET[0],(T*) &gfET_buf[0], 1, (dir_limit*2)*9 ,   0);
    }

    if(mode == 1){

      if(fd.Nv[0] == fd.nv[0]){

        fd0.nv[3]=1;fd0.Nv[3]=1;fd0.iniv[3]=0;
        fd0.nv[2]=fd.nv[3];fd0.Nv[2]=fd.Nv[3];fd0.iniv[2]=fd.iniv[3];
        fd0.nv[1]=fd.nv[2];fd0.Nv[1]=fd.Nv[2];fd0.iniv[1]=fd.iniv[2];

        fd0.nv[0]=fd.nv[1]*fd.nv[0];
        fd0.Nv[0]=fd.Nv[1]*fd.Nv[0];
        fd0.iniv[0]=fd.iniv[1]*fd.nv[0] + fd.iniv[0];

        fd0.set_variable();
        vec_re0 = new Vec_redistribute(fd0, GPU);
        //fd0.check_mem();
        //fd0.print_info();
        //fdp_new0 = new fft_desc_basic();
        //duel_rorate = 1;

      }
      if(fd.Nv[0] != fd.nv[0]){
        int bN = fd.Nmpi/fd.mt;
        /////need to get this number??
        int bi = fd.get_mi_curr(3);

        fd0.nv[3]=bN;fd0.Nv[3]=1;fd0.iniv[3]=bi;

        fd0.nv[2]=fd.nv[3];fd0.Nv[2]=fd.Nv[3];fd0.iniv[2]=fd.iniv[3];

        fd0.nv[1]=fd.nv[2];fd0.nv[0]=fd.nv[1]*fd.nv[0];
        fd0.Nv[1]=fd.nv[2];fd0.Nv[0]=fd.nv[1]*fd.nv[0];
        fd0.iniv[1]=0;fd0.iniv[0]=0;

        fd0.set_variable();
        ////fd0.check_mem();
        ////fd0.print_info();

        vec_re0 = new Vec_redistribute(fd , GPU);
        vec_re1 = new Vec_redistribute(fd0, GPU);

        //fd0 = new fft_desc_basic();
        //fdp_new0 = new fft_desc_basic(geo);
        //duel_rorate = 0;
      }
    }
    /////only map for 3D and 4D

    flag_mem_set = true;

  }

  void clear_mem(){
    if(flag_mem_set == false){return;}
    if(vec_re0 != NULL){delete vec_re0;vec_re0 = NULL;}
    if(vec_re1 != NULL){delete vec_re1;vec_re1 = NULL;}
    ////if(fdp_new0 != NULL){delete fdp_new0;fdp_new0 = NULL;}
    ////if(fdp_new1 != NULL){delete fdp_new1;fdp_new1 = NULL;}
    free_buf(buf, GPU);free_buf(src, GPU);
    Bsize = 0;b0 = -1;civa = -1;bsize = -1;
    N_extra = -1;vol_buf =  0;mode = -2;
    flag_mem_set = false;
  }

  template<typename Ty>
  void reorder(bool forward = true, Ty* input=NULL, int b0_set=-2, int civa_set=-2)
  {
    TIMERC("Rotate_vecs Vec redistribute");
    //set_mem<Ty>(b0_set, civa_set, mode);
    set_mem<Ty>(b0_set, civa_set);
    if(mode == -1){ return;}
    Ty* res = (Ty*)buf;int flag = -1;
    Ty* s0 = NULL;if(input == NULL){s0 = (Ty*) src;}else{s0 = input;}
    if( forward){flag = 0  ;}if(!forward){flag = 100;}

    if(mode == 0){vec_re0->reorder(s0, res, b0, civa, flag);}
    if(mode == 1){
      if(fd.Nv[0] == fd.nv[0]){vec_re0->reorder(s0, res, b0, civa, flag);}
      if(fd.Nv[0] != fd.nv[0]){
        if( forward){
          vec_re0->reorder(s0, res, b0*fd.mt, civa, flag);
          vec_re1->reorder(s0, res, b0, civa, flag);
        }
        if(!forward){
          vec_re1->reorder(s0, res, b0, civa, flag);
          vec_re0->reorder(s0, res, b0*fd.mt, civa, flag);
        }
      }
    }

    s0 = NULL;

  }


  ~Rotate_vecs(){
    clear_mem();
  }

};

}

#endif
