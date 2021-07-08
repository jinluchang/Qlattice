// utils_FFT_redistribute.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_FFT_REDISTRIBUTE
#define UTILS_FFT_REDISTRIBUTE

#pragma once
#include "general_funs.h"
#include "fft_desc.h"

namespace qlat{


////////each MPI NT will have all nx,ny,nz
class FFT_single
{
public:
  long noden;int rank;int Nmpi;
  int nx,ny,nz,nt,vol,Nvol;
  int Nx,Ny,Nz,Nt;
  int mx,my,mz,mt;
  int CPU;

  qlat::vector<int> Nv,nv,mv;
  qlat::vector<qlat::vector<int> > Pos0;
  qlat::vector<qlat::vector<int>  > mi_list;

  qlat::vector<int > orderN;
  fft_desc_basic *fd;


  int biva,civa;
  void* sendV; void* recvV;
  /////void* bufV ;

  qlat::vector<int > currsend;
  qlat::vector<int > currspls;
  qlat::vector<int > currrecv;
  qlat::vector<int > currrpls;

  qlat::vector<int > sendM;
  qlat::vector<int > splsM;
  qlat::vector<int > recvM;
  qlat::vector<int > rplsM;

  qlat::vector<LInt > mapcur_Vtoi;

  /////May need to change?
  qlat::vector<int > secT;

  /////void set_mem(int biva_or,int civa_or, qlat::vector<int > secT_or);
  void set_mem(int biva_or,int civa_or);

  template<typename Ty>
  void call_MPI(int flag);

  template<typename Ty>
  void re_order_recv(int flag);

  template<typename Ty>
  void reorder(Ty *sendbuf,Ty *recvbuf,int biva_or,int civa_or,int fft_flag= 0);
  ////void reorder(Ftype *sendbuf,Ftype *recvbuf,int biva_or,int civa_or,qlat::vector<int > secT_or,int fft_flag= 0);

  //int howmany;
  ///FFTtype_complex *src_dat,*res_dat;
  ///FFTtype_plan puse0;
  ///FFTtype_plan puse1;
  qlat::vector<int > map_mpi_vec;

  std::vector<MPI_Comm > vec_comm_list;

  int flag_set_mem;
  //int flag_set_fft;

  ////int flag_set_fft_force;
  ////void set_fft(int civau=-1);
  ////void do_fft(bool fftdir=true);

  FFT_single(fft_desc_basic &fds, int cpu = 0);
  ~FFT_single();

};


FFT_single::FFT_single(fft_desc_basic &fds, int cpu)
{
  TIMER("Construct FFT_single");
  fd = &fds;
  CPU = cpu;

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
  Pos0 = fds.Pos0;
  mi_list = fds.mi_list;

  biva = 1; civa = 1;

  //////Needed for three point functions, need copy from dev
  secT.resize(fd->Nmpi);
  for(int i=0;i<secT.size();i++){secT[i] = fd->Nt;}

  //////if(secT_or.size()!=fd->Nmpi){print0("secT wrong %8d ! \n", int(secT_or.size()));qassert(false);}
  ////secT = secT_or;
  ////Check same number secT for MPI
  ////long mvol = mx*my*mz;
  qlat::vector<int > check;check.resize(0);
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
  for(int mapi=0;mapi<map_mpi_vec.size();mapi++){map_mpi_vec[mapi] = 0;}


  //for(int icomm=0;icomm<vec_comm_list.size();icomm++)
  //{MPI_Comm_free(&vec_comm_list[icomm]);}
  vec_comm_list.resize(0);

  MPI_Comm vec_comm;
  int color_xyz = fd->init;
  MPI_Comm_split(get_comm() ,color_xyz,rank,&vec_comm);
  {
    int int_tem = -1;
    MPI_Comm_rank(vec_comm,&int_tem);
    map_mpi_vec[fd->rank] = int_tem;
  }
  sum_all_size((int*) (&map_mpi_vec[0]),Nmpi);
  vec_comm_list.push_back(vec_comm);


  flag_set_mem = 0;

  //flag_set_fft = 0;
  //howmany = 0;

  /////flag_set_fft_force = 0;

}

void FFT_single::set_mem(int biva_or,int civa_or)
{
  TIMER("FFT single set mem");

  biva = biva_or;civa = civa_or;

  if(biva<=0 or civa<=0){
    print0("Need at least biv = 1, civ %6d !\n",civa);qassert(false);
  }


  ///map i --> nz*ny*nx   to recieve position
  /////2 is the fatest
  /////mapcur_itoV.resize(biva*nx*ny*nz*Nt);
  /////mapcur_Vtoi.resize(biva*nx*ny*nz*Nt);
  mapcur_Vtoi.resize(nx*ny*nz/Nv[orderN[2]]);
  int Nts = secT[fd->rank];
  for(int tmi=0;tmi<mt;tmi++)
  {
    int t0 = Pos0[mi_list[tmi][0]][3];
    if(t0 == Pos0[rank][3])
    for(int temmi=0;temmi<mi_list[tmi].size();temmi++)
    {
      //int off = Nx*Ny*Nz*Nt*temmi;
      //for(int iv=0;iv<Nx*Ny*Nz*Nt;iv++)
      //#pragma omp parallel for
      for(LInt iv=0;iv<Nx*Ny*Nz;iv++)
      {
        LInt Vi = fd->maporder_Nitoipos[mi_list[tmi][temmi]][iv];

        int umi = temmi;
        LInt offnode = Nx*Ny*Nz*Nts*umi;
        LInt temiv = iv%(Nv[0]*Nv[1]*Nv[2]);
        LInt inode = biva*offnode + temiv;
        if(Vi%Nv[orderN[2]] == 0)
        if(Vi/(nv[0]*nv[1]*nv[2]) == 0)
        {
          LInt temVi = Vi%(nv[0]*nv[1]*nv[2]);
          mapcur_Vtoi[temVi/Nv[orderN[2]]]    = inode/Nv[orderN[2]];
        }

      }
    }
  }

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
    int t0 = Pos0[mi_list[tmi][0]][3]; //Determine the initial time
    if(t0 == Pos0[rank][3]) // If current node initial time is the working initial time
    {
      for(int temmi=0;temmi<mi_list[tmi].size();temmi++) //Determine the communication nodes of this node
      {
        //int umi = mi_list[tmi][temmi];
        int umi = temmi;

        LInt rankglobal = mi_list[tmi][umi];
        LInt ranklocal = map_mpi_vec[rankglobal];

        ////no 2, 2 is inside biva
        currsend[ranklocal] = biva*Nx*Ny*Nz*secT[fd->rank]*civa;
        currspls[ranklocal] = biva*Nx*Ny*Nz*secT[fd->rank]*umi*civa;
        currrecv[ranklocal] = biva*Nx*Ny*Nz*secT[fd->rank]*civa;
        currrpls[ranklocal] = biva*Nx*Ny*Nz*secT[fd->rank]*umi*civa;
      }
    }
  }

  flag_set_mem = 1;
}

template<typename Ty>
void FFT_single::call_MPI(int flag)
{
  TIMER("MPI to single node");
  if(flag_set_mem==0){print0("Buf not set. \n");qassert(false);}

  Ty* src;Ty* res;
  if(flag == 0){res = (Ty*) recvV; src = (Ty*) sendV;}
  if(flag == 1){res = (Ty*) sendV; src = (Ty*) recvV;}

  unsigned int off = sizeof(Ty);Ty atem = 0;
  MPI_Datatype curr = MPI_BYTE;unsigned int M_size = 1;
  get_MPI_type(atem, curr, M_size);
  qassert(off%M_size == 0);off = off/M_size;

  ///if(findN && sizeof(Ty)== 8){curr = MPI_FLOAT ;off = off/sizeof(float)  ;findN=false;}
  ///if(findN && sizeof(Ty)==16){curr = MPI_DOUBLE;off = off/sizeof(double) ;findN=false;}
  ///////print0("Check int %d, long %d \n",sizeof(int), sizeof(long));

  for(int n=0;n<Nmpi/mt;n++)sendM[n] = off*currsend[n];
  for(int n=0;n<Nmpi/mt;n++)recvM[n] = off*currrecv[n];
  for(int n=0;n<Nmpi/mt;n++)splsM[n] = off*currspls[n];
  for(int n=0;n<Nmpi/mt;n++)rplsM[n] = off*currrpls[n];


  //////======Copy data
  int ranklocal = map_mpi_vec[fd->rank];
  //int ranklocal = fd->rank;
  qassert(currsend[ranklocal] == currrecv[ranklocal]);
  if(currsend[ranklocal] != 0){
    copy_data(&res[currrpls[ranklocal]], &src[currspls[ranklocal]], currsend[ranklocal], CPU, true);
    sendM[ranklocal] = 0;
    recvM[ranklocal] = 0;
  }

  {MPI_Alltoallv(src,(int*) &sendM[0],(int*) &splsM[0], curr,
                 res,(int*) &recvM[0],(int*) &rplsM[0], curr, vec_comm_list[0]);}


  //////for(int n=0;n<Nmpi;n++){print0("==rank %d, n %d, send %d \n", rank, n, send[n]);}
 
  //////======Copy data

  ///print0("=====min %d %d %d %d\n",send[0     ] , recv[0     ],spls[0        ], rpls[0        ]);
  ///print0("=====max %d %d %d %d\n",send[Nmpi-1] , recv[Nmpi-1],spls[Nmpi/mt-1], rpls[Nmpi/mt-1]);

  src = NULL;res = NULL;
}

template<typename Ty>
void FFT_single::re_order_recv(int flag)
{
  TIMER("Copy and arrange data");
  //int Nt = fd->Nt;
  ////print0("===nx %d, ny %d, nz %d, Nx %d, Ny %d, Nz %d  \n", nx,ny,nz, Nv[0],Nv[1],Nv[2]);
  size_t svol = fd->nx*fd->ny*fd->nz;
  int Nts = secT[fd->rank];
  Ty* recv = (Ty*) recvV;
  Ty* send = (Ty*) sendV;

  auto& Nv = this->Nv;
  auto& orderN = this->orderN;
  auto& civa = this->civa;
  auto& biva = this->biva;
  auto& mapcur_Vtoi = this->mapcur_Vtoi;

  //copy_data(res, src, currsend, CPU, true);
  /////print0("=====Check Nts %d Nv %d \n", Nts, Nv[orderN[2]]);
  //  #pragma omp parallel for
  size_t LoopN = svol/Nv[orderN[2]];
  if(CPU == 0){
  qacc_for(iv, long(biva*Nts * LoopN), {
    int bi   =  iv/(Nts*LoopN);
    int ti   = (iv%(Nts*LoopN))/LoopN;
    size_t i =  iv%(LoopN);
    size_t Aoff = (bi*Nts + ti)*svol;

    size_t Boff = (bi*Nts + ti)*Nv[0]*Nv[1]*Nv[2];

    size_t off  = (Aoff +      i        *Nv[orderN[2]])*civa;
    size_t offV = (Boff + mapcur_Vtoi[i]*Nv[orderN[2]])*civa;

    ////if(flag==0)copy_data(&send[off ],&recv[offV], Nv[orderN[2]]*civa, CPU, false);
    ////if(flag==1)copy_data(&recv[offV],&send[off ], Nv[orderN[2]]*civa, CPU, false);
    for(size_t ip=0;ip<Nv[orderN[2]]*civa;ip++){
      if(flag==0)send[off  + ip] = recv[offV + ip];
      if(flag==1)recv[offV + ip] = send[off  + ip];
    }
  });
  ////qacc_barrier(dummy);

  }
  if(CPU == 1){
  qthread_for(iv, long(biva*Nts * LoopN), {
    int bi   =  iv/(Nts*LoopN);
    int ti   = (iv%(Nts*LoopN))/LoopN;
    size_t i =  iv%(LoopN);
    size_t Aoff = (bi*Nts + ti)*svol;

    size_t Boff = (bi*Nts + ti)*Nv[0]*Nv[1]*Nv[2];

    size_t off  = (Aoff +      i        *Nv[orderN[2]])*civa;
    size_t offV = (Boff + mapcur_Vtoi[i]*Nv[orderN[2]])*civa;

    for(size_t ip=0;ip<Nv[orderN[2]]*civa;ip++){
      if(flag==0)send[off  + ip] = recv[offV + ip];
      if(flag==1)recv[offV + ip] = send[off  + ip];
    }
  });
  
  }


}


//////buf size  --> biva * (nx*ny*nz/(Nx*Ny*Nz)) * Nx*Ny*Nz * civa * sizeof(Ty)
//////Data will be modified for sendbuf and recvbuf, results on sendbuf
//////vector maps ---> biva,mz,my,mx, (Nxyz), civ --> fd.get_mi_curr()*biva + nbi, (nxyz), civ ;
template<typename Ty>
void FFT_single::reorder(Ty *sendbuf,Ty *recvbuf,int biva_or,int civa_or,int fft_flag)
{
  TIMER("FFT reorder ");
  if(flag_set_mem==0){set_mem(biva_or,civa_or);
    //if(fft_flag%2 != 0)set_fft();
  }
  if(flag_set_mem==1){if(biva != biva_or or civa != civa_or){
    set_mem(biva_or,civa_or);
    //if(howmany != civa_or/2){if(fft_flag%2 != 0)set_fft();}
    }
  }

  //if(fft_flag%2 != 0)if(howmany != civa_or/2)set_fft();

  ////Set the size of civa, civa
  sendV = sendbuf;recvV = recvbuf;
  Ty* recv = (Ty*) recvV;
  Ty* send = (Ty*) sendV;
  ///////int Nt = fd->Nt;
  size_t svol = fd->nx*fd->ny*fd->nz;
  int Nts = secT[fd->rank];

  if(fd->mz*fd->my*fd->mx == 1){
    /////copy_data(recv, send, biva*Nts*svol*civa, CPU, true);
    /////print0("=====Fast mode \n");
    return ;
  }

  if(fft_flag > -2 and fft_flag < 2 )
  {
    call_MPI<Ty >(0);
    /////from recvbuf to sendbuf
    re_order_recv<Ty>(0);

    /////copy_data(recv, send, biva*Nts*svol*civa, CPU, true);
  }

  /////print0("Max in %d %d \n", biva*Nts*svol*civa, INT_MAX);

  ///if(fft_flag == 1){bool fftdir = true;do_fft(fftdir);}
  ///if(fft_flag ==-1){bool fftdir = false;do_fft(fftdir);}

  //////From reorder to original
  if(fft_flag == 100)
  {
    re_order_recv<Ty>(1);
    call_MPI<Ty >(1);
    //memcpy(recv,send, sizeof(Ftype)*biva*Nts*svol*civa);
    //////copy_data(recv, send, biva*Nts*svol*civa, CPU, true);
  }

  send = NULL;recv = NULL;
}

FFT_single::~FFT_single(){

  vec_comm_list.resize(0);

  currsend.resize(0);
  currrecv.resize(0);
  currspls.resize(0);
  currrpls.resize(0);
  mapcur_Vtoi.resize(0);
}


}

#endif
