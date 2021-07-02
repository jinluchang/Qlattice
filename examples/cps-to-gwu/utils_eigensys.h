#ifndef UTILS_EIGENSYS_H
#define UTILS_EIGENSYS_H

#pragma once
#include <qlat/qlat.h>
#include "reduce_V_dev.h"
#include "gammas.h"
#include "fft_desc.h"
#include "utils_Matrix_prod.h"

#define SUMMIT 1

namespace qlat{

struct eigen_ov {
  //layout_minsurface_eo desc;
  //layout_minsurface_eo *desc;
  qlat::fft_desc_basic* fdp;

  EigenM Mvec;     //  nvec --> 2*bfac --> b_size*6
  EigenM Mvec_Sm;  //  Smeared eigensystem 

  EigenM Mvec0;  //  Test eigensystem 
  EigenM Mvec1;  //  Test eigensystem 


  ////int Edouble;
  ////int smear_Eigen;

  /////lattice copied from fd
  int nx;int ny;int nz;int nt;
  int Nx;int Ny;int Nz;int Nt;
  LInt noden;

  ////Eigensystem information
  int n_vec, n0;
  int num_zero;
  int one_minus_halfD;
  double rho;
  double Eeigenerror;
  EigenV eval_self;

  Complexq* ptmp; size_t ptmp_size;
  Complexq* stmp; size_t stmp_size;
  std::vector<Complexq* > Eigenbuf;
  std::vector<Complexq* > Eigendyn;


  ///std::vector< std::vector<int> > ranged;
  /////Group the inner vector for vectorizations
  LInt b_size, bfac;
  LInt bfac_group;

  ///std::vector<int> bi_l;
  ///std::vector<int> bp_l;
  ///std::vector<int> Cur_Xn;


  //////Buffers
  EigenV alpha;
  EigenV alpha_bfac;
  EigenV alpha_list; 
  EigenV eval_list;

  int ncut0,ncut1, ncutgpu, ncutbuf;
  ////std::vector<int > mapN;
  int buffGPU;

  ////EigenV srcM;
  ////EigenV propM;


  //EigenM EProp;
  //EigenM temP;
  //EigenV temj;

  //EigenM temjC;
  //std::vector<EigenV >  alpha_listC;

  ////EigenM MPItem;
  ////EigenM Esource;

  //Eigen::VectorXcd temk;

  ////3pt function
  EigenV alpha_list3pt;
  EigenV alpha_list3ptC;
  EigenV alpha3pt;
  EigenV alpha3pt_ker_low;


  /////Construction and memory allocations
  eigen_ov(qlat::fft_desc_basic &fd,int n_vec_or,long long bsize0=-1,int nV=1);

  void copy_Mvec(int mode = 0, bool random=false);
  void touch_GPU(EigenV &Mres, long long offM=0,long long size=-1, int mode = 1);
  void touch_GPU(EigenM &Mres, long long offM=0,long long size=-1, int mode = 1);
  void untouch_GPU(EigenV &Mres);

  inline Complexq* getEigenP(int ni, size_t xi){
    Complexq* buf = NULL;
    int   chi = xi/(bfac*b_size);
    size_t vi = xi%(bfac*b_size);
    size_t bi = vi/b_size;
    size_t bj = vi%b_size;
    if(ni  < ncutbuf){
      size_t off1 = ((chi*bfac + bi)*ncutbuf + ni )*b_size + bj;
      size_t Lb = bfac_group*ncutbuf*size_t(b_size);
      size_t Li = off1/Lb;
      size_t Lj = off1%Lb;
      /////buf = &(Mvec0[Li][Lj]);
      buf = &(Eigenbuf[Li][Lj]);
    }
    if(ni >= ncutbuf){
      int ncur = ni - ncutbuf;
      #if SUMMIT==1
      int na   = 0;
      #else
      int na   = ncur/ncutgpu;
      #endif

      int nb   = ncur%ncutgpu;

      size_t off1 = (((na*2 + chi)*bfac + bi)*ncutgpu + nb)*b_size + bj;
      size_t Lb = bfac_group*ncutgpu*size_t(b_size);
      size_t Li = off1/Lb;
      size_t Lj = off1%Lb;

      #if SUMMIT==1
      buf = &(Eigendyn[Li][Lj]);
      #else
      buf = &(Mvec0[Li][Lj]);
      #endif
    }
    return buf;
  }

  void setup_bfac(long long bsize0=-1, int nmass=-1);

  void load_eivals(char *enamev,double rho_or,double Eerr=1e-9);
  
  void load_eigen(int icfg, std::string &ov_evecname, io_gwu &io,
    int checknorm = 1, double kappa=0.2,double eigenerror=1e-9);
  void random_eigen();

  void initiallize_mass(std::vector<double> &mass,int Nv=12,int one_minus_halfD_or=1);
  void initiallize_mass();

  ///void seq_L(vector **prop3,double mass_d,double rho,int_vector &map,bool bcast,double_complex* ker_low);
  ///void seq_L(vector **prop3,double mass_d,double rho,int_vector &map,double_complex* ker_low);
  ///void setup_L(vector **prop2);

  ~eigen_ov()
  {
    ////Mvec.resize(0);
    #ifdef QLAT_USE_ACC
    if(ptmp!=NULL)gpuErrchk(cudaFree(ptmp));
    if(stmp!=NULL)gpuErrchk(cudaFree(stmp));
    for(int i=0;i<Eigenbuf.size();i++)if(Eigenbuf[i]!=NULL)gpuErrchk(cudaFree(Eigenbuf[i]));
    for(int i=0;i<Eigendyn.size();i++)if(Eigendyn[i]!=NULL)gpuErrchk(cudaFree(Eigendyn[i]));
    #else
    if(ptmp!=NULL)delete [] ptmp;
    if(stmp!=NULL)delete [] stmp;
    ////if(Eigenbuf!=NULL)delete [] Eigenbuf;
    for(int i=0;i<Eigendyn.size();i++)if(Eigendyn[i]!=NULL)delete [] Eigendyn[i];
    #endif
  }

};

inline Complexq inv_self(const Complexq& lam, double m, double rho,int one_minus_halfD=1)
{
  std::complex<double > tem(lam.real(),lam.imag());
  std::complex<double > v0 = (one_minus_halfD>0)?(1.0-tem/2.0)/(rho*tem+m*(1.0-tem/2.0)):1.0/(rho*tem+m*(1.0-tem/2.0));
  Complexq res(v0.real(),v0.imag());
  return res;
}


inline void resize_EigenM(EigenM& a, size_t n0, size_t n1)
{
  a.resize(0);
  a.resize(n0);
  for(size_t iv=0;iv<n0;iv++)
  {
    a[iv].resize(n1);
    zeroE(a[iv], 1);
  }
}

void eigen_ov::setup_bfac(long long bsize0, int nV)
{
  LInt bcut = 300;
  #ifdef QLAT_USE_ACC
  bcut = 128;
  #endif
  if(bsize0 != -1){bcut = std::abs(bsize0);}
  #ifdef QLAT_USE_ACC
  LInt bfac_or = 6*Nt;
  if(bfac_or < bcut)
  {
    LInt t0 = int(bcut/(6*Nt))+1;
    for(LInt temb=t0;temb>=1;temb--)
    {
      if((noden/Nt)%temb == 0)
      {
        bfac_or = temb*6*Nt;
        break;
      }
    }
  }
  bfac = bfac_or;
  b_size = noden*6/bfac;
  #else
  LInt bsize_or = noden/Nt;
  if(bsize_or > bcut)
  {
    for(LInt temb=bcut;temb<noden/Nt;temb++)
    {
      if((noden/Nt)%temb == 0)
      {
        bsize_or = temb;
        break;
      }
    }
  }
  b_size = bsize_or;
  bfac = noden*6/b_size;
  #endif

  //bfac = bsize_or;
  //b_size = noden*6/bfac;

  /////bfac --> b_size
  //b_size = bsize_or;
  //bfac = noden*6/b_size;

  //b_size = bsize_or;
  //bfac = noden*6/b_size;


  if(bfac%6 !=0){print0("Cannot understand sites on bfac %5ld, bsize %10ld, Total %10ld !\n",
    bfac, b_size, noden*6);abort_r("");}
  if((noden*6)%bfac !=0){print0("Cannot understand sites on node*6/Nt %10ld bfac %10ld!\n",noden*6/Nt,bfac);abort_r("");}
  if(bfac%(Nt*6) !=0){print0("Cannot understand sites on node %10ld bfac %10ld, tsize %5d!\n",noden,bfac,Nt);abort_r("");}

  print0("Lattice sites on node %10ld bfac %10ld bsize %10ld, tsize %5d !\n",noden,bfac,b_size, Nt);
  print0("bfac %10ld, n_vec %10d, b_size %10ld, matrix %10ld x %10ld \n",bfac,n_vec,b_size,2*bfac*n_vec,b_size);


  ///int massN_gpu, bfac_group;
  /////To avoid very large continuous memories on CPU
  bfac_group = 1;

  ////ncutgpu  = 16;
  ////ncutbuf=0;
  /////redefine variables on GPU
  #ifdef QLAT_USE_ACC

  buffGPU    = 2;

  double freeD = 0;double totalD=0;
  size_t freeM = 0;size_t totalM = 0;
  cudaMemGetInfo(&freeM,&totalM);
  freeD = freeM*pow(0.5,30);totalD = totalM*pow(0.5,30);
  struct sysinfo s_info;
  sysinfo(&s_info);

  //int Ns = 12;
  int Ns = nV;if(Ns <=0)Ns = 2;
  ////long long Lat   = noden;
  double memV     = noden*12.0*sizeof(Complexq)*pow(0.5,30);
  double mem_prop = Ns * memV;
  if(totalD < mem_prop + 64*memV){
    print0("===GPU Memory too small, Total %.3e, prop %.3e, 30V %.3e, increase nodes! \n", totalD, mem_prop, 64*memV);
    qassert(false);
  }

  int vfac = ncutgpu;int vini = 8 * vfac;
  int vres = int((totalD*0.90 - mem_prop )/memV); 
  //int gsimple = 32;
  //if(temN >= 32){gsimple = (temN/32)*32;}else{gsimple = 32;}

  if(vres > n_vec ){
    ncutgpu = 0;ncutbuf = n_vec;
    //ncutgpu = n_vec;ncutbuf = 0;
  }
  else{
    if(vres > vini){
      ncutbuf = vres - vini;
      ncutgpu=vini;
      //#if SUMMIT==1
      //ncutgpu=vini;
      //#else
      //ncutgpu=vini;
      //#endif
    }
    else{ncutbuf = 0;ncutgpu = ((vres)/16)*16;if(ncutgpu < 16)ncutgpu = 16;}

    //int temN = int((totalD - 12*memV - mem_prop)/(5*memV));
    //if(ncutgpu > n_vec){ncutgpu = n_vec; ncutbuf = 0;}
    //else{ncutgpu = (temN/16)*16;ncutbuf = 3*ncutgpu;}

    //if(ncutbuf < 16){
    //  ncutbuf = 0; 
    //  temN = int((totalD - 12*memV - mem_prop)/(2*memV));
    //  if(temN >= 32){ncutgpu = (temN/16)*16;}else{ncutgpu = 16;}
    //}

    //temN = int((totalD - 12*memV - mem_prop - 2*ncutgpu*memV)/(memV));
    //if(temN > 16){
    //  ncutbuf=temN;
    //}

  }

  if(ncutgpu == 0 or ncutbuf == 0){bfac_group = 2*bfac;}else{
  bfac_group = 32;
  for(int bini=bfac_group;bini >= 1; bini--){
    size_t tem = get_threads(bini, bfac, 0);
    if(tem != 0 and bfac%tem == 0){bfac_group = tem;break;}
    bfac_group = 1;
  }}

  print0("===CPU free %.3e GB, total %.3e GB; GPU free %.3e GB, total %.3e GB, prop %.3e GB, v %d  %.3e GB, buf %d  %.3e GB, bfac_group %d, E %.3e GB. \n"
          , s_info.freeram*pow(0.5,30),s_info.totalram*pow(0.5,30),freeD,totalD
          , mem_prop,ncutgpu, memV*ncutgpu, ncutbuf, memV*ncutbuf, bfac_group, memV*n_vec);

  #endif

  ////if(bfac < bfac_group or bfac_group == 1){print0("GPU need bfac to be not too small !!!!");qassert(false);}


}


eigen_ov::eigen_ov(qlat::fft_desc_basic &fd,int n_vec_or, long long bsize0, int nV)
{
  ////Default values
  //for(int i=0;i<12;i++) src_tmp0[i]=NULL;
  //desc = &desc_or;
  fdp  = &fd;
  if(fd.order_ch != 0){abort_r("Currently not supported for change fd order.\n ");}

  rho = 1.5;
  Eeigenerror = 1e-11;
  num_zero = 0;
  one_minus_halfD = 1;

  nx = fd.nx;ny = fd.ny;nz = fd.nz;nt = fd.nt;
  Nx = fd.Nx;Ny = fd.Ny;Nz = fd.Nz;Nt = fd.Nt;
  noden = fd.noden;

  n_vec = n_vec_or;
  n0    = 0;
  ptmp = NULL; ptmp_size = 0;
  stmp = NULL; stmp_size = 0;
  Eigenbuf.resize(0);
  Eigendyn.resize(0);


  ncut0      =  20;
  ncut1      = 800;

  ncutgpu    = 8;
  ncutbuf    = 0;

  buffGPU    = 2;
  bfac_group = 1;

  setup_bfac(bsize0, nV);
  fflush_MPI();

  ////smear_Eigen = 0;
  ////bfac = 1;
  ////Default values

  //ranged.resize(4);
  //for(unsigned int i= 0; i< desc->sites_on_node; i++)
  //{
  //  position p = desc->get_position(i,get_node_rank());
  //  ranged[0].push_back(p.x());
  //  ranged[1].push_back(p.y());
  //  ranged[2].push_back(p.z());
  //  ranged[3].push_back(p.t());
  //}
  //for(int i=0;i<4;i++){
  //  std::sort(ranged[i].begin(), ranged[i].end());
  //  ranged[i].erase( unique( ranged[i].begin(), ranged[i].end() ), ranged[i].end() );
  //}
  //Nx=ranged[0].size();Ny=ranged[1].size();Nz=ranged[2].size();Nt=ranged[3].size();
  //nx=desc->nx;ny=desc->ny;nt=desc->nt;nz=desc->nz;
  ////print0("=Need=====Size of dim,  x %5d y %5d z %5d t %5d \n",ranged[0].size(),ranged[1].size(),ranged[2].size(),ranged[3].size());
  //print0("=Need=====Size of dim,  x %5d y %5d z %5d t %5d \n",Nx,Ny,Nz,Nt);

  //temj.resize(b_size);temj.setZero();
  //temP.resize(2*bfac*12,b_size);temP.setZero();

  ////MPItem.resize(bfac,b_size);
  /////Get Postion Mapping
}

void eigen_ov::load_eivals(char *enamev,double rho_or,double Eerr)
{
  rho = rho_or;

  std::vector<double > v,e;
  load_gwu_eigenvalues(v,e, enamev);
  if((n0+n_vec)*2 > int(v.size()))abort_r("Eigen value size too small! ");
  eval_self.resize(n_vec);
  for(int iv=0;iv<n_vec;iv++){
    int off= iv+n0;
    eval_self[iv] = Complexq(v[off*2+0],v[off*2+1]);
    /////////cps base eigen value will not change to complex conjugate
    //////eval_self[iv] = qconj(Complexq(v[off*2+0],v[off*2+1]));
  }

  Eeigenerror = Eerr;
  Ftype rho_tem = rho;
  for(int j=0; j<n_vec; ++j){
    eval_self[j] = eval_self[j]/rho_tem;
    if(abs(eval_self[j]) < Eerr) num_zero += 1;
  }
  print0("num_zero %5d \n",num_zero);

}

inline void checknormM(EigenM &Mvec, Ftype err=1e-3)
{
  TIMER("Check norm Mvec");
  if(Mvec.size() == 0)return;
  LInt noden = Mvec[0].size()/12;
  for(unsigned int iv=0;iv < Mvec.size();iv++)
  {
    Ftype normf = get_norm_vec((Ftype*) &Mvec[iv][0], noden);
    if(fabs(normf - 1.0) > err){
      print0("Eigen vector %5d, norm %.3e wrong. \n",iv, normf);
      abort_r("");
    }
  }
}

inline void rotate_dc(EigenM &Mvec)
{
  TIMER("Rotate eigen d,c");
  if(Mvec.size() == 0)return ;
  if(Mvec[0].size() %12 != 0){abort_r("Mvec size not correct! \n");}
  size_t Nvol = Mvec[0].size()/12;
  Evector tmp;tmp.resize(Nvol*12);
  /////Switch d,c to outter side
  for(unsigned int iv=0;iv<Mvec.size();iv++){
    Complexq* src = (Complexq* ) &(Mvec[iv][0]);
    Complexq* buf = (Complexq* ) &(tmp[0]);
    memcpy(buf,src, Nvol*12*sizeof(Complexq));
    #pragma omp parallel for
    for(size_t isp=0;isp<size_t(Nvol);isp++){
      for(int d=0;d<4;d++)
      for(int c=0;c<3;c++){
        src[(d*3 + c)*Nvol + isp ] = buf[isp*12 + d*3 + c];
      }
    }
  }
}

inline void rotate_cps_to_milc(EigenM &Mvec)
{
  TIMER("Rotate eigen d,c");
  if(Mvec.size() == 0)return ;
  if(Mvec[0].size() %12 != 0){abort_r("Mvec size not correct! \n");}
  ga_matrices_milc  ga_milc;
  size_t Nvol = Mvec[0].size()/12;
  for(unsigned int iv=0;iv<Mvec.size();iv++){
    vecE_gamma(&Mvec[iv][0], ga_milc.ga[0][5],Nvol);
  }
}

void eigen_ov::touch_GPU(EigenV &Mres,long long offM,long long size, int mode )
{
  if(buffGPU==0 or offM <= 0)return;
  #ifdef QLAT_USE_ACC
  int gpu_id = -1;
  cudaGetDevice(&gpu_id);
  long long Total = size;

  if(offM >= Mres.size() or offM <= 0)return;
  if(Total == -1){Total = Mres.size() - offM;}
  if(Total + offM > Mres.size()){Total = Mres.size() - offM;}

  cudaMemAdvise(&Mres[offM], Total*sizeof(Complexq), cudaMemAdviseSetReadMostly, gpu_id);
  if(mode == 1){
  cudaMemPrefetchAsync(&Mres[offM], Total*sizeof(Complexq), gpu_id, cudaStreamLegacy);}

  ////cudaMemPrefetchAsync(&Mres[offM], Total*sizeof(Complexq), gpu_id, cudaStreamLegacy);
  //EigenV tem;tem.resize(Total);
  //qacc_forNB(i, Total, {tem[i] =  Mres[offM + i];});

  ////gpuErrchk(cudaMemcpy( &Mres[offM], &Mres[offM] , Total*sizeof(Complexq), cudaMemcpyHostToDevice));
  ////cudaMemAdvise(&Mres[offM], Total*sizeof(Complexq), cudaMemAdviseSetReadMostly, gpu_id);
  //if(mode == 1){
  //cudaMemPrefetchAsync(&Mres[offM], Total*sizeof(Complexq), gpu_id, cudaStreamLegacy);}
  //#else

  //long offM = mapN[nini*3 + 0]*2*bfac*ncut0*long(b_size);
  //long nloc = mapN[nini*3 + 2];
  //long size = 2*bfac*nloc*long(b_size);
  //if(full){size = Mres.size() - offM;}

  #endif
}

void eigen_ov::touch_GPU(EigenM &Mres,long long offM,long long size, int mode )
{
  if(buffGPU==0 or offM <= 0)return;
  #ifdef QLAT_USE_ACC
  size_t ini = 0;
  for(long i=0;i<Mres.size();i++)
  {
    long cur = Mres[i].size();
    if(size != -1 and ini >= size + offM)break;
    if(size != -1 and ini + cur > size + offM){cur = size + offM - ini;}
    if(ini >= offM){touch_GPU(Mres[i], 0, cur, mode);}
    if(ini < offM and ini+cur > offM){touch_GPU(Mres[i], offM-ini, cur-(offM-ini), mode);}
    ini += Mres[i].size();
  }
  #endif
}


void eigen_ov::untouch_GPU(EigenV &Mres)
{
  #ifdef QLAT_USE_ACC
  int gpu_id = -1;
  cudaGetDevice(&gpu_id);
  long offM = 0;
  size_t Total = Mres.size();
  cudaMemAdvise(&Mres[offM], Total*sizeof(Complexq), cudaMemAdviseUnsetReadMostly, gpu_id);
  #endif
}

void allocate_buf(std::vector<Complexq* > &buf, size_t n0, size_t n1)
{
  TIMER("CUDA Buf mem allocation");
  #ifdef QLAT_USE_ACC
  for(int i=0;i<buf.size();i++)if(buf[i]!=NULL)gpuErrchk(cudaFree(buf[i]));
  #else
  for(int i=0;i<buf.size();i++)if(buf[i]!=NULL)delete [] buf[i];
  #endif

  buf.resize(n0);
  for(int i=0;i<buf.size();i++){
    #ifdef QLAT_USE_ACC
    gpuErrchk(cudaMalloc(&buf[i], n1*sizeof(Complexq)));
    #else
    buf[i] = (Complexq *)malloc(n1*sizeof(Complexq));
    #endif
  }
}

void eigen_ov::copy_Mvec( int mode, bool random )
{
  TIMER("Rotate eigen n_vec");
  EigenM &Msrc = Mvec;
  EigenM &Mres = Mvec0;
  if(Msrc.size() == 0)return ;
  if(Msrc[0].size() %12 != 0){print0("Msrc size not correct! \n");qassert(false);}

  //////size_t Nvol = Msrc[0].size()/12;
  ////int n_vec = Msrc.size();

  #ifndef QLAT_USE_ACC
  ncutgpu = n_vec;
  ncutbuf = 0;
  #endif
  if(ncutgpu > n_vec or ncutgpu == -1){
    ncutgpu = n_vec;ncutbuf=0;
  }

  int nres = n_vec - ncutbuf;
  //mapN.resize(3*nres);
  //for(int i=0;i<nres/ncutgpu + 1;i++)
  //{
  //  for(int j=0;j<ncutgpu;j++)
  //  {
  //    int ni = i*ncutgpu + j;
  //    if(ni < nres)
  //    {
  //      mapN[ni*3+0] = i;
  //      mapN[ni*3+1] = j;
  //      if((i+1)*ncutgpu <= nres){
  //      mapN[ni*3+2] = ncutgpu;
  //      }else{
  //      mapN[ni*3+2] = nres - i*ncutgpu;
  //      }
  //    }
  //  }
  //}
  //////p_vector(mapN);

  long Ng = 0;
  if(ncutgpu!=0){Ng = nres/ncutgpu;if(nres%ncutgpu != 0)Ng += 1;}

  long La = Ng*2*bfac/bfac_group;
  long Lb = bfac_group*ncutgpu*long(b_size);
  resize_EigenM(Mres , La, Lb);

  //////size_t vsize = size_t(ncutbuf) * 2 * bfac * b_size;

  long Lae = 2*bfac/bfac_group;
  long Lbe = bfac_group*ncutbuf*long(b_size);
  if(ncutbuf != 0){allocate_buf(Eigenbuf, Lae, Lbe);}

  #if SUMMIT==1
  long Lax = 2*bfac/bfac_group;
  long Lbx = bfac_group*ncutgpu*long(b_size);
  if(ncutgpu != 0){allocate_buf(Eigendyn, Lax, Lbx);}
  #endif
  
  //#ifdef QLAT_USE_ACC
  //for(int i=0;i<Eigenbuf.size();i++)if(Eigenbuf[i]!=NULL)gpuErrchk(cudaFree(Eigenbuf[i]));
  //for(int i=0;i<Eigendyn.size();i++)if(Eigendyn[i]!=NULL)gpuErrchk(cudaFree(Eigendyn[i]));
  //#else
  //for(int i=0;i<Eigenbuf.size();i++)if(Eigenbuf[i]!=NULL)delete [] Eigenbuf[i];
  //for(int i=0;i<Eigendyn.size();i++)if(Eigendyn[i]!=NULL)delete [] Eigendyn[i];
  //#endif

  //Eigenbuf.resize(Lae);
  //for(int i=0;i<Eigenbuf.size();i++){
  //  #ifdef QLAT_USE_ACC
  //  gpuErrchk(cudaMalloc(&Eigenbuf[i], Lbe*sizeof(Complexq)));
  //  #else
  //  Eigenbuf[i] = (Complexq *)malloc(Lbe*sizeof(Complexq));
  //  #endif
  //}

  ///size_t total = n_vec - ncutbuf;
  ///total = total * 2 * bfac * b_size;
  ///Mres.resize(total);
  //////resize_EigenM(Mres , 1, total);

  print0("Current configuration, n_vec %d, Ng %ld, ncutgpu %d, ncutbuf %d, bfac %ld, b_size %ld, bfac_group %ld \n"
    , n_vec, Ng, ncutgpu, ncutbuf, bfac, b_size, bfac_group);
  fflush_MPI();

  /////#pragma omp parallel for
  if(random){
    if(ncutbuf != 0){
    EigenV tem;tem.resize(Lbe);
    for(int i=0;i<Eigenbuf.size();i++){
      ran_EigenM(tem,1);
      #ifdef QLAT_USE_ACC
      gpuErrchk(cudaMemcpy(  Eigenbuf[i], &tem[0]  , Lbe*sizeof(Complexq),cudaMemcpyHostToDevice));
      #else
      memcpy(      Eigenbuf[i], &tem[0]  , Lbe*sizeof(Complexq));
      #endif
    }
    }
    fflush_MPI();
    if(ncutgpu != 0)ran_EigenM(Mres,1);
  }

  if(!random)
  for(unsigned int ni=0;ni< n_vec;ni++)
  for(unsigned int chi=0;chi<2;chi++)
  for(LInt bi=0;bi<bfac;bi++)
  {
    /////long off1 = ((chi*bfac + bi)*n_vec + ni )*b_size + 0;
    ////////long Vone = 2*bfac*ncutbuf*b_size;
    //long off1 = 0;
    //if(ni >= ncutbuf)
    //{
    //  int na   = mapN[(ni-ncutbuf)*3+0];
    //  int nb   = mapN[(ni-ncutbuf)*3+1];
    //  int nloc = mapN[(ni-ncutbuf)*3+2];
    //  off1 = na*2*bfac*ncutgpu*b_size + ((chi*bfac + bi)*nloc + nb )*b_size + 0;
    //  tmp = &Mres[off1];
    //}else{
    //  off1 = ((chi*bfac + bi)*ncutbuf + ni )*b_size + 0;
    //  tmp = &Eigenbuf[off1];
    //}
    ////Mres[off1] = Msrc[ni][off0];
    /////size_t xi = bi*b_size + 0;
    //tmp = getEigenP(ni, (chi*bfac+bi)*b_size + 0);

    long off0 = (chi*bfac + bi)*b_size + 0;
    Complexq* tmp;
    tmp = getEigenP(ni, off0);
    EigenV& Mtem = Msrc[ni];
    qacc_forNB(bj, long(b_size), {
      //tmp[bj] = Msrc[ni][off0+bj];
      tmp[bj] = Mtem[off0+bj];
    });
  }
  qacc_barrier(dummy);

  //#ifdef QLAT_USE_ACC
  touch_GPU(Mres,0, La*Lb,0);

  //if(buffGPU==1){
  //touch_GPU(Mres);
  //qacc_barrier(dummy);}

  if(buffGPU==2){
  long temM_size = 2*bfac*ncutgpu*long(b_size);
  touch_GPU(Mres, 0, 2*temM_size);
  qacc_barrier(dummy);}


  //#endif
}

void eigen_ov::load_eigen(int icfg, std::string &ov_evecname, io_gwu &io,
  int checknorm, double kappa,double eigenerror)
{
  TIMER("=====Loading Eigen=====");
  char ename[500], enamev[500];
  sprintf(ename,ov_evecname.c_str(),icfg);
  sprintf(enamev,"%s.eigvals", ename);
  print0("Vector File name: %s \n", ename );
  print0("Values File name: %s \n", enamev);

  //////Load eigen values
  double rho_tem = 4 - 1.0/(2*kappa);
  load_eivals(enamev, rho_tem, eigenerror);


  print_mem_info("Before Eigen Memory Allocate");
  resize_EigenM(Mvec, n_vec, 2*6*noden);
  print_mem_info("Eigen Memory Allocate Done");

  //////Load eigen system
  //////shift of evectors n0, defaut 0
  load_gwu_eigen(ename, Mvec, io, n0, n0+n_vec);

  //////Check the norm of the vectors
  if(checknorm == 1){checknormM(Mvec);}

  /////Rotate d,c to outer vectors
  rotate_dc(Mvec);

  copy_Mvec(0 );

  //////rotate_cps_to_milc(Mvec);

  print_mem_info("Eigen Memory Load Done");
}

void eigen_ov::random_eigen()
{
  TIMER("=====Loading Eigen=====");
  eval_self.resize(n_vec);ran_EigenM(eval_self);
  
  //resize_EigenM(Mvec, n_vec, 12);ran_EigenM(Mvec);
  //resize_EigenM(Mvec, n_vec, 2*6*noden);
  //ran_EigenM(Mvec);
  resize_EigenM(Mvec, n_vec, 2*6);
  copy_Mvec(0, true );
  /////size_t total = n_vec;total = total * 2 * bfac * b_size;
  /////Mvec0.resize(total);ran_EigenM(Mvec0);

  print_mem_info("Eigen Memory Load Done");
}

//{
//int one_minus_halfD=1
//initiallize(mass,one_minus_halfD);
//}

//////Nv source number, 12 for prop
void eigen_ov::initiallize_mass(std::vector<double> &mass,int Ns,int one_minus_halfD_or)
{
  TIMER("Set up store memories");
  ////std::vector<double> mass = mass_or;
  int mN = mass.size();
  one_minus_halfD = one_minus_halfD_or;

  //if(alpha.size()     !=2*n_vec*Ns         )alpha.resize(2*n_vec*Ns);
  //if(alpha_bfac.size()!=2*n_vec*Ns*bfac    )alpha_bfac.resize(2*n_vec*Ns * bfac);
  int nlarge = ncutgpu;if(ncutbuf > ncutgpu)nlarge = ncutbuf;
  if(alpha.size()     !=2*nlarge*Ns         )alpha.resize(2*nlarge*Ns);
  if(alpha_bfac.size()!=2*nlarge*Ns*bfac    )alpha_bfac.resize(2*nlarge*Ns * bfac);

  if(alpha_list.size()!=2*nlarge * mN*Ns    )alpha_list.resize(2*nlarge * mN*Ns);
  /////if(srcM.size()      !=2*bfac*Ns*b_size   )srcM.resize(2*bfac * Ns*b_size);
  /////if(propM.size()     !=2*mN*bfac*Ns*b_size)propM.resize(2*mN*bfac * Ns*b_size);

  if(eval_list.size() !=mN*n_vec           )eval_list.resize(mN*n_vec);

  zeroE(alpha_bfac,0,false);zeroE(alpha,0,false);
  zeroE(alpha_list,0,false);
  /////zeroE(alpha_bfac,1);zeroE(alpha,1);

  size_t tmp_Ls = 2*bfac * Ns*b_size;
  size_t tmp_Lp = 2*mN*bfac * Ns*b_size;

  if(stmp_size != tmp_Ls){
    TIMER("CUDA prop mem allocation");
    stmp_size=tmp_Ls;
    #ifdef QLAT_USE_ACC
    if(stmp!=NULL){gpuErrchk(cudaFree(stmp));}
    gpuErrchk( cudaMalloc(&stmp, stmp_size*sizeof(Complexq)));
    #else
    if(stmp!=NULL){delete [] stmp;}
    stmp = (Complexq *)malloc(stmp_size*sizeof(Complexq));
    #endif
  }
  if(ptmp_size != tmp_Lp){
    TIMER("CUDA prop mem allocation");
    ptmp_size=tmp_Lp;
    #ifdef QLAT_USE_ACC
    if(ptmp!=NULL){gpuErrchk(cudaFree(ptmp));}
    gpuErrchk( cudaMalloc(&ptmp, ptmp_size*sizeof(Complexq)));
    #else
    if(ptmp!=NULL){delete [] ptmp;}
    ptmp = (Complexq *)malloc(ptmp_size*sizeof(Complexq));
    #endif
  }


  //zeroE(propM,0,false);

  //zeroE(alpha_list);
  //zeroE(srcM);

  #pragma omp parallel for
  for(int mki=0;mki< mN*n_vec;mki++)
  {
    int mi = mki/n_vec;
    int kn = mki%n_vec;
    eval_list[mi*n_vec + kn] = inv_self(eval_self[kn], mass[mi], rho, one_minus_halfD);
  }

  qacc_barrier(dummy);

  //for(int mi=0;mi<mN;mi++)
  //{
  //  #pragma omp parallel for
  //  for(int kn=0;kn<n_vec;kn++)
  //  {
  //    eval_list[mi*n_vec + kn] = inv_self(eval_self[kn], mass[mi], rho, one_minus_halfD);
  //  }
  //}

  //#ifdef QLAT_USE_ACC
  //alpha_bfac.resize(2*n_vec*Ns * b_size);zeroE(alpha_bfac);
  //#else
  //alpha_bfac.resize(2*n_vec*Ns * bfac);zeroE(alpha_bfac);
  //#endif
  /////resize_EigenM(alpha_list ,2*n_vec, mN*Ns);
  ////resize_EigenM(eval_list, mN, n_vec);

  ////print0("size %10d, %10d, %10d %10d \n",int(eval_list.size()), int(eval_list[0].size()), mN, n_vec);

  //alpha_list.resize(2*n_vec);
  //alpha_list.resize(2*n_vec,mN*Ns);

  //Esource.resize(2*bfac*Ns,b_size);zeroE(Esource);

  //alpha_list3pt.resize(2*Nt*n_vec*Ns);zeroE(alpha_list3pt);
  //alpha3pt.resize(n_vec*Ns*2*Nt);zeroE(alpha);

  //alpha3pt_ker_low.resize(n_vec*Ns*2*nt);zeroE(alpha3pt_ker_low);

}

void eigen_ov::initiallize_mass()
{
  std::vector<double> mass;mass.push_back(0.0);
  initiallize_mass(mass,12,one_minus_halfD);
}


void prop_L_device(eigen_ov& ei,Complexq *src,Complexq *props, int Ns,std::vector<double> &mass,int one_minus_halfD_or=1)
{
  TIMER_FLOPS("==prop_L");

  int mN   = mass.size();
  /////int Ns   = src.size()/(12*noden);

  long long Lat = ei.noden;
  long long vGb = Lat*12;
  int Fcount0 = 3 + 1;
  int Fcount1 = 3 + 1;
  long long Tfloat = ei.n_vec*Ns*mN*vGb*Fcount1 + ei.n_vec*Ns*vGb*Fcount0;
  timer.flops += Tfloat;
  double mem = Lat*12*(ei.n_vec + Ns + Ns*mN)*8.0;
  print0("Memory size %.3e GB, %.3e Gflop \n", 
    mem/(1024.0*1024*1024), Tfloat/(1024.0*1024*1024));
  ///////qlat::get_num_node()


  /////#ifdef QLAT_USE_ACC

  ei.initiallize_mass(mass, Ns, one_minus_halfD_or);
  ei.touch_GPU(ei.eval_list);
  //////auto& Mvec0       = ei.Mvec0;

  auto& alpha       = ei.alpha;
  auto& alpha_bfac  = ei.alpha_bfac;
  auto& alpha_list  = ei.alpha_list;
  auto& eval_list   = ei.eval_list;
  LInt& bfac        = ei.bfac;
  int& n_vec        = ei.n_vec;
  LInt& b_size      = ei.b_size;
  int& ncutgpu      = ei.ncutgpu;
  /////int& ncut0        = ei.ncut0;
  //////int& ncut1        = ei.ncut1;
  int& num_zero     = ei.num_zero;
  int& ncutbuf      = ei.ncutbuf;

  ///if(int(src.size()) != Ns*2*bfac*b_size){abort_r("src size wrong!");}
  ///if(int(props.size())  != mN*Ns*2*bfac*b_size){abort_r("prop size wrong!");}

  ////cudaMemcpy(     ei.stmp, &src[0]  , ei.stmp_size*sizeof(Complexq),cudaMemcpyDeviceToDevice);
  ////cudaMemcpyAsync(ei.ptmp, &props[0], ei.ptmp_size*sizeof(Complexq),cudaMemcpyDeviceToDevice);
  ////cudaMemsetAsync(ei.ptmp, 0, ei.ptmp_size*sizeof(Complexq));

  bool dummy_test = false;
  //bool dummy_test = true;
  //
  /////EigenV tem;tem.resize(ncutgpu*b_size);

  ////#if SUMMIT==1
  ////if(ncutgpu!=0){TIMER("SUMMIT COPY");
  ////long Lax = 2*bfac/ei.bfac_group;
  ////long Lbx = ei.bfac_group*ncutbuf*long(b_size);
  ////for(long a=0;a<Lax;a++){
  ////cudaMemcpy(&ei.Eigendyn[a][0], &ei.Mvec0[0 + a][0], Lbx*sizeof(Complexq),cudaMemcpyHostToDevice);}}
  ////#endif

  int Ng = 0;if(ncutgpu!=0){Ng = ei.n_vec/ncutgpu + 1;}
  int nini = 0;
  for(int ng=0;ng<Ng + 1;ng++)
  {
    ////int nini = ncutgpu*ng;
    if(nini >= ei.n_vec)break;
    /////long Vone = 2*bfac*ncutgpu*b_size;
    int  ncur = 0;
    size_t offM = 0;
    size_t temM_size = 2*bfac*ncutgpu*long(b_size);
    if(nini>=ncutbuf){
      offM = ((nini-ncutbuf)/ncutgpu)*2*bfac*ncutgpu*size_t(b_size);
      //////ei.mapN[(nini-ncutbuf)*3 + 0]*2*bfac*ncutgpu*long(b_size);
      ncur = ncutgpu;
      //////ei.mapN[(nini-ncutbuf)*3 + 2];
    }else{ncur = ncutbuf;}
    //#if SUMMIT==1
    //if(ncutbuf==0)qacc_barrier(dummy);
    //#endif
    #if SUMMIT==1
    if(nini>=ncutbuf)
    {TIMER("SUMMIT COPY");
    long Lax = 2*bfac/ei.bfac_group;
    long Lbx = ei.bfac_group*ncutgpu*long(b_size);
    long noff = offM/Lbx;
    for(long a=0;a<Lax;a++){
    cudaMemcpy(&ei.Eigendyn[a][0], &ei.Mvec0[noff + a][0], Lbx*sizeof(Complexq),cudaMemcpyHostToDevice);}
    qacc_barrier(dummy);}
    #endif

    {
    long m = ncur;
    long n = Ns;
    long w = b_size;

    TIMER_FLOPS("vec reduce");
    long long vGb = 2*bfac*m*n*w;
    int Fcount0 = 3 + 1;  
    timer.flops += vGb*Fcount0;

    //Complexq* rp = &alpha_bfac[(bini + 0)*n_vec*Ns     + 0];
    //Complexq* Ep = &Mvec0[     (bini + 0)*n_vec*b_size + 0];
    //Complexq* sp = &src[      (bini + 0)*Ns*b_size    + 0];

    //matrix_prod_gpu(Ep, sp, rp, m,n, w , bcut, true, dummy_test);

    //print0("current m %d, n %d, w %d , 2*bfac %d .\n", m, n, w, 2*bfac);
    //fflush_MPI();
    ///{
    ///MPI_Barrier(get_comm());
    ///fflush(stdout);}

    //unsigned long bini = 0;
    //Complexq* rp = &alpha_bfac[(bini + 0)*m*n + 0];
    //Complexq* Ep = NULL;
    //if(nini>=ncutbuf){Ep=&Mvec0[offM     +(bini + 0)*m*w + 0];}
    //else{Ep =&ei.Eigenbuf[                (bini + 0)*m*w + 0];}
    //Complexq* sp = &src[       (bini + 0)*n*w + 0];
    //matrix_prod_gpu(Ep, sp, rp, m,n, w , 2*bfac, true, dummy_test);

    //for(long bini=0;bini<2*bfac;bini++)
    //{
    //  Complexq* rp = &alpha_bfac[(bini + 0)*m*n + 0];
    //  Complexq* sp = &src[       (bini + 0)*n*w + 0];

    //  Complexq* Ep = NULL;
    //  if(nini>=ncutbuf){Ep=&Mvec0[offM     +(bini + 0)*m*w + 0];}
    //  else{Ep =&ei.Eigenbuf[                (bini + 0)*m*w + 0];}
    //  matrix_prod_gpu(Ep, sp, rp, m,n, w , 1, true, dummy_test);
    //}

    unsigned long bini = 0;
    int bcut = ei.bfac_group;
    for(int bcg=0;bcg<(2*bfac)/bcut + 1; bcg++)
    {
      if(bini >= 2*bfac)break;
      if(bini + bcut > 2*bfac){bcut = 2*bfac - bini;}
      Complexq* rp = &alpha_bfac[(bini + 0)*m*n + 0];
      ////Complexq* Ep = &Mvec0[offM+(bini + 0)*m*w + 0];
      Complexq* Ep = ei.getEigenP(nini, bini*b_size);
      //Complexq* Ep = &ei.Mvec0[0][0];
      Complexq* sp = &src[       (bini + 0)*n*w + 0];
      ///Complexq* sp = &ei.stmp[       (bini + 0)*n*w + 0];

      ////print0("==m %d, n %d, w %d, bini %d, bcut %d, 2*bfac %d \n", m, n, w, bini, bcut, 2*bfac);
      ////fflush_MPI();
      matrix_prod_gpu(Ep, sp, rp, m,n, w , bcut, true, dummy_test);
      bini += bcut;
    }
    /////Output not conflict, GPU end here
    qacc_barrier(dummy);

    //if(ei.buffGPU == 2){ei.touch_GPU(src);ei.touch_GPU(Mvec0, offM, 2*temM_size);}
    //if(ei.buffGPU == 2){ei.touch_GPU(src);ei.touch_GPU(Mvec0, offM+temM_size, temM_size);}
    //#if SUMMIT==0
    //{TIMER("Load memory");
    //if(ei.buffGPU == 2){ei.touch_GPU(ei.Mvec0, offM + temM_size, temM_size, 2);}
    //if(dummy_test == true)qacc_barrier(dummy);}
    //#endif

    }

    {
    TIMER("Reduce alpha");
    //#pragma omp parallel for
    //for(long coff=0;coff< 2*ncur*Ns;coff++)
    //{
    //  int chi = coff/(ncur*Ns);
    //  long xi = coff%(ncur*Ns);
    //  for(int bi=0;bi<bfac;bi++){
    //    alpha[coff] += alpha_bfac[chi*bfac*ncur*Ns + bi*ncur*Ns + xi];
    //  }
    //}
    ///qthread_for
    qacc_for(coff, long(2*ncur*Ns),{
      int chi = coff/(ncur*Ns);
      long xi = coff%(ncur*Ns);
      for(int bi=0;bi<bfac;bi++){
        alpha[coff] += alpha_bfac[chi*bfac*ncur*Ns + bi*ncur*Ns + xi];
      }
    });
    }

    /////reduce_vec(&alpha_bfac[0], &alpha[0], bfac, 2*n_vec*Ns);

    {
    TIMER("Global sum");
    sum_all_size(reinterpret_cast<Ftype* > (&alpha[0]),2*(2*ncur*Ns));
    }

    {
    TIMER("Get alpha list")
    Complexq Iimag(0.0,1.0);
    Complexq Two2(2.0,0.0);
    //#pragma omp parallel for
    //for(long coff=0;coff< 2*mN*Ns*ncur;coff++)
    //{
    //  int chi   =  coff/(mN*Ns*ncur);
    //  int mi    = (coff/(Ns*ncur  ))%mN;
    //  int is    = (coff/(ncur     ))%Ns;
    //  int kn    = (coff             )%ncur;
    //  Complexq li = eval_list[mi*ei.n_vec + kn + nini];
    //  //unsigned long offA = ((chi*n_vec+kn)*mN+mi)*Ns+is;
    //  unsigned long offA = ((chi*mN+mi)*Ns+is)*ncur + kn;

    //  if(kn+nini >= num_zero){
    //    alpha_list[offA] = Two2*(li.real()*alpha[(chi*ncur+kn)*Ns+is]+Iimag*li.imag()*alpha[((1-chi)*ncur+kn)*Ns+is]);}

    //  if(kn+nini < num_zero){
    //    alpha_list[offA] = li*alpha[(chi*ncur+kn)*Ns+is];}
    //}

    qacc_for(coff, long(2*mN*Ns*ncur),{
      int chi   =  coff/(mN*Ns*ncur);
      int mi    = (coff/(Ns*ncur  ))%mN;
      int is    = (coff/(ncur     ))%Ns;
      int kn    = (coff             )%ncur;
      unsigned long offA = ((chi*mN+mi)*Ns+is)*ncur + kn;

      if(kn + nini >= n_vec){alpha_list[offA] = 0.0;}else{
      //unsigned long offA = ((chi*n_vec+kn)*mN+mi)*Ns+is;
      Complexq li = eval_list[mi*n_vec + kn + nini];
      if(kn+nini >= num_zero){
        alpha_list[offA] = Two2*(li.real()*alpha[(chi*ncur+kn)*Ns+is]+Iimag*li.imag()*alpha[((1-chi)*ncur+kn)*Ns+is]);}
      if(kn+nini  < num_zero){
        alpha_list[offA] = li*alpha[(chi*ncur+kn)*Ns+is];}
      }
    });
    }


    {
    //qacc_for(coff, long(2*mN*Ns*bfac),{});
    long m = mN*Ns;
    long n = b_size;
    long w = ncur;

    //TIMER("vec multi");
    TIMER_FLOPS("vec multi");
    long long vGb = 2*bfac*m*n*w;
    int Fcount0 = 3 + 1;  
    timer.flops += vGb*Fcount0;

    if(ncur<ei.n_vec){zeroE(alpha_bfac,0,false);zeroE(alpha,0,false);}

    for(int coff=0;coff<2*bfac;coff++)
    {
      long chi = coff/bfac;
      long bi  = coff%bfac;

      Complexq* rp = &props[(chi*bfac+bi)*mN*Ns*b_size + 0];
      ////Complexq* rp = &ei.ptmp[(chi*bfac+bi)*mN*Ns*b_size + 0];
      Complexq* ap = &alpha_list[(chi*mN*Ns+0)*w + 0];
      /////Complexq* Ep = &Mvec0[offM + (chi*bfac + bi)*w*b_size + 0];

      Complexq* Ep = ei.getEigenP(nini, coff*b_size);

      //Complexq* Ep = NULL;
      //if(nini>=ncutbuf){Ep=&Mvec0[offM + (chi*bfac + bi)*w*b_size + 0];}
      //else{Ep = &ei.Eigenbuf[            (chi*bfac + bi)*w*b_size + 0];}

      ////Complexq* Ep = &ei.Eigenbuf[     (chi*bfac + bi)*w*b_size + 0];
      ////if(ng!=0){Ep=&Mvec0[offM -Vone + (chi*bfac + bi)*w*b_size + 0];}

      matrix_prod_gpu(ap, Ep, rp, m, n, w ,1, false,  dummy_test, true);
    }
    /////Output not conflict, GPU end here
    qacc_barrier(dummy);
    ///nini += ncur;

    ////qacc_for(isp, long(ei.ptmp_size),{
    ////  props[isp] = ei.ptmp[isp];
    ////});

    /////cudaMemcpy(&props[0], ei.ptmp, ei.ptmp_size*sizeof(Complexq),cudaMemcpyDeviceToDevice);
    }
    /////qacc_barrier(dummy);
    nini += ncur;
  }


  ei.untouch_GPU(eval_list);
  //////#endif

}

//void eigen_ov::prop_L(EigenV &src,EigenV &props,std::vector<double> &mass,int one_minus_halfD_or)
//{
//  TIMER_FLOPS("==prop_L");
//
//  int mN   = mass.size();
//  int Ns   = src.size()/(12*noden);
//
//  long long Lat = noden;
//  long long vGb = Lat*12;
//  int Fcount0 = 3 + 1;
//  int Fcount1 = 3 + 1;
//  long long Tfloat = n_vec*Ns*mN*vGb*Fcount1 + n_vec*Ns*vGb*Fcount0;
//  timer.flops += Tfloat;
//  double mem = Lat*12*(n_vec + Ns + Ns*mN)*8.0;
//  print0("Memory size %.3e GB, %.3e Gflop \n", 
//    mem/(1024.0*1024*1024), Tfloat/(1024.0*1024*1024));
//  ///////qlat::get_num_node()
//
//  ////Eigen::setNbThreads(1);
//  #ifdef QLAT_USE_ACC
//  prop_L_GPU(*this, src, props, mass, one_minus_halfD_or);
//  #else
//
//  //initiallize_mass(mass, Ns, one_minus_halfD_or);
//
//  //if(int(src[0].size()) != 2*bfac*b_size){abort_r("src size wrong!");}
//  //if(int(props.size())  != mN*Ns){abort_r("prop size wrong!");}
//  //if(int(props[0].size()) != 2*bfac*b_size){abort_r("prop vec wrong!");}
//
//  ////EigenV srcM;srcM.resize(2*bfac * Ns*b_size);zeroE(srcM);
//  ////EigenV propM;propM.resize(mN*2*bfac * Ns*b_size);zeroE(propM);
//
//  ////auto& alpha = this -> alpha;
//  ////auto& alpha_bfac = this->alpha_bfac;
//  ////auto& Mvec0 = this->Mvec0;
//  ////auto& alpha_list = this->alpha_list;
//  ////auto& eval_list = this->eval_list;
//  ////int& bfac  = this->bfac;
//  ////int& n_vec = this->n_vec;
//  ////int& b_size = this->b_size;
//  ////int& ncut0 = this->ncut0;
//  ////int& ncut1 = this->ncut1;
//  ////int& num_zero = this->num_zero;
//  ////eigen_ov& f0 = *this;
//
//  ///////Rowmajor, M(3,4) continuous in 4
//  ///////Colmajor, M(3,4) continuous in 3
//  ///////src, Ns, d,c --> vol
//  ///////Copy source to matrix
//
//  //int mode_E = 1;
//
//  //if(mode_E == 1)
//  //{
//  ////resize_EigenM(srcM, 2*bfac * Ns*b_size);
//  ////long Lv = Ns*b_size;
//  //{
//  //TIMER("Copy source");
//  ////qacc_for(coff, long(2*bfac*Ns*b_size),{
//  ////  int chi =  coff/(bfac*Ns*b_size);
//  ////  int bi  = (coff/(Ns*b_size)    )%bfac;
//  ////  int is  = (coff/b_size         )%Ns;
//  ////  int bj  = (coff                )%b_size;
//  ////  srcM[(chi*bfac+bi)*Ns*b_size + is*b_size + bj] = src[is][(chi*bfac+bi)*b_size + bj];
//  ////});
//
//  ////////CPU copy of source vectors
//  //#pragma omp parallel for
//  //for(long coff=0;coff < 2*bfac*Ns; coff++)
//  //{
//  //  long chi =  coff/(bfac*Ns       );
//  //  long bi  = (coff/(Ns       )    )%bfac;
//  //  long is  = (coff                )%Ns;
//
//  //  //srcM[chi*bfac+bi][bj*Ns + is] = src[is][(chi*bfac+bi)*b_size + bj];
//  //  for(int bj=0;bj<b_size;bj++){
//  //    srcM[(chi*bfac+bi)*Ns*b_size  + is*b_size + bj] = src[is][(chi*bfac+bi)*b_size + bj];}
//  //  //srcM[(chi*bfac+bi)*Lv+ is*b_size + bj] = src[is][(chi*bfac+bi)*b_size + bj];
//  //}
//
//  ////////GPU copy of source vectors
//  ////qacc_for(coff, long(2*bfac*Ns),{
//  ////  int chi =  coff/(bfac*Ns       );
//  ////  int bi  = (coff/(Ns       )    )%bfac;
//  ////  int is  = (coff                )%Ns;
//  ////  //srcM[chi*bfac+bi][bj*Ns + is] = src[is][(chi*bfac+bi)*b_size + bj];
//  ////  for(int bj=0;bj<b_size;bj++){
//  ////    srcM[(chi*bfac+bi)*Ns*b_size  + is*b_size + bj] = src[is][(chi*bfac+bi)*b_size + bj];}
//  ////  //srcM[(chi*bfac+bi)*Lv+ is*b_size + bj] = src[is][(chi*bfac+bi)*b_size + bj];
//  ////});
//
//  //}
//
//  //{
//  //TIMER("vec reduce");
//
//  //{
//  //int Ng = n_vec/ncut0 + 1;
//  //qacc_for(coff, long(2*Ng*bfac),{
//  //  long chi =  coff/(Ng*bfac);
//  //  long ng  = (coff/bfac    )%Ng;
//  //  long bi  = (coff         )%bfac;
//
//  //  long ncur = ncut0;
//  //  long nini = ncut0*ng;
//  //  if(nini + ncur > n_vec){ncur = n_vec - nini;}
//  //  if(ncur > 0){
//  //    Complexq* rp = &alpha_bfac[chi*bfac*n_vec*Ns + bi*n_vec*Ns + nini*Ns];
//  //    Complexq* Ep = &Mvec0[(chi*bfac + bi)*n_vec*b_size + nini*b_size];
//  //    Complexq* sp = &srcM[(chi*bfac+bi)*Ns*b_size + 0];
//
//  //    EM re(rp, ncur, Ns); EM Ev(Ep, ncur, b_size); EMC sv(sp, b_size, Ns);
//  //    re += Ev.conjugate() * sv;
//
//  //    //for(int si=0;si<Ns;si++)
//  //    //for(int ni=0;ni<ncur;ni++)
//  //    //for(int bm=0;bm<b_size;bm++)
//  //    //{
//  //    //  rp[ni*Ns + si] += qconj(Ep[ni*b_size+bm]) * sp[si*b_size+bm];
//  //    //}
//
//
//  //  }
//  //});
//  //}
//
//  //qacc_for(coff, long(2*n_vec*Ns),{
//  //  int chi = coff/(n_vec*Ns);
//  //  long xi = coff%(n_vec*Ns);
//  //  for(int bi=0;bi<bfac;bi++){
//  //    alpha[coff] += alpha_bfac[chi*bfac*n_vec*Ns + bi*n_vec*Ns + xi];
//  //  }
//  //});
//
//  ///////reduce_vec(&alpha_bfac[0], &alpha[0], bfac, 2*n_vec*Ns);
//  //}
//  //}
//
//  //{
//  //TIMER("Global sum");
//  //sum_all_size(reinterpret_cast<Ftype* > (&alpha[0]),2*alpha.size());
//  //}
//
//  ////double check_As = 0.0;
//  ////for(int ai = 0;ai<alpha.size();ai++)
//  ////{
//  ////  //int kn = ai/(2*12);
//  ////  //int chi = (ai%(2*12))/12;
//  ////  //int k = ai%12;
//  ////  //print0("alpha, kn %5d, chi %1d, k %3d, real %13.8f imag %13.8f \n", kn, chi, k,  alpha[ai].real(), alpha[ai].imag());
//  ////  check_As += alpha[ai].real();
//  ////  check_As += alpha[ai].imag();
//  ////}
//  ////print0("alpha %13.8f \n",check_As);
//
//  //{
//  //TIMER("Get alpha list")
//  //Complexq Iimag(0.0,1.0);
//  //Complexq Two2(2.0,0.0);
//  //qacc_for(coff, long(2*mN*Ns*n_vec),{
//  //  long chi   =  coff/(mN*Ns*n_vec);
//  //  long mi    = (coff/(Ns*n_vec  ))%mN;
//  //  long is    = (coff/(n_vec     ))%Ns;
//  //  long kn    = (coff             )%n_vec;
//  //  Complexq li = eval_list[mi*n_vec + kn];
//  //  unsigned long offA = ((chi*n_vec+kn)*mN+mi)*Ns+is;
//
//  //  if(kn >= num_zero){
//  //    alpha_list[offA] = Two2*(li.real()*alpha[(chi*n_vec+kn)*Ns+is]+Iimag*li.imag()*alpha[((1-chi)*n_vec+kn)*Ns+is]);}
//
//  //  if(kn < num_zero){
//  //    alpha_list[offA] = li*alpha[(chi*n_vec+kn)*Ns+is];}
//  //});
//  //}
//
//  //{
//  //TIMER("vec multi");
//  ////qacc_for(coff, long(2*mN*Ns*bfac),{});
//
//  //long ncur = ncut1;
//  //long nini = 0;
//  //for(int ng=0;ng<n_vec/ncut1 + 1;ng++)
//  //{
//  //  if(nini == n_vec){break;}
//  //  if(nini + ncur > n_vec){ncur = n_vec - nini;}
//
//  //  qacc_for(coff, long(2*bfac),{
//  //    long chi = coff/bfac;
//  //    long bi  = coff%bfac;
//
//  //    Complexq* rp = &propM[(chi*bfac+bi)*mN*Ns*b_size];
//  //    Complexq* ap = &alpha_list[(chi*n_vec + nini)*mN*Ns];
//  //    Complexq* Ep = &Mvec0[(chi*bfac + bi)*n_vec*b_size + nini*b_size];
//
//  //    EM  re(rp, mN*Ns , b_size);
//  //    EMC al(ap, mN*Ns, ncur);
//  //    EM  Ev(Ep, ncur, b_size);
//  //    re += al * Ev;
//
//  //    //for(int mi=0;mi<mN*Ns;mi++)
//  //    //for(int bi=0;bi<b_size;bi++)
//  //    //for(int ni=0;ni<ncur;ni++)
//  //    //{
//  //    //  rp[mi*b_size + bi] += ap[ni*mN*Ns + mi] * Ep[ni*b_size + bi];
//  //    //}
//
//  //  });
//  //  nini += ncur;
//  //}
//
//  //qacc_for(coff, long(2*mN*Ns*bfac),{
//  //  int chi = coff/(mN*Ns*bfac);
//  //  int off = coff%(mN*Ns*bfac);
//  //  int mi  =  off/(Ns*bfac);
//  //  int si  = (off/bfac    )%(Ns);
//  //  int bi  = (off         )%bfac;
//
//  //  EM r0(&props[(mi*Ns+ si)][(chi*bfac+bi)*b_size + 0], 1, b_size);
//  //  EM r1(&propM[(chi*bfac+bi)*mN*Ns*b_size + (mi*Ns + si)*b_size], 1 , b_size);
//  //  r0 += r1;
//  //});
//
//  ////#pragma omp parallel for
//  ////for(long coff=0; coff < long(2*mN*Ns*bfac); coff++)
//  ////{
//  ////  long chi = coff/(mN*Ns*bfac);
//  ////  long off = coff%(mN*Ns*bfac);
//  ////  long mi  =  off/(Ns*bfac);
//  ////  long si  = (off/bfac    )%(Ns);
//  ////  long bi  = (off         )%bfac;
//
//  ////  EM r0(&props[(mi*Ns+ si)][(chi*bfac+bi)*b_size + 0], 1, b_size);
//  ////  EM r1(&propM[(chi*bfac+bi)*mN*Ns*b_size + (mi*Ns + si)*b_size], 1 , b_size);
//  ////  r0 += r1;
//  ////}
//
//
//  //}
//
//  #endif
//
//  ///Eigen::setNbThreads(0);
//  ////return ;
//}

//void prop_L_GPU(eigen_ov& ei,EigenV &src,EigenV &props,std::vector<double> &mass,int one_minus_halfD_or)
//{
//  return ;
//
//  //#ifdef QLAT_USE_ACC
//
//  //int mN   = mass.size();
//  //int Ns   = src.size()/(12*ei.noden);
//  //ei.initiallize_mass(mass, Ns, one_minus_halfD_or);
//
//  //auto& alpha       = ei.alpha;
//  //auto& alpha_bfac  = ei.alpha_bfac;
//  //auto& Mvec0       = ei.Mvec0;
//  //auto& alpha_list  = ei.alpha_list;
//  //auto& eval_list   = ei.eval_list;
//  //int& bfac         = ei.bfac;
//  ////int& n_vec        = ei.n_vec;
//  //int& b_size       = ei.b_size;
//  //int& ncutgpu      = ei.ncutgpu;
//  ///////int& ncut0        = ei.ncut0;
//  ////////int& ncut1        = ei.ncut1;
//  //int& num_zero     = ei.num_zero;
//
//  //if(int(src.size()) != Ns*2*bfac*b_size){abort_r("src size wrong!");}
//  //if(int(props.size())  != mN*Ns*2*bfac*b_size){abort_r("prop size wrong!");}
//
//  //cudaMemcpy(     ei.stmp, &src[0]  , ei.stmp_size*sizeof(Complexq),cudaMemcpyDeviceToDevice);
//  //cudaMemcpyAsync(ei.ptmp, &props[0], ei.ptmp_size*sizeof(Complexq),cudaMemcpyDeviceToDevice);
//  //////cudaMemsetAsync(ei.ptmp, 0, ei.ptmp_size*sizeof(Complexq));
//
//  //bool dummy_test = false;
//  //////bool dummy_test = true;
//
//  //int Ng = ei.n_vec/ncutgpu + 1;
//  //for(int ng=0;ng<Ng;ng++)
//  //{
//  //  int nini = ncutgpu*ng;
//  //  if(nini >= ei.n_vec)break;
//  //  int  ncur = ei.mapN[nini*3 + 2];
//  //  size_t offM = ei.mapN[nini*3 + 0]*2*bfac*ncutgpu*long(b_size);
//  //  long temM_size = 2*bfac*ncutgpu*long(b_size);
//
//  //  {
//  //  long m = ncur;
//  //  long n = Ns;
//  //  long w = b_size;
//
//  //  TIMER_FLOPS("vec reduce");
//  //  long long vGb = 2*bfac*m*Ns*b_size;
//  //  int Fcount0 = 3 + 1;  
//  //  timer.flops += vGb*Fcount0;
//
//  //  //Complexq* rp = &alpha_bfac[(bini + 0)*n_vec*Ns     + 0];
//  //  //Complexq* Ep = &Mvec0[     (bini + 0)*n_vec*b_size + 0];
//  //  //Complexq* sp = &src[      (bini + 0)*Ns*b_size    + 0];
//
//  //  //matrix_prod_gpu(Ep, sp, rp, m,n, w , bcut, true, dummy_test);
//
//  //  unsigned long bini = 0;
//  //  int bcut = ei.bfac_group;
//  //  for(int bcg=0;bcg<(2*bfac)/bcut + 1; bcg++)
//  //  {
//  //    if(bini >= 2*bfac)break;
//  //    if(bini + bcut > 2*bfac){bcut = 2*bfac - bini;}
//  //    Complexq* rp = &alpha_bfac[(bini + 0)*m*n + 0];
//  //    Complexq* Ep = &Mvec0[offM+(bini + 0)*m*w + 0];
//  //    //Complexq* sp = &src[       (bini + 0)*n*w + 0];
//  //    Complexq* sp = &ei.stmp[       (bini + 0)*n*w + 0];
//
//  //    //print0("==m %d, n %d, w %d, bini %d, bcut %d, 2*bfac %d \n", m, n, w, bini, bcut, 2*bfac);
//  //    //fflush_MPI();
//  //    matrix_prod_gpu(Ep, sp, rp, m,n, w , bcut, true, dummy_test);
//  //    bini += bcut;
//  //  }
//  //  ///////Output not conflict, GPU end here
//  //  qacc_barrier(dummy);
//
//  //  //if(ei.buffGPU == 2){ei.touch_GPU(src);ei.touch_GPU(Mvec0, offM, 2*temM_size);}
//  //  //if(ei.buffGPU == 2){ei.touch_GPU(src);ei.touch_GPU(Mvec0, offM+temM_size, temM_size);}
//  //  if(ei.buffGPU == 2){ei.touch_GPU(Mvec0, offM+temM_size, temM_size, 2);}
//
//  //  }
//
//  //  {
//  //  TIMER("Reduce alpha");
//  //  //#pragma omp parallel for
//  //  //for(long coff=0;coff< 2*ncur*Ns;coff++)
//  //  //{
//  //  //  int chi = coff/(ncur*Ns);
//  //  //  long xi = coff%(ncur*Ns);
//  //  //  for(int bi=0;bi<bfac;bi++){
//  //  //    alpha[coff] += alpha_bfac[chi*bfac*ncur*Ns + bi*ncur*Ns + xi];
//  //  //  }
//  //  //}
//
//  //  qthread_for(coff, long(2*ncur*Ns),{
//  //    int chi = coff/(ncur*Ns);
//  //    long xi = coff%(ncur*Ns);
//  //    for(int bi=0;bi<bfac;bi++){
//  //      alpha[coff] += alpha_bfac[chi*bfac*ncur*Ns + bi*ncur*Ns + xi];
//  //    }
//  //  });
//  //  }
//
//  //  /////reduce_vec(&alpha_bfac[0], &alpha[0], bfac, 2*n_vec*Ns);
//
//  //  {
//  //  TIMER("Global sum");
//  //  sum_all_size(reinterpret_cast<Ftype* > (&alpha[0]),2*alpha.size());
//  //  }
//
//  //  {
//  //  TIMER("Get alpha list")
//  //  Complexq Iimag(0.0,1.0);
//  //  Complexq Two2(2.0,0.0);
//  //  //#pragma omp parallel for
//  //  //for(long coff=0;coff< 2*mN*Ns*ncur;coff++)
//  //  //{
//  //  //  int chi   =  coff/(mN*Ns*ncur);
//  //  //  int mi    = (coff/(Ns*ncur  ))%mN;
//  //  //  int is    = (coff/(ncur     ))%Ns;
//  //  //  int kn    = (coff             )%ncur;
//  //  //  Complexq li = eval_list[mi*ei.n_vec + kn + nini];
//  //  //  //unsigned long offA = ((chi*n_vec+kn)*mN+mi)*Ns+is;
//  //  //  unsigned long offA = ((chi*mN+mi)*Ns+is)*ncur + kn;
//
//  //  //  if(kn+nini >= num_zero){
//  //  //    alpha_list[offA] = Two2*(li.real()*alpha[(chi*ncur+kn)*Ns+is]+Iimag*li.imag()*alpha[((1-chi)*ncur+kn)*Ns+is]);}
//
//  //  //  if(kn+nini < num_zero){
//  //  //    alpha_list[offA] = li*alpha[(chi*ncur+kn)*Ns+is];}
//  //  //}
//
//  //  qthread_for(coff, long(2*mN*Ns*ncur),{
//  //    int chi   =  coff/(mN*Ns*ncur);
//  //    int mi    = (coff/(Ns*ncur  ))%mN;
//  //    int is    = (coff/(ncur     ))%Ns;
//  //    int kn    = (coff             )%ncur;
//  //    Complexq li = eval_list[mi*ei.n_vec + kn + nini];
//  //    //unsigned long offA = ((chi*n_vec+kn)*mN+mi)*Ns+is;
//  //    unsigned long offA = ((chi*mN+mi)*Ns+is)*ncur + kn;
//
//  //    if(kn+nini >= num_zero){
//  //      alpha_list[offA] = Two2*(li.real()*alpha[(chi*ncur+kn)*Ns+is]+Iimag*li.imag()*alpha[((1-chi)*ncur+kn)*Ns+is]);}
//
//  //    if(kn+nini < num_zero){
//  //      alpha_list[offA] = li*alpha[(chi*ncur+kn)*Ns+is];}
//  //  });
//  //  }
//
//
//  //  {
//  //  //qacc_for(coff, long(2*mN*Ns*bfac),{});
//  //  long m = mN*Ns;
//  //  long n = b_size;
//  //  long w = ncur;
//
//  //  //TIMER("vec multi");
//  //  TIMER_FLOPS("vec multi");
//  //  long long vGb = 2*bfac*m*n*w;
//  //  int Fcount0 = 3 + 1;  
//  //  timer.flops += vGb*Fcount0;
//
//  //  if(ncur<ei.n_vec){zeroE(alpha_bfac,0,false);zeroE(alpha,0,false);}
//
//
//  //  for(int coff=0;coff<2*bfac;coff++)
//  //  {
//  //    long chi = coff/bfac;
//  //    long bi  = coff%bfac;
//
//  //    //Complexq* rp = &props[(chi*bfac+bi)*mN*Ns*b_size + 0];
//  //    Complexq* rp = &ei.ptmp[(chi*bfac+bi)*mN*Ns*b_size + 0];
//  //    Complexq* ap = &alpha_list[(chi*mN*Ns+0)*w + 0];
//  //    Complexq* Ep = &Mvec0[offM + (chi*bfac + bi)*w*b_size + 0];
//  //    matrix_prod_gpu(ap, Ep, rp, m, n, w ,1, false,  dummy_test, true);
//  //  }
//  //  /////Output not conflict, GPU end here
//  //  qacc_barrier(dummy);
//  //  ///nini += ncur;
//
//  //  ////qacc_for(isp, long(ei.ptmp_size),{
//  //  ////  props[isp] = ei.ptmp[isp];
//  //  ////});
//
//  //  cudaMemcpy(&props[0], ei.ptmp, ei.ptmp_size*sizeof(Complexq),cudaMemcpyDeviceToDevice);
//  //  }
//  //  /////qacc_barrier(dummy);
//  //}
//
//  ////////ei.untouch_GPU(src);
//
//  //#endif
//}

//LInt eigen_ov::get_threads(LInt thread, LInt Total)
//{
//  for(LInt temb=thread;temb<Total;temb++)
//  {
//    if(Total%temb == 0)
//    {
//      return temb;
//    }
//  }
//  return 1;
//
//  //for(int temb=thread;temb<b_size;temb++)
//  //{
//  //  if(b_size%temb == 0)
//  //  {
//  //    return temb;
//  //  }
//  //}
//  return 32;
//}


}


#endif

