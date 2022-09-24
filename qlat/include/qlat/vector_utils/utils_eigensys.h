// utils_eigensys.h
// Gen Wang
// Sep. 2021

#ifndef UTILS_EIGENSYS_H
#define UTILS_EIGENSYS_H

#pragma once
#include "utils_gammas.h"
#include "utils_fft_desc.h"
#include "utils_Matrix_prod.h"
#include "utils_construction.h"
#include "check_fun.h"
#include "utils_smear_vecs.h"

///#define SUMMIT 0

#define Elocal std::vector<std::vector<Complexq > >
////#define Vlocal qlat::vector_acc<Complexq >
#define Vlocal qlat::vector_gpu<Complexq >
#define EIGENERROR 1e-11

namespace qlat{

struct eigen_ov {
  qlat::fft_desc_basic* fdp;

  Elocal Mvec;     //  nvec --> 2*bfac --> b_size*6
  Elocal Mvec_Sm;  //  Smeared eigensystem 

  bool enable_smearE;

  ////int Edouble;
  ////int smear_Eigen;

  ///move_index mv_civ;

  /////lattice copied from fd
  int nx;int ny;int nz;int nt;
  int Nx;int Ny;int Nz;int Nt;
  LInt noden;

  ////Eigensystem information
  int n_vec;
  int num_zero;
  int one_minus_halfD;
  double rho;
  double Eeigenerror;
  EigenV eval_self;

  /////Group the inner vector for vectorizations
  long b_size, bfac;

  bool gpu_mem_set;
  int ncutgpu, ncutbuf;
  long bfac_group;
  long BFAC_GROUP_CPU;


  int nV_prop;
  double extra_mem_factor;

  ///Vlocal stmp;
  ///Vlocal ptmp;

  std::vector< Vlocal > Eigenbuf;
  std::vector< Vlocal > Eigendyn;
  std::vector< Vlocal > Eigenbuf_Sm;
  std::vector< Vlocal > Eigendyn_Sm;
  int npos_Eigenbuf;int npos_Eigendyn;

  Vlocal alpha;
  Vlocal alpha_bfac;
  Vlocal alpha_list;
  Vlocal eval_list;
  std::vector<double > massL;

  qlat::vector_acc<Complexq > eval_tem;
  /////int ncut0,ncut1;

  /////int buffGPU;

  ////3pt function
  //EigenV alpha_list3pt;
  //EigenV alpha_list3ptC;
  //EigenV alpha3pt;
  //EigenV alpha3pt_ker_low;

  /////Construction and memory allocations
  eigen_ov(qlat::fft_desc_basic &fd,int n_vec_or, long long bsize0=-1, double extra_mem_factor_set = 0.82);

  void copy_evec_to_GPU(int nini);
  template <typename Ty >
  void copy_FieldM_to_Mvec(Ty* src, int ncur, int sm = 0, int dir = 1 , bool data_GPU = false);
  template <typename Ty >
  void copy_to_FieldM(Ty* src, int ncur, int sm = 0, bool data_GPU = false){
    copy_FieldM_to_Mvec(src, ncur, sm, 0, data_GPU);
  }


  template <typename Ty >
  void copy_FieldM_to_Mvec(qlat::FieldM<Ty , 12>& src, int ncur, int sm = 0, int dir = 1);
  template <typename Ty >
  void copy_to_FieldM(qlat::FieldM<Ty , 12>& src, int ncur, int sm = 0){
    copy_FieldM_to_Mvec(src, ncur, sm, 0);
  }
  void load_eigen_Mvec(const std::string& ename, io_vec  &io, int sm = 0, int nini=0, int checknorm = 1);
  void save_eigen_Mvec(const std::string& ename, io_vec  &io, int sm = 0);

  template <typename Tg >
  void smear_eigen(const std::string& Ename_Sm, io_vec  &io,
    const GaugeFieldT<Tg >& gf, const double width, const int step,
    const CoordinateD& mom = CoordinateD(), const bool smear_in_time_dir = false);


  Complexq* getEigenP(int ni, size_t xi, int sm = 0, int mode_initial = 0);

  void setup_bfac(long long bsize0=-1);

  void load_eivals(const std::string& enamev,double rho_or,double Eerr=EIGENERROR, int nini=0);
  
  void load_eigen(const std::string& ov_evecname, io_vec &io,
    int checknorm = 1, double kappa=0.2,double eigenerror=EIGENERROR, int nini=0);

  void random_eigen(int sm = 0, int seed = 1234);

  void print_info()
  {
    print_mem_info();
    double memV     = noden*12.0*sizeof(Complexq)*pow(0.5,30);
    double mem_prop = nV_prop * 12 *memV;
    print0("Lattice sites on node %10ld bfac %10ld bsize %10ld, tsize %5d !\n",noden,bfac,b_size, Nt);
    print0("num_zero %5d \n",num_zero);
    print0("bfac %10ld, n_vec %10d, b_size %10ld, matrix %10ld x %10ld \n",bfac,n_vec,b_size,2*bfac*n_vec,b_size);
    print0("===prop %.3e GB, v %d  %.3e GB, buf %d  %.3e GB, bfac_group %d, E %.3e GB. \n"
            , mem_prop, ncutgpu, memV*ncutgpu, ncutbuf, memV*ncutbuf, int(bfac_group), memV*n_vec);
  }

  /////void checknormM(Ftype err=1e-3);

  void initialize_mass(std::vector<double> &mass,int Nv=12,int one_minus_halfD_or=1, int nprop=1);
  void initialize_mass();

  ///void seq_L(vector **prop3,double mass_d,double rho,int_vector &map,bool bcast,double_complex* ker_low);
  ///void seq_L(vector **prop3,double mass_d,double rho,int_vector &map,double_complex* ker_low);
  ///void setup_L(vector **prop2);

  void setup_gpufac(int nprop=1);
  void allocate_GPU_mem(int nprop=1);
  void clear_GPU_mem();

  ~eigen_ov()
  {
    clear_GPU_mem();
    fdp = NULL;
  }

};

void eigen_ov::setup_bfac(long long bsize0)
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

  //bfac   = 768;
  //b_size = noden*6/bfac;

  if(bfac%6 !=0){print0("Cannot understand sites on bfac %5ld, bsize %10ld, Total %10ld !\n",
    bfac, b_size, noden*6);abort_r("");}
  if((noden*6)%bfac !=0){print0("Cannot understand sites on node*6/Nt %10ld bfac %10ld!\n",noden*6/Nt,bfac);abort_r("");}
  if(bfac%(Nt*6) !=0){print0("Cannot understand sites on node %10ld bfac %10ld, tsize %5d!\n",noden,bfac,Nt);abort_r("");}

}

eigen_ov::eigen_ov(qlat::fft_desc_basic &fd,int n_vec_or, long long bsize0, double extra_mem_factor_set)
{
  fdp  = &fd;
  if(fd.order_ch != 0){abort_r("Currently not supported for change fd order.\n ");}

  rho = 1.5;
  Eeigenerror = EIGENERROR;
  num_zero = 0;
  one_minus_halfD = 1;

  nx = fd.nx;ny = fd.ny;nz = fd.nz;nt = fd.nt;
  Nx = fd.Nx;Ny = fd.Ny;Nz = fd.Nz;Nt = fd.Nt;
  noden = fd.noden;

  n_vec = n_vec_or;
  ////ptmp.resize(0);stmp.resize(0);

  Eigenbuf.resize(0);
  Eigendyn.resize(0);
  Eigenbuf_Sm.resize(0);
  Eigendyn_Sm.resize(0);
  npos_Eigenbuf = -1;npos_Eigendyn = -1;

  enable_smearE = false;

  alpha.resize(0);alpha_bfac.resize(0);
  alpha_list.resize(0);eval_list.resize(0);

  ncutbuf = n_vec;
  ncutgpu = 0;
  bfac_group = -1;
  gpu_mem_set = false;

  /////buffGPU    = 2;
  nV_prop = 1;
  extra_mem_factor = extra_mem_factor_set;

  setup_bfac(bsize0);
  BFAC_GROUP_CPU = 1;
  //BFAC_GROUP_CPU = 2*bfac;
  fflush_MPI();

}

void eigen_ov::setup_gpufac(int nprop)
{
  /////To avoid very large continuous memories on CPU
  bfac_group = 1;
  ////bfac_group = 2*bfac;
  nV_prop    = nprop;

  /////redefine variables on GPU
  #ifdef QLAT_USE_ACC
  double sm_factor = 1.0;
  if(enable_smearE){sm_factor = 0.5;}

  /////buffGPU    = 2;
  size_t freeM = 0;size_t totalM = 0;
  cudaMemGetInfo(&freeM,&totalM);
  //double freeD = 0;
  //freeD = freeM*pow(0.5,30);
  double totalD=0;
  totalD = totalM*pow(0.5,30);
  struct sysinfo s_info;
  sysinfo(&s_info);

  //int Ns = 12;
  int Ns = nV_prop;if(Ns <=0)Ns = 2;
  ////long long Lat   = noden;
  double memV     = noden*12.0*sizeof(Complexq)*pow(0.5,30);
  double mem_prop = Ns * 12.0 * memV;
  if(totalD < mem_prop + 64*memV){
    print0("===GPU Memory too small, Total %.3e, prop %.3e, 30V %.3e, increase nodes! \n", totalD, mem_prop, 64*memV);
    qassert(false);
  }

  int vfac = ncutgpu;int vini = 8 * vfac;
  int vres = int((totalD*extra_mem_factor*sm_factor - mem_prop )/memV); 
  if(fdp->rank != 0){vres=0;};sum_all_size(&vres, 1);
  /////TODO Need global sum and average the final result?
  /////TODO need to test the continuous memory less thant 8GB

  if(vres > n_vec ){
    ncutgpu = 0;ncutbuf = n_vec;
    //ncutgpu = n_vec;ncutbuf = 0;
  }
  else{
    if(vres > vini){
      ncutbuf = vres - vini;
      ncutgpu=vini;
    }
    else{ncutbuf = 0;ncutgpu = ((vres)/16)*16;if(ncutgpu < 16)ncutgpu = 16;}
  }

  if(ncutgpu == 0 or ncutbuf == 0){bfac_group = 2*bfac;}else{
  bfac_group = 32;
  for(int bini=bfac_group;bini >= 1; bini--){
    size_t tem = get_threads(bini, bfac, 0);
    if(tem != 0 and bfac%tem == 0){bfac_group = tem;break;}
    bfac_group = 1;
  }}

  #endif
}

void eigen_ov::allocate_GPU_mem(int nprop)
{
  setup_gpufac(nprop);
  #ifndef QLAT_USE_ACC
  return ;
  #endif

  /////if(ncutbuf <=0 or ncutgpu <= 0){return ;}
  //if(sm == 0){if(Eigenbuf.size() != 0 or Eigendyn.size() != 0){return ;}}
  //if(sm == 1){if(Eigenbuf_Sm.size() != 0 or Eigendyn_Sm.size() != 0){return ;}}

  LInt Lae = 2*bfac/bfac_group;
  LInt Lbe = bfac_group*ncutbuf*long(b_size);
  LInt Lax = 2*bfac/bfac_group;
  LInt Lbx = bfac_group*ncutgpu*long(b_size);

  if(Eigenbuf.size() == Lae and Eigendyn.size() == Lax){return ;}

  { if(ncutbuf != 0){allocate_buf(Eigenbuf, Lae, Lbe);}
    if(ncutgpu != 0){allocate_buf(Eigendyn, Lax, Lbx);}}

  { if(ncutgpu != 0 and enable_smearE == true){allocate_buf(Eigendyn_Sm, Lax, Lbx);}
    if(ncutbuf != 0 and enable_smearE == true){allocate_buf(Eigenbuf_Sm, Lae, Lbe);}}
  gpu_mem_set = true;
}

void eigen_ov::clear_GPU_mem()
{
  for(LInt i=0;i<Eigenbuf.size();i++){Eigenbuf[i].resize(0);}  Eigenbuf.resize(0);
  for(LInt i=0;i<Eigendyn.size();i++){Eigendyn[i].resize(0);}  Eigendyn.resize(0);
  for(LInt i=0;i<Eigenbuf_Sm.size();i++){Eigenbuf_Sm.resize(0);}  Eigenbuf_Sm.resize(0);
  for(LInt i=0;i<Eigendyn_Sm.size();i++){Eigendyn_Sm.resize(0);}  Eigendyn_Sm.resize(0);

  alpha.resize(0);alpha_bfac.resize(0);
  alpha_list.resize(0);eval_list.resize(0);

  /////ptmp.resize(0);stmp.resize(0);

  gpu_mem_set = false;
}

void eigen_ov::load_eivals(const std::string& enamev,double rho_or,double Eerr, int nini)
{
  rho = rho_or;

  std::vector<double > v,e;
  load_gwu_eigenvalues(v,e, enamev.c_str());
  if((nini+n_vec)*2 > int(v.size()))abort_r("Eigen value size too small! ");
  eval_self.resize(n_vec);
  for(int iv=0;iv<n_vec;iv++){
    int off= iv + nini;
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

}

template <typename Ty >
void eigen_ov::copy_FieldM_to_Mvec(Ty* src, int ncur, int sm, int dir , bool data_GPU)
{
  TIMERA("COPY Eigen Vectors Mvec");
  if(ncur >= n_vec){abort_r("Cannot copy to position larger than n_vec ! \n");}
  /////int nread = nb - ba;
  Complexq* s1 = NULL;Ty* s0 = NULL;
  s0 = src;

  int GPU_cpy = 0;
  if(data_GPU){
    if(dir == 1){GPU_cpy = 3;} // from device to host
    if(dir == 0){GPU_cpy = 2;} // from host to device
  }

  //////move d,c to outter loop
  long sizeF = noden;
  if(dir == 1){fdp->mv_civ.dojob(s0, s0, 1, 12, sizeF, 1, 1, data_GPU);}

  ////a factor of 2 by chiral
  LInt total = 2*bfac*b_size;
  for(LInt xini=0;xini < total/b_size; xini++)
  {
    LInt xi = xini*b_size;
    s1 = getEigenP(ncur, xi, sm, 1);
    /////cudaMemcpy(s1, &s0[xi], b_size*sizeof(Complexq), cudaMemcpyHostToDevice);
    if(dir == 1){cpy_data_thread(s1, &s0[xi], b_size, GPU_cpy, false);}
    if(dir == 0){cpy_data_thread(&s0[xi], s1, b_size, GPU_cpy, false);}
  }
  qacc_barrier(dummy);

  if(dir == 0){fdp->mv_civ.dojob(s0, s0, 1, 12, sizeF, 0, 1, data_GPU);}
  s0 = NULL; s1 = NULL;
}


template <typename Ty >
void eigen_ov::copy_FieldM_to_Mvec(qlat::FieldM<Ty , 12>& src, int ncur, int sm, int dir )
{
  TIMERA("COPY Eigen Vectors Mvec");
  if(ncur >= n_vec){abort_r("Cannot copy to position larger than n_vec ! \n");}
  if(dir == 0){if(!src.initialized){
    Geometry geo;fdp->get_geo(geo);
    src.init(geo);
  }}
  Ty* s0 = (Ty*) qlat::get_data(src).data();
  int data_GPU = 0; ////do copies only from CPU here
  copy_FieldM_to_Mvec(s0, ncur, sm, dir, data_GPU);
}

void eigen_ov::copy_evec_to_GPU(int nini)
{
  #ifndef QLAT_USE_ACC
  return ;
  #else

  int mode_dyn = 0;
  if(nini < 0){abort_r("Copy initial negative! \n");}
  if(nini <  ncutbuf){mode_dyn = 0;}
  if(nini >= ncutbuf){mode_dyn = 1;}

  int n0 = 0;int n1 = 0;
  if(mode_dyn == 0){n0=0;n1 = n0 + ncutbuf; }
  if(mode_dyn == 1){
    int maxN = (nini + ncutgpu)/ncutgpu;
    n0=(maxN-1)*ncutgpu ;n1 = n0 + ncutgpu; 
  }

  if(mode_dyn == 0){if(npos_Eigenbuf == n1){return ;}else{npos_Eigenbuf = n1;}}
  if(mode_dyn == 1){if(npos_Eigendyn == n1){return ;}else{npos_Eigendyn = n1;}}

  {
  TIMERB("COPY Eigen Vectors to GPU");
  Complexq* s1;Complexq* s0;
  /////a factor of 2 by chiral
  LInt total = 2*bfac*b_size;
  for(int ncur=n0;ncur<n1;ncur++)
  {
    if(ncur >= n_vec){break;}
    for(LInt xini=0;xini < total/b_size; xini++)
    {
      LInt xi = xini*b_size;
      /////print0("Test num ncur %d, xi %d, bfac %d \n", ncur, xi, bfac);

      s1 = getEigenP(ncur, xi, 0, 0);s0 = getEigenP(ncur, xi, 0, 1);
      //////cudaMemcpy(s1, s0, b_size*sizeof(Complexq), cudaMemcpyHostToDevice);
      cpy_data_thread(s1, s0 , b_size, 2, false);

      if(enable_smearE == true){
        s1 = getEigenP(ncur, xi, 1, 0);s0 = getEigenP(ncur, xi, 1, 1);
        /////cudaMemcpy(s1, s0, b_size*sizeof(Complexq), cudaMemcpyHostToDevice);
        cpy_data_thread(s1, s0 , b_size, 2, false);
      }
    }
  }
  qacc_barrier(dummy);
  }
  #endif

}

Complexq* eigen_ov::getEigenP(int ni, size_t xi, int sm, int mode_initial )
{
  if(sm == 1 and enable_smearE == false){abort_r("Smeared Eigen system not loaded or computed! \n");}
  if(ni >= n_vec){abort_r("Request ni too large! \n");}

  Complexq* buf = NULL;
  Complexq* s0 = NULL;
  int   chi = xi/(bfac*b_size);
  size_t vi = xi%(bfac*b_size);
  size_t bi = vi/b_size;
  size_t bj = vi%b_size;

  /////On CPU use Mvec and Mvec_Sm
  #ifndef QLAT_USE_ACC
  mode_initial = 1;
  #endif

  if(mode_initial == 1)
  {
    size_t off1 = ((chi*bfac + bi)*n_vec + ni )*b_size + bj;
    size_t Lb = BFAC_GROUP_CPU*n_vec*size_t(b_size);
    size_t Li = off1/Lb;
    size_t Lj = off1%Lb;

    if(sm == 0){s0 = Mvec[Li].data();}
    if(sm == 1){s0 = Mvec_Sm[Li].data();}

    buf = &(s0[Lj]);
    return buf;
  }
  ///////only allocate GPU memory when get the vectors on GPU, reuse later
  //////Free the memory when needed, clear_GPU_mem()
  //////copy_evec_to_GPU(ni);
  if(bfac_group <= 0){abort_r("Setup bfac_group first !\n");}

  if(ni  < ncutbuf){
    size_t off1 = ((chi*bfac + bi)*ncutbuf + ni )*b_size + bj;
    size_t Lb = bfac_group*ncutbuf*size_t(b_size);
    size_t Li = off1/Lb;
    size_t Lj = off1%Lb;
    if(sm == 0){s0 = Eigenbuf[Li].data();}
    if(sm == 1){s0 = Eigenbuf_Sm[Li].data();}
    buf = &(s0[Lj]);
  }

  if(ni >= ncutbuf){
    int ncur = ni - ncutbuf;
    int na   = 0;
    int nb   = ncur%ncutgpu;

    size_t off1 = (((na*2 + chi)*bfac + bi)*ncutgpu + nb)*b_size + bj;
    size_t Lb = bfac_group*ncutgpu*size_t(b_size);
    size_t Li = off1/Lb;
    size_t Lj = off1%Lb;

    if(sm == 0){s0 = Eigendyn[Li].data();}
    if(sm == 1){s0 = Eigendyn_Sm[Li].data();}
    buf = &(s0[Lj]);
  }

  s0 = NULL;

  return buf;
}

inline void resize_EigenM(Elocal& a, size_t n0, size_t n1)
{
  a.resize(0);
  a.resize(n0);
  for(size_t iv=0;iv<n0;iv++)
  {
    a[iv].resize(n1);
    ////zeroE(a[iv], 1);
    zero_Ty(a[iv].data(), n1, 0);
  }
}

template <typename Tg >
void eigen_ov::smear_eigen(const std::string& Ename_Sm, io_vec  &io,
  const GaugeFieldT<Tg >& gf, const double width, const int step,
  const CoordinateD& mom, const bool smear_in_time_dir)
{
  TIMER("smear eigen system");
  ////load if exist
  if(Ename_Sm != std::string("NONE") and get_file_size_MPI(Ename_Sm.c_str(), true) != 0){
    load_eigen_Mvec(Ename_Sm, io, 1);
    return ;
  }

  long La = 2*bfac/BFAC_GROUP_CPU;
  long Lb = BFAC_GROUP_CPU*n_vec*long(b_size);
  print_mem_info("Before Eigen Memory Allocate");
  resize_EigenM(Mvec_Sm , La, Lb);enable_smearE = true;
  print_mem_info("Eigen Memory Allocate Done");

  const int each = 12;
  long Ncopy = fdp->Nvol * 12;
  qlat::vector_gpu<Complexq > buf;buf.resize(each * Ncopy);

  move_index mv_idx;
  std::vector<long > job =  job_create(n_vec, each);
  for(LInt ji = 0; ji < job.size()/2 ; ji++)
  {
    int flag = 0;
    ////copy to buf
    for(int iv=0;iv<job[ji*2 + 1];iv++){copy_to_FieldM(&buf[iv*Ncopy], job[ji*2 + 0] + iv,  0, true);}
    flag = 0;mv_idx.dojob(buf.data(), buf.data(), 1, each , fdp->Nvol*12, flag, 1, true);

    smear_propagator_gwu_convension_inner<Complexq, 4,each  , Tg>(buf.data(), gf, width, step, mom, smear_in_time_dir);

    flag = 1;mv_idx.dojob(buf.data(), buf.data(), 1, each , fdp->Nvol*12, flag, 1, true);
    ////copy from buf
    for(int iv=0;iv<job[ji*2 + 1];iv++){copy_FieldM_to_Mvec(&buf[iv*Ncopy], job[ji*2 + 0] + iv,  1, 1, true);}
  }

  ////erase smear
  {
  /////const SmearPlanKey& skey = get_smear_plan_key<Complexq, 4, each>(gf.geo(), smear_in_time_dir);
  ////qlat::clear_all_caches();
  get_smear_plan_cache().clear();
  }

  if(Ename_Sm != std::string("NONE")){
    save_eigen_Mvec(Ename_Sm, io, 1);
  }
}

void eigen_ov::save_eigen_Mvec(const std::string& ename, io_vec  &io, int sm)
{
  if(sm == 1 and enable_smearE == false){print0("Could not save smear eigen without set it up.");return ;}
  const int nini = 0;
  const int ntotal = nini + n_vec;
  const bool read = false;
  inputpara in_write_eigen;
  FILE* file_write  = open_eigensystem_file(ename.c_str(), nini, ntotal, read , io , in_write_eigen , 2);

  int each = io.ionum;
  std::vector<qlat::FieldM<Complexq , 12> > buf;buf.resize(each);
  for(int iv=0;iv<each;iv++){buf[iv].init(io.geop);}

  std::vector<long > job =  job_create(n_vec, each);
  for(LInt ji = 0; ji < job.size()/2 ; ji++)
  {
    for(int iv=0;iv<job[ji*2 + 1];iv++){copy_to_FieldM(buf[iv], job[ji*2 + 0] + iv, sm );}
    /////write to file
    load_eigensystem_vecs(file_write ,   buf, io , in_write_eigen , 0, job[ji*2 + 1]);
  }

  close_eigensystem_file(file_write , io , in_write_eigen );
  print_mem_info("Eigen Memory Write Done");
}

void eigen_ov::load_eigen_Mvec(const std::string& ename, io_vec  &io, int sm, int nini, int checknorm)
{
  int ntotal = nini + n_vec;
  Ftype norm_err  = 1e-3;

  long La = 2*bfac/BFAC_GROUP_CPU;
  long Lb = BFAC_GROUP_CPU*n_vec*long(b_size);
  print_mem_info("Before Eigen Memory Allocate");
  if(sm == 0){resize_EigenM(Mvec    , La, Lb);}
  if(sm == 1){resize_EigenM(Mvec_Sm , La, Lb);enable_smearE = true;}
  print_mem_info("Eigen Memory Allocate Done");

  inputpara in_read_eigen;
  FILE* file_read  = open_eigensystem_file(ename.c_str(), nini, ntotal, true , io , in_read_eigen , 2);

  int each = io.ionum;
  std::vector<qlat::FieldM<Complexq , 12> > buf;buf.resize(each);
  for(int iv=0;iv<each;iv++){buf[iv].init(io.geop);}

  std::vector<long > job =  job_create(n_vec, each);
  for(LInt ji = 0; ji < job.size()/2 ; ji++)
  {
    ////int n0 = nini + job[ji*2 + 0]; int n1 = n0 + job[ji*2 + 1]; 
    /////load from file
    load_eigensystem_vecs(file_read ,   buf, io , in_read_eigen , 0, job[ji*2 + 1]);
    ////copy to Mvec or Mvec_Sm
    for(int iv=0;iv<job[ji*2 + 1];iv++){
      if(checknorm == 1 and sm == 0){
        Ftype* Psrc = (Ftype*) qlat::get_data(buf[iv]).data();
        Ftype normf = get_norm_vec(Psrc, noden);
        if(fabs(normf - 1.0) > norm_err){
          print0("Eigen vector %d, norm %.3e wrong. \n", int(job[ji*2 + 0] + iv), normf);
          abort_r("");
        }
      }

      copy_FieldM_to_Mvec(buf[iv], job[ji*2 + 0] + iv, sm  );
    }
  }

  close_eigensystem_file(file_read , io , in_read_eigen );

  print_mem_info("Eigen Memory Load Done");

  //////mv_civ.free_mem();
}

void eigen_ov::load_eigen(const std::string& ov_evecname, io_vec  &io,
  int checknorm, double kappa,double eigenerror, int nini)
{
  TIMERB("=====Loading Eigen=====");
  char enamev[600];
  ////sprintf(ename, "%s", ov_evecname);
  sprintf(enamev,"%s.eigvals", ov_evecname.c_str());
  print0("Vector File name: %s \n", ov_evecname.c_str() );
  print0("Values File name: %s \n", enamev);
  //////Load eigen values
  double rho_tem = 4 - 1.0/(2*kappa);
  load_eivals(std::string(enamev), rho_tem, eigenerror, nini);

  load_eigen_Mvec(ov_evecname, io, 0 ,nini, checknorm);

}

void eigen_ov::random_eigen(int sm, int seed)
{
  TIMERB("=====Loading random Eigen=====");
  eval_self.resize(n_vec);random_EigenM(eval_self, 0, seed + 10);
  
  long La = 2*bfac/BFAC_GROUP_CPU;
  long Lb = BFAC_GROUP_CPU*n_vec*long(b_size);
  print_mem_info("Before Eigen Memory Allocate");
  if(sm == 0){
    resize_EigenM(Mvec    , La, Lb);
    for(long iv=0;iv<La;iv++)random_Ty(Mvec[iv].data()   , Lb, 0, seed + 20 + iv);
  }
  if(sm == 1){
    resize_EigenM(Mvec_Sm , La, Lb);
    for(long iv=0;iv<La;iv++)random_Ty(Mvec_Sm[iv].data(), Lb, 0, seed + 30 + iv);
    enable_smearE = true;
  }

  print_mem_info("Eigen random Memory Load Done");

}

//////Nv source number, 12 for prop
void eigen_ov::initialize_mass(std::vector<double> &mass, int Ns, int one_minus_halfD_or, int nprop)
{
  TIMERB("Set up store memories");
  ////std::vector<double> mass = mass_or;
  massL = mass;
  int mN = mass.size();
  one_minus_halfD = one_minus_halfD_or;

  if(gpu_mem_set == false or nV_prop != nprop){allocate_GPU_mem(nprop);}

  LInt nlarge = ncutgpu;if(ncutbuf > ncutgpu)nlarge = ncutbuf;

  alpha.resize(2*nlarge*Ns);alpha_bfac.resize(2*nlarge*Ns*bfac);
  alpha_list.resize(2*nlarge * mN*Ns);

  zero_Ty(alpha.data(), alpha.size(), 1, false);
  zero_Ty(alpha_bfac.data(), alpha_bfac.size(), 1, false);
  zero_Ty(alpha_list.data(), alpha_list.size(), 1, false);

  ////alpha.clear(false);alpha_bfac.clear(false);alpha_list.clear(false);

  //size_t tmp_Ls = 2*bfac * Ns*b_size;
  //size_t tmp_Lp = 2*mN*bfac * Ns*b_size;
  //stmp.resize(tmp_Ls);ptmp.resize(tmp_Lp);
  ////stmp.clear(false);ptmp.clear(false);

  if(eval_tem.size() != long(mN*n_vec)){eval_tem.resize(mN*n_vec);}
  #pragma omp parallel for
  for(int mki=0;mki< mN*n_vec;mki++)
  {
    int mi = mki/n_vec;
    int kn = mki%n_vec;
    eval_tem[mi*n_vec + kn] = inv_self(eval_self[kn], mass[mi], rho, one_minus_halfD);
  }

  eval_list.resize( mN*n_vec);
  cpy_data_thread(eval_list.data(), (Complexq*) qlat::get_data(eval_tem).data(), eval_list.size(),  1, false);
  ///eval_list.copy_from(eval_tem, 1, 1);

  /////allocate and copy GPU memory
  copy_evec_to_GPU(0);

  qacc_barrier(dummy);

}

void eigen_ov::initialize_mass()
{
  std::vector<double> mass;mass.push_back(0.0);
  initialize_mass(mass,12,one_minus_halfD);
}

void prop_L_device(eigen_ov& ei,Complexq *src,Complexq *props, int Ns, std::vector<double> &mass, int mode_sm = 0,int one_minus_halfD_or=1)
{
  int mN   = mass.size();

  TIMER_FLOPS("==prop_L");

  long long Lat = ei.noden;
  long long vGb = Lat*12;
  int Fcount0 = 6 + 2;
  int Fcount1 = 6 + 2;
  long long Tfloat = ei.n_vec*Ns*mN*vGb*Fcount1 + ei.n_vec*Ns*vGb*Fcount0;
  timer.flops += Tfloat;
  //double mem = Lat*12*(ei.n_vec + Ns + Ns*mN)*sizeof(Complexq);
  //print0("Memory size %.3e GB, %.3e Gflop \n", 
  //  mem/(1024.0*1024*1024), Tfloat/(1024.0*1024*1024));
  ///////qlat::get_num_node()

  ei.initialize_mass(mass, Ns, one_minus_halfD_or, Ns*mN);
  ////touch_GPU(ei.eval_list, ei.eval_list_size);

  int sm0 = 0; int sm1 = 0;
  if(mode_sm == 0){sm0 = 0; sm1 = 0;}
  if(mode_sm == 1){sm0 = 0; sm1 = 1;}
  if(mode_sm == 2){sm0 = 1; sm1 = 0;}
  if(mode_sm == 3){sm0 = 1; sm1 = 1;}

  Complexq* alpha       = ei.alpha.data();
  Complexq* alpha_bfac  = ei.alpha_bfac.data();
  Complexq* alpha_list  = ei.alpha_list.data();
  Complexq* eval_list   = ei.eval_list.data();

  long& bfac        = ei.bfac;
  long& b_size      = ei.b_size;

  int&  n_vec       = ei.n_vec;
  int& ncutgpu     = ei.ncutgpu;
  int& ncutbuf     = ei.ncutbuf;
  int& num_zero    = ei.num_zero;

  bool dummy_test = false;

  int Ng = 0;if(ncutgpu!=0){Ng = ei.n_vec/ncutgpu + 1;}
  int nini = 0;
  for(int ng=0;ng<Ng + 1;ng++)
  {
    if(nini >= ei.n_vec)break;
    int  ncur = 0;
    if(nini>=ncutbuf){
      ncur = ncutgpu;
    }else{ncur = ncutbuf;}

    ei.copy_evec_to_GPU(nini);

    {
    long m = ncur;
    long n = Ns;
    long w = b_size;

    TIMER_FLOPS("vec reduce");
    long long vGb = 2*bfac*m*n*w;
    int Fcount0   = 6 + 2;  
    timer.flops  += vGb*Fcount0;

    std::vector<long > jobA = job_create(2*bfac, ei.BFAC_GROUP_CPU);
    if(nini > ncutbuf){jobA = job_create(2*bfac, ei.bfac_group);}
    for(LInt jobi=0;jobi < jobA.size()/2; jobi++)
    {
      long bini = jobA[jobi*2 + 0]; long bcut = jobA[jobi*2+1];
      Complexq* rp = &alpha_bfac[(bini + 0)*m*n + 0];
      Complexq* Ep = ei.getEigenP(nini, bini*b_size, sm0);
      Complexq* sp = &src[       (bini + 0)*n*w + 0];
      matrix_prod_gpu(Ep, sp, rp, m,n, w , bcut, true, dummy_test);
    }

    /////Output not conflict, GPU end here
    qacc_barrier(dummy);

    }

    {
    TIMERA("Reduce alpha");
    qacc_for(coff, long(2*ncur*Ns),{
      int chi = coff/(ncur*Ns);
      long xi = coff%(ncur*Ns);
      for(long bi=0;bi<bfac;bi++){
        alpha[coff] += alpha_bfac[chi*bfac*ncur*Ns + bi*ncur*Ns + xi];
      }
    });
    }

    {
    TIMERA("Global sum");
    sum_all_size((Ftype*) (ei.alpha.data()), 2*(2*ncur*Ns), 1);
    }

    ////print_sum(ei.alpha.data(), ei.alpha.size(), "=====alpha", 1);

    {
    TIMERA("Get alpha list")
    Complexq Iimag(0.0,1.0);
    Complexq Two2(2.0,0.0);
    qacc_for(coff, long(2*mN*Ns*ncur),{
      int chi   =  coff/(mN*Ns*ncur);
      int mi    = (coff/(Ns*ncur  ))%mN;
      int is    = (coff/(ncur     ))%Ns;
      int kn    = (coff             )%ncur;
      unsigned long offA = ((chi*mN+mi)*Ns+is)*ncur + kn;

      if(kn + nini >= n_vec){alpha_list[offA] = 0.0;}else{
      Complexq li = eval_list[mi*n_vec + kn + nini];
      if(kn+nini >= num_zero){
        alpha_list[offA] = Two2*(li.real()*alpha[(chi*ncur+kn)*Ns+is]+Iimag*li.imag()*alpha[((1-chi)*ncur+kn)*Ns+is]);}
      if(kn+nini  < num_zero){
        alpha_list[offA] = li*alpha[(chi*ncur+kn)*Ns+is];}
      }
    });
    }


    {
    long m = mN*Ns;
    long n = b_size;
    long w = ncur;

    //TIMER("vec multi");
    TIMER_FLOPS("vec multi");
    long long vGb = 2*bfac*m*n*w;
    int Fcount0   = 2*(3 + 1);  
    timer.flops  += vGb*Fcount0;

    if((nini + ncur) < ei.n_vec){
      ////zero_Ty(ei.alpha_bfac, ei.alpha_bfac_size, 0, false);zero_Ty(ei.alpha, ei.alpha_size, 0, false);
      //ei.alpha_bfac.clear(false);
      //ei.alpha.clear(false);
      zero_Ty(ei.alpha.data(), ei.alpha.size(), 1, false);
      zero_Ty(ei.alpha_bfac.data(), ei.alpha_bfac.size(), 1, false);
    }

    for(long coff=0;coff<2*bfac;coff++)
    {
      long chi = coff/bfac;
      long bi  = coff%bfac;

      Complexq* rp = &props[(chi*bfac+bi)*mN*Ns*b_size + 0];
      ////Complexq* rp = &ei.ptmp[(chi*bfac+bi)*mN*Ns*b_size + 0];
      Complexq* ap = &alpha_list[(chi*mN*Ns+0)*w + 0];

      Complexq* Ep = ei.getEigenP(nini, coff*b_size, sm1);

      matrix_prod_gpu(ap, Ep, rp, m, n, w ,1, false,  dummy_test, true);
    }
    /////Output not conflict, GPU end here
    qacc_barrier(dummy);
    ///nini += ncur;

    }
    nini += ncur;
  }

  ////print_sum(props, 2*bfac*  mN*Ns*b_size, "=====props", 1);

  /////untouch_GPU(ei.eval_list, ei.eval_list_size);
  //////#endif

}

////dir == 0, from src to EigenM res
////dir == 1 from EigenM res to src
template <typename Ty >
void copy_eigen_prop_to_EigenM(Ty* src, EigenM& res, LInt b_size, int nmass, qlat::fft_desc_basic& fd, int dir = 0,int GPU = 1)
{
  if(nmass == 0){return ;}
  if(dir == 0){ini_propE(res, nmass, fd, false);}

  int Ns    = 12*nmass;
  int  NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  LInt total = 6*NTt*Nxyz;
  if(total % b_size != 0){abort_r("eigen system configurations wrong! \n");}
  LInt bfac = total/(b_size);
  LInt each  = Nxyz; if(b_size < Nxyz){each = b_size;}
  LInt group = (2*total)/each;

  //for(unsigned int d1=0;d1<12;d1++)
  //for(int ti=0; ti < NTt; ti++)
  for(int d0=0;d0<Ns;d0++)
  for(LInt gi=0;gi<group;gi++)
  {
    LInt mi = gi*each;

    ////index for res
    LInt d1 =  mi/(NTt*Nxyz);
    LInt ti = (mi/(Nxyz))%NTt;
    LInt vi =  mi%(Nxyz);

    ////index for src
    int chi = mi/(total);
    LInt xi = mi%(total);
    long bi = xi/b_size;
    long bj = xi%b_size;

    Ty* s0       = &src[(chi*bfac+bi)*Ns*b_size  + d0*b_size + bj];
    Complexq* tems = (Complexq*) qlat::get_data(res[(d0*12 + d1)*NTt+ti]).data();
    Complexq* s1 = &tems[vi];
    if(dir == 0){cpy_data_thread(s1, s0, each , GPU, false);}
    if(dir == 1){cpy_data_thread(s0, s1, each , GPU, false);}

    /////useS[(chi*bfac+bi)*Ns*b_size  + is*b_size + bj] = s0[(d0*12 + d1)*NTt+ti][vi];

  }
  qacc_barrier(dummy);

}

/////res in format src 12* sink 12 --> Nt * Nxyz
template <typename Ty >
void FieldM_src_to_FieldM_prop(qlat::FieldM<Ty , 1>& src, qlat::FieldM<Ty , 12*12>& res, int GPU = true, bool dummy = true)
{
  qlat::Geometry& geo = src.geo();

  if(!res.initialized){res.init(geo);}

  //bool do_ini = true;
  //if(res.size() == src.size())if(res[src.size()-1].initialized){do_ini = false;}
  //if(do_ini){res.resize(nV);for(int iv=0;iv<nV;iv++){res[iv].init(geo);}}

  //std::vector<int > nv, Nv, mv;
  //geo_to_nv(geo, nv, Nv,mv);
  long Ncopy = geo.local_volume();

  Ty* s0 = NULL; Ty* s1 = NULL;Ty* st = NULL;
  ///for(int iv=0;iv<nV;iv++)
  s0 = (Ty*) qlat::get_data(src).data();
  st = (Ty*) qlat::get_data(res).data();
  for(unsigned int d0=0;d0<12;d0++)
  {
    //////diagonal elements
    s1 = &st[(d0*12+d0)*Ncopy + 0];
    cpy_data_thread(s1, s0, Ncopy , GPU, false);
  }
  if(dummy)qacc_barrier(dummy);
}

template <typename Ty >
void FieldM_src_to_FieldM_prop(std::vector<qlat::FieldM<Ty , 1> >& src, std::vector<qlat::FieldM<Ty , 12*12> >& res, int GPU = true)
{
  if(src.size() == 0){return ;}
  ////qlat::Geometry& geo = src[0].geo();
  long nV = src.size();
  if(res.size() != src.size()){res.resize(nV);}
  for(int iv=0;iv<nV;iv++)FieldM_src_to_FieldM_prop(src[iv], res[iv], GPU, false);
  qacc_barrier(dummy);

  //bool do_ini = true;if(res.size() == src.size())if(res[src.size()-1].initialized){do_ini = false;}
  //if(do_ini){res.resize(nV);for(int iv=0;iv<nV;iv++){res[iv].init(geo);}}

  ////std::vector<int > nv, Nv, mv;
  ////geo_to_nv(geo, nv, Nv,mv);
  //long Ncopy = geo.local_volume();

  //Ty* s0 = NULL; Ty* s1 = NULL;Ty* st = NULL;
  //for(int iv=0;iv<nV;iv++)
  //{
  //  s0 = (Ty*) qlat::get_data(src[iv]).data();
  //  st = (Ty*) qlat::get_data(res[iv]).data();
  //  for(unsigned int d0=0;d0<12;d0++)
  //  {
  //    //////diagonal elements
  //    s1 = &st[(d0*12+d0)*Ncopy + 0];
  //    cpy_data_thread(s1, s0, Ncopy , GPU, false);
  //  }
  //}
  //qacc_barrier(dummy);

}

////assumed civ == n*12 with n the source indices, 12 the sink indices 
template <typename Ty, int civ >
void copy_eigen_src_to_FieldM(qlat::vector_gpu<Ty >& src, std::vector<qlat::FieldM<Ty , civ> >& res, LInt b_size, qlat::fft_desc_basic& fd, int dir = 0, int GPU = 1, bool rotate = false)
{
  if(civ%12 != 0){abort_r("FieldM type not supported!\n");}
  unsigned int nV = 0;int cfac = civ/12;

  int  NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  LInt sizeF = NTt*Nxyz;
  LInt total = 6*sizeF;
  if(total % b_size != 0){abort_r("eigen system configurations wrong! \n");}

  if(dir == 0){
    long dsize = src.size();
    if(dsize%(2*total) != 0){abort_r("src size wrong!\n");};
    nV  = dsize/(2*total);
    if(nV%(cfac) != 0){abort_r("res civ wrong!\n");}
    unsigned int ntem = nV/cfac;

    bool do_ini = true;if(res.size() == ntem)if(res[ntem-1].initialized){do_ini = false;}
    if(do_ini){
      //////print0("initial Fprop. \n");
      Geometry geo;fd.get_geo(geo);
      res.resize(0);res.resize(ntem);
      for(LInt iv=0;iv<res.size();iv++){res[iv].init(geo);}}
  }
  if(dir == 1){
    nV = res.size() * cfac;
    src.resize(nV * 2*total);
  }

  /////rotate FieldM, from Vol->civ to civ->Vol
  if(dir == 1 and rotate == true){
    for(LInt iv=0;iv<res.size();iv++){
      Ty* s0 = (Ty*) qlat::get_data(res[iv]).data();
      fd.mv_civ.dojob(s0, s0, 1, civ, sizeF, 1, 1, GPU);
    }
  }

  LInt bfac = total/(b_size);
  LInt each  = Nxyz; if(b_size < Nxyz){each = b_size;}
  LInt group = (2*total)/each;

  Ty* psrc       = src.data();
  Ty* s0 = NULL; Ty* s1 = NULL;Ty* st = NULL;
  for(unsigned int d0=0;d0<nV;d0++)
  for(LInt gi=0;gi<group;gi++)
  {
    LInt mi = gi*each;

    ////index for res
    LInt d1 =  mi/(NTt*Nxyz);
    LInt ti = (mi/(Nxyz))%NTt;
    LInt vi =  mi%(Nxyz);
    int d0a = d0/cfac;
    int d0b = d0%cfac;

    ////index for src
    int chi = mi/(total);
    LInt xi = mi%(total);
    long bi = xi/b_size;
    long bj = xi%b_size;

    s0 = &psrc[(chi*bfac+bi)*nV*b_size  + d0*b_size + bj];
    st = (Complexq*) qlat::get_data(res[d0a]).data();
    s1 = &st[((d0b*12 + d1)*NTt+ti)*Nxyz + vi];

    if(dir == 0){cpy_data_thread(s1, s0, each , GPU, false);}
    if(dir == 1){cpy_data_thread(s0, s1, each , GPU, false);}
  }

  qacc_barrier(dummy);

  if(dir == 0 and rotate == true){
    for(LInt iv=0;iv<res.size();iv++){
      Ty* s0 = (Ty*) qlat::get_data(res[iv]).data();
      fd.mv_civ.dojob(s0, s0, 1, civ, sizeF, 0, 1, GPU);
    }
  }
}

template <typename Ty, int civ >
void copy_FieldM_to_eigen_src(std::vector<qlat::FieldM<Ty , civ> >& src, qlat::vector_gpu<Ty >& res, LInt b_size, qlat::fft_desc_basic& fd, int GPU = 1, bool rotate = false)
{
  copy_eigen_src_to_FieldM(res, src, b_size, 1,fd, GPU, rotate);
}



}

#undef  Elocal
#undef  EIGENERROR
////#undef  Vlocal

#endif

