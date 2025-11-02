// utils_eigen_ov.sh
// Gen Wang
// Sep. 2021

#ifndef UTILS_EIGEN_OV_H
#define UTILS_EIGEN_OV_H

#pragma once
#include "utils_gammas.h"
#include "utils_fft_desc.h"
#include "utils_Matrix_prod.h"
#include "utils_check_fun.h"
#include "utils_sector_funs.h"
#include "utils_smear_vecs.h"
#include "utils_props_type.h"
#include "utils_Smear_vecs_2link.h"
#include "utils_io_vec.h"

///#define SUMMIT 0

////#define Vlocal qlat::vector<Complexq >
#define Elocal std::vector<qlat::vector<Complexq > >
#define Vlocal qlat::vector_gpu<Complexq >
#define EIGENERROR 1e-11

namespace qlat{

struct eigen_ov {
  Geometry geo;

  Elocal Mvec;     //  nvec --> 2*bfac --> b_size*6
  Elocal Mvec_Sm;  //  Smeared eigensystem 

  bool enable_smearE;

  ////int Edouble;
  ////int smear_Eigen;

  ///move_index mv_civ;

  /////lattice copied from fd
  Int nx;int ny;int nz;int nt;
  Int Nx;int Ny;int Nz;int Nt;
  LInt noden;

  ////Eigensystem information
  Int n_vec;
  Int num_zero;
  Int one_minus_halfD;
  double rho;
  double Eeigenerror;
  EigenV eval_self;

  /////Group the inner vector for vectorizations
  Long b_size, bfac;

  bool gpu_mem_set;
  Int ncutgpu, ncutbuf;
  Long bfac_group;
  Long BFAC_GROUP_CPU;


  Int nV_prop;
  double extra_mem_factor;

  ///Vlocal stmp;
  ///Vlocal ptmp;

  std::vector< Vlocal > Eigenbuf;
  std::vector< Vlocal > Eigendyn;
  std::vector< Vlocal > Eigenbuf_Sm;
  std::vector< Vlocal > Eigendyn_Sm;
  Int npos_Eigenbuf;int npos_Eigendyn;

  Vlocal alpha;
  Vlocal alpha_buf;
  Vlocal alpha_bfac;
  Vlocal alpha_list;
  Vlocal eval_list;
  std::vector<double > massL;

  qlat::vector<Complexq > eval_tem;
  /////int ncut0,ncut1;

  // buf for pointers
  vector<Complexq* > rpL;
  vector<Complexq* > EpL;
  vector<Complexq* > spL;

  Int Ns_buf;
  Int nprop_buf;

  ////3pt function
  //EigenV alpha_list3pt;
  //EigenV alpha_list3ptC;
  //EigenV alpha3pt;
  //EigenV alpha3pt_ker_low;

  /////Construction and memory allocations
  eigen_ov(const Geometry& geo_,Int n_vec_or, Long bsize0=-1, double extra_mem_factor_set = 0.82);

  void copy_evec_to_GPU(Int nini);
  template <typename Ty >
  void copy_FieldM_to_Mvec(Ty* src, Int ncur, Int sm = 0, Int dir = 1 , bool data_GPU = false);
  template <typename Ty >
  void copy_to_FieldM(Ty* src, Int ncur, Int sm = 0, bool data_GPU = false){
    copy_FieldM_to_Mvec(src, ncur, sm, 0, data_GPU);
  }


  template <typename Ty >
  void copy_FieldM_to_Mvec(qlat::FieldM<Ty , 12>& src, Int ncur, Int sm = 0, Int dir = 1);
  template <typename Ty >
  void copy_to_FieldM(qlat::FieldM<Ty , 12>& src, Int ncur, Int sm = 0){
    copy_FieldM_to_Mvec(src, ncur, sm, 0);
  }
  void load_eigen_Mvec(const std::string& ename, Int sm = 0, Int nini=0, Int checknorm = 1);

  template <typename Tg >
  void load_eigen_Mvec_smear(const std::string& ename, 
    const GaugeFieldT<Tg >& gf,
    std::vector< qlat::GaugeFieldT<Ftype > >& gfL,
    std::vector<qlat::vector_gpu<qlat::ComplexT<Ftype > > >& propS,
    Int nini = 0, Int checknorm = 0,  
    const double src_width  = 0.0, const Int src_step  = 0,
    const double sink_width = 0.0, const Int sink_step = 0,
    const CoordinateD& src_mom  = CoordinateD(),
    const CoordinateD& sink_mom = CoordinateD(),
    const bool src_smear_in_time_dir  = false,
    const bool sink_smear_in_time_dir = false
  );
  void save_eigen_Mvec(const std::string& ename, Int sm = 0, Int save_type = 2);

  template <typename Tg >
  void smear_eigen(const std::string& Ename_Sm,
    const GaugeFieldT<Tg >& gf, const double width, const Int step,
    const CoordinateD& mom = CoordinateD(), const bool smear_in_time_dir = false);


  Complexq* getEigenP(Int ni, size_t xi, Int sm = 0, Int mode_initial = 0);

  void setup_bfac(Long bsize0=-1);

  void load_eivals(const std::string& enamev,double rho_or,double Eerr=EIGENERROR, Int nini=0);
  
  void load_eigen(const std::string& ov_evecname,
    Int checknorm = 1, double kappa=0.2,double eigenerror=EIGENERROR, Int nini=0);

  void random_eigen(Int sm = 0, Int seed = 1234);

  void print_info()
  {
    print_mem_info();
    double memV     = noden*12.0*sizeof(Complexq)*pow(0.5,30);
    double mem_prop = nV_prop * 12 *memV;
    qmessage("Lattice sites on node %10ld bfac %10ld bsize %10ld, tsize %5d !\n",noden,bfac,b_size, Nt);
    qmessage("num_zero %5d \n",num_zero);
    qmessage("bfac %10ld, n_vec %10d, b_size %10ld, matrix %10ld x %10ld \n",bfac,n_vec,b_size,2*bfac*n_vec,b_size);
    qmessage("===prop %.3e GB, v %d  %.3e GB, buf %d  %.3e GB, bfac_group %d, E %.3e GB, smear_buf %1d. \n"
            , mem_prop, ncutgpu, memV*ncutgpu, ncutbuf, memV*ncutbuf, int(bfac_group), memV*n_vec, int(enable_smearE));
  }

  /////void checknormM(Ftype err=1e-3);
  void initialize_mass(const std::vector<double>& mass,Int Nv=12,Int one_minus_halfD_or=1, Int nprop=1);
  void initialize_mass();
  void alpha_clear();

  ///void seq_L(vector **prop3,double mass_d,double rho,int_vector &map,bool bcast,double_complex* ker_low);
  ///void seq_L(vector **prop3,double mass_d,double rho,int_vector &map,double_complex* ker_low);
  ///void setup_L(vector **prop2);

  void setup_gpufac(Int nprop=1);
  void allocate_GPU_mem(Int nprop=1);
  void clear_GPU_mem(Int cpu_also = 0);
  void print_norm_Mvec();

  ~eigen_ov()
  {
    clear_GPU_mem(1);
  }

};

void eigen_ov::setup_bfac(Long bsize0)
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

  if(bfac%6 !=0){qmessage("Cannot understand sites on bfac %5ld, bsize %10ld, Total %10ld !\n",
    bfac, b_size, noden*6);abort_r("");}
  if((noden*6)%bfac !=0){qmessage("Cannot understand sites on node*6/Nt %10ld bfac %10ld!\n",noden*6/Nt,bfac);abort_r("");}
  if(bfac%(Nt*6) !=0){qmessage("Cannot understand sites on node %10ld bfac %10ld, tsize %5d!\n",noden,bfac,Nt);abort_r("");}

}

eigen_ov::eigen_ov(const Geometry& geo_,Int n_vec_or, Long bsize0, double extra_mem_factor_set)
{
  geo = geo_;
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  ///if(fd.order_ch != 0){abort_r("Currently not supported for change fd order.\n ");}

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

  alpha.resize(0);
  alpha_buf.resize(0);
  alpha_bfac.resize(0);
  alpha_list.resize(0);eval_list.resize(0);

  ncutbuf = n_vec;
  ncutgpu = 0;
  bfac_group = -1;
  gpu_mem_set = false;

  nV_prop = 1;
  extra_mem_factor = extra_mem_factor_set;

  setup_bfac(bsize0);
  BFAC_GROUP_CPU = 1;
  //BFAC_GROUP_CPU = 2*bfac;

  Ns_buf = 0;
  nprop_buf = 0;

  fflush_MPI();

}

void eigen_ov::setup_gpufac(Int nprop)
{
  TIMERA("setup_gpufac");
  /////To avoid very large continuous memories on CPU
  bfac_group = 1;
  ////bfac_group = 2*bfac;
  nV_prop    = nprop;

  /////redefine variables on GPU
  #ifdef QLAT_USE_ACC
  double sm_factor = 1.0;
  if(enable_smearE){sm_factor = 0.5;}

  size_t freeM = 0;size_t totalM = 0;
  qacc_ErrCheck(qacc_MemGetInfo(&freeM,&totalM));
  //double freeD = 0;
  //freeD = freeM*pow(0.5,30);
  double totalD=0;
  totalD = totalM*pow(0.5,30);

  //int Ns = 12;
  Int Ns = nV_prop;if(Ns <=0)Ns = 2;
  ////long long Lat   = noden;
  double memV     = noden*12.0*sizeof(Complexq)*pow(0.5,30);
  double mem_prop = Ns * memV;////Ns have a factor of 12 already
  if(totalD < mem_prop + 6*12*memV){
    qmessage("===GPU Memory too small, Total %.3e, prop %5d %.3e, 6 prop %.3e, increase nodes! \n", totalD, Ns, mem_prop, 6*12*memV);
    Qassert(false);
  }

  Int vfac = ncutgpu;int vini = 8 * vfac;
  Int vres = int((totalD*extra_mem_factor*sm_factor - mem_prop )/memV); 
  if(qlat::get_id_node() != 0){vres=0;};sum_all_size(&vres, 1);
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
  for(Int bini=bfac_group;bini >= 1; bini--){
    size_t tem = get_threads(bini, bfac, 0);
    if(tem != 0 and bfac%tem == 0){bfac_group = tem;break;}
    bfac_group = 1;
  }}

  #endif
}

void eigen_ov::allocate_GPU_mem(Int nprop)
{
  setup_gpufac(nprop);
  #ifndef QLAT_USE_ACC
  return ;
  #endif

  /////if(ncutbuf <=0 or ncutgpu <= 0){return ;}
  //if(sm == 0){if(Eigenbuf.size() != 0 or Eigendyn.size() != 0){return ;}}
  //if(sm == 1){if(Eigenbuf_Sm.size() != 0 or Eigendyn_Sm.size() != 0){return ;}}

  LInt Lae = 2*bfac/bfac_group;
  LInt Lbe = bfac_group*ncutbuf*Long(b_size);
  LInt Lax = 2*bfac/bfac_group;
  LInt Lbx = bfac_group*ncutgpu*Long(b_size);

  if(Eigenbuf.size() == Lae and Eigendyn.size() == Lax){return ;}

  {
  TIMERA("allocate_GPU_mem");
  {
    if(ncutbuf != 0){allocate_buf(Eigenbuf, Lae, Lbe);}
    if(ncutgpu != 0){allocate_buf(Eigendyn, Lax, Lbx);}
  }

  { if(ncutgpu != 0 and enable_smearE == true){allocate_buf(Eigendyn_Sm, Lax, Lbx);}
    if(ncutbuf != 0 and enable_smearE == true){allocate_buf(Eigenbuf_Sm, Lae, Lbe);}}
  gpu_mem_set = true;
  //resize restore buffer positions
  npos_Eigenbuf = -1;npos_Eigendyn = -1;
  }
}

void eigen_ov::clear_GPU_mem(Int cpu_also)
{
  for(LInt i=0;i<Eigenbuf.size();i++){Eigenbuf[i].resize(0);}  Eigenbuf.resize(0);
  for(LInt i=0;i<Eigendyn.size();i++){Eigendyn[i].resize(0);}  Eigendyn.resize(0);
  for(LInt i=0;i<Eigenbuf_Sm.size();i++){Eigenbuf_Sm.resize(0);}  Eigenbuf_Sm.resize(0);
  for(LInt i=0;i<Eigendyn_Sm.size();i++){Eigendyn_Sm.resize(0);}  Eigendyn_Sm.resize(0);

  alpha.resize(0);
  alpha_buf.resize(0);
  alpha_bfac.resize(0);
  alpha_list.resize(0);eval_list.resize(0);

  /////ptmp.resize(0);stmp.resize(0);
  if(cpu_also == 1)
  {
    Mvec.resize(0);
    Mvec_Sm.resize(0);
  }

  gpu_mem_set = false;
}

void eigen_ov::load_eivals(const std::string& enamev,double rho_or,double Eerr, Int nini)
{
  rho = rho_or;

  std::vector<double > v,e;
  load_txt_eigenvalues(v,e, enamev.c_str());
  if((nini+n_vec)*2 > int(v.size()))abort_r("Eigen value size too small! ");
  eval_self.resize(n_vec);
  for(Int iv=0;iv<n_vec;iv++){
    Int off= iv + nini;
    eval_self[iv] = Complexq(v[off*2+0],v[off*2+1]);
    /////////cps base eigen value will not change to complex conjugate
    //////eval_self[iv] = qconj(Complexq(v[off*2+0],v[off*2+1]));
  }

  num_zero = 0;// initialize the number of zeros
  Eeigenerror = Eerr;
  Ftype rho_tem = rho;
  for(Int j=0; j<n_vec; ++j){
    eval_self[j] = eval_self[j]/rho_tem;
    if(std::sqrt(qnorm(eval_self[j])) < Eerr) num_zero += 1; // this will square the values
    //if(abs(eval_self[j]) < Eerr) num_zero += 1;
    //qmessage("n %3d, abs %.8e, norm %.8e, zero %3d \n", j, abs(eval_self[j]), qnorm(eval_self[j]), num_zero);
  }

}

template <typename Ty >
void eigen_ov::copy_FieldM_to_Mvec(Ty* src, Int ncur, Int sm, Int dir , bool data_GPU)
{
  TIMERA("COPY Eigen Vectors Mvec");
  if(ncur >= n_vec){abort_r("Cannot copy to position larger than n_vec ! \n");}
  /////int nread = nb - ba;
  Complexq* s1 = NULL;Ty* s0 = NULL;
  s0 = src;
  move_index mv_civ;

  Int GPU_cpy = 0;
  if(data_GPU){
    if(dir == 1){GPU_cpy = 3;} // from device to host
    if(dir == 0){GPU_cpy = 2;} // from host to device
  }

  //////move d,c to outter loop
  Long sizeF = noden;
  if(dir == 1){mv_civ.dojob(s0, s0, 1, 12, sizeF, 1, 1, data_GPU);}

  ////a factor of 2 by chiral
  LInt total = 2*bfac*b_size;
  for(LInt xini=0;xini < total/b_size; xini++)
  {
    LInt xi = xini*b_size;
    s1 = getEigenP(ncur, xi, sm, 1);
    /////qacc_Memcpy(s1, &s0[xi], b_size*sizeof(Complexq), qacc_MemcpyHostToDevice);
    if(dir == 1){cpy_data_thread(s1, &s0[xi], b_size, GPU_cpy, QFALSE);}
    if(dir == 0){cpy_data_thread(&s0[xi], s1, b_size, GPU_cpy, QFALSE);}
  }
  qacc_barrier(dummy);

  if(dir == 0){mv_civ.dojob(s0, s0, 1, 12, sizeF, 0, 1, data_GPU);}
  s0 = NULL; s1 = NULL;
}


template <typename Ty >
void eigen_ov::copy_FieldM_to_Mvec(qlat::FieldM<Ty , 12>& src, Int ncur, Int sm, Int dir )
{
  TIMERA("COPY Eigen Vectors Mvec");
  if(ncur >= n_vec){abort_r("Cannot copy to position larger than n_vec ! \n");}
  if(dir == 0){if(!src.initialized){
    src.init(geo);
  }}
  Ty* s0 = (Ty*) qlat::get_data(src).data();
  Int data_GPU = 0; ////do copies only from CPU here
  copy_FieldM_to_Mvec(s0, ncur, sm, dir, data_GPU);
}

void eigen_ov::copy_evec_to_GPU(Int nini)
{
  (void) nini;
  #ifndef QLAT_USE_ACC
  return ;
  #else

  Int mode_dyn = 0;
  if(nini < 0){abort_r("Copy initial negative! \n");}
  if(nini <  ncutbuf){mode_dyn = 0;}
  if(nini >= ncutbuf){mode_dyn = 1;}

  Int n0 = 0;int n1 = 0;
  if(mode_dyn == 0){n0=0;n1 = n0 + ncutbuf; }
  if(mode_dyn == 1){
    Int maxN = (nini + ncutgpu)/ncutgpu;
    n0=(maxN-1)*ncutgpu ;n1 = n0 + ncutgpu; 
  }

  if(mode_dyn == 0){if(npos_Eigenbuf == n1){return ;}else{npos_Eigenbuf = n1;}}
  if(mode_dyn == 1){if(npos_Eigendyn == n1){return ;}else{npos_Eigendyn = n1;}}

  {
  TIMERB("COPY Eigen Vectors to GPU");
  Complexq* s1;Complexq* s0;
  /////a factor of 2 by chiral
  LInt total = 2*bfac*b_size;
  for(Int ncur=n0;ncur<n1;ncur++)
  {
    if(ncur >= n_vec){break;}
    for(LInt xini=0;xini < total/b_size; xini++)
    {
      LInt xi = xini*b_size;
      /////qmessage("Test num ncur %d, xi %d, bfac %d \n", ncur, xi, bfac);

      s1 = getEigenP(ncur, xi, 0, 0);s0 = getEigenP(ncur, xi, 0, 1);
      //////qacc_Memcpy(s1, s0, b_size*sizeof(Complexq), qacc_MemcpyHostToDevice);
      cpy_data_thread(s1, s0 , b_size, 2, QFALSE);

      if(enable_smearE == true){
        s1 = getEigenP(ncur, xi, 1, 0);s0 = getEigenP(ncur, xi, 1, 1);
        /////qacc_Memcpy(s1, s0, b_size*sizeof(Complexq), qacc_MemcpyHostToDevice);
        cpy_data_thread(s1, s0 , b_size, 2, QFALSE);
      }
    }
  }
  qacc_barrier(dummy);
  }
  #endif

}

Complexq* eigen_ov::getEigenP(Int ni, size_t xi, Int sm, Int mode_initial )
{
  if(sm == 1 and enable_smearE == false){abort_r("Smeared Eigen system not loaded or computed! \n");}
  if(ni >= n_vec){abort_r("Request ni too large! \n");}

  Complexq* buf = NULL;
  Complexq* s0 = NULL;
  Int   chi = xi/(bfac*b_size);
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
  if(bfac_group <= 0){abort_r("Setup bfac_group first !\n");}
  // copy_evec_to_GPU(ni); // may need to check whether it's done or not ...

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
    Int ncur = ni - ncutbuf;
    Int na   = 0;
    Int nb   = ncur%ncutgpu;

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
void eigen_ov::smear_eigen(const std::string& Ename_Sm,
  const GaugeFieldT<Tg >& gf, const double width, const Int step,
  const CoordinateD& mom, const bool smear_in_time_dir)
{
  TIMER("smear eigen system");
  ////load if exist
  if(Ename_Sm != std::string("NONE") and get_file_size_MPI(Ename_Sm.c_str(), true) != 0){
    load_eigen_Mvec(Ename_Sm, 1);
    return ;
  }

  Long La = 2*bfac/BFAC_GROUP_CPU;
  Long Lb = BFAC_GROUP_CPU*n_vec*Long(b_size);
  print_mem_info("Before Eigen Memory Allocate");
  resize_EigenM(Mvec_Sm , La, Lb);enable_smearE = true;
  print_mem_info("Eigen Memory Allocate Done");

  const Int each = 12;
  Long Ncopy = geo.local_volume() * 12;
  qlat::vector_gpu<Complexq > buf;buf.resize(each * Ncopy);

  move_index mv_idx;
  std::vector<Long > job =  job_create(n_vec, each);
  for(LInt ji = 0; ji < job.size()/2 ; ji++)
  {
    Int flag = 0;
    ////copy to buf
    for(Int iv=0;iv<job[ji*2 + 1];iv++){copy_to_FieldM(&buf[iv*Ncopy], job[ji*2 + 0] + iv,  0, true);}
    flag = 0;mv_idx.dojob(buf.data(), buf.data(), 1, each , geo.local_volume()*12, flag, 1, true);

    smear_propagator_gwu_convension_inner<Complexq, 4,each  , Tg>(buf.data(), gf, width, step, mom, smear_in_time_dir);

    flag = 1;mv_idx.dojob(buf.data(), buf.data(), 1, each , geo.local_volume()*12, flag, 1, true);
    ////copy from buf
    for(Int iv=0;iv<job[ji*2 + 1];iv++){copy_FieldM_to_Mvec(&buf[iv*Ncopy], job[ji*2 + 0] + iv,  1, 1, true);}
  }

  ////erase smear
  {
  /////const SmearPlanKey& skey = get_smear_plan_key<Complexq, 4, each>(gf.geo(), smear_in_time_dir);
  ////qlat::clear_all_caches();
  get_smear_plan_cache().clear();
  }

  if(Ename_Sm != std::string("NONE")){
    save_eigen_Mvec(Ename_Sm, 1);
  }
}

void eigen_ov::save_eigen_Mvec(const std::string& ename, Int sm, Int save_type)
{
  if(sm == 1 and enable_smearE == false){qmessage("Could not save smear eigen without set it up.");return ;}
  io_vec& io_use = get_io_vec_plan(geo);

  const Int nini = 0;
  const Int ntotal = nini + n_vec;
  const bool read = false;
  inputpara in_write_eigen;
  FILE* file_write  = open_eigensystem_file(ename.c_str(), nini, ntotal, read , io_use , in_write_eigen , save_type);

  Int each = io_use.ionum;
  std::vector<qlat::FieldM<Complexq , 12> > buf;buf.resize(each);
  for(Int iv=0;iv<each;iv++){buf[iv].init(io_use.geo());}

  std::vector<Long > job =  job_create(n_vec, each);
  for(LInt ji = 0; ji < job.size()/2 ; ji++)
  {
    for(Int iv=0;iv<job[ji*2 + 1];iv++){copy_to_FieldM(buf[iv], job[ji*2 + 0] + iv, sm );}
    /////write to file
    load_eigensystem_vecs(file_write ,   buf, io_use , in_write_eigen , 0, job[ji*2 + 1]);
  }

  close_eigensystem_file(file_write , io_use , in_write_eigen );
  print_mem_info("Eigen Memory Write Done");
}

void eigen_ov::print_norm_Mvec()
{
  qlat::FieldM<Complexq , 12> buf;
  buf.init(geo);

  const Int Nsm = enable_smearE ? 2 : 1 ;
  for(Int sm = 0; sm < Nsm; sm++)
  for(Int ni = 0; ni < n_vec; ni++)
  {
    copy_to_FieldM(buf, ni, sm);
    Ftype* Psrc = (Ftype*) qlat::get_data(buf).data();
    Ftype normf = get_norm_vec(Psrc, noden);
    qmessage("Eigen vector %8d, sm %1d, norm 1.0 + %+.8e flag .", int(ni), sm, normf - 1.0);
  }
}

// some strange float eigensystem name with .s appendix
inline std::string get_eigen_name_string(const std::string& ename)
{
  Int found_file = 0;
  std::string name_read = ename;
  if(get_file_size_MPI(name_read, true) > 0){
    found_file = 1;
  };
  // check other file name
  if(found_file == 0){
    name_read = name_read + ".s";
    if(get_file_size_MPI(name_read, true) > 0){
      found_file = 1;
    }
    if(found_file == 0){
      qmessage("eigen not found %s \n", ename.c_str());
      Qassert(false);
    }
  }
  return name_read;
}

template <typename Tg >
void eigen_ov::load_eigen_Mvec_smear(const std::string& ename, 
    const GaugeFieldT<Tg >& gf,
    std::vector< qlat::GaugeFieldT<Ftype > >& gfL,
    std::vector<qlat::vector_gpu<qlat::ComplexT<Ftype > > >& propS,
    Int nini, Int checknorm, 
    const double src_width , const Int src_step ,
    const double sink_width, const Int sink_step,
    const CoordinateD& src_mom ,
    const CoordinateD& sink_mom,
    const bool src_smear_in_time_dir ,
    const bool sink_smear_in_time_dir
  )
{
  TIMER("load_eigen_Mvec_smear");
  std::string name_read = get_eigen_name_string(ename);
  Int ntotal = nini + n_vec;
  Ftype norm_err  = 1e-3;
  std::string val = get_env(std::string("q_eigen_norm_err"));
  if(val != ""){norm_err = stringtodouble(val);}
  Int print_norms = 0; 
  val = get_env(std::string("q_eigen_print_norm"));
  if(val != ""){print_norms = stringtonum(val);}

  io_vec& io_use = get_io_vec_plan(geo);
  const Int each_io = io_use.ionum;

  Long La = 2*bfac/BFAC_GROUP_CPU;
  Long Lb = BFAC_GROUP_CPU*n_vec*Long(b_size);
  print_mem_info("Before Eigen Memory Allocate");

  // setup smear tags
  resize_EigenM(Mvec    , La, Lb);
  resize_EigenM(Mvec_Sm , La, Lb);enable_smearE = true;

  print_mem_info("Eigen Memory Allocate Done");

  inputpara in_read_eigen;
  FILE* file_read  = open_eigensystem_file(name_read.c_str(), nini, ntotal, true , io_use , in_read_eigen , 2);

  std::vector<double > widthL = {sink_width, src_width};
  std::vector<Int    > stepL  = {sink_step, src_step};
  std::vector<CoordinateD > momL;momL.resize(2);
  momL[0] = sink_mom;momL[1] = src_mom;
  std::vector<bool > smear_in_time_dirL = {sink_smear_in_time_dir, src_smear_in_time_dir};

  {
  TIMER("load eigen from DISC");
  std::vector<qlat::FieldM<Complexq , 12> > buf ;buf.resize(each_io);
  for(Int iv=0;iv<each_io;iv++){buf[iv].init(io_use.geo());}
  std::vector<Long > job =  job_create(n_vec, each_io);
  for(LInt ji = 0; ji < job.size()/2 ; ji++)
  {
    ////int n0 = nini + job[ji*2 + 0]; Int n1 = n0 + job[ji*2 + 1]; 
    /////load from file
    load_eigensystem_vecs(file_read ,   buf, io_use , in_read_eigen , 0, job[ji*2 + 1]);
    ////copy to Mvec or Mvec_Sm
    for(Int iv=0;iv<job[ji*2 + 1];iv++){
      if(checknorm == 1){
        Ftype* Psrc = (Ftype*) qlat::get_data(buf[iv]).data();
        Ftype normf = get_norm_vec(Psrc, noden);
        if(print_norms == 1){
        qmessage("Eigen vector %8d, norm 1.0 + %+.8e flag . \n", int(job[ji*2 + 0] + iv + nini), normf - 1.0);}
        if(fabs(normf - 1.0) > norm_err){
          qmessage("Eigen vector %d, norm 1.0 + %.8e wrong. \n", int(job[ji*2 + 0] + iv + nini), normf - 1.0);
          abort_r("");
        }
      }

      copy_FieldM_to_Mvec(buf[iv], job[ji*2 + 0] + iv, 0  );
    }

  }
  }
  close_eigensystem_file(file_read , io_use , in_read_eigen );
  clear_io_vec_cache();

  {
  TIMER("smear eigen for src and sink");
  const Int each_sm = 12;
  const Long Ncopy = geo.local_volume() * 12;
  std::vector<qlat::vector_gpu<Complexq > > bufL;bufL.resize(2);
  std::vector<Long > job =  job_create(n_vec, each_sm);

  move_index mv_idx;
  Int flag = 0;
  for(LInt ji = 0; ji < job.size()/2 ; ji++)
  {
    ///do smearings and copy to correct place
    for(Int sm = 0 ; sm < 2; sm++)
    {
      if(Long(bufL[sm].size()) != Long(each_sm * Ncopy)){bufL[sm].resize(each_sm * Ncopy);}
      ////copy to buf
      for(Int iv=0;iv<job[ji*2 + 1];iv++){copy_to_FieldM(&bufL[sm][iv*Ncopy], job[ji*2 + 0] + iv,  0, true);}
      if(stepL[sm] == 0){continue;}

      Complexq* bufP = bufL[sm].data();
      flag = 0;mv_idx.dojob(bufP, bufP, 1, each_sm , geo.local_volume()*12, flag, 1, true);

      if(gfL.size() == 0){
        smear_propagator_gwu_convension_inner<Complexq, 4, each_sm  , Tg>(bufP, gf, 
          widthL[sm], stepL[sm], momL[sm], smear_in_time_dirL[sm]);
      }else{
        Qassert(smear_in_time_dirL[sm] == 0);//not transferred to 2link smearings and may not work
        smear_propagator_gwu_convension_2shift_modi<Complexq, Ftype, 4, each_sm>(bufP, geo, gfL, 
          widthL[sm], stepL[sm], propS, -1, 0, momL[sm]);
      }
      flag = 1;mv_idx.dojob(bufP, bufP, 1, each_sm , geo.local_volume()*12, flag, 1, true);
    }
    for(Int sm = 0 ; sm < 2; sm++)
    {
      for(Int iv=0;iv<job[ji*2 + 1];iv++){copy_FieldM_to_Mvec(&bufL[sm][iv*Ncopy], job[ji*2 + 0] + iv,  sm, 1, true);}
    }
  }
  }

  print_mem_info("Eigen Memory Load Done");
  //////mv_civ.free_mem();
  get_smear_plan_cache().clear();
}

void eigen_ov::load_eigen_Mvec(const std::string& ename, Int sm, Int nini, Int checknorm)
{
  std::string name_read = get_eigen_name_string(ename);

  Int ntotal = nini + n_vec;
  Ftype norm_err  = 1e-3;
  std::string val = get_env(std::string("q_eigen_norm_err"));
  if(val != ""){norm_err = stringtodouble(val);}
  Int print_norms = 0; 
  val = get_env(std::string("q_eigen_print_norm"));
  if(val != ""){print_norms = stringtonum(val);}

  io_vec& io_use = get_io_vec_plan(geo);

  Long La = 2*bfac/BFAC_GROUP_CPU;
  Long Lb = BFAC_GROUP_CPU*n_vec*Long(b_size);
  print_mem_info("Before Eigen Memory Allocate");
  if(sm == 0){resize_EigenM(Mvec    , La, Lb);}
  if(sm == 1){resize_EigenM(Mvec_Sm , La, Lb);enable_smearE = true;}
  print_mem_info("Eigen Memory Allocate Done");

  inputpara in_read_eigen;
  FILE* file_read  = open_eigensystem_file(name_read.c_str(), nini, ntotal, true , io_use , in_read_eigen , 2);

  Int each = io_use.ionum;
  std::vector<qlat::FieldM<Complexq , 12> > buf;buf.resize(each);
  for(Int iv=0;iv<each;iv++){buf[iv].init(io_use.geo());}

  std::vector<Long > job =  job_create(n_vec, each);
  for(LInt ji = 0; ji < job.size()/2 ; ji++)
  {
    ////int n0 = nini + job[ji*2 + 0]; Int n1 = n0 + job[ji*2 + 1]; 
    /////load from file
    load_eigensystem_vecs(file_read ,   buf, io_use , in_read_eigen , 0, job[ji*2 + 1]);
    ////copy to Mvec or Mvec_Sm
    for(Int iv=0;iv<job[ji*2 + 1];iv++){
      if(checknorm == 1 and sm == 0){
        Ftype* Psrc = (Ftype*) qlat::get_data(buf[iv]).data();
        Ftype normf = get_norm_vec(Psrc, noden);
        if(print_norms == 1){
        qmessage("Eigen vector %8d, norm 1.0 + %+.8e flag . \n", int(job[ji*2 + 0] + iv + nini), normf - 1.0);}
        if(fabs(normf - 1.0) > norm_err){
          qmessage("Eigen vector %d, norm 1.0 + %.8e wrong. \n", int(job[ji*2 + 0] + iv + nini), normf - 1.0);
          abort_r("");
        }
      }

      copy_FieldM_to_Mvec(buf[iv], job[ji*2 + 0] + iv, sm  );
    }
  }

  close_eigensystem_file(file_read , io_use , in_read_eigen );

  print_mem_info("Eigen Memory Load Done");

  //////mv_civ.free_mem();
}

void eigen_ov::load_eigen(const std::string& ov_evecname,
  Int checknorm, double kappa,double eigenerror, Int nini)
{
  TIMERB("=====Loading Eigen=====");

  std::string enamev = ssprintf("%s.eigvals", ov_evecname.c_str());
  qmessage("Vector File name: %s \n", ov_evecname.c_str() );
  qmessage("Values File name: %s \n", enamev.c_str());
  //////Load eigen values
  double rho_tem = 4 - 1.0/(2*kappa);
  load_eivals(enamev, rho_tem, eigenerror, nini);

  load_eigen_Mvec(ov_evecname, 0 ,nini, checknorm);

}

void eigen_ov::random_eigen(Int sm, Int seed)
{
  TIMERB("=====Loading random Eigen=====");
  eval_self.resize(n_vec);random_EigenM(eval_self, 0, seed + 10);
  
  Long La = 2*bfac/BFAC_GROUP_CPU;
  Long Lb = BFAC_GROUP_CPU*n_vec*Long(b_size);
  print_mem_info("Before Eigen Memory Allocate");
  if(sm == 0){
    resize_EigenM(Mvec    , La, Lb);
    for(Long iv=0;iv<La;iv++)random_Ty(Mvec[iv].data()   , Lb, 0, seed + 20 + iv);
  }
  if(sm == 1){
    resize_EigenM(Mvec_Sm , La, Lb);
    for(Long iv=0;iv<La;iv++)random_Ty(Mvec_Sm[iv].data(), Lb, 0, seed + 30 + iv);
    enable_smearE = true;
  }

  // set up bfac and zero mass 
  initialize_mass();
  print_mem_info("Eigen random Memory Load Done");

}

void eigen_ov::alpha_clear(){
  // clear alpha bufs which is needed before matrix prod
  zero_Ty(alpha.data(), alpha.size(), 1, QFALSE);
  zero_Ty(alpha_buf.data(), alpha_buf.size(), 1, QFALSE);
  zero_Ty(alpha_bfac.data(), alpha_bfac.size(), 1, QFALSE);
  zero_Ty(alpha_list.data(), alpha_list.size(), 1, QFALSE);
  qacc_barrier(dummy);
}

//////Nv source number, 12 for prop
void eigen_ov::initialize_mass(const std::vector<double>& mass, Int Ns, Int one_minus_halfD_or, Int nprop)
{
  Int need_update_size = 0;

  if(mass.size() != massL.size()){need_update_size = 1;}
  else{
    if(mass.size() == 0){need_update_size = 1;}
    for(unsigned int mi=0;mi<mass.size();mi++){
      if(mass[mi] != massL[mi]){need_update_size = 1;}
    }
    if(Ns != Ns_buf){need_update_size = 1;}
    if(one_minus_halfD_or != one_minus_halfD){need_update_size = 1;}
    if(nprop != nprop_buf){need_update_size = 1;}
  }
  // need_update_size = 1;

  if(need_update_size == 1){
    TIMERB("Set up store memories");
    massL = mass;
    one_minus_halfD = one_minus_halfD_or;
    Ns_buf = Ns;
    nprop_buf = nprop;

    ////std::vector<double> mass = mass_or;
    const Int mN = mass.size();

    // less the memory needed, allocate only large props needed
    #ifdef QLAT_USE_ACC
    if(gpu_mem_set == false or nprop > nV_prop + 12){
      qmessage("===mem %1d, n %5d, %5d, pos %5d %5d \n", int(gpu_mem_set), int(nprop), int(nV_prop), int(npos_Eigenbuf), int(npos_Eigendyn));
      allocate_GPU_mem(nprop);
    }
    #endif

    LInt nlarge = ncutgpu;if(ncutbuf > ncutgpu)nlarge = ncutbuf;

    alpha.resize(2*nlarge*Ns );
    alpha_buf.resize(alpha.size());
    alpha_bfac.resize(2*nlarge*Ns*bfac);
    alpha_list.resize(2*nlarge * mN*Ns);

    ////alpha.clear(false);alpha_bfac.clear(false);alpha_list.clear(false);

    //size_t tmp_Ls = 2*bfac * Ns*b_size;
    //size_t tmp_Lp = 2*mN*bfac * Ns*b_size;
    //stmp.resize(tmp_Ls);ptmp.resize(tmp_Lp);
    ////stmp.clear(false);ptmp.clear(false);

    if(eval_tem.size() != Long(mN*n_vec)){eval_tem.resize(mN*n_vec);}
    #pragma omp parallel for
    for(Int mki=0;mki< mN*n_vec;mki++)
    {
      Int mi = mki/n_vec;
      Int kn = mki%n_vec;
      eval_tem[mi*n_vec + kn] = inv_self(eval_self[kn], mass[mi], rho, one_minus_halfD);
    }

    if(Long(eval_list.size()) != Long( mN*n_vec)){eval_list.resize( mN*n_vec);}
    cpy_data_thread(eval_list.data(), (Complexq*) qlat::get_data(eval_tem).data(), eval_list.size(),  1, QFALSE);

    /////allocate and copy GPU memory, not needed,
    // only when use will need copy
    //copy_evec_to_GPU(0);

    const Long Nbfac = 2 * bfac;
    if(rpL.size() != Nbfac)rpL.resize(Nbfac);
    if(EpL.size() != Nbfac)EpL.resize(Nbfac);
    if(spL.size() != Nbfac)spL.resize(Nbfac);
  }

  alpha_clear();
}

void eigen_ov::initialize_mass()
{
  std::vector<double> mass;mass.push_back(0.0);
  initialize_mass(mass,12,one_minus_halfD);
}

/*
  conj_prop : 
    false --> from src to sink prop, 
    true from sink to src prop ? with gamma5's without sink gamma5 since Ns is not know
  tsrcL :
    construct only from src time t_src, other time-slice will be ignored
  result format :
    (2, bfac, mN, Ns, b_size)
    Nb_eigen_prod = bfac / ( 6 * Nt);
    ### Qassert(bfac % ( 6 * Nt) == 0);
    bfac --> [6, Nt, Nb_eigen_prod]
*/
void prop_L_device(eigen_ov& ei,Complexq *src, Complexq *props, Int Ns, const std::vector<double> &mass, Int mode_sm = 0, 
  const std::vector<Int >& tsrcL = std::vector<Int>(), const bool conj_prop = false, Int one_minus_halfD_or=1)
{
  TIMER_FLOPS("==prop_L");

  const Int mN   = mass.size();
  long long Lat = ei.noden;
  long long vGb = Lat*12;
  Int Fcount0 = 6 + 2;
  Int Fcount1 = 6 + 2;
  long long Tfloat = ei.n_vec*Ns*mN*vGb*Fcount1 + ei.n_vec*Ns*vGb*Fcount0;
  timer.flops += Tfloat;
  //double mem = Lat*12*(ei.n_vec + Ns + Ns*mN)*sizeof(Complexq);
  //qmessage("Memory size %.3e GB, %.3e Gflop \n", 
  //  mem/(1024.0*1024*1024), Tfloat/(1024.0*1024*1024));
  ///////qlat::get_num_node()

  ei.initialize_mass(mass, Ns, one_minus_halfD_or, Ns*mN);
  ////touch_GPU(ei.eval_list, ei.eval_list_size);

  Int sm0 = 0; Int sm1 = 0;
  if(mode_sm == 0){sm0 = 0; sm1 = 0;}//pt - pt
  if(mode_sm == 1){sm0 = 0; sm1 = 1;}//pt - sm
  if(mode_sm == 2){sm0 = 1; sm1 = 0;}//sm - pt
  if(mode_sm == 3){sm0 = 1; sm1 = 1;}//sm - sm

  Complexq* alpha       = ei.alpha.data();
  Complexq* alpha_buf   = ei.alpha_buf.data();
  Complexq* alpha_bfac  = ei.alpha_bfac.data();
  Complexq* alpha_list  = ei.alpha_list.data();
  Complexq* eval_list   = ei.eval_list.data();
  vector<Complexq* >& rpL = ei.rpL;
  vector<Complexq* >& EpL = ei.EpL;
  vector<Complexq* >& spL = ei.spL;

  const Geometry& geo = ei.geo;
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  const Long& bfac        = ei.bfac;
  const Long& b_size      = ei.b_size;

  const Int&  n_vec      = ei.n_vec;
  const Int& ncutgpu     = ei.ncutgpu;
  const Int& ncutbuf     = ei.ncutbuf;
  const Int& num_zero    = ei.num_zero;
  std::vector<bool > t_multi;t_multi.resize(fd.nt);
  Long Nb_eigen_prod = 0;
  if(tsrcL.size() != 0){
    Qassert(bfac % ( 6 * fd.Nt) == 0);
    Nb_eigen_prod = bfac / ( 6 * fd.Nt);
  }

  for(Int ti=0;ti<fd.nt;ti++){
    t_multi[ti] = true;
    if(tsrcL.size() != 0){
      // only do finit time slice
      bool find = false;
      for(unsigned int si=0;si<tsrcL.size();si++)
      {
        if(ti == tsrcL[si]){find = true;break;}
      }
      if(!find){t_multi[ti] = false;}
    }
  }

  // prop conj if no-zero
  if(conj_prop){
    Complexq* r = props;
    const Long Nd = 2 * bfac * mN * Ns * b_size;
    qacc_for(isp, Nd, {
      r[isp] = qconj(r[isp]);
    })
  }

  ////QBOOL dummy_test = QFALSE;

  Int Ng = 0;if(ncutgpu!=0){Ng = ei.n_vec/ncutgpu + 1;}
  Int nini = 0;
  /////each group have ncutgpu of vectors
  for(Int ng=0;ng<Ng + 1;ng++)
  {
    if(nini >= ei.n_vec)break;
    Int  ncur = 0;
    if(nini>=ncutbuf){
      ncur = ncutgpu;
    }else{ncur = ncutbuf;}

    ei.copy_evec_to_GPU(nini);

    {
    const Long m = ncur;
    const Long n = Ns;
    const Long w = b_size;

    TIMERA("prop low vec reduce");
    //TIMER_FLOPS("vec reduce");
    //long long vGb = 2*bfac*m*n*w;
    //int Fcount0   = 6 + 2;  
    //timer.flops  += vGb*Fcount0;

    Long count_b = 0;
    for(Long bini=0;bini<2*bfac;bini++){
      bool add_meas = true;
      if(Nb_eigen_prod != 0)
      {
        //const Long bini  = chi*bfac + (bsc*Nt + ti)*Nb_eigen_prod + bi;
        const Int ti =  ((bini % bfac) / Nb_eigen_prod) % fd.Nt;
        const Int t0 =  ti + fd.init;// current time slice
        if(!t_multi[t0]){
          add_meas = false;
        }
      }

      if(add_meas){
        rpL[count_b] = &alpha_bfac[(bini + 0)*m*n + 0];
        EpL[count_b] = ei.getEigenP(nini, bini*b_size, sm0);
        spL[count_b] = &src[       (bini + 0)*n*w + 0];
        count_b += 1;
      }
    }
    //if(Nb_eigen_prod == 0){Qassert(count_b == 2*bfac);}
    if(count_b > 0){
      /*
        !conj_prop normal mode with ev conj
         conj_prop then all without conj and conj at coef state
      */
      matrix_prodP(EpL.data(), spL.data(), rpL.data(), m,n, w , count_b, !conj_prop);
    }
    qacc_barrier(dummy);

    //std::vector<Long > jobA = job_create(2*bfac, ei.BFAC_GROUP_CPU);
    //if(nini > ncutbuf){jobA = job_create(2*bfac, ei.bfac_group);} //// vector size groups
    //for(LInt jobi=0;jobi < jobA.size()/2; jobi++)
    //{
    //  Long bini = jobA[jobi*2 + 0]; Long bcut = jobA[jobi*2+1];
    //  Complexq* rp = &alpha_bfac[(bini + 0)*m*n + 0];
    //  Complexq* Ep = ei.getEigenP(nini, bini*b_size, sm0);
    //  Complexq* sp = &src[       (bini + 0)*n*w + 0];
    //  //matrix_prod(Ep, sp, rp, m,n, w , bcut, true, dummy_test);
    //  matrix_prod(Ep, sp, rp, m,n, w , bcut, true);
    //}

    /////Output not conflict, GPU end here

    }

    {
    TIMERA("Reduce alpha");
    qacc_for(coff, Long(2*ncur*Ns),{
      Int chi = coff/(ncur*Ns);
      Long xi = coff%(ncur*Ns);
      if(!conj_prop){
        for(Long bi=0;bi<bfac;bi++){
          alpha[coff] += alpha_bfac[chi*bfac*ncur*Ns + bi*ncur*Ns + xi];
        }
      }else{
        const double factor_chi = chi ? -1.0 : 1.0;
        for(Long bi=0;bi<bfac;bi++){
          alpha[coff] += qconj(alpha_bfac[chi*bfac*ncur*Ns + bi*ncur*Ns + xi]) * factor_chi;
        }
      }
    });
    }

    {
    TIMERA("Global sum");
    sum_all_size((Ftype*) (ei.alpha.data()), (Ftype*) (ei.alpha_buf.data()), 2*(2*ncur*Ns), 1);
    }

    ////print_sum(ei.alpha.data(), ei.alpha.size(), "=====alpha", 1);

    {
    TIMERA("Get alpha list")
    const Complexq Iimag(0.0,1.0);
    const Complexq Two2(2.0,0.0);
    qacc_for(coff, Long(2*mN*Ns*ncur),{
      const Int chi   =  coff/(mN*Ns*ncur);
      const Int mi    = (coff/(Ns*ncur  ))%mN;
      const Int is    = (coff/(ncur     ))%Ns;
      const Int kn    = (coff             )%ncur;
      const long long offA = ((chi*mN+mi)*Ns+is)*ncur + kn;
      const Complexq a0 = alpha_buf[(chi*ncur+kn)*Ns+is];
      const Complexq a1 = alpha_buf[((1-chi)*ncur+kn)*Ns+is];

      if(kn + nini >= n_vec){
        alpha_list[offA] = 0.0;
      }else{
        Complexq li = eval_list[mi*n_vec + kn + nini];
        if(kn+nini >= num_zero){
          alpha_list[offA] = Two2*(li.real()*a0+Iimag*li.imag()*a1);
        }
        if(kn+nini  < num_zero){
          alpha_list[offA] = li*a0;
        }
      }
      if(conj_prop and chi){
        //const double factor_chi = chi ? -1.0 : 1.0;
        alpha_list[offA] *= -1.0;
      }
    });
    }


    {
    const Long m = mN*Ns;
    const Long n = b_size;
    const Long w = ncur;

    //TIMER("vec multi");
    TIMER_FLOPS("vec multi");
    long long vGb = 2*bfac*m*n*w;
    Int Fcount0   = 2*(3 + 1);  
    timer.flops  += vGb*Fcount0;

    if((nini + ncur) < ei.n_vec){
      ////zero_Ty(ei.alpha_bfac, ei.alpha_bfac_size, 0, false);zero_Ty(ei.alpha, ei.alpha_size, 0, false);
      //ei.alpha_bfac.clear(false);
      //ei.alpha.clear(false);
      zero_Ty(ei.alpha.data(), ei.alpha.size(), 1, QFALSE);
      zero_Ty(ei.alpha_bfac.data(), ei.alpha_bfac.size(), 1, QFALSE);
    }

    for(Long coff=0;coff<2*bfac;coff++){
      Long chi = coff/bfac;
      Long bi  = coff%bfac;
      rpL[coff] = &props[(chi*bfac+bi)*mN*Ns*b_size + 0];
      EpL[coff] = &alpha_list[(chi*mN*Ns+0)*w + 0];
      spL[coff] = ei.getEigenP(nini, coff*b_size, sm1);
    }
    matrix_prodP(EpL.data(), spL.data(), rpL.data(), m, n, w , 2*bfac, false, true);
    //for(Long coff=0;coff<2*bfac;coff++)
    //{
    //  Long chi = coff/bfac;
    //  Long bi  = coff%bfac;

    //  Complexq* rp = &props[(chi*bfac+bi)*mN*Ns*b_size + 0];
    //  ////Complexq* rp = &ei.ptmp[(chi*bfac+bi)*mN*Ns*b_size + 0];
    //  Complexq* ap = &alpha_list[(chi*mN*Ns+0)*w + 0];

    //  Complexq* Ep = ei.getEigenP(nini, coff*b_size, sm1);

    //  //matrix_prod(ap, Ep, rp, m, n, w ,1, false,  dummy_test, true);
    //  matrix_prod(ap, Ep, rp, m, n, w ,1, false, true);
    //}
    /////Output not conflict, GPU end here
    qacc_barrier(dummy);
    ///nini += ncur;

    }
    nini += ncur;
  }

  // prop conj 
  if(conj_prop){
    // (chi*bfac+bi)*mN*Ns*b_size + 0
    Complexq* r = props;
    const Long Nd = 2 * bfac * mN * Ns * b_size;
    qacc_for(isp, Nd, {
      r[isp] = qconj(r[isp]);
    })
  }

  //ei.print_info();
  ////print_sum(props, 2*bfac*  mN*Ns*b_size, "=====props", 1);
  /////untouch_GPU(ei.eval_list, ei.eval_list_size);
  //////#endif
}

/*
  Get low mode of propagtors or remove the low modes
  inputs propagators with 12 x 12, only one mass 
  Ngroup to shrink buffer usage
  Ty can be different than Complexq
  Always under Complexq type for eigensystem
  res and src can be the same if MassL is one
*/
template <typename Ty>
void ov_prop_L(std::vector<FieldG<Ty > >& res, std::vector<FieldG<Ty >>& src, eigen_ov& ei, 
  vector_gpu<Complexq>& buf_eig, const std::vector<double >& massL, Int mode_sm = 0, 
  const Int c_add = 0, const Int Ngroup_ = -1, const std::vector<Int >& tsrcL = std::vector<Int>(),
   const bool conj_prop = false, Int one_minus_halfD_or=1, const std::vector<Int > sec_info = std::vector<Int>()){

  // includes src vector 12
  const Int Ns = src.size();
  const Int Nmass = massL.size();
  Qassert(Ns > 0 and Nmass > 0);
  const Geometry& geo = src[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  const LInt total = 12 * fd.Nvol ;
  const Int Ngroup = Ngroup_ == -1 ? Ns : Ngroup_;
  const Int Ndc = 144;
  Qassert(c_add == 0 or c_add == -1 or c_add == 1);
  Qassert(Ngroup >= 1);

  for(LInt si=0;si<src.size();si++){
    Qassert(src[si].initialized and src[si].mem_order == QLAT_OUTTER and src[si].multiplicity == Ndc);
  }
  if(res.size() != src.size() * Nmass){res.resize(src.size() * Nmass);}
  for(LInt si=0;si<res.size();si++){
    if(!res[si].initialized){
      //res[si].init(geo, Ndc, QMGPU, QLAT_OUTTER);
      res[si].init_size(src[0]);
    }
    Qassert(res[si].initialized and res[si].mem_order == QLAT_OUTTER and res[si].multiplicity == Ndc);
    //clear
    //if(c_add == 0){
    //  set_zero(res[si]);
    //}
  }

  const Int Nbuf_T = (tsrcL.size() >= 2 and sec_info.size() == 2) ? 2 : 1;
  Int src_dT  = 0;
  Int src_ini = 0;
  vector<Int > src_time;
  vector<Int> src_t_order;
  if(Nbuf_T == 2){
    src_dT  = sec_info[0];
    src_ini = sec_info[1];
    vector<Int > map_sec;
    get_map_sec(map_sec, src_ini, src_dT, fd.nt);
    get_src_times(src_time, src_t_order, map_sec, src_ini, src_dT);
  }

  std::vector<std::vector<Int > > tbufL;tbufL.resize(tsrcL.size());
  if(Nbuf_T == 2){
    for(unsigned int ti=0;ti<tsrcL.size();ti++){
      tbufL[ti].resize(1);
      tbufL[ti][0] = tsrcL[ti];
    }
  }else{
    tbufL.resize(1);
    tbufL[0].resize(tsrcL.size());
    for(unsigned int ti=0;ti<tsrcL.size();ti++){
      tbufL[0][ti] = tsrcL[ti];
    }
  }

  // check whether time slice is for each sources
  if(Nbuf_T == 2){
    Qassert(fd.nt % src_dT == 0 and int(tsrcL.size()) == fd.nt / src_dT);
    Qassert(ei.bfac % ( 6 * fd.Nt) == 0);
    for(unsigned int ti=0;ti<tsrcL.size();ti++){
      const Int t0 = tsrcL[ti];
      Qassert( (t0 - src_ini + fd.nt) % src_dT == 0 );
    }
  }

  // continus memeory allocations
  buf_eig.resizeL(Ngroup * 12 * (1 + Nbuf_T * Nmass) * total);
  Complexq* buf_src  = &buf_eig[0 ];// size Ngroup * 12 * (1) * total
  Complexq* buf_res  = &buf_eig[Ngroup * 12 * (1 + 0) * total];// size Ngroup * 12 * (0 + Nmass) * total
  Complexq* buf_res1 = NULL;
  if(Nbuf_T == 2){buf_res1 = &buf_eig[Ngroup * 12 * (1 + Nmass) * total];}// size Ngroup * 12 * (0 + Nmass) * total
  const Int b_size = ei.b_size;

  // pointers for fields
  std::vector<FieldG<Ty > > srcP;
  std::vector<FieldG<Ty > > resP;
  std::vector<Long > jobA = job_create(Ns, Ngroup);
  for(LInt jobi=0;jobi < jobA.size()/2; jobi++)
  {
    Long bini = jobA[jobi*2 + 0]; Long bcut = jobA[jobi*2+1];
    const Int Nk = bcut * 12;
    if(Long(srcP.size()) != bcut){srcP.resize(0);srcP.resize(bcut);}
    if(Long(resP.size()) != Nmass * bcut){resP.resize(0);resP.resize(Nmass * bcut);}
    for(Int bi=0;bi<bcut;bi++){
      const Int bj = bini + bi;
      srcP[bi].set_pointer(src[bj]);
      for(Int mi=0;mi<Nmass;mi++){
        resP[mi*bcut + bi].set_pointer(res[mi*Ns + bj]);
      }
    }

    copy_bsize_prop_to_FieldP<Complexq, FieldG<Ty >, 0>(srcP, buf_src, Nk, b_size, fd, 1, false, 1);
    //buf_res.set_zero();
    zero_Ty(buf_res, bcut * 12 * Nmass * total, 1);

    /*
      actual low mode props
      mass --> Nsrc
    */
    //qmessage("Nt %5d \n", int(tbufL[0].size()));
    prop_L_device(ei, buf_src, buf_res , Nk, massL, mode_sm, tbufL[0], conj_prop, one_minus_halfD_or);
    if(Nbuf_T == 2){
      for(unsigned int ti=1;ti<tbufL.size();ti++){
        zero_Ty(buf_res1, bcut * 12 * Nmass * total, 1);
        prop_L_device(ei, buf_src, buf_res1, Nk, massL, mode_sm, tbufL[ti], conj_prop, one_minus_halfD_or);
        /* 
          copy results
          (2, bfac, mN, Ns, b_size)
          bfac --> [6, Nt, Nb_eigen_prod]
        */
        const Long bfac  = ei.bfac;
        const Long mN    = massL.size();
        const Long Ns    = Nk;
        const Long b_size = ei.b_size;
        const Long Nt    = fd.Nt;
        const Long Nb_eigen_prod = bfac / ( 6 * Nt );
        const Long Ndata = 2 * bfac * mN * Ns * b_size;
        const Long node_tini = fd.init;
        //const Long Nsrc_perT   = fd.nt / src_dT;
        //Qassert(Nsrc_perT == Long(tbufL.size()));
        const Int curr_src_time = tbufL[ti][0];

        qacc_for(isp, Ndata, {
          const Int bi = (isp / (mN * Ns * b_size) ) % bfac;
          const Int t0 = ( bi / Nb_eigen_prod ) % Nt + node_tini;
          //const Int seci = map_sec[t0];
          // map src number, first shit and each src have two sectors, then treat final sections
          //const Int src_num = (( seci + 1 ) / 2 ) % Nsrc_perT;
          //if(src_num == int(ti))
          if(src_time[t0] == curr_src_time)
          {
            buf_res[isp] = buf_res1[isp];
          }
        });
      }
    }
    //Ty norm0 = buf_res.norm2();
    //qmessage("Test %+.8e \n", norm0.real());

    // templates for equal, add, subtract
    if(c_add ==  0)copy_bsize_prop_to_FieldP<Complexq, FieldG<Ty >, 0>(resP, buf_res, Nk*Nmass, b_size, fd, 1, false, 0);
    if(c_add ==  1)copy_bsize_prop_to_FieldP<Complexq, FieldG<Ty >, 1>(resP, buf_res, Nk*Nmass, b_size, fd, 1, false, 0);
    if(c_add == -1)copy_bsize_prop_to_FieldP<Complexq, FieldG<Ty >,-1>(resP, buf_res, Nk*Nmass, b_size, fd, 1, false, 0);
    //for(LInt bi=0;bi<resP.size();bi++)
    //{
    //  double n0 = norm_FieldG(resP[bi]);
    //  qmessage("Test %+.8e \n", n0);
    //}
  }

}


}

#undef  Elocal
#undef  EIGENERROR
#undef  Vlocal

#endif

