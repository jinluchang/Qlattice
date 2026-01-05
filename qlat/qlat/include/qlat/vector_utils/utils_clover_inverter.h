#ifndef UTILS_CLOVER_INVERTER_H
#define UTILS_CLOVER_INVERTER_H

#pragma once

#include <quda.h>
#include <split_grid.h>
#include <tune_quda.h>
#include <deflation.h>
#include <invert_quda.h>

#include <cstdlib>
#include "utils_float_type.h"
#include "quda_para.h"
#include "general_funs.h"
#include "utils_io_vec.h"
#include "utils_eo_copies.h"

/////use inv_param.dirac_order = QUDA_INTERNAL_DIRAC_ORDER; when we want to work on GPU memeories

//  an overall sign for gamma5 compare to cps base
#define QUDA_GAMMA QUDA_DEGRAND_ROSSI_GAMMA_BASIS
// for internal quda dslash
//#define QUDA_GAMMA QUDA_UKQCD_GAMMA_BASIS

namespace qlat
{  //

struct quda_clover_inverter {
  Long V;
  Int X[4];
  QudaGaugeParam  gauge_param;
  qlat::vector<qlat::ComplexD > quda_gf_default;
  QudaInvertParam   inv_param;
  quda::SolverParam solverParam;
  //quda::Solver *solve_cg;
  //bool CG_reset;

  box<Geometry> geo;
  qlat::vector<Long > map_index;
  Int solve_mode ;
  /////QudaInvertParam df_param;

  QudaEigParam    eig_param;
  Int nvec;
  double inv_residue;
  double inv_time;
  Int    inv_iter;
  double inv_gflops;
  Int add_high;
  Int num_src;
  Int num_src_inv;
  Int quda_verbos;
  Int prec_type_check;

  quda::DiracMatrix* mat;
  //quda::DiracMatrix* mat_pc;
  //quda::DiracMatrix* mat_Mdag;
  //quda::DiracMatrix* mat_MMdag;
  //quda::DiracMatrix* mat_MdagM;

  //quda::DiracM* m_cg;
  //quda::DiracM* mSloppy;
  //quda::DiracM* mPre;
  //quda::DiracM* mEig;

  double mass_mat;
  //double mass_eig;
  double mass_value;

  ////eigen related
  quda::Dirac *dirac;
  quda::Dirac *dirac_pc;

  quda::Dirac *dirac_cg;
  quda::Dirac *dSloppy;
  quda::Dirac *dPre;
  quda::Dirac *dEig;

  Int spinor_site_size;
  qlat::vector<Long > mapI;// field copy index buffers

  std::vector<signed char > QUDA_clover;
  std::vector<signed char > QUDA_clover_inv;

  quda::ColorSpinorParam cs_cpu;
  quda::ColorSpinorParam cs_gpu;
  quda::ColorSpinorParam cs_gpuH;
  quda::ColorSpinorParam cs_gpuD;
  quda::ColorSpinorParam cs_gpuF;

  //=== split gauge buffers
  std::vector<Int> split_grid;
  Int split_mrh;
  QudaInvertParam   split_param;
  QudaInvertParam   split_param_copy;
  quda::ColorSpinorParam cs_param_split;
  quda::ColorSpinorParam devbuf_param;
  std::vector<std::vector<quda::ColorSpinorField> > split_fields;
  std::vector<vector_gpu<ComplexD > > split_fields_buffer;
  std::vector<quda::ColorSpinorField>    dev_buf;
  quda::CommKey split_key ;
  quda::Dirac *split_dirac;
  bool split_setup;
  bool split_comm_buff;
  Int num_sub_partition;
  inline void rotate_quda_split(const bool swap_gauge = true);
  inline void rotate_to_quda_split(){
    if(split_comm_buff == true){return ;}
    else{rotate_quda_split(true);}
  }
  inline void rotate_from_quda_split(){
    if(split_comm_buff == false){return ;}
    else{rotate_quda_split(true);}
  }
  inline void copy_to_buf_split(std::vector<quda::ColorSpinorField>& src, const Int i);
  inline void copy_from_buf_split(std::vector<quda::ColorSpinorField>& res, const Int i);
  inline void set_quda_split(const std::vector<Int >& grid, const Int num_src);
  inline void free_quda_split();
  //===

  //quda::ColorSpinorParam cs_copy;

  //quda::ColorSpinorField *csrc, *cres;
  //quda::ColorSpinorField *cpu_src, *cpu_res;
  quda::ColorSpinorField *gsrc, *gres;
  quda::ColorSpinorField *gquda0, *gquda1;
  quda::ColorSpinorField *gsrcH, *gresH;
  ///quda::ColorSpinorField *gsrcD, *gresD;
  quda::ColorSpinorField *gtmp1D, *gtmp2D, *gtmp3D, *gtmp4D, *gtmp5D;
  quda::ColorSpinorField *gtmp1F, *gtmp2F;
  quda::ColorSpinorField *gadd; ////for low mode addition
  //quda::ColorSpinorField *gtmp_invG0;

  quda::ColorSpinorField *ctmp0, *ctmp1, *ctmp2;
  quda::ColorSpinorField *gtmp0, *gtmp1, *gtmp2;
  //quda::ColorSpinorField *gsrcF, *gresF; ////format to match gpu solver
  ///bool singleE;

  bool io_rotate_bfac;

  ////0 for wilson, clover, 1 for stagger
  Int fermion_type;
  Int check_residue;
  bool clover_alloc;
  bool clover_setup;

  quda_clover_inverter(const Geometry& geo_, QudaTboundary t_boundary, Int num_src_=1);

  inline void free_mem();
  inline void setup_link(qlat::ComplexD* quda_gf);

  inline void setup_clover(const double kappa, const double clover_csw);
  inline void set_dirac_mrh(const Int mrh);

  template<typename Ty>
  inline void do_inv(Ty* res, Ty* src, const double kappa, const double err = 1e-10, const Int niter = 10000);

  /*
    copy to res under CPS basis
    from gquda0/gquda1 under internal quda gamma basis (QUDA_UKQCD_GAMMA_BASIS)
  */
  template<typename Ty>
  void copy_field_to_qlat(Ty* res, const int field_num = 0, const int dir = 1);

  template<typename Ty>
  void gamma_ukqcd_to_qlat(vector_cs<Ty >& eigen);

  template<typename Ty>
  void copy_field_from_qlat(Ty* src, const int field_num = 0)
  {
    copy_field_to_qlat(src, field_num, 0);
  }

  void print_plaq();

  //inline void setup_inv_kappa(const double kappa);
  inline void clear_mat();
  inline void free_csfield(const Int mode = 0);
  inline void alloc_csfield_cpu();
  inline void alloc_csfield_gpu();

  //inline void setup_inv_param_prec(Int prec_type = 0, bool force_reload = false);
  inline void setup_gauge_param(QudaTboundary t_boundary);

  inline void random_src(const Int seed);

  inline void save_prop(const void* srcP, const char* filename);

  ~quda_clover_inverter();

};

inline void quda_clover_inverter::setup_gauge_param(QudaTboundary t_boundary)
{
  TIMER("setup_gauge_param");
  for (Int mu = 0; mu < 4; mu++) {gauge_param.X[mu] = X[mu];}

  ////tuning flag for not high mode

  gauge_param.type = QUDA_WILSON_LINKS; //// or QUDA_SU3_LINKS
  // Slowest changing to fastest changing: even-odd, mu, x_cb_4d, row, column,
  // complex See the code later in this file to see the conversion between
  // Grid inde and Quda index.
  gauge_param.gauge_order = QUDA_MILC_GAUGE_ORDER;

  gauge_param.t_boundary = t_boundary; ////fermion_t_boundary
  //gauge_param.t_boundary = QUDA_PERIODIC_T; ////fermion_t_boundary
  //QUDA_ANTI_PERIODIC_T

  ////fine tune the precisions if needed
  gauge_param.cpu_prec               = QUDA_DOUBLE_PRECISION;
  gauge_param.cuda_prec              = QUDA_DOUBLE_PRECISION;

  gauge_param.cuda_prec_sloppy       = QUDA_SINGLE_PRECISION;
  gauge_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
  gauge_param.cuda_prec_eigensolver  = QUDA_SINGLE_PRECISION;

  //gauge_param.cuda_prec_sloppy       = QUDA_DOUBLE_PRECISION;
  //gauge_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
  //gauge_param.cuda_prec_eigensolver  = QUDA_DOUBLE_PRECISION;

  //gauge_param.reconstruct                   = QUDA_RECONSTRUCT_NO;  ////gauge_param.reconstruct = link_recon;
  //gauge_param.reconstruct_sloppy            = QUDA_RECONSTRUCT_NO;
  //gauge_param.reconstruct_precondition      = QUDA_RECONSTRUCT_NO;
  //gauge_param.reconstruct_eigensolver       = QUDA_RECONSTRUCT_NO;
  //gauge_param.reconstruct_refinement_sloppy = QUDA_RECONSTRUCT_NO;  /////link_recon

  gauge_param.reconstruct                   = QUDA_RECONSTRUCT_12;  ////gauge_param.reconstruct = link_recon;
  gauge_param.reconstruct_sloppy            = QUDA_RECONSTRUCT_12;
  gauge_param.location = QUDA_CPU_FIELD_LOCATION;
  gauge_param.reconstruct_precondition      = QUDA_RECONSTRUCT_NO;
  gauge_param.reconstruct_eigensolver       = QUDA_RECONSTRUCT_NO;
  gauge_param.reconstruct_refinement_sloppy = QUDA_RECONSTRUCT_NO;  /////link_recon

  gauge_param.anisotropy = 1.0;

  // For HISQ, this must always be set to 1.0, since the tadpole
  // correction is baked into the coefficients for the first fattening.
  // The tadpole doesn't mean anything for the second fattening
  // since the input fields are unitarized.
  gauge_param.tadpole_coeff = 1.0;

  gauge_param.ga_pad = 0;
  gauge_param.mom_ga_pad = 0;

  /////gauge fix paras
  gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;

  //////gauge_param.struct_size = sizeof(gauge_param);

  gauge_param.scale = 1.0;
  gauge_param.staggered_phase_type = QUDA_STAGGERED_PHASE_MILC;

  gauge_param.anisotropy = 1.0; /////anisotropy

  Int pad_size = 0;
  // For multi-GPU, ga_pad must be large enough to store a time-slice
  Int x_face_size = gauge_param.X[1] * gauge_param.X[2] * gauge_param.X[3] / 2;
  Int y_face_size = gauge_param.X[0] * gauge_param.X[2] * gauge_param.X[3] / 2;
  Int z_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[3] / 2;
  Int t_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[2] / 2;
  pad_size = std::max({x_face_size, y_face_size, z_face_size, t_face_size});
  gauge_param.ga_pad = pad_size;
  gauge_param.struct_size = sizeof(gauge_param);
  ////===END of gauge_param
}

quda_clover_inverter::quda_clover_inverter(const Geometry& geo_, QudaTboundary t_boundary, Int num_src_)
{
  TIMER("quda_clover_inverter_constuctor");
  /////set up gauge parameters
  geo.set(geo_);
  V = geo().local_volume();
  Qassert(num_src_ > 0);
  num_src = num_src_;
  num_src_inv = num_src;
  prec_type_check = -2;

  for (Int mu = 0; mu < 4; mu++) {X[mu] = geo().node_site[mu];}
  ////===Start of gauge_param
  gauge_param = newQudaGaugeParam();
  ////quda_gf_default = NULL;
  inv_param   = newQudaInvertParam();

  setup_gauge_param(t_boundary);

  add_high = 1;
  solve_mode = 0;

  clover_alloc = false;
  clover_setup = false;
  //csrc = NULL; cres = NULL;
  gsrc = NULL; gres = NULL;
  gquda0 = NULL;gquda1 = NULL;
  gsrcH= NULL; gresH= NULL;
  //gsrcD= NULL; gresD= NULL;
  gtmp1D = NULL;gtmp2D = NULL;gtmp3D = NULL;
  gtmp4D = NULL;gtmp5D = NULL;
  gtmp1F = NULL;gtmp2F = NULL;
  gadd = NULL;
  //gtmp_invG0 = NULL;

  ctmp0 = NULL; ctmp1 = NULL; ctmp2 = NULL;
  gtmp0 = NULL; gtmp1 = NULL; gtmp2 = NULL;

  mat = NULL;
  //mat_pc = NULL;
  //mat_Mdag  = NULL;
  //mat_MMdag = NULL;
  //mat_MdagM = NULL;

  //m_cg = NULL;
  //mSloppy = NULL;
  //mPre = NULL;
  //mEig = NULL;

  //CG_reset = true;
  //solve_cg = NULL;

  dirac = NULL;
  dirac_cg = NULL;
  dirac_pc = NULL;
  dSloppy = NULL;
  dPre    = NULL;
  dEig    = NULL;

  split_dirac = NULL;
  split_comm_buff = false;

  nvec = 0;
  check_residue = 0;

  spinor_site_size = 0;
  io_rotate_bfac = false;

  mass_mat = -1000000;
  ///mass_eig = -1000000;
  mass_value  = -100000;
  inv_time = 0.0;
  inv_iter = 0;
  inv_gflops = 0.0;

  quda_verbos = 0;
  std::string val = qlat::get_env(std::string("qlat_quda_verbos"));
  if(val != ""){quda_verbos = stringtonum(val);}
}

inline void quda_clover_inverter::setup_link(qlat::ComplexD* quda_gf)
{
  TIMER("setup_link");
  /////load gauge to quda GPU default position
  freeGaugeQuda();
  loadGaugeQuda((void *) quda_gf, &gauge_param);
  //if(apply_stag_phase == 1 and gauge_with_phase == false){
  //  applyGaugeFieldScaling_long((qlat::ComplexD*) quda_gf, V/2, &gauge_param, QUDA_STAGGERED_DSLASH);
  //  loadGaugeQuda((void *) quda_gf, &gauge_param);
  //  gauge_with_phase = true;
  //}
  //else{
  //  //print_plaq();
  //}
  ////quda_gf_default = (void *) quda_gf; //required to reload gauge with prec
  if(quda_gf != quda_gf_default.data()){
    quda_gf_default.resize(V *3*3*4 );
    cpy_data_thread(&quda_gf_default[0], quda_gf, V*3*3*4, false);
  }
}

inline void quda_clover_inverter::print_plaq()
{
  ///// Compute plaquette as a sanity check from interal GPU gauge
  double plaq[3];
  plaqQuda(plaq);
  QudaGaugeObservableParam param = newQudaGaugeObservableParam();
  param.compute_plaquette = QUDA_BOOLEAN_TRUE;
  param.compute_qcharge = QUDA_BOOLEAN_TRUE;
  gaugeObservablesQuda(&param);
  if(quda::comm_rank_global()== 0)printfQuda("Computed plaquette is %.8e (spatial = %.8e, temporal = %.8e), topological charge = %.8e \n",
      plaq[0], plaq[1], plaq[2],  param.qcharge);
}

////mode = 0 all field, mode = 1 cpu field , mode = 2 gpu field
inline void quda_clover_inverter::free_csfield(const Int mode)
{
  TIMER("free_csfield");
  if(mode == 0 or mode == 1){
  //if(csrc  != NULL){delete csrc;csrc=NULL;}if(cres != NULL){delete cres;cres=NULL;}
  if(ctmp0 != NULL){delete ctmp0;ctmp0=NULL;}
  if(ctmp1 != NULL){delete ctmp1;ctmp1=NULL;}
  if(ctmp2 != NULL){delete ctmp2;ctmp2=NULL;}
  }

  if(mode == 0 or mode == 2){
  //if(gsrcD != NULL){delete gsrcD;gsrcD=NULL;}if(gresD!= NULL){delete gresD;gresD=NULL;}
  if(gsrc  != NULL){delete gsrc ;gsrc =NULL;}if(gres != NULL){delete gres ;gres =NULL;}
  if(gsrcH != NULL){delete gsrcH;gsrcH=NULL;}if(gresH!= NULL){delete gresH;gresH=NULL;}
  if(gtmp0 != NULL){delete gtmp0;gtmp0=NULL;}

  if(gquda0  != NULL){delete gquda0 ;gquda0 =NULL;}if(gquda1 != NULL){delete gquda1 ;gquda1 =NULL;}

  if(gtmp1 != NULL){delete gtmp1;gtmp1=NULL;}
  if(gtmp2 != NULL){delete gtmp2;gtmp2=NULL;}
  if(gadd  != NULL){delete gadd;gadd=NULL;}
  }

  if(mode == 0 or mode == 3){
  if(gtmp1D!= NULL){delete gtmp1D;gtmp1D=NULL;}
  if(gtmp2D!= NULL){delete gtmp2D;gtmp2D=NULL;}
  if(gtmp3D!= NULL){delete gtmp3D;gtmp3D=NULL;}
  if(gtmp4D!= NULL){delete gtmp4D;gtmp4D=NULL;}
  if(gtmp5D!= NULL){delete gtmp5D;gtmp5D=NULL;}
  if(gtmp1F!= NULL){delete gtmp1F;gtmp1F=NULL;}
  if(gtmp2F!= NULL){delete gtmp2F;gtmp2F=NULL;}}

  //if(gtmp_invG0!= NULL){delete gtmp_invG0;gtmp_invG0=NULL;}
}

inline void quda_clover_inverter::alloc_csfield_gpu()
{
  TIMER("alloc_csfield_gpu");
  free_csfield(2);
  quda::ColorSpinorParam cs_tmp(cs_gpu);
  cs_tmp.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);
  cs_tmp.gammaBasis = QUDA_GAMMA;
  gsrc  = quda::ColorSpinorField::Create(cs_tmp);
  gres  = quda::ColorSpinorField::Create(cs_tmp);
  gtmp0 = quda::ColorSpinorField::Create(cs_tmp);
  gtmp1 = quda::ColorSpinorField::Create(cs_tmp);
  gtmp2 = quda::ColorSpinorField::Create(cs_tmp);

  quda::ColorSpinorParam cs_gpu_quda(cs_gpu);
  cs_gpu_quda.gammaBasis = QUDA_UKQCD_GAMMA_BASIS;// default quda dslash basis...
  gquda0 = quda::ColorSpinorField::Create(cs_gpu_quda);
  gquda1 = quda::ColorSpinorField::Create(cs_gpu_quda);

  cs_gpuH.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);
  cs_gpuH.gammaBasis = QUDA_GAMMA;
  //cs_gpuH.pc_type = QUDA_4D_PC;
  gsrcH = quda::ColorSpinorField::Create(cs_gpuH); 
  gresH = quda::ColorSpinorField::Create(cs_gpuH);
  //cs_tmp.setPrecision(QUDA_DOUBLE_PRECISION, QUDA_DOUBLE_PRECISION, true);
  //gsrcD = quda::ColorSpinorField::Create(cs_tmp);
  //gresD = quda::ColorSpinorField::Create(cs_tmp);
}

inline void quda_clover_inverter::alloc_csfield_cpu()
{
  TIMER("alloc_csfield_cpu");
  free_csfield(1);
  //bool pc_solution = (inv_param.solution_type == QUDA_MATPC_SOLUTION) ||
  //  (inv_param.solution_type == QUDA_MATPCDAG_MATPC_SOLUTION);
  bool pc_solution = false;
  void* temV = NULL;
  quda::GaugeField cpuGauge(gauge_param);
  quda::ColorSpinorParam cpuParam_tem(temV, inv_param, cpuGauge.X(), pc_solution, inv_param.input_location);
  cs_cpu = quda::ColorSpinorParam(cpuParam_tem);
  //cs_cpu.setPrecision(inv_param.cpu_prec);
  cs_cpu.setPrecision(QUDA_DOUBLE_PRECISION); ////double for all cpu cs field
  cs_cpu.create = QUDA_ZERO_FIELD_CREATE;
  cs_cpu.location = QUDA_CPU_FIELD_LOCATION;
  cs_cpu.is_composite  = true;
  cs_cpu.is_component  = false;
  cs_cpu.composite_dim = num_src;

  cs_gpu = quda::ColorSpinorParam(cs_cpu);cs_gpu.location = QUDA_CUDA_FIELD_LOCATION;
  cs_gpu.create = QUDA_ZERO_FIELD_CREATE;
  //cs_gpu.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);
  cs_gpu.setPrecision(QUDA_DOUBLE_PRECISION, QUDA_DOUBLE_PRECISION, true); ////double for all cpu cs field
  ////cs_gpu.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec_eigensolver, true);

  cs_gpu.siteSubset = QUDA_FULL_SITE_SUBSET;
  //if (inv_param.solution_type == QUDA_MAT_SOLUTION || inv_param.solution_type == QUDA_MATDAG_MAT_SOLUTION) {
  //  cs_gpu.siteSubset = QUDA_FULL_SITE_SUBSET;
  //} else {
  //  cs_gpu.siteSubset = QUDA_PARITY_SITE_SUBSET;
  //  cs_gpu.x[0] /= 2;
  //}

  //cs_gpu.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec_eigensolver, true);
  //cs_gpu.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);
  //cs_gpu.setPrecision(inv_param.cuda_prec_eigensolver, inv_param.cuda_prec_eigensolver, true);
  //if(cs_gpu.nSpin != 1)
  {
    cs_gpu.gammaBasis = QUDA_GAMMA;
    cs_cpu.gammaBasis = QUDA_GAMMA;
    //cs_gpu.pc_type = QUDA_4D_PC;
    //cs_cpu.pc_type = QUDA_4D_PC;
  }

  //csrc  = quda::ColorSpinorField::Create(cs_cpu);
  //cres  = quda::ColorSpinorField::Create(cs_cpu);
  ctmp0 = quda::ColorSpinorField::Create(cs_cpu);
  ctmp1 = quda::ColorSpinorField::Create(cs_cpu);
  ctmp2 = quda::ColorSpinorField::Create(cs_cpu);

  //bool single_file=true;
  if(fermion_type == 0){pc_solution = false;}
  if(fermion_type == 1){pc_solution = true ;}
  quda::ColorSpinorParam gpuParam_tem(temV, inv_param, cpuGauge.X(), pc_solution, QUDA_CUDA_FIELD_LOCATION);

  cs_gpuH = gpuParam_tem;
  cs_gpuH.create = QUDA_ZERO_FIELD_CREATE;
  cs_gpuH.is_composite  = true;
  cs_gpuH.is_component  = false;
  cs_gpuH.composite_dim = num_src;


  cs_gpuF= quda::ColorSpinorParam(gpuParam_tem);
  cs_gpuF.create = QUDA_ZERO_FIELD_CREATE;
  cs_gpuF.setPrecision(QUDA_SINGLE_PRECISION, QUDA_SINGLE_PRECISION, true);
  cs_gpuF.is_composite  = false;
  cs_gpuF.is_component  = false;

  ////cs_gpuF.setPrecision(QUDA_DOUBLE_PRECISION, QUDA_DOUBLE_PRECISION, true);

  cs_gpuD= quda::ColorSpinorParam(gpuParam_tem);
  cs_gpuD.is_composite  = false;
  cs_gpuD.is_component  = false;
  cs_gpuD.create = QUDA_ZERO_FIELD_CREATE;
  cs_gpuD.setPrecision(QUDA_DOUBLE_PRECISION, QUDA_DOUBLE_PRECISION, true);
}

/////double prec prop save
inline void quda_clover_inverter::save_prop(const void* srcP, const char* filename)
{
  TIMER("quda save_prop");
  const Int n0 = 1;
  std::vector<qlat::FieldM<qlat::ComplexD , 3> > prop;prop.resize(n0);
  for(Int n = 0; n < prop.size(); n++){prop[n].init(geo());}

  qlat::ComplexD* src = (qlat::ComplexD*) srcP;
  Long Nvol = geo().local_volume() * n0 ;

  for(Int n=0;n<n0;n++){
    quda_cf_to_qlat_cf_test(prop[n], &src[n*Nvol]);
  }   

  std::string VECS_TYPE("STAGGERED_Prop");
  std::string INFO_LIST = ssprintf("mass %.8f", inv_param.mass);

  bool single_file = false;
  bool read = false;
  qlat::load_qlat_noisesT(filename, prop, read, single_file, VECS_TYPE, std::string("NONE"), 0, n0, false);
}

inline void quda_clover_inverter::free_mem(){
  TIMER("quda free_mem");
  V = 0;
  for(Int i=0;i<4;i++){X[i] = 0;}

  free_csfield();
  if(clover_setup == true){freeCloverQuda();}

  QUDA_clover.resize(0);
  QUDA_clover_inv.resize(0);

  nvec = 0;

  clear_mat();
}

template<typename Ty>
void quda_clover_inverter::gamma_ukqcd_to_qlat(vector_cs<Ty >& eigen){
  TIMER("gamma_ukqcd_to_qlat");
  for(int vi=0;vi<eigen.nvec;vi++){
    eigen.copy_to(  (ComplexT<double>*) (*gquda0).data(), vi);
    *gsrc = *gquda0;
    quda_cf_to_qlat_cf_pointer((ComplexT<double>*) (*gres).data(), (ComplexT<double>*) (*gsrc).data(), 12, geo(), mapI );
    eigen.copy_from((ComplexT<double>*) (*gres).data(), vi);
  }
}

template<typename Ty>
void quda_clover_inverter::copy_field_to_qlat(Ty* res, const int field_num, const int dir )
{
  TIMER("copy_field_to_qlat");
  if(field_num != 0 and field_num != 1){Qassert(false);}
  if(dir != 0 and dir != 1){Qassert(false);}
  if(dir == 1){
    if(field_num == 0){*gsrc = *gquda0;}
    if(field_num == 1){*gsrc = *gquda1;}
    quda_cf_to_qlat_cf_pointer(res, (ComplexT<double>*) (*gsrc).data(), 12, geo(), mapI );
  }

  if(dir == 0){
    qlat_cf_to_quda_cf_pointer((ComplexT<double>*) (*gsrc).data(), res, 12, geo(), mapI );
    if(field_num == 0){*gquda0 = *gsrc;}
    if(field_num == 1){*gquda1 = *gsrc;}
  }
}

inline void quda_clover_inverter::setup_clover(const double kappa, const double clover_csw)
{
  if(clover_setup  == false or kappa != inv_param.kappa or inv_param.clover_csw != clover_csw)
  {
    if(clover_setup == true){freeCloverQuda();}
    fermion_type = 0;
    spinor_site_size = 12;

    /////===Start of Inv parameters
    if(clover_csw != 0){
    inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH;
    }else{
    inv_param.dslash_type = QUDA_WILSON_DSLASH;
    }

    /////double kappa      = kappa;
    double anisotropy = gauge_param.anisotropy;
    inv_param.kappa = kappa;
    inv_param.mass = 0.5 / inv_param.kappa - (1.0 + 3.0 / anisotropy);

    // Use 3D or 4D laplace
    //===inv_param.laplace3D = laplace3D;
    inv_param.Ls = 1;

    inv_param.cpu_prec                      = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec                     = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec_sloppy              = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_refinement_sloppy   = QUDA_SINGLE_PRECISION;

    inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
    inv_param.gamma_basis = QUDA_GAMMA;

    //inv_param.dirac_order = QUDA_DIRAC_ORDER;
    inv_param.dirac_order = QUDA_INTERNAL_DIRAC_ORDER;

    inv_param.clover_cpu_prec               = QUDA_DOUBLE_PRECISION;
    inv_param.clover_cuda_prec              = QUDA_DOUBLE_PRECISION;

    inv_param.clover_cuda_prec_sloppy       = QUDA_SINGLE_PRECISION;
    inv_param.clover_cuda_prec_refinement_sloppy = QUDA_SINGLE_PRECISION;

    ////related to eigensystem
    inv_param.cuda_prec_eigensolver         = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_precondition        = QUDA_SINGLE_PRECISION;
    inv_param.clover_cuda_prec_precondition = QUDA_SINGLE_PRECISION;
    inv_param.clover_cuda_prec_eigensolver  = QUDA_SINGLE_PRECISION;

    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
    // Use kappa * csw or supplied clover_coeff
    bool compute_clover_trlog = false;
    //bool compute_clover_trlog = true;
    inv_param.clover_csw = clover_csw;
    inv_param.clover_coeff = inv_param.clover_csw * inv_param.kappa;
    /////===unknow para
    inv_param.compute_clover_trlog = compute_clover_trlog ? 1 : 0;

    // General parameter setup
    inv_param.inv_type = QUDA_CG_INVERTER;

    // solution_type specifies *what* system is to be solved.
    // solve_type specifies *how* the system is to be solved.
    //
    // We have the following four cases (plus preconditioned variants):
    //
    // solution_type    solve_type    Effect
    // -------------    ----------    ------
    // MAT              DIRECT        Solve Ax=b
    // MATDAG_MAT       DIRECT        Solve A^dag y = b, followed by Ax=y
    // MAT              NORMOP        Solve (A^dag A) x = (A^dag b)
    // MATDAG_MAT       NORMOP        Solve (A^dag A) x = b
    // MAT              NORMERR       Solve (A A^dag) y = b, then x = A^dag y
    //
    // We generally require that the solution_type and solve_type
    // preconditioning match.  As an exception, the unpreconditioned MAT
    // solution_type may be used with any solve_type, including
    // DIRECT_PC and NORMOP_PC.  In these cases, preparation of the
    // preconditioned source and reconstruction of the full solution are
    // taken care of by Dirac::prepare() and Dirac::reconstruct(),
    // respectively.

    // Eventually we want the unpreconditioned solution.
    inv_param.solution_type = QUDA_MAT_SOLUTION;

    // NORMOP_PC means preconditioned normal operator MdagM
    //inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
    inv_param.solve_type = QUDA_NORMOP_SOLVE;

    // There might be a performance difference.
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;

    inv_param.dagger = QUDA_DAG_NO;

    //inv_param.mass_normalization   = QUDA_KAPPA_NORMALIZATION;
    inv_param.mass_normalization   = QUDA_MASS_NORMALIZATION;
    inv_param.solver_normalization =  QUDA_SOURCE_NORMALIZATION;
    //inv_param.solver_normalization = QUDA_DEFAULT_NORMALIZATION;
    //inv_param.solver_normalization = QUDA_KAPPA_NORMALIZATION;

    ///unknown
    Int gcrNkrylov = 10;
    QudaCABasis ca_basis = QUDA_POWER_BASIS;
    double ca_lambda_min = 0.0;
    double ca_lambda_max = -1.0;

    inv_param.pipeline = 1;
    inv_param.Nsteps = 2;
    inv_param.gcrNkrylov = gcrNkrylov;
    inv_param.ca_basis = ca_basis;
    inv_param.ca_lambda_min = ca_lambda_min;
    inv_param.ca_lambda_max = ca_lambda_max;

    inv_param.tol     = 1e-10;
    inv_param.maxiter = 100000;
    inv_param.tol_restart = 0.0005;
    //if (tol_hq == 0 && tol == 0) {
    //  errorQuda("qudaInvert: requesting zero residual\n");
    //  exit(1);
    //}

    // require both L2 relative and heavy quark residual to determine convergence
    inv_param.residual_type = static_cast<QudaResidualType_s>(0);
    inv_param.residual_type = static_cast<QudaResidualType_s>(inv_param.residual_type | QUDA_L2_RELATIVE_RESIDUAL);

    inv_param.tol_hq = 1e-5; // specify a tolerance for the residual for heavy quark residual
    //inv_param.tol_hq = 1e-10; // specify a tolerance for the residual for heavy quark residual

    //// Offsets used only by multi-shift solver
    //// These should be set in the application code. We set the them here by way of
    //// example
    //inv_param.num_offset = multishift;
    //for (Int i = 0; i < inv_param.num_offset; i++) inv_param.offset[i] = 0.06 + i * i * 0.1;
    //// these can be set individually
    //for (Int i = 0; i < inv_param.num_offset; i++) {
    //  inv_param.tol_offset[i] = inv_param.tol;
    //  inv_param.tol_hq_offset[i] = inv_param.tol_hq;
    //}

    // This is for Quda's sophisticated reliable update. 0.1 should be good.
    inv_param.reliable_delta = 0.1;
    inv_param.use_alternative_reliable = true;
    inv_param.use_sloppy_partial_accumulator = 0;
    inv_param.solution_accumulator_pipeline = 1;
    inv_param.max_res_increase = 1;

    // domain decomposition preconditioner parameters
    // inv_param.inv_type_precondition = precon_type;

    //inv_param.schwarz_type = precon_schwarz_type;
    //inv_param.precondition_cycle = precon_schwarz_cycle;
    //inv_param.tol_precondition = tol_precondition;
    //inv_param.maxiter_precondition = maxiter_precondition;
    //inv_param.omega = 1.0;

    inv_param.input_location  = QUDA_CUDA_FIELD_LOCATION;
    inv_param.output_location = QUDA_CUDA_FIELD_LOCATION;

    //inv_param.input_location  = QUDA_CUDA_FIELD_LOCATION;
    //inv_param.output_location = QUDA_CUDA_FIELD_LOCATION;

    // QUDA_DEBUG_VERBOSE is too nasty.

    inv_param.extlib_type = QUDA_EIGEN_EXTLIB;

    // Whether or not to use native BLAS LAPACK
    //inv_param.native_blas_lapack = (native_blas_lapack ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE);
    inv_param.native_blas_lapack = QUDA_BOOLEAN_TRUE;

    //// Whether or not use fused kernels for Mobius
    ////inv_param.use_mobius_fused_kernel = use_mobius_fused_kernel ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
    //inv_param.use_mobius_fused_kernel = QUDA_BOOLEAN_FALSE;


    ///////===check this variable
    std::array<int, 4> grid_partition = {1, 1, 1, 1};
    inv_param.split_grid[0] = grid_partition[0];
    inv_param.split_grid[1] = grid_partition[1];
    inv_param.split_grid[2] = grid_partition[2];
    inv_param.split_grid[3] = grid_partition[3];

    inv_param.struct_size = sizeof(inv_param);

    ////===END of Inv parameters

    {
    //////===operator define
    //void *clover = nullptr;
    //void *clover_inv = nullptr;

    Int clover_site_size           = 72; // real numbers per block-diagonal clover matrix
    ////size_t host_clover_data_type_size = (cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
    Int host_clover_data_type_size = sizeof(quda::Complex)/2;
    ////size_t host_spinor_data_type_size = (cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
    Int host_spinor_data_type_size = sizeof(quda::Complex)/2;

    QUDA_clover.resize(    V * clover_site_size * host_clover_data_type_size);
    QUDA_clover_inv.resize(V * clover_site_size * host_spinor_data_type_size);

    //constructHostCloverField((void*) quda_clover.data(), (void*) quda_clover_inv.data(), inv_param);
    /////===host

    bool compute_clover = true;
    if(inv_param.clover_csw == 0){compute_clover = false;}
    //double norm = 0.00; // clover components are random numbers in the range (-norm, norm)
    //double diag = 1.0;  // constant added to the diagonal
    ////if (!compute_clover) constructQudaCloverField((void*) quda_clover.data(), norm, diag, inv_param.clover_cpu_prec, V);

    //loadCloverQuda((void*) quda_clover.data(), (void*)  quda_clover_inv.data(), &inv_param);
    //if(in.clover_csw == 0){
    //inv_param.compute_clover_inverse = 0;
    //inv_param.return_clover_inverse  = 0;
    //}
    /////===host
    ///// Load the clover terms to the device
    //loadCloverQuda((void*) quda_clover.data(), (void*)  quda_clover_inv.data(), &inv_param);
    //loadCloverQuda((void*) quda_clover.data(), (void*)  quda_clover_inv.data(), &inv_param);
    inv_param.compute_clover = compute_clover;
    if(compute_clover){
      inv_param.return_clover = 1;
      inv_param.compute_clover_inverse = 1;
      inv_param.return_clover_inverse  = 1;

      loadCloverQuda((void*) QUDA_clover.data(), (void*)  QUDA_clover_inv.data(), &inv_param);
    }else{
      inv_param.compute_clover_inverse = 0;
      inv_param.return_clover_inverse  = 0;
    }
    }
    clover_setup = true;
    if(clover_alloc == false){
      alloc_csfield_cpu();
      alloc_csfield_gpu();
      clover_alloc = true;
    }

    Int verbos = quda_verbos;
    if(verbos <= -1)
    {
      inv_param.verbosity   = QUDA_SILENT;
      setVerbosity(QUDA_SILENT);
    }

    if(verbos == 0)
    {
      inv_param.verbosity   = QUDA_SUMMARIZE;
    }
    if(verbos == 1)
    {
      setVerbosity(QUDA_VERBOSE);
      inv_param.verbosity   = QUDA_VERBOSE;
    }

    if(verbos == 2)
    {
      setVerbosity(QUDA_DEBUG_VERBOSE);
      inv_param.verbosity   = QUDA_VERBOSE;
    }
    if(dirac != NULL){delete dirac;dirac = NULL;}
    const bool pc_solve = false;
    quda::DiracParam diracParam;
    setDiracParam(diracParam, &inv_param, pc_solve);
    dirac     = quda::Dirac::create(diracParam);
  }
}

inline void quda_clover_inverter::set_dirac_mrh(const Int mrh)
{
  if(dirac != NULL){delete dirac;dirac = NULL;}
  inv_param.num_src = mrh;
  const bool pc_solve = false;
  quda::DiracParam diracParam;
  setDiracParam(diracParam, &inv_param, pc_solve);
  dirac     = quda::Dirac::create(diracParam);
}

inline void quda_clover_inverter::free_quda_split(){
  split_fields.resize(0);
     dev_buf.resize(0);
  split_fields_buffer.resize(0);
  split_setup = false;
  if(split_dirac != NULL){
    // swap current gauges with the split ones if split ones are not null
    swapGaugeSplit(true);
    // free the split gauges and swap the split gauges with the buffered gauges, if split gauges are not      null
    swapGaugeSplit(false);
    delete split_dirac;split_dirac = NULL;
  }
}

inline void quda_clover_inverter::clear_mat()
{
  TIMER("clear_mat");
  if(dirac != NULL){delete dirac;dirac=NULL;}
  if(dirac_cg != NULL){delete dirac_cg;dirac_cg=NULL;}
  if(dirac_pc != NULL){delete dirac_pc;dirac_pc=NULL;}
  //if(split_dirac != NULL){delete split_dirac;split_dirac = NULL;}
  free_quda_split();
  if(dSloppy != NULL){delete dSloppy;dSloppy=NULL;}
  if(dPre != NULL){delete dPre;dPre=NULL;}
  if(dEig != NULL){delete dEig;dEig=NULL;}

  if(mat       != NULL){delete mat      ;  mat       = NULL;}
  //if(mat_pc    != NULL){delete mat_pc   ;  mat_pc    = NULL;}
  //if(mat_Mdag  != NULL){delete mat_Mdag ;  mat_Mdag  = NULL;}
  //if(mat_MMdag != NULL){delete mat_MMdag;  mat_MMdag = NULL;}
  //if(mat_MdagM != NULL){delete mat_MdagM;  mat_MdagM = NULL;}

  //if(m_cg != NULL){delete m_cg;  m_cg = NULL;}
  //if(mSloppy != NULL){delete mSloppy;  mSloppy = NULL;}
  //if(mPre != NULL){delete mPre;  mPre = NULL;}
  //if(mEig != NULL){delete mEig;  mEig = NULL;}
}

inline void quda_clover_inverter::random_src(const Int seed)
{
  for(Int i=0;i<int(seed)%20 + 20*quda::comm_rank_global();i++){quda::comm_drand();}
  Int random_mode = 1;
  random_Ty((qlat::ComplexD*) gsrc->data(), gsrc->Volume() * spinor_site_size, 1, seed + 111111, random_mode);
  quda::blas::ax(0.05, *gsrc);
  quda::Complex* tmp = (quda::Complex*) (gsrc->data());
  Long totalN = gsrc->Volume() * spinor_site_size;
}

template<typename Ty>
inline void quda_clover_inverter::do_inv(Ty* res, Ty* src, const double kappa, const double err, const Int niter )
{
  TIMER_FLOPS("QUDA CG");
  timeval tm0,tm1;gettimeofday(&tm0, NULL);gettimeofday(&tm1, NULL);

  //setup_inv_param_prec(prec_type); ////restore precisions
  inv_param.tol = err;
  ///if(err < 1e-6 and err > 1e-10){inv_param.tol_restart    = err*5e+3;}
  inv_param.maxiter = niter;
  setup_clover(kappa, inv_param.clover_csw);

  Qassert((void*) (*gsrc).data() != (void*) src);
  qlat_cf_to_quda_cf_test((*gsrc), src, geo(), map_index);

  //ctmp0 = quda::ColorSpinorField::Create(cs_cpu);
  //ctmp1 = quda::ColorSpinorField::Create(cs_cpu);
  //qudaMemcpy((void*)ctmp0->data(), (void*)(*gsrc).data(), ctmp0->Volume() * spinor_site_size * sizeof(quda::Complex), qudaMemcpyDeviceToHost);

  Qassert(fermion_type == 0);

  //*gsrcH = *ctmp0;

  //quda::blas::zero(*gresH);
  //*gsrcH = *gsrc;////quda base rotation
  //invertQuda((void*) (*gresH).data(), (void*) (*gsrcH).data(), &inv_param);
  //*gres  = *gresH;

  quda::blas::zero(*gres);
  invertQuda((void*) (*gres).data(), (void*) (*gsrc).data(), &inv_param);

  //*ctmp1 = *gresH;

  //qudaMemcpy((void*)(*gres).data(), (void*)ctmp1->data(), ctmp0->Volume() * spinor_site_size * sizeof(quda::Complex), qudaMemcpyHostToDevice);
  //for(Int si=0;si<(*gsrc).CompositeDim();si++)
  //{
  //  invertQuda((void*) (*gres).Component(si).data(), (void*) (*gsrc).Component(si).data(), &inv_param);
  //}
  //invertQuda(res, src, &inv_param);

  Qassert((void*) (*gres).data() != (void*) res);
  quda_cf_to_qlat_cf_test(res, (*gres), geo(), map_index);

  gettimeofday(&tm1, NULL);double time0 = tm1.tv_sec - tm0.tv_sec;time0 += (tm1.tv_usec - tm0.tv_usec)/1000000.0;

  inv_param.secs += 1e-25;
  if(quda_verbos >= -1)
  if(quda::comm_rank_global() == 0){
    qmessage("Done: %8d iter / %.6f secs = %.3f Gflops, Cost %.3f Gflops, %.6f secs \n",
          inv_param.iter, inv_param.secs, inv_param.gflops / inv_param.secs, inv_param.gflops, time0);
  }
  inv_time   = time0;
  inv_iter   = inv_param.iter;
  inv_gflops = inv_param.gflops / (inv_param.secs * qlat::get_num_node());
  timer.flops +=  inv_param.gflops * 1024 * 1024 * 1024 / qlat::get_num_node();
}

inline void quda_clover_inverter::rotate_quda_split(const bool swap_gauge){
  TIMER("rotate_quda_split");
  Qassert(split_setup);
  if(split_comm_buff == true){
    quda::push_communicator(quda::default_comm_key);
    split_comm_buff = false;
  }else{
    quda::push_communicator(split_key);
    split_comm_buff = true;
  }
  //if(split){
  //  quda::push_communicator(quda::default_comm_key);
  //}
  //else{
  //  quda::push_communicator(split_key);
  //}
  if(swap_gauge){swapGaugeSplit(true);}
  updateR();
  quda::comm_barrier();
  //
  //QudaInvertParam   split_param;
  //QudaInvertParam   split_param_copy;
  //quda::ColorSpinorParam cs_param_split;
  //quda::ColorSpinorParam devbuf_param;
  //std::vector<quda::ColorSpinorField> _collect_b;
  //std::vector<quda::ColorSpinorField> _collect_x;
  //std::vector<quda::ColorSpinorField>    dev_buf;
  //quda::CommKey split_key ;
  //quda::Dirac *split_dirac;
}

inline void quda_clover_inverter::set_quda_split(const std::vector<Int >& grid, const Int num_src){
  TIMER("set_quda_split");
  free_quda_split();// in case have previous buffers
  split_grid = grid;
  split_mrh  = num_src;
  //std::vector<std::string > grid_s = stringtolist(split_info);
  //std::vector<Int > grid;grid.resize(4);
  //for(int i=0;i<4;i++){grid[i] = stringtonum(grid_s[i]);}
  Qassert(grid.size() == 4);
  split_setup = true;
  //const Long V = geo.local_volume();
  //
  const Geometry& geo_ = geo();
  qlat::vector<Int> Nv,nv,mv;
  geo_to_nv(geo_, nv, Nv, mv);
  for(int i=0;i<4;i++){Qassert( mv[i] >= grid[i] and mv[i] % grid[i] == 0);}
  //
  split_param = inv_param;
  //qlat::vector<Int> Nv,nv,mv;
  //geo_to_nv(geo, nv, Nv, mv);
  //
  split_param.split_grid[0] = grid[0];
  split_param.split_grid[1] = grid[1];
  split_param.split_grid[2] = grid[2];
  split_param.split_grid[3] = grid[3];
  split_key = {grid[0], grid[1], grid[2], grid[3]};
  num_sub_partition = quda::product(split_key);
  //qmessage("Quda splited MRH num_src %5d sub %5d \n", int(num_src), int(num_sub_partition));
  Qassert(num_src >= num_sub_partition and num_src % num_sub_partition == 0);
  const Int dslash_mrh = num_src / num_sub_partition;
  split_param.num_src = dslash_mrh * num_sub_partition;
  split_param.num_src_per_sub_partition = dslash_mrh;
  //winv.inv_param.num_src = dslash_mrh;
  //
  //{
  //quda::ColorSpinorParam cs_param_ref(QL[0][0]);
  //cs_param_ref.create = QUDA_REFERENCE_FIELD_CREATE;
  //std::vector<quda::ColorSpinorField > s;
  //std::vector<quda::ColorSpinorField > r;
  //s.resize(dslash_mrh);
  //r.resize(dslash_mrh);
  //for(int bi=0;bi<in.bfac;bi++)
  //{
  //  TIMER_FLOPS("MRH dslash");
  //  //const double flops0 = quda::Tunable::flops_global();
  //  for(int j=0;j<num_sub_partition;j++){
  //    for(int m=0;m<dslash_mrh;m++){
  //      cs_param_ref.v = QL[0][j*dslash_mrh+m].data();
  //      s[m] = quda::ColorSpinorField(cs_param_ref);
  //      cs_param_ref.v = QL[1][j*dslash_mrh+m].data();
  //      r[m] = quda::ColorSpinorField(cs_param_ref);
  //    }
  //    winv.dirac->M(r, s);
  //    qacc_barrier(dummy);
  //  }
  //  //winv.dirac->M(QL[1], QL[0]);
  //  //const double flops1 = quda::Tunable::flops_global();
  //  //qmessage("quda flops %.3f per site \n", (flops1 - flops0)/(geo.local_volume()));
  //  //timer.flops += (flops1 - flops0) / qlat::get_num_node();
  //  timer.flops += geo.local_volume() * 1368 * dslash_mrh * num_sub_partition;
  //}
  //}
  UpdateSplitGauge(&split_param, false, false, split_param.split_grid); // split the gauges into split_key form
  quda::comm_barrier();
  //
  cs_param_split = quda::ColorSpinorParam(*gquda0);
  cs_param_split.setPrecision(split_param.cuda_prec, split_param.cuda_prec, true); // Native format
  cs_param_split.location = QUDA_CUDA_FIELD_LOCATION;                    // Device side
  cs_param_split.create = QUDA_REFERENCE_FIELD_CREATE;
  //
  // Expand the geometry for the collected fields
  split_fields_buffer.resize(5);
  split_fields.resize(split_fields_buffer.size());
  size_t V_tot = geo().local_volume() * 12;// srh clover size
  for (int d = 0; d < quda::CommKey::n_dim; d++) {
    cs_param_split.x[d] *= split_key[d]; 
    V_tot *= split_key[d];
  }
  for(unsigned int j=0;j<split_fields.size();j++){
    split_fields_buffer[j].resize(split_param.num_src_per_sub_partition * V_tot);
    split_fields[j].resize(split_param.num_src_per_sub_partition);
    for(Int i=0;i<split_param.num_src_per_sub_partition;i++){
      cs_param_split.v   = &split_fields_buffer[j][i* V_tot];
      split_fields[j][i] = quda::ColorSpinorField(cs_param_split);
    }
  }
  //
  // We will use these dev_buf fields to download (if needed) and convert
  // external fields into internal foramt
  devbuf_param = quda::ColorSpinorParam(*gquda0);
  devbuf_param.location = cs_param_split.location;                   // same location as collected (for copyOffset)
  devbuf_param.setPrecision(split_param.cuda_prec, split_param.cuda_prec, true); // Native format
  //std::vector<quda::ColorSpinorField> dev_buf(num_sub_partition, devbuf_param);
  dev_buf.resize(num_sub_partition);
  for(Int i=0;i<num_sub_partition;i++){
    dev_buf[i] = quda::ColorSpinorField(devbuf_param);
  }
  //
  //
  //for (int n = 0; n < split_param.num_src_per_sub_partition; n++) {
  //  // Download and change to Native Order and split
  //  for (int j = 0; j < num_sub_partition; j++) dev_buf[j].copy(QL[0][n * num_sub_partition + j]);
  //  quda::split_field(_collect_b[n], {dev_buf.begin(), dev_buf.end()}, split_key, QUDA_4D_PC);
  //}
  //
  rotate_quda_split(false);
  //quda::push_communicator(split_key);
  //updateR();
  //quda::comm_barrier();
  //
  split_param_copy = split_param;
  split_param_copy.input_location  = cs_param_split.location;
  split_param_copy.output_location = cs_param_split.location;
  //
  // Important: Don't use accessors for external formats any more
  // Since input fields are in Native order now
  split_param_copy.dirac_order = QUDA_INTERNAL_DIRAC_ORDER;
  //
  // We need to set the cpu_prec in the param_copy, because the op() passed in
  // to us will try to create wrappers to the pointers we pass in. They expect
  // the input spinors to be on the host, and will use param_copy.cpu_prec to set
  // the precision. We want to avoid the situation, where the internal prec and the
  // cpu_prec are somehow different.
  split_param_copy.cpu_prec = split_fields[0][0].Precision();
  //
  //quda::Dirac *dirac;
  // create dirac under splitted communicator
  if(split_dirac != NULL){delete split_dirac;split_dirac = NULL;}
  const bool pc_solve = false;
  quda::DiracParam diracParam;
  setDiracParam(diracParam, &split_param_copy, pc_solve);
  split_dirac     = quda::Dirac::create(diracParam);
  //{
  //  // Make a copy of the params we can mess with
  //  // Set solver input/output param location
  //  //
  //  // Do the solves
  //  //std::vector<void *> x_raw(param->num_src_per_sub_partition);
  //  //std::vector<void *> b_raw(param->num_src_per_sub_partition);
  //  //for (auto i = 0u; i < x_raw.size(); i++) x_raw[i] = _collect_x[i].data();
  //  //for (auto i = 0u; i < b_raw.size(); i++) b_raw[i] = _collect_b[i].data();
  //  //op(x_raw, b_raw, param_copy, args...);
  //  for(int bi=0;bi<in.bfac;bi++)
  //  {
  //    TIMER_FLOPS("SRH splitted dslash");
  //    //const double flops0 = quda::Tunable::flops_global();
  //    dirac->M(_collect_x, _collect_b);
  //    //const double flops1 = quda::Tunable::flops_global();
  //    //qmessage("quda flops %.3f per site \n", (flops1 - flops0)/(geo.local_volume()));
  //    qacc_barrier(dummy);
  //    //timer.flops += (flops1 - flops0) / qlat::get_num_node();
  //    timer.flops += geo.local_volume() * 1368 * dslash_mrh * num_sub_partition;
  //  }
  //}
  //
  // back to the default communicator, now join the param entries
  //auto split_rank = quda::comm_rank();
  // rotate back to original communicator
  rotate_quda_split(true);
  //quda::push_communicator(quda::default_comm_key);
  //updateR();
  //quda::comm_barrier();
  //quda::joinInvertParam(split_param, param_copy, split_key, split_rank);
  //
  // Join spinors: _h_x are aliases to host pointers in 'external order: QDP++, QDP-JIT, etc'
  //for (int n = 0; n < split_param.num_src_per_sub_partition; n++) {
  //  // join fields
  //  quda::join_field({dev_buf.begin(), dev_buf.end()}, _collect_x[n], split_key, QUDA_4D_PC);
  //  //
  //  // export to desired location and layout
  //  for (int j = 0; j < num_sub_partition; j++) QL[2][n * num_sub_partition + j].copy(dev_buf[j]);
  //}
  //
  /* 
    swap to original link and delete buffer
    always keep the splitted gauge
  */
  //swapGaugeSplit(true);
}

inline void quda_clover_inverter::copy_to_buf_split(std::vector<quda::ColorSpinorField>& src, const Int i){
  TIMER("copy_to_buf_split");
  rotate_from_quda_split();// swap to default comm to split fields
  Qassert(src.size() == split_param.num_src_per_sub_partition * num_sub_partition);
  Qassert(i < split_fields.size());
  for (int n = 0; n < split_param.num_src_per_sub_partition; n++) {
    // Download and change to Native Order and split
    for (int j = 0; j < num_sub_partition; j++) dev_buf[j].copy(src[n * num_sub_partition + j]);
    quda::split_field_device(split_fields[i][n], {dev_buf.begin(), dev_buf.end()}, split_key, QUDA_4D_PC);
  }
  rotate_to_quda_split();// swap to splitted comms
}

inline void quda_clover_inverter::copy_from_buf_split(std::vector<quda::ColorSpinorField>& res, const Int i){
  TIMER("copy_from_buf_split");
  rotate_from_quda_split();
  Qassert(res.size() == split_param.num_src_per_sub_partition * num_sub_partition);
  Qassert(i < split_fields.size());
  for (int n = 0; n < split_param.num_src_per_sub_partition; n++) {
    quda::join_field_device({dev_buf.begin(), dev_buf.end()}, split_fields[i][n], split_key, QUDA_4D_PC);
    //
    // export to desired location and layout
    for (int j = 0; j < num_sub_partition; j++) res[n * num_sub_partition + j].copy(dev_buf[j]);
  }
}

quda_clover_inverter::~quda_clover_inverter()
{
  free_mem();
}

// wilson dslash flops
inline double get_wilson_full_flops(const Geometry& geo)
{
  int nv[4];
  for(int i=0;i<4;i++){
    //nv[i] = geo.node_site[i] * geo.geon.size_node[i];
    nv[i] = geo.node_site[i];
  }
  const double fac = 1368;// need confirm from Quda
  double flops     = nv[0]  * nv[1] * nv[2] * nv[3] * fac;
  return flops;
}

template <typename Float> void constructCloverField(Float *res, double norm, double diag, Long V)
{

  Float c = 2.0 * norm / RAND_MAX;

  for (Int i = 0; i < V; i++) {
    for (Int j = 0; j < 72; j++) { res[i * 72 + j] = c * rand() - norm; }

    // impose clover symmetry on each chiral block
    for (Int ch = 0; ch < 2; ch++) {
      res[i * 72 + 3 + 36 * ch] = -res[i * 72 + 0 + 36 * ch];
      res[i * 72 + 4 + 36 * ch] = -res[i * 72 + 1 + 36 * ch];
      res[i * 72 + 5 + 36 * ch] = -res[i * 72 + 2 + 36 * ch];
      res[i * 72 + 30 + 36 * ch] = -res[i * 72 + 6 + 36 * ch];
      res[i * 72 + 31 + 36 * ch] = -res[i * 72 + 7 + 36 * ch];
      res[i * 72 + 32 + 36 * ch] = -res[i * 72 + 8 + 36 * ch];
      res[i * 72 + 33 + 36 * ch] = -res[i * 72 + 9 + 36 * ch];
      res[i * 72 + 34 + 36 * ch] = -res[i * 72 + 16 + 36 * ch];
      res[i * 72 + 35 + 36 * ch] = -res[i * 72 + 17 + 36 * ch];
    }

    for (Int j = 0; j < 6; j++) {
      res[i * 72 + j] += diag;
      res[i * 72 + j + 36] += diag;
    }
  }
}

void constructQudaCloverField(void *clover, double norm, double diag, QudaPrecision precision, Long V)
{
  if (precision == QUDA_DOUBLE_PRECISION)
    constructCloverField((double *)clover, norm, diag, V);
  else
    constructCloverField((float *)clover, norm, diag,  V);
}

template<typename Td>
void get_clover_prop(quda_clover_inverter& qinv, qlat::FermionField4dT<Td >& src, qlat::FermionField4dT<Td >& prop,
       const double kappa, const double err, const Int niter)
{
  qlat::ComplexT<Td >* srcP = (qlat::ComplexT<Td >*) qlat::get_data(src).data();

  const Geometry& geo = src.geo();
  if(!prop.initialized){prop.init(geo);} ////allocate mem for prop
  qlat::ComplexT<Td >* propP = (qlat::ComplexT<Td >*) qlat::get_data(prop).data();

  qinv.do_inv(propP, srcP, kappa, err, niter);
}

template<typename Td>
void get_clover_prop(quda_clover_inverter& qinv, qlat::Propagator4dT<Td >& src, qlat::Propagator4dT<Td >& prop,
      std::vector<qlat::FermionField4dT<Td > >& b0, std::vector<qlat::FermionField4dT<Td > >& b1,
       const double kappa, const double err, const Int niter)
{
  const Geometry& geo = src.geo();
  if(!prop.initialized){prop.init(geo);} ////allocate mem for prop

  if(b0.size() != 12){b0.resize(12);}
  if(b1.size() != 12){b1.resize(12);}
  for(Int i=0;i<12;i++){if(!b0[i].initialized){b0[i].init(geo);}}

  prop4d_to_Fermion(b0, src);
  for(Int i=0;i<12;i++){
    get_clover_prop(qinv, b0[i], b1[i], kappa, err, niter);
  }
  Fermion_to_prop4d(prop, b1);
}

template<typename Td>
void get_clover_prop(quda_clover_inverter& qinv, qlat::Propagator4dT<Td >& src, qlat::Propagator4dT<Td >& prop,
       const double kappa, const double err, const Int niter)
{
  std::vector<qlat::FermionField4dT<Td > > b0;
  std::vector<qlat::FermionField4dT<Td > > b1;
  get_clover_prop(qinv, src, prop, b0, b1, kappa, err, niter);
}

}  // namespace qlat

#endif


