#ifndef UTILS_QUDA_INVERTER_H
#define UTILS_QUDA_INVERTER_H

#pragma once

#include <quda.h>
////#include <gauge_force_quda.h>
#include <tune_quda.h>
#include <deflation.h>
#include <invert_quda.h>

// for gauge fixing
// #include <unitarization_links.h>
#include <gauge_tools.h>

#include <cstdlib>
#include "utils_float_type.h"
#include "quda_para.h"
#include "general_funs.h"
#include "utils_io_vec.h"
#include "utils_eo_copies.h"
#include "utils_Gaugefield_tools.h"
#include "utils_eigen_load.h"

//static quda::TimeProfile profileEigensolve("eigensolveQuda");
//static quda::TimeProfile profileInvertC("prefilesolveQudaC");

namespace qlat
{  //

//#define DIMCG        10
//#define DIMCG_DOUBLE 5

struct quda_inverter {
  Long V;
  Int X[4];
  QudaGaugeParam  gauge_param;
  ////void* quda_gf_default;
  qlat::vector<qlat::ComplexD > quda_gf_default;
  QudaInvertParam inv_param;
  //quda::SolverParam solverParam;
  //quda::Solver *solve_cg;
  //bool CG_reset;

  box<Geometry> geoB;
  qlat::vector<Long > map_index;
  qlat::FieldM<int8_t, 1> eo;//buffer for eo signs
  Int solve_mode ;
  /////QudaInvertParam df_param;

  QudaEigParam    eig_param;
  //int nvec;
  double inv_residue;
  double inv_time;
  Int    inv_iter;
  double inv_gflops;
  //int add_high;
  //int num_src;
  //int num_src_inv;
  Int quda_verbos;
  Int prec_type_check;

  ///std::vector<std::vector<quda::Complex > > evecs;

  //quda::EigenSolver *eig_solveK;
  //quda::EigenSolver *eig_solveF;
  quda::DiracMatrix* mat;
  quda::DiracMatrix* mat_E;
  quda::DiracMatrix* mat_pc;
  quda::DiracMatrix* mat_Mdag;
  quda::DiracMatrix* mat_MMdag;
  quda::DiracMatrix* mat_MdagM;

  //quda::DiracM* m_cg;
  //quda::DiracM* mSloppy;
  //quda::DiracM* mPre;
  //quda::DiracM* mEig;

  //void* df_preconditioner;

  double mass_mat;
  //double mass_eig;
  double mass_value;

  ////eigen related
  quda::Dirac *dirac;
  quda::Dirac *dirac_pc;
  quda::Dirac *dirac_pc_odd;

  //quda::Dirac *dirac_cg;
  //quda::Dirac *dSloppy;
  //quda::Dirac *dPre;
  //quda::Dirac *dEig;

  //std::vector<quda::Complex > evals_ZERO;
  //std::vector<double        > evals_ERR ;

  //std::vector<quda::ColorSpinorField > kSpace;
  //std::vector<quda::ColorSpinorField > fSpace;
  //////std::vector<quda::ColorSpinorField > kSpace_cpu;
  //std::vector<quda::Complex > evalsK;
  //std::vector<quda::Complex > evalsF;

  std::vector<quda::ColorSpinorField> out_multishift_gpu;
  std::vector<quda::ColorSpinorField> src_multishift_gpu;

  std::vector<quda::ColorSpinorField > ZSpace;
  std::vector<quda::Complex > evalsZ;
  ///std::vector<quda::ColorSpinorField *> evecs_;
  Int spinor_site_size;
  Int max_src_qinv;

  //quda::ColorSpinorParam cs_cpu;
  quda::ColorSpinorParam cs_gpu;
  quda::ColorSpinorParam cs_gpuH;

  //quda::ColorSpinorField *csrc, *cres;
  //quda::ColorSpinorField *cpu_src, *cpu_res;
  quda::ColorSpinorField *gsrc, *gres;
  quda::ColorSpinorField *gsrcH, *gresH;
  ///quda::ColorSpinorField *gsrcD, *gresD;


  // should be used only for buffers of MRH inversions
  qlat::vector_gpu<ComplexT<double> > buf_inv;

  // used for eigensystem construction of double / single
  quda::ColorSpinorParam cs_gpuD;
  quda::ColorSpinorParam cs_gpuF;
  quda::ColorSpinorField *gtmp1D, *gtmp2D, *gtmp3D, *gtmp4D, *gtmp5D;
  quda::ColorSpinorField *gtmp1F, *gtmp2F;
  quda::ColorSpinorField *gadd; ////for low mode addition
  //quda::ColorSpinorField *gtmp_invG0;

  //quda::ColorSpinorField *ctmp0, *ctmp1, *ctmp2;
  quda::ColorSpinorField *gtmp0, *gtmp1, *gtmp2;
  //quda::ColorSpinorField *gsrcF, *gresF; ////format to match gpu solver
  ///bool singleE;

  bool io_rotate_bfac;

  ////0 for wilson, clover, 1 for stagger
  Int fermion_type;
  bool gauge_with_phase;
  //int use_eigen_pc;

  /*  TODO add redeflate for QUDA
      default precision float, double
      eigen related pointers and prec
  */

  std::vector<void* > eigenL;
  std::vector<Int   > eigen_precL;
  bool eigen_with_nvec;
  std::vector<Int > eigen_with_nvecL;
  //void* eigen0;//float  eigen
  //void* eigen1;//double eigen for high prec

  Int check_residue;
  bool clover_alloc;

  //int Nshift;
  //int shift_masses;
  //std::vector<quda::ColorSpinorField> out_multishift_gpu;

  std::vector<std::vector<quda::Complex > > quda_clover;
  std::vector<quda::Complex > quda_clover_inv;

  quda_inverter(const Geometry& geo_, QudaTboundary t_boundary);

  inline void free_mem();
  inline void setup_link(qlat::ComplexD* quda_gf, const Int apply_stag_phase = 0, const Int force_phase = 0);
  inline void gauge_fix(qlat::ComplexD* quda_gf, Int gf_gauge_dir = 3, Int gf_maxiter = 10000, double gf_tolerance = 1e-6, Int fix_type = 0);

  //default e0 single, e1 double
  template<typename T1, typename T2>
  inline void setup_eigen(eigen_cs<T1>* e0, eigen_cs<T2>* e1);

  // even-even eigensystem, could be used with dslash to get full prop
  // mode == 0, deflate, mode == 1 get low prop
  // buf_prec = 0 : input doubleC, buf_prec = 1 : input singleC
  inline void deflate_Ty(qlat::vector<void* >& Pres, qlat::vector<void* >& Psrc, double mass, Int buf_prec, Int mode = 0, Int clear = 1, const bool data_from_quda = true);

  // src need to be unchanged if res != src
  // res could be the same as src
  // use gtmp0
  template<typename Ty, typename Tk>
  inline void get_low_prop(qlat::vector<Ty* >& res, qlat::vector<Ty* >& src, qlat::vector<Tk* >& buf, Int qlat_format = 1);

  template<typename Ty>
  inline void get_inversion_bufs(qlat::vector<Ty* >& res, const Int nvecs, const Int halfV = 0){
    LInt Vd = geoB().local_volume() * 3 * sizeof(Ty) / sizeof(qlat::ComplexT<double>);
    if(halfV == 1){Vd = Vd / 2;}
    Qassert(Vd / (geoB().local_volume() * 3) <= 64);//at most 64 vectors in buf
    const LInt totalD = Vd * nvecs;
    buf_inv.resizeL(totalD);
    res.resize(nvecs);
    for(Int iv=0;iv<nvecs;iv++)
    {
      res[iv] = (Ty*) &buf_inv[iv*Vd];
    }
  }

  // low mode from gsrc to gres
  // use gtmp1
  inline void prepare_low_prop();

  //inline void setup_clover(const double clover_csw);
  //initialize inv parameters, allocate initial cpu and gpu buffers
  inline void setup_stagger_inv();

  //inline void update_eigen_mass(const double mass, bool force = false);

  //inline void setup_inc_eigencg(const Int n_ev, const Int n_kr, const Int n_conv, const Int df_grid, const double tol, const double inc_tol, const double tol_restart, const Int restart_n, const Int pipeline, const Int inv_type = 1);

  //template<typename Ty>
  //inline void do_inv(Ty* res, Ty* src, const double mass, const double err = 1e-10, const Int niter = 10000 , const Int prec_type = 0);
  //inline void deflate(quda::ColorSpinorField &sol, const quda::ColorSpinorField& src,
  //    std::vector<quda::ColorSpinorField> &evecs, const std::vector<quda::Complex > &evals, bool accumulate = false);

  //inline void deflate(std::vector<quda::ColorSpinorField* > &sol, const std::vector<quda::ColorSpinorField* > &src,
  //    std::vector<quda::ColorSpinorField> &evecs, const std::vector<quda::Complex > &evals, bool accumulate = false);

  inline void callMultiSrcQuda(qlat::vector<void* >& res, qlat::vector<void* >& src, Int max_src = -1);

  void print_plaq();

  //inline void setup_inv_mass(const double mass);
  // setup masses and Dirac
  inline void setup_mat_mass(const double mass, const bool force_do = false);
  inline void clear_mat();

  inline void free_csfield(const Int mode = 0);

  // ctmp0, ctmp1, ctmp2, 
  inline void alloc_csfield_initial();

  // set up ColorSpinorParam from inv
  // gsrc, gres, gtmp0, gtmp1, gtmp2, cs_gpuH, gsrcH, gresH
  inline void alloc_csfield_gpu();
  //need check residual CG

  //inline void save_evecs(const char* filename, const bool read = false, const Int ndouble = -1, const Int split_save = 1 );
  //inline void save_evecsF(const char* filename, const bool read = false);
  //inline void check_residualF();

  // ghost function to set precision for inv and gauge, not loaded at all, need setup_inv_param_prec to force load gauge again or use setup_link
  inline void setup_inv_param_prec_type(Int prec_type = 0);
  // automatic change presion based on current gauge, qinv buffers, update dslash
  inline void setup_inv_param_prec(Int prec_type = 0, bool force_reload = false);
  inline void setup_gauge_param(QudaTboundary t_boundary);

  inline void random_src(const Int seed);

  //inline void reconstruct_full(const double mass = -1);
  //inline void save_prop(const void* srcP, const char* filename);

  // create solverParam needed by invertQuda_COPY_single
  //inline void setup_CG();
  //inline void clear_CG(){
  //  return ;
  //  //if(solve_cg != NULL){delete solve_cg;solve_cg = NULL;}
  //}
  //inline void invertQuda_COPY(quda::ColorSpinorField& x, quda::ColorSpinorField& b, Int solve_mode_ = 0);
  // the COPY could use single precison as input
  //inline void invertQuda_COPY(quda::ColorSpinorField& res, quda::ColorSpinorField& src, Int solve_mode_ = 0);
  inline void invertQuda_COPY_single(quda::ColorSpinorField& res, quda::ColorSpinorField& src);

  inline void get_param_ref(quda::ColorSpinorParam& param, const bool is_double ){
    param   = cs_gpu;
    param.is_composite  = false;
    param.is_component  = false;
    if(is_double){
      param.setPrecision(QUDA_DOUBLE_PRECISION, QUDA_DOUBLE_PRECISION, true);
    }else{
      param.setPrecision(QUDA_SINGLE_PRECISION, QUDA_SINGLE_PRECISION, true);
    }
    param.create = QUDA_REFERENCE_FIELD_CREATE;
  }

  template <class Ty>
  inline void qlat_cf_to_quda_cf_P(void* res, Ty* src, qlat::vector<Long >& map, quda::ColorSpinorParam& param){
    param.v = (void*) res;
    quda::ColorSpinorField* Qvec = NULL;
    Qvec = new quda::ColorSpinorField(param);
    qlat_cf_to_quda_cf(*Qvec, src, geoB(), map);
    delete Qvec;
  }

  template <class Ty>
  inline void quda_cf_to_qlat_cf_cf_P(Ty* res, void* src, qlat::vector<Long >& map, quda::ColorSpinorParam& param){
    param.v = (void*) src;
    quda::ColorSpinorField* Qvec = NULL;
    Qvec = new quda::ColorSpinorField(param);
    quda_cf_to_qlat_cf(res, *Qvec, geoB(), map);
    delete Qvec;
  }

  /*
    transform from quda double prec order to single prec order 
    prec = 0 double, prec = 1 single
  */
  inline void quda_to_single_order(vector<void* >& src, const Int prec, const Int dir = 0){
    TIMER("quda_to_single_order");
    if(prec == 0){
      return ;
    }
    const size_t Ndata = size_t(V) * 3 / 2;
    //
    qlat::vector_gpu<int8_t >& quda_buf = get_vector_gpu_plan<int8_t >(0, "quda_field_copy_buffers", 1);
    quda_buf.resizeL(Ndata * sizeof(ComplexT<double>) / sizeof(int8_t));
    //
    quda::ColorSpinorParam pd = cs_gpuH;
    quda::ColorSpinorParam pf = cs_gpuH;
    pd.setPrecision(QUDA_DOUBLE_PRECISION, QUDA_DOUBLE_PRECISION, true);
    pf.setPrecision(QUDA_SINGLE_PRECISION, QUDA_SINGLE_PRECISION, true);
    //
    pd.create = QUDA_REFERENCE_FIELD_CREATE;
    pf.create = QUDA_REFERENCE_FIELD_CREATE;
    //
    quda::ColorSpinorField* Qd = NULL;
    quda::ColorSpinorField* Qf = NULL;
    pd.v = (void*) quda_buf.data();
    Qd = new quda::ColorSpinorField(pd);
    //
    for(Long iv=0;iv<src.size();iv++)
    {
      //display_mem_type(src[iv]);
      pf.v = (void*) src[iv];
      Qf = new quda::ColorSpinorField(pf);
      //
      if(dir == 0){
        cpy_GPU((ComplexT<double>* ) quda_buf.data(), (ComplexT<float >* ) src[iv], Ndata);
        *Qf = *Qd;
      }
      if(dir == 1){
        *Qd = *Qf;
        cpy_GPU((ComplexT<float >* ) src[iv], (ComplexT<double>* ) quda_buf.data(), Ndata);
      }
      delete Qf;
    }
    delete Qd;
  }

  /*
    transform from quda single prec order to double prec order 
  */
  inline void quda_to_double_order(vector<void* >& src, const Int prec){
    quda_to_single_order(src, prec, 1);
  }

  ~quda_inverter();
};

quda_inverter::quda_inverter(const Geometry& geo_, QudaTboundary t_boundary)
{
  TIMER("quda_inverter_constuctor");
  /////set up gauge parameters
  geoB.set(geo_);
  V = geoB().local_volume();
  //Qassert(num_src_ > 0);
  //num_src = num_src_;
  //num_src_inv = num_src;
  prec_type_check = -2;

  for (Int mu = 0; mu < 4; mu++) {
    X[mu] = geoB().node_site[mu];
    Qassert(geoB().node_site[mu] % 2 == 0);// needed for eo inverter
  }
  ////===Start of gauge_param
  gauge_param = newQudaGaugeParam();
  ////quda_gf_default = NULL;
  inv_param   = newQudaInvertParam();

  setup_gauge_param(t_boundary);

  //add_high = 1;
  //solve_mode = 0;

  clover_alloc = false;
  //csrc = NULL; cres = NULL;
  gsrc = NULL; gres = NULL;
  gsrcH= NULL; gresH= NULL;
  //gsrcD= NULL; gresD= NULL;
  gtmp1D = NULL;gtmp2D = NULL;gtmp3D = NULL;
  gtmp4D = NULL;gtmp5D = NULL;
  gtmp1F = NULL;gtmp2F = NULL;
  gadd = NULL;
  //gtmp_invG0 = NULL;

  //ctmp0 = NULL; ctmp1 = NULL; ctmp2 = NULL;
  gtmp0 = NULL; gtmp1 = NULL; gtmp2 = NULL;

  //eig_solveK = NULL;
  //eig_solveF = NULL;
  mat = NULL;
  mat_pc = NULL;
  mat_Mdag  = NULL;
  mat_MMdag = NULL;
  mat_E = NULL;
  mat_MdagM = NULL;

  //m_cg = NULL;
  //mSloppy = NULL;
  //mPre = NULL;
  //mEig = NULL;
  //CG_reset = true;

  //solve_cg = NULL;
  dirac = NULL;
  //dirac_cg = NULL;
  dirac_pc     = NULL;
  dirac_pc_odd = NULL;
  //dSloppy = NULL;
  //dPre    = NULL;
  //dEig    = NULL;

  eigen_with_nvec = false;

  //df_preconditioner = NULL;
  //nvec = 0;
  //use_eigen_pc = 0;
  check_residue = 0;

  spinor_site_size =  0;
  max_src_qinv     =  qlat::get_env_long_default(std::string("qlat_quda_stag_mrh"), 9);
  gauge_with_phase = false;
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

  //default stagger
  setup_stagger_inv();
  ////default mass zero
  //setup_mat_mass(0);
  ////default precision
  //setup_inv_param_prec(0);
}

inline void quda_inverter::setup_gauge_param(QudaTboundary t_boundary)
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

  // improve performance a little
  gauge_param.reconstruct                   = QUDA_RECONSTRUCT_13;  ////gauge_param.reconstruct = link_recon;
  gauge_param.reconstruct_sloppy            = QUDA_RECONSTRUCT_13;
  gauge_param.location = QUDA_CPU_FIELD_LOCATION;

  //gauge_param.reconstruct                   = QUDA_RECONSTRUCT_NO;  ////gauge_param.reconstruct = link_recon;
  //gauge_param.reconstruct_sloppy            = QUDA_RECONSTRUCT_NO;

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

template<typename T1, typename T2>
inline void quda_inverter::setup_eigen(eigen_cs<T1>* e0, eigen_cs<T2>* e1)
{
  eigenL.resize(2);
  eigen_precL.resize(2);
  eigen_with_nvecL.resize(2);
  //double precision eigen first
  Qassert(sizeof(T1) >= sizeof(T2));

  eigenL[0] = e0;
  eigenL[1] = e1;
  eigen_precL[0] = 1 - get_data_type_is_double<T1>();
  eigen_precL[1] = 1 - get_data_type_is_double<T2>();

  for(Int i=0;i<2;i++){
    Int total_n = 0;
    void* Ea    = eigenL[i];
    Int   prec  = eigen_precL[i];
    if(Ea != NULL){
      if(prec == 0){
        eigen_cs<qlat::ComplexT<double  >>* eigen = (eigen_cs<qlat::ComplexT<double  >>*) Ea;
        total_n += eigen->get_nvec();
      }
      if(prec == 1){
        eigen_cs<qlat::ComplexT<float   >>* eigen = (eigen_cs<qlat::ComplexT<float   >>*) Ea;
        total_n += eigen->get_nvec();
      }
    }
    eigen_with_nvecL[i] = total_n;
    if(total_n > 0){
      eigen_with_nvec = true;
    }
  }
}

inline void quda_inverter::setup_link(qlat::ComplexD* quda_gf, const Int apply_stag_phase, const Int force_phase)
{
  TIMER("setup_link");
  /////load gauge to quda GPU default position
  freeGaugeQuda();
  //{
  //qlat::ComplexD res =  vec_norm2(quda_gf, quda_gf, geo().local_volume() * 4 * 9, QMGPU, 64);
  //qmessage("gauge norm %.8e %.8e \n",  res.real(), res.imag());
  //}
  if(apply_stag_phase == 1 and (gauge_with_phase == false or force_phase == 1)){
    applyGaugeFieldScaling_long((qlat::ComplexD*) quda_gf, V/2, &gauge_param, QUDA_STAGGERED_DSLASH);
    loadGaugeQuda((void *) quda_gf, &gauge_param);
    gauge_with_phase = true;
  }
  else{
    loadGaugeQuda((void *) quda_gf, &gauge_param);
    //print_plaq();
  }
  ////quda_gf_default = (void *) quda_gf; //required to reload gauge with prec
  if(quda_gf != quda_gf_default.data()){
    quda_gf_default.resize(V *3*3*4 );
    cpy_data_thread(&quda_gf_default[0], quda_gf, V*3*3*4, false);
  }
  if(gauge_with_phase == false){errorQuda("Quda stagger need link phases! \n");}
  //print_plaq();
}

inline void quda_inverter::print_plaq()
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
inline void quda_inverter::free_csfield(const Int mode)
{
  TIMER("free_csfield");
  Qassert(mode == 0 or mode == 2 or mode == 3);
  //if(mode == 0 or mode == 1){
  ////if(csrc  != NULL){delete csrc;csrc=NULL;}if(cres != NULL){delete cres;cres=NULL;}
  //if(ctmp0 != NULL){delete ctmp0;ctmp0=NULL;}
  //if(ctmp1 != NULL){delete ctmp1;ctmp1=NULL;}
  //if(ctmp2 != NULL){delete ctmp2;ctmp2=NULL;}
  //}

  if(mode == 0 or mode == 2){
  //if(gsrcD != NULL){delete gsrcD;gsrcD=NULL;}if(gresD!= NULL){delete gresD;gresD=NULL;}
  if(gsrc  != NULL){delete gsrc ;gsrc =NULL;}if(gres != NULL){delete gres ;gres =NULL;}
  if(gsrcH != NULL){delete gsrcH;gsrcH=NULL;}if(gresH!= NULL){delete gresH;gresH=NULL;}

  if(gtmp0 != NULL){delete gtmp0;gtmp0=NULL;}
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

inline void quda_inverter::alloc_csfield_gpu()
{
  TIMER("alloc_csfield_gpu");
  free_csfield(2);
  //quda::ColorSpinorParam cs_tmp(cs_gpu);
  cs_gpu.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);
  cs_gpu.is_composite  = false;
  cs_gpu.is_component  = false;
  cs_gpu.composite_dim = 1;

  //gsrc  = quda::ColorSpinorField::Create(cs_gpu);
  gsrc  = quda::ColorSpinorField::Create(cs_gpu);
  gres  = quda::ColorSpinorField::Create(cs_gpu);
  gtmp0 = quda::ColorSpinorField::Create(cs_gpu);
  gtmp1 = quda::ColorSpinorField::Create(cs_gpu);
  gtmp2 = quda::ColorSpinorField::Create(cs_gpu);

  cs_gpuH.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);
  //cs_gpuH.setPrecision(QUDA_DOUBLE_PRECISION, inv_param.cuda_prec, true);
  cs_gpuH.is_composite  = false;
  cs_gpuH.is_component  = false;
  cs_gpuH.composite_dim = 1;
  //if (cs_gpuH.pc_type != QUDA_5D_PC && cs_gpuH.pc_type != QUDA_4D_PC){
  //  errorQuda("self Unexpected pc_type %d", cs_gpuH.pc_type);  
  //}

  gsrcH = quda::ColorSpinorField::Create(cs_gpuH); 
  gresH = quda::ColorSpinorField::Create(cs_gpuH);

  //cs_tmp.setPrecision(QUDA_DOUBLE_PRECISION, QUDA_DOUBLE_PRECISION, true);
  //gsrcD = quda::ColorSpinorField::Create(cs_tmp);
  //gresD = quda::ColorSpinorField::Create(cs_tmp);
}

inline void quda_inverter::alloc_csfield_initial()
{
  TIMER("alloc_csfield_initial");
  //free_csfield(1);
  //bool pc_solution = (inv_param.solution_type == QUDA_MATPC_SOLUTION) ||
  //  (inv_param.solution_type == QUDA_MATPCDAG_MATPC_SOLUTION);
  bool pc_solution = false;
  void* temV = NULL;
  quda::GaugeField gpuGauge(gauge_param);
  quda::ColorSpinorParam gpuParam0(temV, inv_param, gpuGauge.X(), pc_solution, QUDA_CUDA_FIELD_LOCATION);
  //cs_gpu.pc_type = QUDA_4D_PC;
  cs_gpu = quda::ColorSpinorParam(gpuParam0);
  ////cs_gpu.setPrecision(QUDA_DOUBLE_PRECISION); ////double for all cpu cs field
  cs_gpu.create = QUDA_ZERO_FIELD_CREATE;
  cs_gpu.location = QUDA_CUDA_FIELD_LOCATION;
  cs_gpu.is_composite  = false;
  cs_gpu.is_component  = false;
  cs_gpu.composite_dim = 1;

  //cs_gpu = quda::ColorSpinorParam(cs_cpu);cs_gpu.location = QUDA_CUDA_FIELD_LOCATION;
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
  if(cs_gpu.nSpin != 1) cs_gpu.gammaBasis = QUDA_UKQCD_GAMMA_BASIS;

  //csrc  = quda::ColorSpinorField::Create(cs_cpu);
  //cres  = quda::ColorSpinorField::Create(cs_cpu);
  //ctmp0 = quda::ColorSpinorField::Create(cs_cpu);
  //ctmp1 = quda::ColorSpinorField::Create(cs_cpu);
  //ctmp2 = quda::ColorSpinorField::Create(cs_cpu);

  //bool single_file=true;
  //if(fermion_type == 0){pc_solution = false;}
  Qassert(fermion_type == 1);
  //if(fermion_type == 1){}
  pc_solution = true ;
  quda::ColorSpinorParam gpuParam1(temV, inv_param, gpuGauge.X(), pc_solution, QUDA_CUDA_FIELD_LOCATION);
  cs_gpuH = quda::ColorSpinorParam(gpuParam1);
  cs_gpuH.create = QUDA_ZERO_FIELD_CREATE;
  cs_gpuH.is_composite  = false;
  cs_gpuH.is_component  = false;
  cs_gpuH.composite_dim = 1;

  cs_gpuF= quda::ColorSpinorParam(gpuParam1);
  cs_gpuF.create = QUDA_ZERO_FIELD_CREATE;
  cs_gpuF.setPrecision(QUDA_SINGLE_PRECISION, QUDA_SINGLE_PRECISION, true);
  cs_gpuF.is_composite  = false;
  cs_gpuF.is_component  = false;

  ////cs_gpuF.setPrecision(QUDA_DOUBLE_PRECISION, QUDA_DOUBLE_PRECISION, true);

  cs_gpuD= quda::ColorSpinorParam(gpuParam1);
  cs_gpuD.is_composite  = false;
  cs_gpuD.is_component  = false;
  cs_gpuD.create = QUDA_ZERO_FIELD_CREATE;
  cs_gpuD.setPrecision(QUDA_DOUBLE_PRECISION, QUDA_DOUBLE_PRECISION, true);

  //quda::ColorSpinorParam cs_gpu_tem = cs_gpu;
  //cs_gpu_tem.setPrecision(QUDA_SINGLE_PRECISION, QUDA_SINGLE_PRECISION, true);
  //gsrcF = quda::ColorSpinorField::Create(cs_gpu_tem);
  //////=====START construct Staggered color spin parameters
  //cs_cpu.nColor = 3;
  //cs_cpu.nSpin = 1;
  //cs_cpu.nDim = 5;
  //for (Int d = 0; d < 4; d++) cs_cpu.x[d] = gauge_param.X[d];
  //bool pc = isPCSolution(inv_param.solution_type);
  //if (pc) cs_cpu.x[0] /= 2;
  //cs_cpu.x[4] = 1;
  //cs_cpu.pc_type = QUDA_4D_PC;
  //cs_cpu.siteSubset = pc ? QUDA_PARITY_SITE_SUBSET : QUDA_FULL_SITE_SUBSET;
  //// Lattice vector data properties
  //cs_cpu.pad = 0;
  //cs_cpu.siteOrder = QUDA_EVEN_ODD_SITE_ORDER;
  //cs_cpu.fieldOrder = QUDA_SPACE_SPIN_COLOR_FIELD_ORDER;
  //cs_cpu.gammaBasis = inv_param.gamma_basis;
  //cs_cpu.create = QUDA_ZERO_FIELD_CREATE;
  //cs_cpu.location = QUDA_CPU_FIELD_LOCATION;
  //////=====END construct Staggered color spin parameters


}

/////double prec prop save
//inline void quda_inverter::save_prop(const void* srcP, const char* filename)
//{
//  TIMER("quda save_prop");
//  ////qlat::Coordinate total_site;
//  ////qlat::Coordinate node_site = qlat::get_size_node();
//  ////for(Int d=0;d<4;d++){total_site[d] = X[d] * node_site[d];}
//  ////Geometry geo;geo.init(total_site, 1);
//
//  const Int n0 = 1;
//  std::vector<qlat::FieldM<qlat::ComplexD , 3> > prop;prop.resize(n0);
//  for(Int n = 0; n < prop.size(); n++){prop[n].init(geo);}
//
//  qlat::ComplexD* src = (qlat::ComplexD*) srcP;
//  Long Nvol = geo.local_volume() * n0 ;
//
//  for(Int n=0;n<n0;n++){
//    quda_cf_to_qlat_cf(prop[n], &src[n*Nvol]);
//  }   
//
//  std::string VECS_TYPE("STAGGERED_Prop");
//  std::string INFO_LIST = ssprintf("mass %.8f", inv_param.mass);
//
//  bool single_file = false;
//  bool read = false;
//  qlat::load_qlat_noisesT(filename, prop, read, single_file, VECS_TYPE, std::string("NONE"), 0, n0, false);
//}

//inline void quda_inverter::save_evecsF(const char* filename, const bool read)
//{
//  TIMERB("save_evecsF");
//  if(nvec <= 0 ){return ;}
//  std::string filename0, filename1;
//  filename0 = ssprintf("%s.full", filename);
//  double mre = 0.0;
//
//  quda::ColorSpinorField& c0 = (*ctmp0).Component(0);
//  //quda::ColorSpinorParam cs0 = quda::ColorSpinorParam(cs_gpu);
//  //cs0.is_composite  = false;cs0.is_component  = false;
//  //quda::ColorSpinorParam cs0 = quda::ColorSpinorParam(cs_gpu);
//  //cs0.is_composite  = false;cs0.is_component  = false;
//  ////quda::ColorSpinorField& c1 = (*ctmp1).Component(0);
//  ///quda::ColorSpinorField c0 = quda::ColorSpinorField(cs0);
//
//  Int n0 = 0;
//  if(read == true)
//  {
//    inputpara in0;////inputpara in1;
//    in0.load_para(filename0.c_str(), false);
//    n0 = in0.nvec;
//
//    ///read eigen value masses
//    std::vector<std::string > mL0 = stringtolist(in0.INFO_LIST);
//    double mi0 = stringtodouble(mL0[1]);
//    mre= mi0;
//  }else{n0 = ZSpace.size();}
//
//  std::vector<qlat::FieldM<qlat::ComplexD , 3> > eigD;eigD.resize(n0);
//
//  qlat::Coordinate total_site;
//  qlat::Coordinate node_site = qlat::get_size_node();
//  for(Int d=0;d<4;d++){total_site[d] = X[d] * node_site[d];}
//  qlat::Geometry geo;geo.init(total_site);
//
//  for(unsigned int n = 0; n < eigD.size(); n++){eigD[n].init(geo);}
//
//  if(read == true)
//  {
//    //for (Int i = 0; i < ZSpace.size(); i++){delete ZSpace[i];}ZSpace.resize(0);
//    //for(Int n=0;n<n0;n++)ZSpace.push_back(quda::ColorSpinorField::Create(cs_gpu));
//    ZSpace.resize(n0);
//    for(Int n=0;n<n0;n++)ZSpace[n] = quda::ColorSpinorField(cs_gpu);
//  }
//
//  if(read == false)
//  for(Int n=0;n<n0;n++){
//    (c0) = (ZSpace[n]);
//    quda_cf_to_qlat_cf(eigD[n], (qlat::ComplexD*) c0.data());
//  }   
//
//
//  std::string VECS_TYPE("STAGGERED_Eigensystem");
//  std::string INFO_LIST = ssprintf("mass %.8f", mass_value);
//
//  bool single_file = false;
//  qlat::load_qlat_noisesT(filename0.c_str(), eigD, read, single_file, VECS_TYPE, INFO_LIST, 0, n0, io_rotate_bfac);
//
//  if(read == true )
//  for(Int n=0;n<n0;n++){
//    qlat_cf_to_quda_cf((qlat::ComplexD*) c0.data(),  eigD[n]);
//    (ZSpace[n]) = (c0);
//  }   
//
//  filename1 = ssprintf("%s.evals", filename0.c_str());
//  std::vector<double > values, errors;
//  values.resize(2*n0); errors.resize(n0);
//  if(read == false){
//    Qassert(evalsZ.size() != evals_ZERO.size());
//    Qassert(evals_ZERO.size() != 0 and evals_ERR.size() != 0);
//    for(Int n=0;n<nvec;n++){
//      values[n*2+0] = evals_ZERO[n].real();
//      values[n*2+1] = evals_ZERO[n].imag();
//      errors[n] = evals_ERR[n]; 
//    }
//    save_txt_eigenvalues(values, errors, filename1.c_str(), "Staggered Fermions");
//  }
//
//  if(read == true){
//    load_txt_eigenvalues(values, errors, filename1.c_str());Qassert(n0 <= values.size());
//    evals_ZERO.resize(n0);evals_ERR.resize(n0);
//    /////for(Int n=  0;n<n0;n++){evals_ZERO[n] = quda::Complex(values[n*2+0], values[n*2+1]);}    
//    for(Int n=  0;n<n0;n++){evals_ZERO[n] = quda::Complex(values[n*2+0] - 4.0*mre*mre, values[n*2+1]);}    
//    for(Int n=  0;n<n0;n++){evals_ERR[ n] = errors[n];}    
//
//    evalsZ.resize(n0);evalsK.resize(n0, 0.0);evalsF.resize(0, 0.0);
//    update_eigen_mass(0.0, true);
//  }
//}

//inline void quda_inverter::save_evecs(const char* filename, const bool read, const Int ndouble, const Int split_save )
//{
//  TIMERB("save_evecs");
//  if(nvec <= 0){return ;}
//  ////Qassert(nvec == kSpace.size());Qassert(nvec == evals.size());
//  //bool singleE = false;
//  //if(compute == true ){singleE = false;}
//  //if(compute == false){singleE = true;}
//
//  std::string filename0, filename1;
//  filename0 = ssprintf("%s", filename);
//  filename1 = ssprintf("%s.single", filename);
//  double mre = 0.0;
//  quda::ColorSpinorField& c0 = (*ctmp0).Component(0);
//
//  Int nsave = nvec;
//  if(fermion_type == 1){
//    nsave = nvec/2;
//    if(nvec%2 != 0){qmessage("Even Odd savings, nvec  2 not 0, nvec %d \n", nvec);Qassert(false);}
//  }
//
//  //bool single_file=true;
//  bool single_file=false;
//  std::string VECS_TYPE("STAGGERED_Eigensystem");
//  //char infoL[500];
//  //ssprintf(infoL,"mass %.8f", mass_value);
//  //std::string INFO_LIST(infoL);
//  std::string INFO_LIST = ssprintf("mass %.8f", mass_value);
//
//  Int n0 = nsave;int n1 = 0;
//  if(read == false){
//    if(split_save == 0){n0=nsave;n1=0;}
//    if(split_save == 1){
//    if(ndouble < 0 or ndouble >= nvec){
//      n0 = nsave;n1=0;}
//    else{
//      n0 = ndouble;if(fermion_type == 1){n0 = ndouble/2;}
//      n1 = nsave - n0;}
//    }
//  }else{
//    ////ndouble no meaning here
//    //inputpara in0;inputpara in1;
//    //in0.load_para(filename0, false);
//    //n0 = in0.nvec;
//    //if(nsave <= n0){n0 = nsave;n1 = 0;}
//    //else{
//    //  in1.load_para(filename1, false);n1 = in1.nvec;
//    //  Qassert(nsave <= n0+n1);
//    //  n1 = nsave - n0  ;
//    //}
//
//    inputpara in0;inputpara in1;
//    if(split_save == 0){
//      in0.load_para(filename0.c_str(), false);
//      n0 = in0.nvec;
//      Qassert(n0 >= nsave);
//      if(ndouble == -1 or ndouble >= nvec){
//        n0 = nsave;n1=0;
//      }
//      else{
//        n0 = ndouble;if(fermion_type == 1){n0 = ndouble/2;}
//        n1 = nsave - n0;
//      }
//
//      ///read eigen value masses
//      std::vector<std::string > mL0 = stringtolist(in0.INFO_LIST);
//      double mi0 = stringtodouble(mL0[1]);
//      mre = mi0;
//    }
//
//    if(split_save == 1){
//      ////ignore ndouble
//      in0.load_para(filename0.c_str(), false);
//      in1.load_para(filename1.c_str(), false);
//      Qassert(in0.nvec + in1.nvec >= nsave);
//      if(nsave <= in0.nvec){n0 = nsave;n1 = 0;}
//      if(nsave  > in0.nvec){n0 = in0.nvec;n1 = nsave -n0 ;}
//
//      ///read eigen value masses
//      std::vector<std::string > mL0 = stringtolist(in0.INFO_LIST);
//      std::vector<std::string > mL1 = stringtolist(in1.INFO_LIST);
//      double mi0 = stringtodouble(mL0[1]);
//      double mi1 = stringtodouble(mL1[1]);
//      Qassert(mi0 == mi1);
//      mre = mi0;
//    }
//
//  }
//
//  Int ns0 = n0;int ns1 = n1;
//  if(nsave == nvec/2){ns0 = n0*2;ns1 = n1*2;}
//
//  qlat::Coordinate total_site;
//  qlat::Coordinate node_site = qlat::get_size_node();
//  for(Int d=0;d<4;d++){total_site[d] = X[d] * node_site[d];}
//  qlat::Geometry geo;geo.init(total_site);
//
//  //std::vector<qlat::FieldM<qlat::ComplexF, 3> > eig;eig.resize(nsave);
//  Long Nvol = geo.local_volume() * spinor_site_size;
//  std::vector<qlat::FieldM<qlat::ComplexD , 3> > eigD;eigD.resize(n0);
//  //////if(read == false){eigD.resize(n0 + n1);}
//  for(Int n = 0; n < eigD.size(); n++){eigD[n].init(geo);}
//
//  std::vector<qlat::FieldM<qlat::ComplexF, 3> > eigF;eigF.resize(n1);
//  for (Int n = 0; n < eigF.size(); n++){eigF[n].init(geo);}
//
//
//  if(read == false){
//    for(Int n = 0; n < n0; n++){
//      if(fermion_type == 0){c0 = (kSpace[n]);}
//      if(fermion_type == 1){
//        c0.Even() = (kSpace[2*n+0]);
//        c0.Odd()  = (kSpace[2*n+1]);
//      }
//      //qlat::ComplexD* q = (qlat::ComplexD*) qlat::get_data(buf).data();
//      //cpy_data_thread(q, (qlat::ComplexD*) ctmp0->data(), Nvol, 0);
//      quda_cf_to_qlat_cf(eigD[n], (qlat::ComplexD*) c0.data());
//    }
//
//    ////Maybe needed
//    //eig_solveK->orthonormalizeMGS(fSpace, ns1);
//    for(Int n = 0; n < n1; n++){
//      if(fermion_type == 0){c0 = (kSpace[n + n0]);}
//      if(fermion_type == 1){
//        c0.Even() = (kSpace[2*(n+n0)+0]);
//        c0.Odd()  = (kSpace[2*(n+n0)+1]);
//      }
//      //qlat::ComplexD* q = (qlat::ComplexD*) qlat::get_data(buf).data();
//      //cpy_data_thread(q, (qlat::ComplexD*) ctmp0->data(), Nvol, 0);
//      quda_cf_to_qlat_cf(eigF[n], (qlat::ComplexD*) c0.data());
//    }
//  }
//
//  if(read == false){
//    if(n0!=0){single_file = false;qlat::load_qlat_noisesT(filename0.c_str(), eigD, read, single_file, VECS_TYPE, INFO_LIST, 0, n0, false);}
//    if(n1!=0){single_file = true ;qlat::load_qlat_noisesT(filename1.c_str(), eigF, read, single_file, VECS_TYPE, INFO_LIST, 0, n1, false);}
//
//  }
//  if(read == true){
//    if(n0!=0){qlat::load_qlat_noisesT(filename0.c_str(), eigD, read, single_file, VECS_TYPE, INFO_LIST, 0,   n0, false);}
//
//    if(split_save == 0)
//    if(n1!=0){qlat::load_qlat_noisesT(filename0.c_str(), eigF, read, single_file, VECS_TYPE, INFO_LIST,n0,n0+n1, false);}
//    if(split_save == 1)
//    if(n1!=0){qlat::load_qlat_noisesT(filename1.c_str(), eigF, read, single_file, VECS_TYPE, INFO_LIST, 0,   n1, false);}
//  }
//
//  /////===load to gpu
//  if(read == true){
//  //for (Int i = 0; i < kSpace.size(); i++) delete kSpace[i];kSpace.resize(0);
//  //for (Int i = 0; i < fSpace.size(); i++) delete fSpace[i];fSpace.resize(0);
//  kSpace.resize(ns0);
//  fSpace.resize(ns1);
//
//  eig_param.n_ev_deflate = ns0;eig_solveK = quda::EigenSolver::create(&eig_param, *mat_E);
//  eig_param.n_ev_deflate = ns1;eig_solveF = quda::EigenSolver::create(&eig_param, *mat_E);
//  quda::ColorSpinorParam csD = quda::ColorSpinorParam(cs_gpuD);
//  quda::ColorSpinorParam csF = quda::ColorSpinorParam(cs_gpuF);
//  csD.is_composite  = false;csD.is_component  = false;
//  csF.is_composite  = false;csF.is_component  = false;
//
//
//  for(Int n = 0; n < ns0; n++){kSpace[n] = quda::ColorSpinorField(csD);}
//  for(Int n = 0; n < ns1; n++){fSpace[n] = quda::ColorSpinorField(csF);}
//  for(Int n = 0; n < n0; n++){
//    qlat_cf_to_quda_cf((qlat::ComplexD*) c0.data(), eigD[n]);
//    if(fermion_type == 0){(kSpace[n]) = c0;}
//    if(fermion_type == 1){
//      (kSpace[2*n+0]) = c0.Even();
//      (kSpace[2*n+1]) = c0.Odd() ;
//    }
//  }
//
//  for(Int n = 0; n < n1; n++){
//    qlat_cf_to_quda_cf((qlat::ComplexD*) c0.data(), eigF[n]);
//    if(fermion_type == 0){(fSpace[n]) = c0;}
//    if(fermion_type == 1){
//      (fSpace[2*n+0]) = c0.Even();
//      (fSpace[2*n+1]) = c0.Odd() ;
//    }
//  }
//  }
//  //eig_solveF->orthonormalizeMGS(fSpace, ns1);
//
//  std::string fileE = ssprintf("%s.evals", filename);
//
//  std::vector<double > values, errors;
//  values.resize(2*nvec); errors.resize(nvec);
//  if(read == false){
//    Qassert(evals_ZERO.size() != 0 and evals_ERR.size() != 0);
//    for(Int n=0;n<nvec;n++){
//      values[n*2+0] = evals_ZERO[n].real();
//      values[n*2+1] = evals_ZERO[n].imag();
//      errors[n] = evals_ERR[n]; 
//    }
//    save_txt_eigenvalues(values, errors, fileE.c_str(), "Staggered Fermions");
//  }
//
//  if(read == true){
//    load_txt_eigenvalues(values, errors, fileE.c_str());Qassert(ns0 + ns1 <= values.size()/2);
//
//    evals_ZERO.resize(ns0+ns1);evals_ERR.resize(ns0+ns1);
//    for(Int n=  0;n<ns0+ns1;n++){evals_ZERO[n] = quda::Complex(values[n*2+0] - 4.0*mre*mre, values[n*2+1]);}    
//
//    for(Int n=  0;n<ns0+ns1;n++){evals_ERR[ n] = errors[n];}    
//
//    evalsK.resize(ns0, 0.0);evalsF.resize(ns1, 0.0);
//    update_eigen_mass(0.0, true);
//  }
//
//}

inline void quda_inverter::free_mem(){
  TIMER("quda free_mem");
  buf_inv.resize(0);
  V = 0;
  for(Int i=0;i<4;i++){X[i] = 0;}

  //for (Int i = 0; i < evecs_.size(); i++) delete evecs_[i];
  //evecs_.resize(0);
  //evecs.resize(0);
  //kSpace.resize(0);
  //fSpace.resize(0);
  //ZSpace.resize(0);
  ////kSpace_cpu.resize(0);
  //for (Int i = 0; i < kSpace.size(); i++) delete kSpace[i];kSpace.resize(0);
  //for (Int i = 0; i < fSpace.size(); i++) delete fSpace[i];fSpace.resize(0);
  //for (Int i = 0; i < ZSpace.size(); i++) delete ZSpace[i];ZSpace.resize(0);
  //for (Int i = 0; i < kSpace_cpu.size(); i++) delete kSpace_cpu[i];kSpace_cpu.resize(0);
  //evalsK.resize(0);
  //evalsF.resize(0);
  //evalsZ.resize(0);
  free_csfield();

  out_multishift_gpu.resize(0);
  src_multishift_gpu.resize(0);

  //if(eig_solveF != NULL){delete eig_solveF;eig_solveF=NULL;}
  //if(eig_solveK != NULL){delete eig_solveK;eig_solveK=NULL;}
  //nvec = 0;

  //if(df_preconditioner != NULL){destroyDeflationQuda(df_preconditioner);df_preconditioner=NULL;}
  //clear_CG();
  clear_mat();
}

//inline void quda_inverter::setup_clover(const double clover_csw)
//{
//  fermion_type = 0;
//  if(fermion_type == 0){spinor_site_size = 12;}
//
//  /////===Start of Inv parameters
//  inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH;
//
//  /////double kappa      = kappa;
//  double anisotropy = gauge_param.anisotropy;
//  inv_param.kappa = 0.150;
//  inv_param.mass = 0.5 / inv_param.kappa - (1.0 + 3.0 / anisotropy);
//
//  // Use 3D or 4D laplace
//  //===inv_param.laplace3D = laplace3D;
//  inv_param.Ls = 1;
//
//  inv_param.cpu_prec                      = QUDA_DOUBLE_PRECISION;
//  inv_param.cuda_prec                     = QUDA_DOUBLE_PRECISION;
//  inv_param.cuda_prec_sloppy              = QUDA_SINGLE_PRECISION;
//  inv_param.cuda_prec_refinement_sloppy   = QUDA_SINGLE_PRECISION;
//
//
//  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
//  inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
//  inv_param.dirac_order = QUDA_DIRAC_ORDER;
//
//
//  inv_param.clover_cpu_prec               = QUDA_DOUBLE_PRECISION;
//  inv_param.clover_cuda_prec              = QUDA_DOUBLE_PRECISION;
//
//  inv_param.clover_cuda_prec_sloppy       = QUDA_SINGLE_PRECISION;
//  inv_param.clover_cuda_prec_refinement_sloppy = QUDA_SINGLE_PRECISION;
//
//  ////related to eigensystem
//  inv_param.cuda_prec_eigensolver         = QUDA_SINGLE_PRECISION;
//  inv_param.cuda_prec_precondition        = QUDA_SINGLE_PRECISION;
//  inv_param.clover_cuda_prec_precondition = QUDA_SINGLE_PRECISION;
//  inv_param.clover_cuda_prec_eigensolver  = QUDA_SINGLE_PRECISION;
//
//  inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
//  // Use kappa * csw or supplied clover_coeff
//  bool compute_clover_trlog = false;
//  //bool compute_clover_trlog = true;
//  inv_param.clover_csw = clover_csw;
//  inv_param.clover_coeff = inv_param.clover_csw * inv_param.kappa;
//  /////===unknow para
//  inv_param.compute_clover_trlog = compute_clover_trlog ? 1 : 0;
//
//  // General parameter setup
//  inv_param.inv_type = QUDA_CG_INVERTER;
//
//  // solution_type specifies *what* system is to be solved.
//  // solve_type specifies *how* the system is to be solved.
//  //
//  // We have the following four cases (plus preconditioned variants):
//  //
//  // solution_type    solve_type    Effect
//  // -------------    ----------    ------
//  // MAT              DIRECT        Solve Ax=b
//  // MATDAG_MAT       DIRECT        Solve A^dag y = b, followed by Ax=y
//  // MAT              NORMOP        Solve (A^dag A) x = (A^dag b)
//  // MATDAG_MAT       NORMOP        Solve (A^dag A) x = b
//  // MAT              NORMERR       Solve (A A^dag) y = b, then x = A^dag y
//  //
//  // We generally require that the solution_type and solve_type
//  // preconditioning match.  As an exception, the unpreconditioned MAT
//  // solution_type may be used with any solve_type, including
//  // DIRECT_PC and NORMOP_PC.  In these cases, preparation of the
//  // preconditioned source and reconstruction of the full solution are
//  // taken care of by Dirac::prepare() and Dirac::reconstruct(),
//  // respectively.
//
//  // Eventually we want the unpreconditioned solution.
//  inv_param.solution_type = QUDA_MAT_SOLUTION;
//
//  // NORMOP_PC means preconditioned normal operator MdagM
//  //inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
//  inv_param.solve_type = QUDA_NORMOP_SOLVE;
//
//  // There might be a performance difference.
//  inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
//
//  inv_param.dagger = QUDA_DAG_NO;
//
//  //inv_param.mass_normalization   = QUDA_KAPPA_NORMALIZATION;
//  inv_param.mass_normalization   = QUDA_MASS_NORMALIZATION;
//  inv_param.solver_normalization =  QUDA_SOURCE_NORMALIZATION;
//  //inv_param.solver_normalization = QUDA_DEFAULT_NORMALIZATION;
//  //inv_param.solver_normalization = QUDA_KAPPA_NORMALIZATION;
//
//  ///unknown
//  Int gcrNkrylov = 10;
//  QudaCABasis ca_basis = QUDA_POWER_BASIS;
//  double ca_lambda_min = 0.0;
//  double ca_lambda_max = -1.0;
//
//  inv_param.pipeline = 1;
//  inv_param.Nsteps = 2;
//  inv_param.gcrNkrylov = gcrNkrylov;
//  inv_param.ca_basis = ca_basis;
//  inv_param.ca_lambda_min = ca_lambda_min;
//  inv_param.ca_lambda_max = ca_lambda_max;
//
//  inv_param.tol     = 1e-10;
//  inv_param.maxiter = 100000;
//  inv_param.tol_restart = 0.0005;
//  //if (tol_hq == 0 && tol == 0) {
//  //  errorQuda("qudaInvert: requesting zero residual\n");
//  //  exit(1);
//  //}
//
//  // require both L2 relative and heavy quark residual to determine convergence
//  inv_param.residual_type = static_cast<QudaResidualType_s>(0);
//  inv_param.residual_type = static_cast<QudaResidualType_s>(inv_param.residual_type | QUDA_L2_RELATIVE_RESIDUAL);
//
//  inv_param.tol_hq = 1e-5; // specify a tolerance for the residual for heavy quark residual
//  //inv_param.tol_hq = 1e-10; // specify a tolerance for the residual for heavy quark residual
//
//  //// Offsets used only by multi-shift solver
//  //// These should be set in the application code. We set the them here by way of
//  //// example
//  //inv_param.num_offset = multishift;
//  //for (Int i = 0; i < inv_param.num_offset; i++) inv_param.offset[i] = 0.06 + i * i * 0.1;
//  //// these can be set individually
//  //for (Int i = 0; i < inv_param.num_offset; i++) {
//  //  inv_param.tol_offset[i] = inv_param.tol;
//  //  inv_param.tol_hq_offset[i] = inv_param.tol_hq;
//  //}
//
//  // This is for Quda's sophisticated reliable update. 0.1 should be good.
//  inv_param.reliable_delta = 0.1;
//  inv_param.use_alternative_reliable = true;
//  inv_param.use_sloppy_partial_accumulator = 0;
//  inv_param.solution_accumulator_pipeline = 1;
//  inv_param.max_res_increase = 1;
//
//  // domain decomposition preconditioner parameters
//  // inv_param.inv_type_precondition = precon_type;
//
//  //inv_param.schwarz_type = precon_schwarz_type;
//  //inv_param.precondition_cycle = precon_schwarz_cycle;
//  //inv_param.tol_precondition = tol_precondition;
//  //inv_param.maxiter_precondition = maxiter_precondition;
//  //inv_param.omega = 1.0;
//
//  //inv_param.input_location  = QUDA_CPU_FIELD_LOCATION;
//  //inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
//
//  inv_param.input_location  = QUDA_CUDA_FIELD_LOCATION;
//  inv_param.output_location = QUDA_CUDA_FIELD_LOCATION;
//
//  // QUDA_DEBUG_VERBOSE is too nasty.
//
//  inv_param.extlib_type = QUDA_EIGEN_EXTLIB;
//
//  // Whether or not to use native BLAS LAPACK
//  //inv_param.native_blas_lapack = (native_blas_lapack ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE);
//  inv_param.native_blas_lapack = QUDA_BOOLEAN_TRUE;
//
//  //// Whether or not use fused kernels for Mobius
//  ////inv_param.use_mobius_fused_kernel = use_mobius_fused_kernel ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
//  //inv_param.use_mobius_fused_kernel = QUDA_BOOLEAN_FALSE;
//
//
//  ///////===check this variable
//  std::array<int, 4> grid_partition = {1, 1, 1, 1};
//  inv_param.split_grid[0] = grid_partition[0];
//  inv_param.split_grid[1] = grid_partition[1];
//  inv_param.split_grid[2] = grid_partition[2];
//  inv_param.split_grid[3] = grid_partition[3];
//
//  inv_param.struct_size = sizeof(inv_param);
//
//  ////===END of Inv parameters
//
//  ///int do_inv = 1;
//  {
//  //////===operator define
//  //void *clover = nullptr;
//  //void *clover_inv = nullptr;
//
//  Int clover_site_size           = 72; // real numbers per block-diagonal clover matrix
//  ////size_t host_clover_data_type_size = (cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
//  Int host_clover_data_type_size = sizeof(quda::Complex)/2;
//  ////size_t host_spinor_data_type_size = (cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
//  Int host_spinor_data_type_size = sizeof(quda::Complex)/2;
//
//  quda_clover.resize(    V * clover_site_size * host_clover_data_type_size);
//  quda_clover_inv.resize(V * clover_site_size * host_spinor_data_type_size);
//
//  //constructHostCloverField((void*) quda_clover.data(), (void*) quda_clover_inv.data(), inv_param);
//  /////===host
//
//  bool compute_clover = true;
//  //if(in.clover_csw == 0){compute_clover = false;}
//  //double norm = 0.00; // clover components are random numbers in the range (-norm, norm)
//  //double diag = 1.0;  // constant added to the diagonal
//  ////if (!compute_clover) constructQudaCloverField((void*) quda_clover.data(), norm, diag, inv_param.clover_cpu_prec, V);
//  inv_param.compute_clover = compute_clover;
//  if (compute_clover) inv_param.return_clover = 1;
//  inv_param.compute_clover_inverse = 1;
//  inv_param.return_clover_inverse  = 1;
//
//  //loadCloverQuda((void*) quda_clover.data(), (void*)  quda_clover_inv.data(), &inv_param);
//  //if(in.clover_csw == 0){
//  //inv_param.compute_clover_inverse = 0;
//  //inv_param.return_clover_inverse  = 0;
//  //}
//  /////===host
//  ///// Load the clover terms to the device
//  //loadCloverQuda((void*) quda_clover.data(), (void*)  quda_clover_inv.data(), &inv_param);
//  //loadCloverQuda((void*) quda_clover.data(), (void*)  quda_clover_inv.data(), &inv_param);
//  loadCloverQuda((void*) quda_clover.data(), (void*)  quda_clover_inv.data(), &inv_param);
//  }
//
//  alloc_csfield_cpu();
//  alloc_csfield_gpu();
//
//}

inline void quda_inverter::setup_stagger_inv()
{
  TIMER("setup_stagger_inv");
  fermion_type = 1;
  if(fermion_type == 1){spinor_site_size = 3 ;}
  //if(gauge_with_phase == false){errorQuda("Quda stagger need link phases! \n");}

  /////===Start of Inv parameters
  inv_param.dslash_type = QUDA_STAGGERED_DSLASH;
  ////inv_param.mass  = fermion_mass;
  //inv_param.mass  = 0.2;
  inv_param.mass  = 0.2;

  // Use 3D or 4D laplace
  //===inv_param.laplace3D = laplace3D;
  inv_param.Ls = 1;

  inv_param.cpu_prec                      = QUDA_DOUBLE_PRECISION;
  inv_param.cuda_prec                     = QUDA_DOUBLE_PRECISION;
  //inv_param.cuda_prec_precondition   = QUDA_DOUBLE_PRECISION;

  inv_param.cuda_prec_sloppy              = QUDA_SINGLE_PRECISION;
  inv_param.cuda_prec_refinement_sloppy   = QUDA_SINGLE_PRECISION;
  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
  inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  ////inv_param.gamma_basis = QUDA_UKQCD_GAMMA_BASIS;
  //inv_param.dirac_order = QUDA_DIRAC_ORDER;
  inv_param.dirac_order = QUDA_INTERNAL_DIRAC_ORDER;

  ////related to eigensystem
  inv_param.cuda_prec_eigensolver         = QUDA_SINGLE_PRECISION;
  inv_param.cuda_prec_precondition        = QUDA_SINGLE_PRECISION;

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
  ////inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
  inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
  //inv_param.solve_type = QUDA_NORMOP_SOLVE;
  //if(in.solver_type == 0){inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;}
  //if(in.solver_type == 1){inv_param.solve_type = QUDA_NORMERR_PC_SOLVE;}
  //if(in.solver_type == 2){inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;}


  // There might be a performance difference.
  inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;

  inv_param.dagger = QUDA_DAG_NO;

  //inv_param.mass_normalization   = QUDA_KAPPA_NORMALIZATION;
  inv_param.mass_normalization   = QUDA_MASS_NORMALIZATION;
  inv_param.solver_normalization =  QUDA_SOURCE_NORMALIZATION;
  //inv_param.solver_normalization = QUDA_DEFAULT_NORMALIZATION;
  //inv_param.solver_normalization = QUDA_DEFAULT_NORMALIZATION;

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

  //inv_param.input_location  = QUDA_CPU_FIELD_LOCATION;
  //inv_param.output_location = QUDA_CPU_FIELD_LOCATION;

  inv_param.input_location  = QUDA_CUDA_FIELD_LOCATION;
  inv_param.output_location = QUDA_CUDA_FIELD_LOCATION;

  // QUDA_DEBUG_VERBOSE is too nasty.

  inv_param.extlib_type = QUDA_EIGEN_EXTLIB;

  // Whether or not to use native BLAS LAPACK
  //inv_param.native_blas_lapack = (native_blas_lapack ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE);
  inv_param.native_blas_lapack = QUDA_BOOLEAN_TRUE;
  ///inv_param.native_blas_lapack = QUDA_BOOLEAN_FALSE;

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

  alloc_csfield_initial();
  alloc_csfield_gpu();
}

////0, double to single, 1 single to half, 10 double to double, 11 single to single, 12 half to half
////may still cost time even using the comparisons, TODO need to investigate the issues
inline void quda_inverter::setup_inv_param_prec_type(Int prec_type)
{
  TIMER("setup_inv_param_prec_type");
  inv_param.cpu_prec                      = QUDA_DOUBLE_PRECISION;
  if(prec_type == 0)
  {
    inv_param.cpu_prec                      = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec                     = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec_sloppy              = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_refinement_sloppy   = QUDA_SINGLE_PRECISION;
  }

  if(prec_type == 1)
  {
    inv_param.cpu_prec                      = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec                     = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_sloppy              = QUDA_HALF_PRECISION;
    inv_param.cuda_prec_refinement_sloppy   = QUDA_HALF_PRECISION;
  }

  if(prec_type == 2)
  {
    inv_param.cpu_prec                      = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec                     = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec_sloppy              = QUDA_HALF_PRECISION;
    inv_param.cuda_prec_refinement_sloppy   = QUDA_HALF_PRECISION;
  }

  
  if(prec_type == 10)
  {
    inv_param.cpu_prec                      = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec                     = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec_sloppy              = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec_refinement_sloppy   = QUDA_DOUBLE_PRECISION;
  }
 
  if(prec_type == 11)
  {
    inv_param.cpu_prec                      = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec                     = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_sloppy              = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_refinement_sloppy   = QUDA_SINGLE_PRECISION;
  }

  Qassert(prec_type != 12);
  //if(prec_type == 12)
  //{
  //  inv_param.cpu_prec                      = QUDA_SINGLE_PRECISION;
  //  inv_param.cuda_prec                     = QUDA_HALF_PRECISION;
  //  inv_param.cuda_prec_sloppy              = QUDA_HALF_PRECISION;
  //  inv_param.cuda_prec_refinement_sloppy   = QUDA_HALF_PRECISION;
  //}

  gauge_param.cuda_prec              = inv_param.cuda_prec;
  gauge_param.cuda_prec_sloppy       = inv_param.cuda_prec_sloppy;
  gauge_param.struct_size = sizeof(gauge_param);
  prec_type_check = -2;////need reload if this function is called

}

////0, double to single, 1 single to half, 10 double to double, 11 single to single, 12 half to half
////may still cost time even using the comparisons, TODO need to investigate the issues
inline void quda_inverter::setup_inv_param_prec(Int prec_type, bool force_reload)
{
  if(prec_type_check == prec_type and force_reload == false){return ;}
  TIMER("setup_inv_param_prec");
  if(prec_type_check == -2){force_reload = true;}///force_reload if setup_inv_param_prec_type is called
  inv_param.cpu_prec                      = QUDA_DOUBLE_PRECISION;
  if(prec_type == 0)
  {
    inv_param.cuda_prec                     = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec_sloppy              = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_refinement_sloppy   = QUDA_SINGLE_PRECISION;
  }

  if(prec_type == 1)
  {
    inv_param.cuda_prec                     = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_sloppy              = QUDA_HALF_PRECISION;
    inv_param.cuda_prec_refinement_sloppy   = QUDA_HALF_PRECISION;
  }

  if(prec_type == 2)
  {
    inv_param.cuda_prec                     = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec_sloppy              = QUDA_HALF_PRECISION;
    inv_param.cuda_prec_refinement_sloppy   = QUDA_HALF_PRECISION;
  }

  
  if(prec_type == 10)
  {
    inv_param.cuda_prec                     = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec_sloppy              = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec_refinement_sloppy   = QUDA_DOUBLE_PRECISION;
  }
 
  if(prec_type == 11)
  {
    inv_param.cuda_prec                     = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_sloppy              = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_refinement_sloppy   = QUDA_SINGLE_PRECISION;
  }

  Qassert(prec_type != 12);
  //if(prec_type == 12)
  //{
  //  inv_param.cuda_prec                     = QUDA_HALF_PRECISION;
  //  inv_param.cuda_prec_sloppy              = QUDA_HALF_PRECISION;
  //  inv_param.cuda_prec_refinement_sloppy   = QUDA_HALF_PRECISION;
  //}

  if( gauge_param.cuda_prec != inv_param.cuda_prec or gauge_param.cuda_prec_sloppy != inv_param.cuda_prec_sloppy or force_reload)
  {
    Qassert(quda_gf_default.size() != 0);
    gauge_param.cuda_prec              = inv_param.cuda_prec;
    gauge_param.cuda_prec_sloppy       = inv_param.cuda_prec_sloppy;
    gauge_param.struct_size = sizeof(gauge_param);
    //gauge_param.cuda_prec_precondition = inv_param.cuda_prec_sloppy;
    //gauge_param.cuda_prec_eigensolver  = inv_param.cuda_prec_sloppy;

    //////gauge_param.cpu_prec               = QUDA_DOUBLE_PRECISION;
    //////gauge_param.cuda_prec              = QUDA_DOUBLE_PRECISION;

    //setup_gauge_param(gauge_param.t_boundary);
    //Qassert(fermion_type == 1);
    //freeGaugeQuda();
    setup_link(&quda_gf_default[0], 0);
    alloc_csfield_gpu();
    //CG_reset = true;
    setup_mat_mass(mass_mat, true);
  }
  prec_type_check = prec_type;
}

inline void quda_inverter::clear_mat()
{
  TIMER("clear_mat");
  if(dirac != NULL){delete dirac;dirac=NULL;}
  //if(dirac_cg != NULL){delete dirac_cg;dirac_cg=NULL;}
  if(dirac_pc != NULL){delete dirac_pc;dirac_pc=NULL;}
  if(dirac_pc_odd != NULL){delete dirac_pc_odd;dirac_pc_odd=NULL;}
  //if(dSloppy != NULL){delete dSloppy;dSloppy=NULL;}
  //if(dPre != NULL){delete dPre;dPre=NULL;}
  //if(dEig != NULL){delete dEig;dEig=NULL;}

  if(mat       != NULL){delete mat      ;  mat       = NULL;}
  if(mat_pc    != NULL){delete mat_pc   ;  mat_pc    = NULL;}
  if(mat_Mdag  != NULL){delete mat_Mdag ;  mat_Mdag  = NULL;}
  if(mat_MMdag != NULL){delete mat_MMdag;  mat_MMdag = NULL;}
  if(mat_MdagM != NULL){delete mat_MdagM;  mat_MdagM = NULL;}

  //if(m_cg != NULL){delete m_cg;  m_cg = NULL;}
  //if(mSloppy != NULL){delete mSloppy;  mSloppy = NULL;}
  //if(mPre != NULL){delete mPre;  mPre = NULL;}
  //if(mEig != NULL){delete mEig;  mEig = NULL;}

  mat_E = NULL;
  //clear_CG();
}

// equal mass to inverter and set verbose
//inline void quda_inverter::setup_inv_mass(const double mass)
//{
//  TIMER("setup_inv_mass");
//  Int verbos = quda_verbos;
//  Qassert(fermion_type == 1);//only stagger
//  if(verbos <= -1)
//  {
//    setVerbosity(QUDA_SILENT);
//    inv_param.verbosity   = QUDA_SILENT;
//  }
//
//  if(verbos == 0)
//  {
//    setVerbosity(QUDA_SUMMARIZE);
//    inv_param.verbosity   = QUDA_SUMMARIZE;
//  }
//  if(verbos == 1)
//  {
//    setVerbosity(QUDA_VERBOSE);
//    inv_param.verbosity   = QUDA_VERBOSE;
//  }
//
//  if(verbos == 2)
//  {
//    setVerbosity(QUDA_DEBUG_VERBOSE);
//    inv_param.verbosity   = QUDA_VERBOSE;
//  }
//
//  //if(fermion_type == 0)
//  //{
//  //  //inv_param.kappa = 0.130;
//  //  //inv_param.mass = 0.5 / kappa - (1.0 + 3.0 / anisotropy);
//  //  double anisotropy = gauge_param.anisotropy;
//  //  inv_param.mass = mass;
//  //  inv_param.kappa = 0.5 / (mass + (1.0 + 3.0 / anisotropy));
//  //  ////double clover_csw setup with clover
//  //  inv_param.clover_coeff = inv_param.clover_csw * inv_param.kappa;
//  //  ////if(quda::comm_rank()== 0)printfQuda("Kappa = %.8f Mass = %.8f\n", inv_param.kappa, inv_param.mass);
//  //  return ;
//  //}
//
//  if(fermion_type == 1){inv_param.mass = mass;return ;}
//
//  errorQuda("Fermion type not found!\n");
//
//}

// set up mass and Dirac
inline void quda_inverter::setup_mat_mass(const double mass, const bool force_do)
{
  TIMER("setup_mat_mass");
  Int flag_do = 0;
  if(mat == NULL){flag_do = 1;}
  if(std::fabs(mass_mat - mass) > 1e-15){flag_do = 1;}
  if(force_do == true){flag_do = 1;}
  if(flag_do == 0){return ;}

  ///setup verbose
  Int verbos = quda_verbos;
  Qassert(fermion_type == 1);//only stagger
  if(verbos <= -1)
  {
    setVerbosity(QUDA_SILENT);
    inv_param.verbosity   = QUDA_SILENT;
  }

  if(verbos == 0)
  {
    setVerbosity(QUDA_SUMMARIZE);
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

  //if(fermion_type == 0)
  //{
  //  //inv_param.kappa = 0.130;
  //  //inv_param.mass = 0.5 / kappa - (1.0 + 3.0 / anisotropy);
  //  double anisotropy = gauge_param.anisotropy;
  //  inv_param.mass = mass;
  //  inv_param.kappa = 0.5 / (mass + (1.0 + 3.0 / anisotropy));
  //  ////double clover_csw setup with clover
  //  inv_param.clover_coeff = inv_param.clover_csw * inv_param.kappa;
  //  ////if(quda::comm_rank()== 0)printfQuda("Kappa = %.8f Mass = %.8f\n", inv_param.kappa, inv_param.mass);
  //  return ;
  //}

  if(fermion_type == 1){inv_param.mass = mass;}
  //errorQuda("Fermion type not found!\n");
  //setup_inv_mass(mass);

  mass_mat = mass;

  /////qmessage("requested precision %d\n",inv_param.cuda_prec);
  //createDiracWithEig(dirac, dSloppy, dPre, dEig, inv_param, pc_solve);
  //createDirac(dirac, dSloppy, dPre, inv_param, pc_solve);

  //QudaSolveType;
  //QudaSolutionType
  clear_mat();
  //bool pc_solve = (inv_param.solve_type == QUDA_DIRECT_PC_SOLVE) ||
  //  (inv_param.solve_type == QUDA_NORMOP_PC_SOLVE) || (inv_param.solve_type == QUDA_NORMERR_PC_SOLVE);
  bool pc_solve = false;

  quda::DiracParam diracParam;
  setDiracParam(diracParam, &inv_param, pc_solve);

  //diracParam.tmp1 = gtmp1;
  //diracParam.tmp2 = gtmp2;

  dirac    = quda::Dirac::create(diracParam);

  inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
  setDiracParam(diracParam, &inv_param, true);
  dirac_pc_odd = quda::Dirac::create(diracParam);

  inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
  setDiracParam(diracParam, &inv_param, true);
  dirac_pc = quda::Dirac::create(diracParam);

  mat        = new quda::DiracM(*dirac);
  mat_pc     = new quda::DiracM(*dirac_pc);
  mat_MdagM  = new quda::DiracMdagM(dirac);
  mat_MMdag  = new quda::DiracMMdag(*dirac);
  mat_Mdag   = new quda::DiracMdag(*dirac);
  //CG_reset = true;

  ///CG related matrix
  //{
  //  //bool pc_solve = (inv_param.solve_type == QUDA_DIRECT_PC_SOLVE) ||
  //  //  (inv_param.solve_type == QUDA_NORMOP_PC_SOLVE) || (inv_param.solve_type == QUDA_NORMERR_PC_SOLVE);
  //  //createDiracWithEig(dirac_cg, dSloppy, dPre, dEig, inv_param, pc_solve);
  //  //m_cg       = new quda::DiracM(*dirac_cg);
  //  //mSloppy    = new quda::DiracM(*dSloppy);
  //  //mPre       = new quda::DiracM(*dPre);
  //  //mEig       = new quda::DiracM(*dEig);
  //}

  //DiracParam diracParam;
  //setDiracParam(diracParam, &inv_param, pc_solve);
  //dirac  = Dirac::create(diracParam);
}

//inline void quda_inverter::setup_inc_eigencg(const Int n_ev, const Int n_kr, const Int n_conv, const Int df_grid, const double tol, const double inc_tol, const double tol_restart, const Int restart_n, const Int pipeline, const Int inv_type)
//{
//  /////// inv_param.solve_type = QUDA_NORMOP_SOLVE;
//  /////// inv_param.solve_type = QUDA_DIRECT_SOLVE;
//  ///inv_param.inv_type = QUDA_CG_INVERTER;
//  if(inv_type == 0){inv_param.inv_type = QUDA_EIGCG_INVERTER    ;inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;}
//  //if(inv_type == 0){inv_param.inv_type = QUDA_EIGCG_INVERTER    ;inv_param.solve_type = QUDA_NORMOP_SOLVE;}
//  //if(inv_type == 1){inv_param.inv_type = QUDA_INC_EIGCG_INVERTER;inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;}
//  if(inv_type == 1){inv_param.inv_type = QUDA_INC_EIGCG_INVERTER;inv_param.solve_type = QUDA_NORMOP_SOLVE;}
//  if(inv_type == 2){inv_param.inv_type = QUDA_INC_EIGCG_INVERTER;inv_param.solve_type = QUDA_DIRECT_SOLVE;}
//  if(inv_type == 3){inv_param.inv_type = QUDA_INC_EIGCG_INVERTER;inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;}
//  if(inv_type == 4){inv_param.inv_type = QUDA_INC_EIGCG_INVERTER;inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;}
//  if(inv_type == 5){inv_param.inv_type = QUDA_INC_EIGCG_INVERTER;inv_param.solve_type = QUDA_NORMERR_PC_SOLVE;}
//  //if(inv_type == 2){inv_param.inv_type = QUDA_GMRESDR_INVERTER ;inv_param.solve_type = QUDA_NORMOP_SOLVE;}
//  //if(inv_type == 2){inv_param.inv_type = QUDA_GMRESDR_INVERTER ;inv_param.solve_type = QUDA_NORMOP_SOLVE;}
//  //if(inv_type == 3){inv_param.inv_type = QUDA_FGMRESDR_INVERTER;inv_param.solve_type = QUDA_NORMOP_SOLVE;}
//  //inv_type_precondition == QUDA_CG_INVERTER;
//
//  //if(fermion_type == 0)inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
//  //if(fermion_type == 1)inv_param.solve_type = QUDA_NORMOP_SOLVE;
//  
//  inv_param.use_init_guess = QUDA_USE_INIT_GUESS_NO;
//  inv_param.inc_tol= inc_tol;
//  if(pipeline ==  0){inv_param.pipeline = n_ev;}
//  else{              inv_param.pipeline = pipeline;}
//  
//
//  ////additional inv_param
//  inv_param.rhs_idx = 0;
//  inv_param.n_ev            = n_ev;
//  inv_param.max_search_dim  = n_kr;
//  inv_param.deflation_grid  = df_grid;
//  //inv_param.deflation_grid  = 16;
//  //inv_param.deflation_grid  = n_ev;
//  inv_param.tol_restart     = tol_restart;
//  //inv_param.tol_restart     = inc_tol*1e-2;
//  ////inv_param.tol_restart     = 5e+3*1e-7;
//
//  inv_param.eigcg_max_restarts = restart_n;
//  inv_param.max_restart_num    = restart_n;
//
//  inv_param.eigenval_tol = tol;
//  //inv_param.eigenval_tol = inc_tol*1e2;
//  //inv_param.inc_tol = inv_param.inc_tol;
//  //inv_param.eigenval_tol = 1e-1;
//  //inv_param.eigenval_tol = tol*1e2;
//
//  // domain decomposition preconditioner parameters
//  inv_param.schwarz_type = QUDA_INVALID_SCHWARZ;
//  inv_param.accelerator_type_precondition = QUDA_INVALID_ACCELERATOR;
//  inv_param.precondition_cycle = 1;
//  ////inv_param.tol_precondition = tol;
//  inv_param.tol_precondition = 1e-2;
//  inv_param.maxiter_precondition = 4;
//  inv_param.omega = 1.0;
//
//  inv_param.gcrNkrylov = 6;
//  inv_param.inv_type_precondition = QUDA_INVALID_INVERTER;
//  //inv_param.inv_type_precondition = QUDA_CG_INVERTER;
//  //inv_param.inv_type_precondition = QUDA_MG_INVERTER;
//  //inv_param.inv_type_precondition = QUDA_SD_INVERTER;
//  inv_param.extlib_type = QUDA_EIGEN_EXTLIB;
//
//  if(inc_tol > 1e-7){
//    inv_param.cuda_prec_eigensolver         = QUDA_SINGLE_PRECISION;
//    inv_param.cuda_prec_precondition        = QUDA_SINGLE_PRECISION;
//    inv_param.clover_cuda_prec_precondition = QUDA_SINGLE_PRECISION;
//    inv_param.clover_cuda_prec_eigensolver  = QUDA_SINGLE_PRECISION;
//
//    inv_param.cuda_prec_sloppy              = QUDA_SINGLE_PRECISION;
//    inv_param.clover_cuda_prec_sloppy       = QUDA_SINGLE_PRECISION;
//
//    gauge_param.cuda_prec_sloppy            = QUDA_SINGLE_PRECISION;
//    gauge_param.cuda_prec_precondition      = QUDA_SINGLE_PRECISION;
//    gauge_param.cuda_prec_eigensolver       = QUDA_SINGLE_PRECISION;
//
//    //inv_param.cuda_prec_refinement_sloppy   = QUDA_SINGLE_PRECISION;
//    //inv_param.clover_cuda_prec_refinement_sloppy = QUDA_SINGLE_PRECISION;
//  }else{
//    inv_param.cuda_prec_eigensolver         = QUDA_DOUBLE_PRECISION;
//    inv_param.cuda_prec_precondition        = QUDA_DOUBLE_PRECISION;
//    inv_param.clover_cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
//    inv_param.clover_cuda_prec_eigensolver  = QUDA_DOUBLE_PRECISION;
//    gauge_param.cuda_prec_precondition      = QUDA_DOUBLE_PRECISION;
//    gauge_param.cuda_prec_eigensolver       = QUDA_DOUBLE_PRECISION;
//
//
//    inv_param.cuda_prec_sloppy              = QUDA_DOUBLE_PRECISION;
//    inv_param.clover_cuda_prec_sloppy       = QUDA_DOUBLE_PRECISION;
//    gauge_param.cuda_prec_sloppy            = QUDA_DOUBLE_PRECISION;
//
//    //inv_param.cuda_prec_sloppy              = QUDA_SINGLE_PRECISION;
//    //inv_param.clover_cuda_prec_sloppy       = QUDA_SINGLE_PRECISION;
//    //gauge_param.cuda_prec_sloppy            = QUDA_SINGLE_PRECISION;
//
//    //inv_param.cuda_prec_refinement_sloppy   = QUDA_DOUBLE_PRECISION;
//    //inv_param.clover_cuda_prec_refinement_sloppy = QUDA_DOUBLE_PRECISION;
//
//    //inv_param.cuda_prec_sloppy              = QUDA_SINGLE_PRECISION;
//    //inv_param.cuda_prec_refinement_sloppy   = QUDA_SINGLE_PRECISION;
//    //inv_param.clover_cuda_prec_sloppy       = QUDA_SINGLE_PRECISION;
//    //inv_param.clover_cuda_prec_refinement_sloppy = QUDA_SINGLE_PRECISION;
//
//  }
//
//
//  inv_param.struct_size = sizeof(inv_param);
//
//  ////clear_mat();
//  alloc_csfield_cpu();
//  alloc_csfield_gpu();
//
//  //cs_cpu.print();
//  //cs_gpu.print();
//
//  eig_param = newQudaEigParam();
//  eig_param.preserve_deflation = QUDA_BOOLEAN_TRUE;
//
//  eig_param.invert_param = &inv_param;
//
//  eig_param.import_vectors = QUDA_BOOLEAN_FALSE;
//  eig_param.run_verify = QUDA_BOOLEAN_FALSE;
//
//  strcpy(eig_param.vec_infile, "");
//  strcpy(eig_param.vec_outfile, "");
//  eig_param.io_parity_inflate = QUDA_BOOLEAN_FALSE;
//
//  eig_param.extlib_type = QUDA_EIGEN_EXTLIB;
//
//  /////TODO mixed precision
//  //eig_param.cuda_prec_ritz =  QUDA_SINGLE_PRECISION;
//  if(inc_tol > 1e-7){
//  eig_param.cuda_prec_ritz =  QUDA_SINGLE_PRECISION;
//  inv_param.cuda_prec_ritz =  QUDA_SINGLE_PRECISION;
//  }else{
//  eig_param.cuda_prec_ritz =  QUDA_DOUBLE_PRECISION;
//  inv_param.cuda_prec_ritz =  QUDA_DOUBLE_PRECISION;}
//  //eig_param.cuda_prec_ritz =  QUDA_DOUBLE_PRECISION;
//  //if(tol > 1e-6){eig_param.cuda_prec_ritz = QUDA_SINGLE_PRECISION;}
//  //else{eig_param.cuda_prec_ritz =  QUDA_DOUBLE_PRECISION;}
//  //eig_param.cuda_prec_ritz = QUDA_DOUBLE_PRECISION;
//  eig_param.location       = QUDA_CUDA_FIELD_LOCATION;
//  eig_param.mem_type_ritz  = QUDA_MEMORY_DEVICE;
//
//  eig_param.n_ev   = n_ev;
//  eig_param.n_kr   = n_kr;
//  if(n_conv == 0){eig_param.n_conv = n_conv;}
//  else{           eig_param.n_conv = n_kr + n_ev;}
//  eig_param.tol    = tol;
//
//  eig_param.max_restarts     = inv_param.eigcg_max_restarts;
//
//  eig_param.n_ev_deflate = eig_param.n_conv;
//  eig_param.require_convergence = QUDA_BOOLEAN_TRUE;
//  ////eig_param.require_convergence = QUDA_BOOLEAN_FALSE;
//
//  eig_param.nk  = eig_param.n_ev;
//  eig_param.np  = eig_param.n_ev * inv_param.deflation_grid;
//
//  if(df_preconditioner != NULL){destroyDeflationQuda(df_preconditioner);df_preconditioner=NULL;}
//  df_preconditioner  = newDeflationQuda(&eig_param);
//  inv_param.deflation_op   = df_preconditioner;
//
//}

//inline void quda_inverter::update_eigen_mass(const double mass, bool force)
//{
//  TIMER("update_eigen_mass");
//  if(mass_value != mass or force)
//  {
//    Qassert(evals_ZERO.size() == evalsK.size() + evalsF.size());
//    Int n0  = evalsK.size();
//    Int n1  = evalsF.size();
//    Int nvec = n0+n1;
//    ///evalsZ.resize(nvec);
//    for(Int ni=0;ni<nvec;ni++)
//    {
//      quda::Complex tmp = quda::Complex(evals_ZERO[ni].real() + 4.0 * mass * mass, 0.0);
//      if(ni <  n0){evalsK[ni     ] = tmp;}
//      if(ni >= n0){evalsF[ni - n0] = tmp;}
//      ///evalsZ[ni] = tmp;
//    }
//  }
//  mass_value = mass;
//}

inline void quda_inverter::random_src(const Int seed)
{
  for(Int i=0;i<int(seed)%20 + 20*quda::comm_rank_global();i++){quda::comm_drand();}
  Int random_mode = 1;
  random_Ty((qlat::ComplexD*) gsrc->data(), gsrc->Volume() * spinor_site_size, 1, seed + 111111, random_mode);
  quda::blas::ax(0.05, *gsrc);
  quda::Complex* tmp = (quda::Complex*) (gsrc->data());
  Long totalN = gsrc->Volume() * spinor_site_size;
  //for(Long i=0;i<totalN/200 + 1;i++)
  //for(Long i=0;i< 50;i++)
  //{
  //  Long vi = Long(quda::comm_drand()*totalN);
  //  double ri = quda::comm_drand()*totalN;
  //  tmp[vi] = ri * quda::Complex(quda::comm_drand(), quda::comm_drand());
  //}
}

/////reconstruct and push_back to Fspace with double prec
//inline void quda_inverter::reconstruct_full(const double mass)
//{
//  TIMER("reconstruct_full");
//  std::vector<quda::Complex  >& vals = evalsZ;
//  if(use_eigen_pc == 0 or nvec <= 0){
//    ////for (Int i = 0; i < ZSpace.size(); i++){delete ZSpace[i];}
//    ZSpace.resize(0); evalsZ.resize(0);
//    return;
//  }
//  Qassert(kSpace.size() + fSpace.size() == nvec);
//  const Int off = kSpace.size();
//  quda::Complex Im = quda::Complex(0.0, 1.0);
//  setup_mat_mass(mass, true);
//  update_eigen_mass(mass, true);
//  double m = mass;
//  if(mass == -1){m = inv_param.mass;}
//
//  //quda::DiracParam diracParam;
//  //inv_param.mass = 0.0;
//  //setDiracParam(diracParam, &inv_param, false);
//  //quda::Dirac *dirac0;
//  //dirac0    = quda::Dirac::create(diracParam);
//
//  vals.resize(2 * nvec);
//
//  quda::ColorSpinorField* src = NULL;
//  quda::Complex           lab = quda::Complex(0.0, 0.0);
//  quda::ColorSpinorField* Czero;
//  quda::ColorSpinorParam cs0 = quda::ColorSpinorParam(cs_gpu);
//  quda::ColorSpinorParam csD = quda::ColorSpinorParam(cs_gpuD);
//  cs0.is_composite  = false;cs0.is_component  = false;
//  csD.is_composite  = false;csD.is_component  = false;
//
//  Czero = quda::ColorSpinorField::Create(cs0);
//  quda::blas::zero(*Czero);
//  src = quda::ColorSpinorField::Create(csD);
//
//  ZSpace.resize(2 * nvec);
//  for(Int n = 0; n < nvec; n++)
//  { 
//    if(n <  off){(*src) = (kSpace[n]      ); lab = evalsK[n].real();}
//    if(n >= off){(*src) = (fSpace[n - off]); lab = evalsF[n-off].real();}
//    lab = std::sqrt(lab.real()/1.0 - 4.0*m*m);
//
//    ////ZSpace.push_back(quda::ColorSpinorField::Create(cs_gpu));
//    ////ZSpace.push_back(quda::ColorSpinorField::Create(cs_gpu));
//    ZSpace[n*2 + 0] = quda::ColorSpinorField(cs0);
//    ZSpace[n*2 + 1] = quda::ColorSpinorField(cs0);
//    quda::blas::zero((ZSpace[n*2 + 0]));
//    quda::blas::zero((ZSpace[n*2 + 1]));
//
//    //quda::blas::caxpy(-1.0*Im*lab, *src, (ZSpace[n*2 + 0]).Even());
//    //quda::blas::caxpy(+1.0*Im*lab, *src, (ZSpace[n*2 + 1]).Even());
//    ///////Doe v_e
//    //dirac->Dslash((ZSpace[n*2 + 0]).Odd(), *src, QUDA_ODD_PARITY);
//    //(ZSpace[n*2 + 1]).Odd() = (ZSpace[n*2 + 0]).Odd();
//    ////quda::blas::ax(-1.0/(std::sqrt(2)*lab.real()), ZSpace[n*2 + 1].Odd());
//    //quda::blas::ax(-1.0/(std::sqrt(2)*lab.real()), ZSpace[n*2 + 0]);
//    //quda::blas::ax(-1.0/(std::sqrt(2)*lab.real()), ZSpace[n*2 + 1]);
//
//    //quda::blas::caxpy(-1.0*Im*lab, *src, (ZSpace[n*2 + 0]).Even());
//    //quda::blas::caxpy(+1.0*Im*lab, *src, (ZSpace[n*2 + 1]).Even());
//    /////Doe v_e
//
//    dirac->Dslash((ZSpace[n*2 + 0]).Even(), *src, QUDA_ODD_PARITY);
//    quda::blas::caxpy(-1.0*Im/(std::sqrt(2)*lab), (ZSpace[n*2 + 0]).Even(), (ZSpace[n*2 + 0]).Odd());
//    quda::blas::caxpy(+1.0*Im/(std::sqrt(2)*lab), (ZSpace[n*2 + 0]).Even(), (ZSpace[n*2 + 1]).Odd());
//
//    (ZSpace[n*2 + 0]).Even() =  (*src);
//    quda::blas::ax(-1.0/(std::sqrt(2)), ZSpace[n*2 + 0].Even());
//
//    (ZSpace[n*2 + 1]).Even() =  (ZSpace[n*2 + 0]).Even();
//
//    ///(ZSpace[n*2 + 1]).Odd() = (ZSpace[n*2 + 0]).Odd();
//    //quda::blas::cax(+1.0*Im/(std::sqrt(2)*lab), (ZSpace[n*2 + 0]).Odd());
//    //quda::blas::cax(-1.0*Im/(std::sqrt(2)*lab), (ZSpace[n*2 + 1]).Odd());
//    //quda::blas::ax(-1.0/(std::sqrt(2)*lab.real()), ZSpace[n*2 + 1].Odd());
//    //quda::blas::ax(-1.0/(std::sqrt(2)*lab.real()), ZSpace[n*2 + 0]);
//    //quda::blas::ax(-1.0/(std::sqrt(2)*lab.real()), ZSpace[n*2 + 1]);
//
//
//    /////m + b |\lambda|
//    vals[n*2 + 0] = 2*m + Im * lab; 
//    vals[n*2 + 1] = 2*m - Im * lab; 
//  }
//
//  delete Czero;
//  delete src;
//
//  //inv_param.mass = m;
//  //delete dirac0;
//
//  check_residualF();
//}

//inline void quda_inverter::check_residualF()
//{
//  TIMER("check_residualF");
//  std::vector<quda::Complex  >& vals = evalsZ;
//  quda::Complex n_unit(-1.0, 0.0);
//  /////qmessage("mass check %.8e \n", mass_mat);
//
//  for(Int n=0;n<nvec*2;n++)
//  {   
//    quda::ColorSpinorField& gs = (*gsrc).Component(0);
//    quda::ColorSpinorField& gr = (*gres).Component(0);
//    quda::ColorSpinorField& g1 = (*gtmp1).Component(0);
//
//    gs = ZSpace[n];
//    gr = ZSpace[n];
//    quda::blas::zero(g1);
//    //(*mat)(*gtmp1, *gres, *gtmp0, *gtmp2);
//    (*mat)( g1, gr );
//
//    quda::Complex srcn = quda::blas::norm2( gs );
//    quda::Complex resn = quda::blas::norm2( g1 );
//    quda::Complex evals  = quda::blas::cDotProduct( gs, g1 ) / quda::blas::norm2(gs);
//    /////quda::Complex factor = sqrt(quda::blas::norm2( g1 )) / sqrt(quda::blas::norm2( gs));
//    quda::blas::caxpby(evals, gs , n_unit, g1 ); ///tmp = ev*src + n * gtmp
//    quda::Complex residual = sqrt(quda::blas::norm2( g1 ) / srcn);
//
//    quda::blas::zero( g1 );
//    //(*mat)(*gtmp1, *gres, *gtmp0, *gtmp2);
//    (*mat)( g1, gr );
//    quda::Complex even_v  = quda::blas::cDotProduct(gs.Even(), g1.Even()) / quda::blas::norm2(gs.Even());
//    quda::blas::caxpby(even_v, gs.Even() , n_unit, g1.Even() );
//    quda::Complex res_e = sqrt(quda::blas::norm2(g1.Even())/ quda::blas::norm2(gs.Even()));
//
//    quda::blas::zero(g1);
//    //(*mat)(*gtmp1, *gres, *gtmp0, *gtmp2);
//    (*mat)( g1, gr );
//    quda::Complex odd_v  = quda::blas::cDotProduct(gs.Odd(), g1.Odd()) / quda::blas::norm2(gs.Odd());
//    quda::blas::caxpby(odd_v, gs.Odd() , n_unit, g1.Odd() );
//    quda::Complex res_o = sqrt(quda::blas::norm2(g1.Odd())/ quda::blas::norm2(gs.Odd()));
//
//    if(quda::comm_rank_global()== 0){
//      //printf("===vec %5d, v %+.1e %+.1e, e %+.1e %+.1e, residual %.3e, %.3e %.3e, %+.1e %+.1e %+.1e %+.1e \n", n, 
//      //    residual.real(), vals[n].real(), vals[n].imag(), evals.real(), evals.imag(), res_e.real(), res_o.real(),
//      //    even_v.real(), even_v.imag(), odd_v.real(), odd_v.imag() );
//      printf("===vec %5d, norm %+.2e v %+.1e %+.1e, e %+.1e %+.1e, residual %+.3e, %+.3e %+.3e \n", n, 1 - std::sqrt(srcn.real()),
//          vals[n].real(), vals[n].imag(), evals.real()-vals[n].real(), evals.imag()-vals[n].imag(),
//          residual.real(), res_e.real(), res_o.real() );
//    }
//  }
//}

inline void quda_inverter::prepare_low_prop()
{
  if(inv_param.cuda_prec == QUDA_DOUBLE_PRECISION){
    qlat::vector<qlat::ComplexT<double>* > src;
    qlat::vector<qlat::ComplexT<double>* > res;
    qlat::vector<qlat::ComplexT<double>* > buf;
    src.resize(1);
    res.resize(1);
    buf.resize(1);
    src[0] = (qlat::ComplexT<double>*) gsrc->data();
    res[0] = (qlat::ComplexT<double>*) gres->data();
    buf[0] = (qlat::ComplexT<double>*) gtmp1->data();
    get_low_prop(res, src, buf, 0);
  }
  if(inv_param.cuda_prec == QUDA_SINGLE_PRECISION){
    qlat::vector<qlat::ComplexT<float>* > src;
    qlat::vector<qlat::ComplexT<float>* > res;
    qlat::vector<qlat::ComplexT<float>* > buf;
    src.resize(1);
    res.resize(1);
    buf.resize(1);
    src[0] = (qlat::ComplexT<float>*) gsrc->data();
    res[0] = (qlat::ComplexT<float>*) gres->data();
    buf[0] = (qlat::ComplexT<float>*) gtmp1->data();
    get_low_prop(res, src, buf, 0);
  }
}

////from gsrc to gtmp2
////TODO move dirac0 to buffers
//inline void quda_inverter::prepare_low_prop(Int mode )
//{
//  TIMER("prepare_low_prop");
//
//  /////std::string val = get_env(std::string("q_low_mode"));
//  /////if(val != ""){mode = stringtonum(val);}
//  
//  if(use_eigen_pc == 0 or nvec <= 0){for(Int di=0;di<num_src_inv;di++){quda::blas::zero((*gres).Component(di));}return;}
//  double m = inv_param.mass;
//  update_eigen_mass(m, true);
//
//  quda::DiracParam diracParam;
//  inv_param.mass = 0.0;
//  setDiracParam(diracParam, &inv_param, false);
//  quda::Dirac *dirac0;
//  dirac0    = quda::Dirac::create(diracParam);
//  ////qmessage("mass %.8e \n", m);
//
//  quda::Complex n_unit(1.0, 0.0);
//  if(mode == 0){
//    ////mode 0 write
//    Int n_defl = 0;
//    //*gtmp1 = *gsrc;
//    /////Doe v_e
//    for(Int di=0;di<num_src_inv;di++)
//    {
//      quda::blas::zero((*gres).Component(di));
//      quda::blas::zero((*gtmp1).Component(di));
//      dirac0->Dslash((*gtmp1).Component(di).Odd(), (*gsrc).Component(di).Odd(), QUDA_EVEN_PARITY);
//      quda::blas::ax(-1.0, (*gtmp1).Component(di).Odd());
//      (*gtmp1).Component(di).Even() = (*gsrc).Component(di).Even();
//
//      if(gtmp1D == NULL){gtmp1D = new quda::ColorSpinorField(cs_gpuD);}
//      if(gtmp2D == NULL){gtmp2D = new quda::ColorSpinorField(cs_gpuD);}
//      if(gtmp1F == NULL){gtmp1F = new quda::ColorSpinorField(cs_gpuF);}
//      if(gtmp2F == NULL){gtmp2F = new quda::ColorSpinorField(cs_gpuF);}
//      //if(gtmp2D == NULL){gtmp2D = new quda::ColorSpinorField(cs_gpuD);}
//      //if(gtmp2F == NULL){gtmp2F = new quda::ColorSpinorField(cs_gpuF);}
//
//      for(Int kf=0;kf<2;kf++)
//      {
//        std::vector<quda::ColorSpinorField > eig_vecs;
//        if(kf==0){n_defl = kSpace.size();eig_vecs.resize(n_defl);for(Int i=0;i<n_defl;i++){eig_vecs[i]=kSpace[i].create_alias();}}
//        if(kf==1){n_defl = fSpace.size();eig_vecs.resize(n_defl);for(Int i=0;i<n_defl;i++){eig_vecs[i]=fSpace[i].create_alias();}}
//
//        if(n_defl == 0){continue ;}
//        std::vector<quda::Complex> A(n_defl);
//        std::vector<quda::Complex> B(n_defl);
//
//        // 1. Take block inner product: (V_i)^dag * vec = A_i
//        std::vector<quda::ColorSpinorField *> buf;
//        std::vector<quda::ColorSpinorField> buf_;buf_.resize(1);
//        ////if(kf == 0){src_.push_back(&(gtmp1->Even()));}
//        if(kf == 0){buf.push_back(gtmp1D);}
//        if(kf == 1){buf.push_back(gtmp1F);}
//        buf_[0] = (*buf[0]).create_alias();
//        *buf[0] = gtmp1->Component(di).Even();
//
//        quda::blas::block::cDotProduct(A, eig_vecs, buf_);
//        ////src_[0] = &(gtmp1->Odd());
//        *buf[0] = gtmp1->Component(di).Odd();
//        quda::blas::block::cDotProduct(B, eig_vecs, buf_);
//        quda::Complex ev;
//        for (Int i = 0; i < n_defl; i++) {
//          if(kf==0){ev = evalsK[i];}
//          if(kf==1){ev = evalsF[i];}
//          quda::Complex l = std::sqrt(ev.real()/1.0 - 4.0*m*m);
//          //quda::Complex l = std::sqrt(evals_ZERO[i].real());
//          quda::Complex ai = A[i];
//          quda::Complex bi = B[i];
//          A[i] = ( 2*m * ai - bi ) / (4*m*m + l*l);
//          //B[i] = (-1.0 * ai - (2.0 * m * bi)/(l*l) ) / (4.0*m*m + l*l);
//          B[i] = (-1.0 * ai - (2.0 * m * bi)/(l*l) ) / (4.0*m*m + l*l);
//          /////-1.0 * ai/(4.0*m*m + l*l) - 2*m * bi/(l*l * (4*m*m + l*l));
//        }
//
//        ////std::vector<quda::ColorSpinorField *> sol;
//        //sol.push_back(&(gres->Even()));
//        quda::blas::zero(*buf[0]);
//        quda::blas::block::caxpy(A, eig_vecs, buf_);
//        /////gres->Even() = gres->Even() + *buf[0];
//        gtmp2->Component(di).Even() = *buf[0];
//
//        quda::blas::caxpby(n_unit, gtmp2->Component(di).Even() , n_unit, gres->Component(di).Even());
//        ////sol[0] = &(gres->Odd());
//
//        quda::blas::zero(*buf[0]);
//        quda::blas::block::caxpy(B, eig_vecs, buf_);
//        gtmp2->Component(di).Odd()  = *buf[0];
//        quda::blas::caxpby(n_unit, gtmp2->Component(di).Odd()  , n_unit, gres->Component(di).Odd());
//        ////gres->Odd()  = gres->Odd() + *buf[0];
//      }
//
//      quda::blas::zero((*gtmp1).Component(di).Odd());
//      dirac0->Dslash((*gtmp1).Component(di).Odd(), (*gres).Component(di).Odd(), QUDA_ODD_PARITY);
//      //dirac->Dslash((*gtmp1).Component(di).Odd(), (*gres).Component(di).Odd(), QUDA_EVEN_PARITY);
//      quda::blas::ax(-1.0, (*gtmp1).Component(di).Odd());
//      (*gres).Component(di).Odd() = (*gtmp1).Component(di).Odd();
//    }
//
//    ////////===check sections uncomment for test 
//    //quda::Complex s0_e = quda::blas::norm2(gsrc->Even());
//    //quda::Complex s0_o = quda::blas::norm2(gsrc->Odd());
//    //quda::Complex n0_e = quda::blas::norm2(gtmp2->Even());
//    //quda::Complex n0_o = quda::blas::norm2(gtmp2->Odd());
//    //quda::Complex n1_e = quda::blas::norm2(gres->Even());
//    //quda::Complex n1_o = quda::blas::norm2(gres->Odd());
//
//    //quda::Complex m_unit( 1.0, 0.0);
//    //quda::Complex n_unit(-1.0, 0.0);
//    //(*gtmp0) =(*gtmp2);
//    //quda::blas::caxpby(m_unit, (*gres).Even() , n_unit, (*gtmp0).Even());
//    //quda::Complex res_e = quda::blas::norm2((*gtmp0).Even());
//    //(*gtmp0) =(*gtmp2);
//    //quda::blas::caxpby(m_unit, (*gres).Odd() , n_unit, (*gtmp0).Odd());
//    //quda::Complex res_o = quda::blas::norm2((*gtmp0).Odd() );
//
//    //if(quda::comm_rank_global()== 0)printf("===Even %.3e, Odd %.3e , s %.1e %.1e , n %.1e %.1e %.1e %.1e \n", res_e.real(), res_o.real(), 
//    //    s0_e.real(), s0_o.real(), n0_e.real(), n0_o.real(), n1_e.real(), n1_o.real());
//  }
//
//  /////eigensystem reconstructed
//  if(mode == 1){
//    if(ZSpace.size() != nvec*2){Qassert(false);}
//    ////mode 1 write
//    ///quda::ColorSpinorField& res = *gtmp2;
//    for(Int di=0;di<num_src_inv;di++)
//    {
//    quda::ColorSpinorField& gs = (*gsrc).Component(di);
//    quda::ColorSpinorField& g2 = (*gtmp2).Component(di);
//    quda::blas::zero(g2);
//
//    // Perform Sum_i V_i * (L_i)^{-1} * (V_i)^dag * vec = vec_defl
//    Int n_defl = ZSpace.size();
//    //std::vector<quda::ColorSpinorField *> eig_vecs;
//    //eig_vecs.reserve(n_defl);
//    //for (Int i = 0; i < n_defl; i++) eig_vecs.push_back(&ZSpace[i]);
//
//    // 1. Take block inner product: (V_i)^dag * vec = A_i
//    std::vector<quda::Complex> s(n_defl * 1);
//    std::vector<quda::ColorSpinorField > src;src.resize(1);
//    src[0] = gs.create_alias();
//    //std::vector<quda::ColorSpinorField *> src_;src_.push_back(&gs);
//    //quda::cvector_ref<quda::ColorSpinorField > SRC_ = {src_.begin(), src_.begin() + 1};
//    //SRC_[0] = *src_[0];
//    ///= const_cast<decltype(src) &>(src);
//    quda::blas::block::cDotProduct(s, ZSpace, src);
//
//    // 2. Perform block caxpy: V_i * (L_i)^{-1} * A_i
//    for (Int i = 0; i < n_defl; i++) { s[i] /= evalsZ[i]; }
//
//    // 3. Accumulate sum vec_defl = Sum_i V_i * (L_i)^{-1} * A_i
//    //if (!accumulate)
//    //  for (auto &x : sol) blas::zero(*x);
//    std::vector<quda::ColorSpinorField > sol;sol.resize(1);
//    sol[0] = g2.create_alias();
//    //std::vector<quda::ColorSpinorField *> sol;sol.push_back(&g2);
//    //quda::cvector_ref<quda::ColorSpinorField > SOL_(sol.begin(), sol.begin()+1);
//    quda::blas::block::caxpy(s, ZSpace, sol);
//    (*gres).Component(di) = g2;
//
//    }
//  }
//
//
//  inv_param.mass = m;
//  delete dirac0;
//  quda::saveTuneCache();
//}

//inline void quda_inverter::deflate(quda::ColorSpinorField &sol, const quda::ColorSpinorField &src, std::vector<quda::ColorSpinorField> &evecs, const std::vector<quda::Complex> &evals, bool accumulate)
//{
//  // FIXME add support for mixed-precison dot product to avoid this copy
//  //if (src.Precision() != evecs[0]->Precision() && !tmp1) {
//  //  quda::ColorSpinorParam param(*evecs[0]);
//  //  tmp1 = new quda::ColorSpinorField(param);
//  //}
//  //if (sol.Precision() != evecs[0]->Precision() && !tmp2) {
//  //  quda::ColorSpinorParam param(*evecs[0]);
//  //  tmp2 = new quda::ColorSpinorField(param);
//  //}
//  quda::ColorSpinorField *tmp1, *tmp2;
//  if( src.Precision() != evecs[0].Precision()){
//    quda::ColorSpinorParam param(evecs[0]);
//    if(evecs[0].Precision() == QUDA_DOUBLE_PRECISION){
//      if(gtmp1D == NULL){gtmp1D = new quda::ColorSpinorField(param);}
//      tmp1 = gtmp1D;
//    }
//    if(evecs[0].Precision() == QUDA_SINGLE_PRECISION){
//      if(gtmp1F == NULL){gtmp1F = new quda::ColorSpinorField(param);}
//      tmp1 = gtmp1F;
//    }
//  }
//  if( sol.Precision() != evecs[0].Precision()){
//    quda::ColorSpinorParam param(evecs[0]);
//    if(evecs[0].Precision() == QUDA_DOUBLE_PRECISION){
//      if(gtmp2D == NULL){gtmp2D = new quda::ColorSpinorField(param);}
//      tmp2 = gtmp2D;
//    }
//    if(evecs[0].Precision() == QUDA_SINGLE_PRECISION){
//      if(gtmp2F == NULL){gtmp2F = new quda::ColorSpinorField(param);}
//      tmp2 = gtmp2F;
//    }
//  }
//  quda::ColorSpinorField *src_tmp = src.Precision() != evecs[0].Precision() ? tmp1 : const_cast<quda::ColorSpinorField *>(&src);
//  quda::ColorSpinorField *sol_tmp = sol.Precision() != evecs[0].Precision() ? tmp2 : const_cast<quda::ColorSpinorField *>(&sol);
//  quda::blas::copy(*src_tmp, src); // no-op if these alias
//
//  std::vector<quda::ColorSpinorField *> src_ {src_tmp};
//  std::vector<quda::ColorSpinorField *> sol_ {sol_tmp};
//  //quda::cvector_ref<quda::ColorSpinorField> src_;src_.resize(1);
//  //quda::cvector_ref<quda::ColorSpinorField > sol_;sol_.resize(1);
//  //src_[0] = *src_tmp;
//  //sol_[0] = *sol_tmp;
//  //std::vector<std::reference_wrapper<quda::ColorSpinorField >> src_;
//  //std::vector<std::reference_wrapper<quda::ColorSpinorField >> sol_;
//  //src_.resize(1);src_[0] = *src_tmp;
//  //sol_.resize(1);sol_[0] = *sol_tmp;
//  deflate(sol_, src_, evecs, evals, accumulate);
//  quda::blas::copy(sol, *sol_tmp); // no-op if these alias
//}

/*
  If double inputs, do it directly
  If single inputs, transform the data orderings with QUDA and then transform back 
*/
inline void quda_inverter::deflate_Ty(qlat::vector<void* >& Pres, qlat::vector<void* >& Psrc, double mass, Int buf_prec , Int mode, Int clear, const bool data_from_quda)
{
  //Qassert(get_data_type<Ty>() == ComplexD_TYPE or get_data_type<Ty>() == ComplexF_TYPE);
  Qassert(buf_prec == 0 or buf_prec == 1);
  Qassert(eigenL.size() == eigen_precL.size());
  const Long Ndata = 3 * geoB().local_volume() / 2; // half volume vectors
  if(eigen_with_nvec == false){
    for(unsigned int vi=0;vi<Pres.size();vi++){
      if(buf_prec == 0){zero_Ty((qlat::ComplexT<double>*) Pres[vi], Ndata, 1, QFALSE);}
      if(buf_prec == 1){zero_Ty((qlat::ComplexT< float>*) Pres[vi], Ndata, 1, QFALSE);}
    }
    qacc_barrier(dummy);
    return ;
  }
  qlat::ComplexT<double> rD;
  qlat::ComplexT< float> rF;
  std::vector<Int> clearL;clearL.resize(eigenL.size());
  Int found_eig = 0;
  for(unsigned int i=0;i<eigenL.size();i++)
  {
    clearL[i] = 0;
    if(eigen_with_nvecL[i] != 0 and found_eig == 0 and clear == 1){
      clearL[i] = 1;
      found_eig = 1;
    }
  }

  //std::vector<void* > eL;eL.resize(2);
  //std::vector<Int   > pL;pL.resize(2);
  //eL[0] = eigen0;eL[1] = eigen1;
  //pL[0] = eigen_prec0;pL[1] = eigen_prec1;

  //// check norm
  //{
  //  for(unsigned int vi=0;vi<Pres.size();vi++){
  //    double d0 = 0.0;
  //    double d1 = 0.0;
  //    if(buf_prec == 0){d0 = norm_vec((qlat::ComplexT<double>*) Psrc[vi], Ndata);}
  //    if(buf_prec == 1){d0 = norm_vec((qlat::ComplexT<float >*) Psrc[vi], Ndata);}
  //    if(buf_prec == 0){d1 = norm_vec((qlat::ComplexT<double>*) Pres[vi], Ndata);}
  //    if(buf_prec == 1){d1 = norm_vec((qlat::ComplexT<float >*) Pres[vi], Ndata);}
  //    qmessage("==vec eigens %5d norm %.8e %.8e \n", int(vi), d0, d1);
  //  }
  //}

  if(data_from_quda){
    quda_to_double_order(Psrc, buf_prec);
    if(Psrc[0] != Pres[0] and clear == 0){
      quda_to_double_order(Pres, buf_prec);
    }
  }
  for(unsigned int i=0;i<eigenL.size();i++)
  {
    void* Ea        = eigenL[i];
    Int  eigen_prec = eigen_precL[i];
    if(Ea != NULL){
      if(eigen_prec == 0){
        eigen_cs<qlat::ComplexT<double  >>* eigen = (eigen_cs<qlat::ComplexT<double  >>*) Ea;
        if(buf_prec == 0){eigen->deflate(Pres, Psrc, rD, mass, mode, true, 0, -1, 1, clearL[i]);}
        if(buf_prec == 1){eigen->deflate(Pres, Psrc, rF, mass, mode, true, 0, -1, 1, clearL[i]);}
      }

      if(eigen_prec == 1){
        eigen_cs<qlat::ComplexT<float  > >* eigen = (eigen_cs<qlat::ComplexT<float  > >*) Ea;
        if(buf_prec == 0){eigen->deflate(Pres, Psrc, rD, mass, mode, true, 0, -1, 1, clearL[i]);}
        if(buf_prec == 1){eigen->deflate(Pres, Psrc, rF, mass, mode, true, 0, -1, 1, clearL[i]);}
      }
    }
  }
  if(data_from_quda){
    quda_to_single_order(Psrc, buf_prec);
    // sometime the pointer is the same
    if(Psrc[0] != Pres[0]){
      quda_to_single_order(Pres, buf_prec);
    }
  }

  //// check norm
  //{
  //  for(unsigned int vi=0;vi<Pres.size();vi++){
  //    double d0 = 0.0;
  //    double d1 = 0.0;
  //    if(buf_prec == 0){d0 = norm_vec((qlat::ComplexT<double>*) Psrc[vi], Ndata);}
  //    if(buf_prec == 1){d0 = norm_vec((qlat::ComplexT<float >*) Psrc[vi], Ndata);}
  //    if(buf_prec == 0){d1 = norm_vec((qlat::ComplexT<double>*) Pres[vi], Ndata);}
  //    if(buf_prec == 1){d1 = norm_vec((qlat::ComplexT<float >*) Pres[vi], Ndata);}
  //    qmessage("==vec eigenr %5d norm %.8e %.8e \n", int(vi), d0, d1);
  //  }
  //}
}

/* 
  copy src to buf with rotation or not
  change buf with dslash
  eigen low mode
  change res with dslash
  rotate res if needed
*/
template<typename Ty, typename Tk>
inline void quda_inverter::get_low_prop(qlat::vector<Ty* >& res, qlat::vector<Ty* >& src, qlat::vector<Tk* >& buf, Int qlat_format)
{
  // double input, single buf
  if(get_data_type_is_double<Ty >() == 1 and get_data_type_is_double<Tk >() == 0){
    Qassert(false);
  }

  //// odd dslash
  quda::ColorSpinorParam param   = cs_gpu;
  get_param_ref(param, get_data_type_is_double<Ty >());
  const Int buf_prec =  get_data_type_is_double<Ty >() ? 0 : 1;

  const Int nsrc = src.size();
  Qassert(nsrc == res.size() and nsrc == buf.size());
  qlat::vector<Ty* > cbuf;cbuf.resize(nsrc);
  for(Int iv=0;iv<nsrc;iv++){
    cbuf[iv] = (Ty*) buf[iv];
  }

  const Long V = geoB().local_volume();
  //if no eigen vectors return 0 vectors with each size V*3
  if(eigen_with_nvec == false){
    for(Int vi=0;vi<nsrc;vi++){
      zero_Ty(res[vi], V*3, 1, QFALSE);
    }
    qacc_barrier(dummy);
    return ;
  }

  // const Int DIM = 3;
  Ty** pB =  res.data();
  Ty** pR = cbuf.data();
  if(qlat_format == 0){
    pB = cbuf.data();
    pR =  res.data();
  }

  // src copy to buf
  for(Int vi=0;vi<nsrc;vi++){
    if(qlat_format == 1){
      qlat_cf_to_quda_cf_P(cbuf[vi], src[vi], map_index, param);
      //{
      //  param.v = (void*) cbuf[vi];
      //  quda::ColorSpinorField* Qvec = NULL;
      //  Qvec = new quda::ColorSpinorField(param);
      //  qlat_cf_to_quda_cf(*Qvec, src[vi], geoB(), map_index);
      //  //qlat_cf_to_quda_cf(cbuf[vi], src[vi], DIM, geoB(), map_index);
      //  delete Qvec;
      //}
      //
      if(pB[vi]!=cbuf[vi]){
        cpy_GPU(pB[vi], cbuf[vi], V*3);
      }
    }else{
      cpy_GPU(pB[vi], src[vi], V*3);
    }
  }


  //if(get_data_type_is_double<Ty >()){
  //param.is_composite  = false;
  //param.is_component  = false;
  ////data double to quda double
  //if(get_data_type_is_double<Ty >()){
  //  param.setPrecision(QUDA_DOUBLE_PRECISION, QUDA_DOUBLE_PRECISION, true);
  //  buf_prec = 0;
  //}else{
  //  param.setPrecision(QUDA_SINGLE_PRECISION, QUDA_SINGLE_PRECISION, true);
  //  buf_prec = 1;
  //}
  //param.create = QUDA_REFERENCE_FIELD_CREATE;

  //if(inv_param.cuda_prec == QUDA_DOUBLE_PRECISION){buf_prec = 0;}
  //if(inv_param.cuda_prec == QUDA_SINGLE_PRECISION){buf_prec = 1;}

  quda::ColorSpinorField& gbuf = (*gtmp0);

  std::vector<quda::ColorSpinorField* > Qvec;Qvec.resize(nsrc * 2);
  // odd dlash and replaced
  for(Int vi=0;vi<nsrc;vi++){
    param.v = (void*) pB[vi];
    Qvec[vi] = new quda::ColorSpinorField(param);
    gbuf.Odd() = Qvec[vi]->Odd();
    dirac->Dslash(gbuf.Even(), gbuf.Odd(), QUDA_EVEN_PARITY);
    Qvec[vi]->Odd() = gbuf.Even();

    param.v = (void*) pR[vi];
    Qvec[nsrc + vi] = new quda::ColorSpinorField(param);
  }

  // copy pointers
  qlat::vector<void* > Psrc;
  qlat::vector<void* > Pres;
  Psrc.resize(2*nsrc);
  Pres.resize(2*nsrc);
  for(Int vi=0;vi<nsrc;vi++)
  {
    Psrc[vi*2 + 0] = Qvec[vi]->Even().data();
    Psrc[vi*2 + 1] = Qvec[vi]->Odd().data();
    Pres[vi*2 + 0] = Qvec[nsrc + vi]->Even().data();
    Pres[vi*2 + 1] = Qvec[nsrc + vi]->Odd().data();
  }

  // actual low mode jobs
  deflate_Ty(Pres, Psrc, mass_mat, buf_prec, 1);

  //dslash
  for(Int vi=0;vi<nsrc;vi++){
    gbuf.Odd()  = Qvec[nsrc + vi]->Odd();
    dirac->Dslash(gbuf.Even(), gbuf.Odd() , QUDA_ODD_PARITY);
    Qvec[nsrc + vi]->Odd() = gbuf.Even();
  }

  // buf copy to res with rotations
  if(qlat_format == 1)
  for(Int vi=0;vi<nsrc;vi++){
    quda_cf_to_qlat_cf_cf_P(pB[vi], pR[vi], map_index, param);
    //quda_cf_to_qlat_cf(pB[vi], pR[vi], DIM, geoB(), map_index);
  }

  for(Int vi=0;vi<2*nsrc;vi++){
    delete Qvec[vi];
  }
}

//inline void quda_inverter::deflate(std::vector<quda::ColorSpinorField* > &sol, const std::vector<quda::ColorSpinorField* > &src,
//                            std::vector<quda::ColorSpinorField> &evecs, const std::vector<quda::Complex> &evals, bool accumulate)
//{
//  TIMER("quda deflate");
//  Int n_defl = evecs.size();
//  if( n_defl == 0){ return; }
//
//  if (getVerbosity() >= QUDA_VERBOSE){printfQuda("Deflating %d vectors\n", n_defl);}
//
//  return ;
//
//  std::vector<quda::ColorSpinorField > src_;
//  std::vector<quda::ColorSpinorField > sol_;
//  const Int Nsrc = src.size();
//  Qassert(src.size() == sol.size());
//  src_.resize(Nsrc);
//  sol_.resize(Nsrc);
//  for(Int i=0;i<Nsrc;i++)
//  {
//    src_[i] = (*src[i]).create_alias();
//    sol_[i] = (*sol[i]).create_alias();
//  }
//
//  // Perform Sum_i V_i * (L_i)^{-1} * (V_i)^dag * vec = vec_defl
//  // for all i computed eigenvectors and values.
//
//  // Pointers to the required Krylov space vectors,
//  // no extra memory is allocated.
//  //std::vector<quda::ColorSpinorField *> eig_vecs;
//  //eig_vecs.reserve(n_defl);
//  //for (Int i = 0; i < n_defl; i++) eig_vecs.push_back(&evecs[i]);
//
//  // 1. Take block inner product: (V_i)^dag * vec = A_i
//  std::vector<quda::Complex> s(n_defl * src_.size());
//  //std::vector<quda::ColorSpinorField *> src_ = const_cast<decltype(src) &>(src);
//  //quda::blas::cDotProduct(s, evecs, src);
//  quda::blas::block::cDotProduct(s, evecs, src_);
//
//  // 2. Perform block caxpy: V_i * (L_i)^{-1} * A_i
//  for (Int i = 0; i < n_defl; i++) { s[i] /= evals[i].real(); }
//
//  // 3. Accumulate sum vec_defl = Sum_i V_i * (L_i)^{-1} * A_i
//  if (!accumulate)
//    for (auto &x : sol_) quda::blas::zero(x);
//  //quda::blas::caxpy(s, eig_vecs, sol);
//  //quda::blas::caxpy(s, evecs, sol);
//  quda::blas::block::caxpy(s, evecs, sol_);
//
//  // Save Deflation tuning
//  quda::saveTuneCache();
//}


//inline void quda_inverter::setup_CG()
//{
//  //solverParam.maxiter = inv_param.maxiter;
//  //solverParam.tol = inv_param.tol;
//
//  //if(CG_reset == false){return;}
//  //{
//  //TIMERA("setup_CG");
//  ////clear_CG();
//  //QudaInvertParam& param = inv_param;
//  //// It was probably a bad design decision to encode whether the system is even/odd preconditioned (PC) in
//  //// solve_type and solution_type, rather than in separate members of QudaInvertParam.  We're stuck with it
//  //// for now, though, so here we factorize everything for convenience.
//  ////bool pc_solution = (param.solution_type == QUDA_MATPC_SOLUTION) ||
//  ////  (param.solution_type == QUDA_MATPCDAG_MATPC_SOLUTION);
//  ////bool pc_solve = (param.solve_type == QUDA_DIRECT_PC_SOLVE) ||
//  ////  (param.solve_type == QUDA_NORMOP_PC_SOLVE) || (param.solve_type == QUDA_NORMERR_PC_SOLVE);
//  ////bool direct_solve = (param.solve_type == QUDA_DIRECT_SOLVE) ||
//  ////  (param.solve_type == QUDA_DIRECT_PC_SOLVE);
//
//  //solverParam = quda::SolverParam(param);
//  //// solve_cg = quda::Solver::create(solverParam, *m_cg, *mSloppy, *mPre, *mEig);
//  //}
//  //if (getVerbosity() >= QUDA_VERBOSE) { printfQuda("Solution = %g\n", blas::norm2(x)); }
//  //CG_reset = false;
//}

inline void quda_inverter::callMultiSrcQuda(qlat::vector<void* >& res, qlat::vector<void* >& src, Int max_src)
{
  TIMER_FLOPS("invertQuda Multi");
  Qassert(src.size() == res.size());
  Int max_src_inv = -1;
  if(max_src_qinv != -1){
    max_src_inv = max_src_qinv;
  }
  if(max_src != -1){
    max_src_inv = max_src;
  }
  const Int num_src = src.size();
  if(num_src <= max_src_inv){
    inv_param.num_src = num_src;
    invertMultiSrcQuda(res.data(), src.data(), &inv_param);
  }else{
    std::vector<Long > jobA = job_create(num_src, max_src_inv);
    Int njobs = jobA.size() / 2;
    Qassert(max_src_inv >= 4 and njobs >= 2);
    if(jobA[(njobs - 1)*2 + 1] == 1 )
    {
      jobA[(njobs - 1)*2 + 1 ] += 2;
      jobA[(njobs - 2)*2 + 1 ] -= 2;

      jobA[(njobs - 1)*2 + 0 ] -= 2;
    }

    std::vector<void* > res_new;
    std::vector<void* > src_new;
    for(LInt jobi=0;jobi < jobA.size()/2; jobi++)
    {
      Long bini = jobA[jobi*2 + 0];
      Long bcut = jobA[jobi*2 + 1];
      src_new.resize(bcut);
      res_new.resize(bcut);
      for(Long i=0;i<bcut;i++)
      {
        src_new[i] = src[bini + i];
        res_new[i] = res[bini + i];
      }
      //qmessage("ini %d, cut %d, tot %d \n", int(bini), int(bcut), int(num_src));
      inv_param.num_src = bcut;
      invertMultiSrcQuda(res_new.data(), src_new.data(), &inv_param);
    }
  }

  timer.flops += inv_param.gflops * 1024 * 1024 * 1024 / qlat::get_num_node();
}

inline void quda_inverter::invertQuda_COPY_single(quda::ColorSpinorField& res, quda::ColorSpinorField& src)
{
  TIMERA("Quda invertQuda_COPY_single");
  Qassert(inv_param.solution_type == QUDA_MATPC_SOLUTION);
  inv_param.cpu_prec = src.Precision();
  //Qassert(src.Precision() == inv_param.cpu_prec);
  qlat::vector<void* > sI;
  qlat::vector<void* > rI;
  sI.resize(1);
  rI.resize(1);
  sI[0] = src.data();
  rI[0] = res.data();
  
  callMultiSrcQuda(rI, sI);
  //invertMultiSrcQuda(rI.data(), sI.data(), &inv_param);
  //if(src.Precision() == QUDA_DOUBLE_PRECISION)
  //{
  //}

  //extra copy needed
  //if(src.Precision() == QUDA_SINGLE_PRECISION)
  //{
  //  (*gsrcH).Component(0) = src;
  //  invertQuda((*gresH).Component(0).data(), (*gsrcH).Component(0).data(), &inv_param);
  //  res = (*gresH).Component(0);
  //}

  ////solverParam.secs = 0;
  ////solverParam.gflops = 0;
  ////solverParam.iter = 0;
  //solve_mode = solve_mode_;
  //QudaInvertParam& param = inv_param;
  ////param.secs = 0;
  ////param.gflops = 0;
  //param.iter = 0;
  //solverParam.iter = 0;
  //auto profile = quda::pushProfile(profileInvertC, param.secs, param.gflops);
  ////quda::profilerStart(__func__);
  ////quda::enable_profiler = true;
  //quda::device::profile::start();

  //bool pc_solution = (param.solution_type == QUDA_MATPC_SOLUTION) ||
  //  (param.solution_type == QUDA_MATPCDAG_MATPC_SOLUTION);
  //bool pc_solve = (param.solve_type == QUDA_DIRECT_PC_SOLVE) ||
  //  (param.solve_type == QUDA_NORMOP_PC_SOLVE) || (param.solve_type == QUDA_NORMERR_PC_SOLVE);
  //bool direct_solve = (param.solve_type == QUDA_DIRECT_SOLVE) ||
  //  (param.solve_type == QUDA_DIRECT_PC_SOLVE);

  //double nb = 1.0;
  //{
  //quda::ColorSpinorField& x = res;
  //quda::ColorSpinorField& b = src;
  ////if(use_eigen_pc == 0){quda::blas::zero(x);}
  //if(param.use_init_guess != QUDA_USE_INIT_GUESS_YES){quda::blas::zero(x);}
  //nb = quda::blas::norm2(b);
  //if(nb == 0.0){errorQuda("Source has zero norm");}

  //if (getVerbosity() >= QUDA_VERBOSE) {
  //  printfQuda("Source: %g\n", nb);
  //  if (param.use_init_guess == QUDA_USE_INIT_GUESS_YES) { printfQuda("Initial guess: %g\n", quda::blas::norm2(x)); }
  //}

  //// rescale the source and solution vectors to help prevent the onset of underflow
  //if (param.solver_normalization == QUDA_SOURCE_NORMALIZATION) {
  //  quda::blas::ax(1.0 / sqrt(nb), b);
  //  quda::blas::ax(1.0 / sqrt(nb), x);
  //}

  //if (getVerbosity() >= QUDA_VERBOSE) {
  //  double nin  = quda::blas::norm2(b);
  //  double nout = quda::blas::norm2(x);
  //  if(qlat::get_id_node() == 0){
  //  printfQuda("Prepared source = %g\n", nin);
  //  printfQuda("Prepared solution = %g\n", nout);
  //  }
  //}

  //quda::massRescale(b, param, false);

  //if (getVerbosity() >= QUDA_VERBOSE) {
  //  double nin = quda::blas::norm2(b);
  //  if(qlat::get_id_node() == 0){
  //  printfQuda("Prepared source post mass rescale = %g\n", nin);
  //  }
  //}
  //if (pc_solution && !pc_solve) {
  //  errorQuda("Preconditioned (PC) solution_type requires a PC solve_type");
  //}
  //}

  //if (direct_solve) {
  //  if(solve_mode == 0){
  //    solverParam.num_src = 1;
  //    quda::Solver* solve_cg = quda::Solver::create(solverParam, *m_cg, *mSloppy, *mPre, *mEig);
  //    (*solve_cg)(res, src);
  //    delete solve_cg;
  //  }else{errorQuda("Solver not added");}
  //}else{errorQuda("Solver not added");}

  //{
  //quda::ColorSpinorField& x = res;
  //////quda::ColorSpinorField& b = src;

  //if (getVerbosity() >= QUDA_VERBOSE) { printfQuda("Solution = %g\n", quda::blas::norm2(x)); }
  //////dirac.reconstruct(x, b, param.solution_type);

  //if (param.solver_normalization == QUDA_SOURCE_NORMALIZATION) {
  //  // rescale the solution
  //  quda::blas::ax(sqrt(nb), x);
  //  //quda::blas::ax(sqrt(nb), b);
  //}

  //if (getVerbosity() >= QUDA_VERBOSE) {
  //  printfQuda("Reconstructed solution: %g\n", quda::blas::norm2(x));
  //}
  //}

  //solverParam.updateInvertParam(param);

  ////inv_param.secs = solverParam.secs;
  ////inv_param.gflops = solverParam.gflops;
  ////quda::comm_allreduce_sum(param.gflops);
  ////inv_param.iter = solverParam.iter;
  //// cache is written out even if a Long benchmarking job gets interrupted
  //quda::saveTuneCache();
  ////quda::profilerStop(__func__);
  //quda::device::profile::stop();
  ////quda::enable_profiler = false;
}

//inline void quda_inverter::invertQuda_COPY(quda::ColorSpinorField& res, quda::ColorSpinorField& src, Int solve_mode_)
//{
//  TIMERA("Quda invertQuda_COPY");
//  //solverParam.secs = 0;
//  //solverParam.gflops = 0;
//  //solverParam.iter = 0;
//  solve_mode = solve_mode_;
//  QudaInvertParam& param = inv_param;
//  //param.secs = 0;
//  //param.gflops = 0;
//  param.iter = 0;
//  solverParam.iter = 0;
//  auto profile = quda::pushProfile(profileInvertC, param.secs, param.gflops);
//  //quda::profilerStart(__func__);
//  ///quda::enable_profiler = true;
//  quda::device::profile::start();
//
//
//  // It was probably a bad design decision to encode whether the system is even/odd preconditioned (PC) in
//  // solve_type and solution_type, rather than in separate members of QudaInvertParam.  We're stuck with it
//  // for now, though, so here we factorize everything for convenience.
//  bool pc_solution = (param.solution_type == QUDA_MATPC_SOLUTION) ||
//    (param.solution_type == QUDA_MATPCDAG_MATPC_SOLUTION);
//  bool pc_solve = (param.solve_type == QUDA_DIRECT_PC_SOLVE) ||
//    (param.solve_type == QUDA_NORMOP_PC_SOLVE) || (param.solve_type == QUDA_NORMERR_PC_SOLVE);
//  bool direct_solve = (param.solve_type == QUDA_DIRECT_SOLVE) ||
//    (param.solve_type == QUDA_DIRECT_PC_SOLVE);
//
//  //quda::Dirac &dirac = *dirac_cg;
//  //quda::Dirac &diracSloppy = *dSloppy;
//  //quda::Dirac &diracPre = *dPre;
//  //quda::Dirac &diracEig = *dEig;
//
//  //quda::ColorSpinorField *in = nullptr;
//  //quda::ColorSpinorField *out = nullptr;
//
//  ////src b, result x
//  //const bool no_clear_guess = (param.use_init_guess == QUDA_USE_INIT_GUESS_YES && !param.chrono_use_resident);
//  /////->Even(), gres->Odd()
//  std::vector<double > nbL;nbL.resize(num_src_inv);
//  for(Int di=0;di<num_src_inv;di++)
//  {
//    quda::ColorSpinorField& x = res.Component(di);
//    quda::ColorSpinorField& b = src.Component(di);
//    if(param.use_init_guess != QUDA_USE_INIT_GUESS_YES){quda::blas::zero(x);}
//    ////if(use_eigen_pc == 0){quda::blas::zero(x);}
//    const double nb  = quda::blas::norm2(b);
//    nbL[di] = nb;
//    if(nb == 0.0){errorQuda("Source has zero norm");}
//
//    //qmessage("==normI di %3d, %+.8e \n", di, nb);
//
//    if (getVerbosity() >= QUDA_VERBOSE) {
//      printfQuda("Source: %g\n", nb);
//      if (param.use_init_guess == QUDA_USE_INIT_GUESS_YES) { printfQuda("Initial guess: %g\n", quda::blas::norm2(x)); }
//    }
//
//    // rescale the source and solution vectors to help prevent the onset of underflow
//    if (param.solver_normalization == QUDA_SOURCE_NORMALIZATION) {
//      quda::blas::ax(1.0 / sqrt(nb), b);
//      quda::blas::ax(1.0 / sqrt(nb), x);
//    }
//
//    if (getVerbosity() >= QUDA_VERBOSE) {
//      double nin  = quda::blas::norm2(b);
//      double nout = quda::blas::norm2(x);
//      if(qlat::get_id_node() == 0){
//      printfQuda("Prepared source = %g\n", nin);
//      printfQuda("Prepared solution = %g\n", nout);
//      }
//    }
//
//    quda::massRescale(b, param, false);
//
//    if (getVerbosity() >= QUDA_VERBOSE) {
//      double nin = quda::blas::norm2(b);
//      if(qlat::get_id_node() == 0){
//      printfQuda("Prepared source post mass rescale = %g\n", nin);
//      }
//    }
//  }
//
//  //dirac.prepare(in, out, x, b, param.solution_type);
//
//  // solution_type specifies *what* system is to be solved.
//  // solve_type specifies *how* the system is to be solved.
//  //
//  // We have the following four cases (plus preconditioned variants):
//  //
//  // solution_type    solve_type    Effect
//  // -------------    ----------    ------
//  // MAT              DIRECT        Solve Ax=b
//  // MATDAG_MAT       DIRECT        Solve A^dag y = b, followed by Ax=y
//  // MAT              NORMOP        Solve (A^dag A) x = (A^dag b)
//  // MATDAG_MAT       NORMOP        Solve (A^dag A) x = b
//  // MAT              NORMERR       Solve (A A^dag) y = b, then x = A^dag y
//  //
//  // We generally require that the solution_type and solve_type
//  // preconditioning match.  As an exception, the unpreconditioned MAT
//  // solution_type may be used with any solve_type, including
//  // DIRECT_PC and NORMOP_PC.  In these cases, preparation of the
//  // preconditioned source and reconstruction of the full solution are
//  // taken care of by Dirac::prepare() and Dirac::reconstruct(),
//  // respectively.
//
//  if (pc_solution && !pc_solve) {
//    errorQuda("Preconditioned (PC) solution_type requires a PC solve_type");
//  }
//
//  //param->use_sloppy_partial_accumulator = 1;
//  Qassert(src.IsComposite());
//  Qassert(res.IsComposite());
//  const Int dim = src.CompositeDim();
//  Qassert(num_src_inv <= dim);
//  if (direct_solve) {
//    if(solve_mode == 0){
//      solverParam.num_src = 1;
//      //quda::ColorSpinorField r = res.Component(0);
//      for(Int di=0;di<num_src_inv;di++)
//      {
//        quda::Solver* solve_cg = quda::Solver::create(solverParam, *m_cg, *mSloppy, *mPre, *mEig);
//        //(*solve_cg)(r, src.Component(di));
//        (*solve_cg)(res.Component(di), src.Component(di));
//        delete solve_cg;
//        // TODO some strange issue for reusing this solve_cg, will zero res[0]....
//        //res.Component(di) = r;//some strange make zeros... for the initial Component
//
//        //double a = quda::blas::norm2(src.Component(di));
//        //double b = quda::blas::norm2(res.Component(di));
//        //qmessage("==normJ di %3d, %+.8e %+.8e \n", di, a, b);
//      }
//
//      //for(Int di=0;di<num_src_inv;di++)
//      //{
//      //  double a = quda::blas::norm2(src.Component(di));
//      //  double b = quda::blas::norm2(res.Component(di));
//      //  qmessage("==normZ di %3d, %+.8e %+.8e \n", di, a, b);
//      //}
//
//    }
//    //if(solve_mode == 1){
//    //  solverParam.num_src = num_src_inv;
//    //  (*solve_cg).blocksolve(res, src);
//    //}
//    //if(solve_mode == 2){
//    //  solverParam.num_src = num_src_inv;
//    //  (*solve_cg).solve(res, src);
//    //}
//
//  }else{errorQuda("Solver not added");}
//
//
//  for(Int di=0;di<num_src_inv;di++)
//  {
//
//    //{
//    //double a = quda::blas::norm2(src.Component(di));
//    //double b = quda::blas::norm2(res.Component(di));
//    //qmessage("==normM di %3d, %+.8e %+.8e, nb %+.8e \n", di, a, b, nbL[di]);
//    //}
//
//    quda::ColorSpinorField& x = res.Component(di);
//    //quda::ColorSpinorField& b = src.Component(di);
//    const double nb = nbL[di];
//
//    if (getVerbosity() >= QUDA_VERBOSE) { printfQuda("Solution = %g\n", quda::blas::norm2(x)); }
//    ////dirac.reconstruct(x, b, param.solution_type);
//
//    if (param.solver_normalization == QUDA_SOURCE_NORMALIZATION) {
//      // rescale the solution
//      quda::blas::ax(sqrt(nb), x);
//      //quda::blas::ax(sqrt(nb), b);
//    }
//
//    if (getVerbosity() >= QUDA_VERBOSE) {
//      printfQuda("Reconstructed solution: %g\n", quda::blas::norm2(x));
//    }
//
//    //{
//    //double a = quda::blas::norm2(src.Component(di));
//    //double b = quda::blas::norm2(res.Component(di));
//    //qmessage("==normK di %3d, %+.8e %+.8e, nb %+.8e \n", di, a, b, nb);
//    //}
//
//  }
//  solverParam.updateInvertParam(param);
//  //////the flops here may be wrong
//  //inv_param.secs = solverParam.secs;
//  //inv_param.gflops = solverParam.gflops;
//  //quda::comm_allreduce_sum(inv_param.gflops);
//  //inv_param.iter = solverParam.iter;
//  // cache is written out even if a Long benchmarking job gets interrupted
//  quda::saveTuneCache();
//  quda::device::profile::stop();
//}

//template<typename Ty>
//inline void quda_inverter::do_inv(Ty* res, Ty* src, const double mass, const double err, const Int niter , const Int prec_type )
//{
//  TIMER_FLOPS("QUDA CG");
//  timeval tm0,tm1;gettimeofday(&tm0, NULL);gettimeofday(&tm1, NULL);
//
//  setup_inv_param_prec(prec_type); ////restore precisions
//  inv_param.tol = err;
//  ///if(err < 1e-6 and err > 1e-10){inv_param.tol_restart    = err*5e+3;}
//  inv_param.maxiter = niter;
//  //solverParam = quda::SolverParam(inv_param);
//  //maxiter, tol, 
//  //if(fermion_type == 0){abort_r("Not suppported!\n");setup_inv_mass(mass);}
//  //if(fermion_type == 1){
//  //  setup_mat_mass(mass);
//  //}
//  setup_mat_mass(mass);
//  //setup_CG();
//
//  //if((void*)csrc->data() != src)qudaMemcpy((void*)csrc->data(), src, csrc->Volume() * spinor_site_size * sizeof(quda::Complex), qudaMemcpyHostToHost);
//  Qassert((void*) (*gsrc).data() != (void*) src);
//  qlat_cf_to_quda_cf((*gsrc), src, geo, map_index);
//  quda::blas::zero(*gres);
//
//  if(eigen_with_nvec == true){
//    inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;
//  }else{
//    inv_param.use_init_guess = QUDA_USE_INIT_GUESS_NO ;
//  }
//
//  //update_eigen_mass(mass, true);
//  //if(use_eigen_pc == 1){
//
//  //  if(fermion_type == 1){
//  //    if(add_high == 0){
//  //      //(*gsrc) = (*csrc);
//  //      prepare_low_prop();
//  //      quda_cf_to_qlat_cf(res, (*gres), geo, map_index);
//  //      //(*cres) = (*gres);
//  //      //if((void*)cres->data() != res){qudaMemcpy(res, (void*)cres->data(), 
//  //      //        cres->Volume() * spinor_site_size * sizeof(quda::Complex), qudaMemcpyHostToHost);}
//  //      gettimeofday(&tm1, NULL);double time0 = tm1.tv_sec - tm0.tv_sec;time0 += (tm1.tv_usec - tm0.tv_usec)/1000000.0;
//  //      if(quda::comm_rank_global() == 0)printfQuda("prepare low Done:  %.6f secs \n", time0);
//  //      return ;
//  //    }
//
//  //    for(Int di=0;di<num_src_inv;di++)
//  //    {
//  //      quda::ColorSpinorField srcP;
//  //      quda::ColorSpinorField sol0;
//  //      quda::ColorSpinorField sol1;
//
//  //      dirac_pc->prepare(srcP, sol0, (*gtmp0).Component(di), (*gsrc).Component(di), QUDA_MAT_SOLUTION);
//  //      dirac_pc->prepare(srcP, sol1, (*gres ).Component(di), (*gsrc).Component(di), QUDA_MAT_SOLUTION);
//
//  //      quda::blas::zero(sol0);quda::blas::zero(sol1);
//  //      if(kSpace.size()!=0){deflate((sol0), (srcP), kSpace, evalsK, false);}
//  //      if(fSpace.size()!=0){deflate((sol1), (srcP), fSpace, evalsF, false);}
//
//  //      quda::Complex lambda  = quda::Complex(1.0, 0.0);
//  //      if(kSpace.size()!=0){lambda  = lambda - quda::blas::cDotProduct(sol0, sol1) / quda::blas::norm2(sol0);}
//  //      quda::blas::caxpy(lambda, sol0, sol1);
//  //    }
//
//  //  }
//  //  //*cres = *gres;
//
//  //}else
//  {
//    if(fermion_type == 1){
//      //for(Int di=0;di<num_src_inv;di++)
//      {
//        quda::ColorSpinorField srcP;
//        quda::ColorSpinorField solP;
//
//        dirac_pc->prepare(srcP, solP, (*gres), (*gsrc), QUDA_MAT_SOLUTION);
//      }
//    }
//  }
//
//  //if(fermion_type == 0){invertQuda((void*)(cres->Even()).data(), (void*)(csrc->Even()).data(), &inv_param);}
//  if(fermion_type == 1)
//  {
//    inv_param.solution_type = QUDA_MATPC_SOLUTION;
//    inv_param.cpu_prec = (*gsrc).Precision();
//    //Qassert((*gsrc).Component(0).Precision() == inv_param.cpu_prec);
//
//    //if(add_high >= 1)
//    {
//      //invertQuda((void*)((*ctmp1).Component(0).Even()).data(), (void*)((*ctmp1).Component(0).Odd()).data(), &inv_param);
//      //for(Int di=0;di<num_src_inv;di++)
//      {
//        //(*gresH).Component(di) = (*gres).Component(di).Even();
//        //(*gsrcH).Component(di) = (*gres).Component(di).Odd();
//        //double a =  quda::blas::norm2((*gsrcH).Component(di));
//        //double b =  quda::blas::norm2((*gresH).Component(di));
//        //qmessage("==srcA %3d, norm2 %+.8e \n", di, a);
//        //qmessage("==resA %3d, norm2 %+.8e \n", di, b);
//
//        //invertQuda((void*)(*gres).Component(di).Even().data(), (void*)(*gres).Component(di).Odd().data(), &inv_param);
//        //invertQuda((void*)(*gresH).Component(di).data(), (void*)(*gsrcH).Component(di).data(), &inv_param);
//        invertQuda((*gres).Even().data(), (*gres).Odd().data(), &inv_param);
//      }
//
//      ////debug print norms
//      //{
//      //  qlat::vector_gpu<qlat::ComplexD > res_tmp;
//      //  res_tmp.resize(gsrc->Volume() * spinor_site_size);
//      //  quda_cf_to_qlat_cf(res_tmp.data(), (*gres).Component(0), geo, map_index);
//      //  qlat::ComplexD normC = res_tmp.norm2();
//      //  double na = quda::blas::norm2((*gres).Component(0).Even());
//      //  double nb = quda::blas::norm2((*gres).Component(0).Odd());
//      //  qmessage("check norm ori, %.20f, %.20f, %.20f \n", normC.real(), na, nb);
//      //  //double nb = quda::blas::norm2(b);
//      //  //qmessage("check norm ori, ");
//      //  //res_tmp.print_norm2();
//      //}
//
//      //invertQuda_COPY(*gresH, *gsrcH, solve_mode);
//
//      //for(Int di=0;di<num_src_inv;di++)
//      //{
//      //  (*gres).Component(di).Even() = (*gresH).Component(di);
//      //  (*gres).Component(di).Odd()  = (*gsrcH).Component(di);
//
//      //  //double a =  quda::blas::norm2((*gsrcH).Component(di));
//      //  //double b =  quda::blas::norm2((*gresH).Component(di));
//      //  //qmessage("==srcB %3d, norm2 %+.8e \n", di, a);
//      //  //qmessage("==resB %3d, norm2 %+.8e \n", di, b);
//      //}
//
//      //{
//      //  qlat::vector_gpu<qlat::ComplexD > res_tmp;
//      //  res_tmp.resize(gsrc->Volume() * spinor_site_size);
//      //  quda_cf_to_qlat_cf(res_tmp.data(), (*gres).Component(0), geo, map_index);
//      //  double na = quda::blas::norm2((*gres).Component(0).Even());
//      //  double nb = quda::blas::norm2((*gres).Component(0).Odd());
//      //  qlat::ComplexD normC = res_tmp.norm2();
//      //  qmessage("check norm ori, %.20f, %.20f, %.20f \n", normC.real(), na, nb);
//      //  //qmessage("check norm ori, ");
//      //  //res_tmp.print_norm2();
//      //}
//
//
//      //quda::ColorSpinorParam cs_cpu1 = quda::ColorSpinorParam(*gsrcH);
//      //cs_cpu1.location = QUDA_CPU_FIELD_LOCATION;
//      //cs_cpu1.is_composite  = false;
//      //cs_cpu1.is_component  = false;
//
//      //quda::ColorSpinorField csrc(cs_cpu1);
//      //quda::ColorSpinorField cres(cs_cpu1);
//
//      //const size_t Nd = (*gsrcH).Component(0).Volume() * spinor_site_size * sizeof(quda::Complex);
//      //qudaMemcpy((void*)cres.data(), (*gresH).Component(0).data(), Nd, qudaMemcpyDeviceToHost);
//      //qudaMemcpy((void*)csrc.data(), (*gsrcH).Component(0).data(), Nd, qudaMemcpyDeviceToHost);
//      ////csrc = (*gsrcH).Component(0);
//      ////cres = (*gresH).Component(0); 
//      //invertQuda((void*)(cres).data(), (void*)(csrc).data(), &inv_param);
//      //qudaMemcpy((*gresH).Component(0).data(), (void*)cres.data(), Nd, qudaMemcpyHostToDevice);
//      //qudaMemcpy((*gsrcH).Component(0).data(), (void*)csrc.data(), Nd, qudaMemcpyHostToDevice);
//
//      ////(*gsrcH).Component(0) = csrc;
//      ////(*gresH).Component(0) = cres; 
//
//    }
//
//    //*gres = *cres;
//    //*gres = *ctmp1;
//    //for(Int di=0;di<num_src_inv;di++)
//    {
//      //(*gres).Component(di) = (*ctmp1).Component(di);
//      dirac_pc->reconstruct((*gres), (*gsrc), QUDA_MAT_SOLUTION);
//      //double a =  quda::blas::norm2((*gsrc).Component(di));
//      //double b =  quda::blas::norm2((*gres).Component(di));
//      //qmessage("==srcC %3d, norm2 %+.8e \n", di, a);
//      //qmessage("==resC %3d, norm2 %+.8e \n", di, b);
//    }
//
//    //*cres = *gres;
//    inv_param.solution_type = QUDA_MAT_SOLUTION;
//  }
//
//  if(check_residue == 1)
//  { 
//    /////TODO not working for multi GPU
//    /////===check residue
//    quda::Complex n_unit(-1.0, 0.0);
//    
//    quda::ColorSpinorParam cs_gpu1 = quda::ColorSpinorParam(*gres);
//    cs_gpu1.is_composite  = false;
//    cs_gpu1.is_component  = false;
//
//    quda::ColorSpinorField buf(cs_gpu1);
//    //for(Int di=0;di<num_src_inv;di++)
//    {
//      quda::ColorSpinorField& src = (*gsrc);
//      quda::ColorSpinorField& res = (*gres);
//      quda::blas::zero(buf);
//      (*mat)((buf), res);
//      qacc_barrier(dummy);
//
//      quda::Complex temp = quda::blas::norm2(res);
//      ////if(quda::comm_rank_global()== 0)printf("===result %.8e \n", temp.real());
//
//      if(fermion_type == 0){quda::blas::ax(1.0/(2.0*inv_param.kappa), buf);}
//
//      quda::Complex evals  = quda::blas::cDotProduct(src, buf) / quda::blas::norm2(src);
//      quda::Complex factor = sqrt(quda::blas::norm2(buf)) / sqrt(quda::blas::norm2(src));
//      quda::blas::caxpby(evals, src , n_unit, buf);
//      //if(fermion_type == 1){
//      //}
//      quda::Complex residual = sqrt(quda::blas::norm2(buf)/ quda::blas::norm2(src));
//      quda::Complex nor_e1 = sqrt( quda::blas::norm2((res).Even()));
//      quda::Complex nor_o1 = sqrt( quda::blas::norm2((res).Odd() ));
//
//      quda::Complex nor_e = sqrt( quda::blas::norm2((src).Even()));
//      quda::Complex nor_o = sqrt( quda::blas::norm2((src).Odd() ));
//      quda::Complex res_e = sqrt(quda::blas::norm2((buf).Even())/ quda::blas::norm2((src)));
//      quda::Complex res_o = sqrt(quda::blas::norm2((buf).Odd())/ quda::blas::norm2((src)));
//
//      if(quda::comm_rank_global()== 0)printf("===solution residual %.3e, factor %.3e, sol norm %.8e, e %.8e , o %.8e \n", 
//            residual.real(), factor.real(), temp.real(), res_e.real(), res_o.real());
//
//      //if(quda::comm_rank_global()== 0)printf("===solution residual %.3e, factor %.3e, sol norm %.8e, e %.8e %.3e %.3e, o %.8e %.3e %.3e \n", 
//      //      residual.real(), factor.real(), temp.real(), res_e.real(), nor_e1.real(), nor_e.real(), res_o.real(), nor_o1.real(), nor_o.real());
//      inv_residue = residual.real();
//    }
//
//    //qacc_barrier(dummy);
//    //fflush_MPI();
//    ////std::vector<quda::Complex> sump(3);
//    ////quda::comm_allreduce_sum(sump);
//    /////===check residue
//  }
//  //if((void*)cres->data() != res){qudaMemcpy(res, (void*)cres->data(), cres->Volume() * spinor_site_size * sizeof(quda::Complex), qudaMemcpyHostToHost);}
//  Qassert((void*) (*gres).data() != (void*) res);
//  //*ctmp0 = *gres;
//  //const size_t Nd = (*gsrc).Component(0).Volume() * spinor_site_size * sizeof(quda::Complex);
//  //qlat::vector<qlat::ComplexD > buf;buf.resize(Nd/sizeof(quda::Complex));
//  //for(Int di=0;di<num_src_inv;di++){
//  //  quda_cf_to_qlat_cf((qlat::ComplexD*) buf.data(), (qlat::ComplexD*) ctmp0->Component(di).data(), geo, 3);
//  //  qudaMemcpy(&res[di*Nd/sizeof(quda::Complex)], (void*) buf.data(), Nd, qudaMemcpyDeviceToDevice);
//  //}
//  quda_cf_to_qlat_cf(res, (*gres), geo, map_index);
//
//  ////{errorQuda("QUDA may have bugs right now for the code on CUDA FIELD below\n");}
//
//  gettimeofday(&tm1, NULL);double time0 = tm1.tv_sec - tm0.tv_sec;time0 += (tm1.tv_usec - tm0.tv_usec)/1000000.0;
//
//  //int verbos = 0;
//  inv_param.secs += 1e-25;
//  //std::string val = qlat::get_env(std::string("qlat_quda_verbos"));
//  //if(val != ""){verbos = stringtonum(val);}
//  if(quda_verbos >= 0 or quda_verbos == -2)
//  {
//    qmessage("Done: %8d iter / %.6f secs = %.3f Gflops %.3f Gflops/GPU, Cost %.3f Gflops \n",
//          inv_param.iter, inv_param.secs, inv_param.gflops / inv_param.secs,
//          inv_param.gflops / (inv_param.secs * qlat::get_num_node()), inv_param.gflops);
//  }
//  inv_time   = time0;
//  inv_iter   = inv_param.iter;
//  inv_gflops = inv_param.gflops / (inv_param.secs * qlat::get_num_node());
//  timer.flops +=  inv_param.gflops * 1024 * 1024 * 1024 / qlat::get_num_node();
//
//  ////quda::saveTuneCache();
//}

quda_inverter::~quda_inverter()
{
  free_mem();
}

/*
  buf will be overwritten 
  prop will be ignored, the precondtioning with eigensystem is used
  buf should be double for QUDA ?
  buf could be high precision but use only Ty
  low_only : 0, do full props, 1, only low mode prop, -1 use test inverter
  inv_even_even : flag to do even-even / odd-odd inversions, odd-odd have no eigen support
    may not need this, odd source will be only a mass rescale
*/
template<typename Ty>
void get_staggered_prop_group(quda_inverter& qinv, qlat::vector<Ty* >& src, qlat::vector<Ty* >& prop,
    const double mass, const double err, const Int niter, Int low_only = 0, const Int prec_type = 0, 
    const bool inv_even_even = true)
{
  TIMER("get_staggered_prop_group");
  //timeval tm0,tm1;gettimeofday(&tm0, NULL);gettimeofday(&tm1, NULL);

  const Int nsrc = src.size();
  Qassert( nsrc == Long(prop.size()) );
  if(nsrc == 0){return ;}
  qlat::vector<Long >& map = qinv.map_index;
  const Geometry& geo = qinv.geoB();
  Qassert(low_only == 0 or low_only == 1 or low_only == -1);///0 for full, 1 for low only, -1 other ways to solve
  const Int restart_cg = 0;


  qlat::vector<Ty* > buf;
  qinv.get_inversion_bufs(buf, nsrc, 0);

  //auto& Ebuf = eigen.Ebuf;
  //auto& Sbuf = eigen.Sbuf;
  const Int DIM = 3;
  const LInt Nd = qinv.geoB().local_volume() * DIM;

  if(qinv.eigen_with_nvec == true){
    qinv.inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;
  }else{
    qinv.inv_param.use_init_guess = QUDA_USE_INIT_GUESS_NO ;
  }
  qinv.setup_inv_param_prec(prec_type); ////restore precisions
  qinv.inv_param.tol = err;
  qinv.inv_param.maxiter = niter;

  qinv.inv_param.iter = 0;
  qinv.inv_param.secs = 0;
  qinv.inv_param.gflops = 0;

  qinv.setup_mat_mass(mass);

  quda::Dirac *dirac_pc = qinv.dirac_pc;
  if(!inv_even_even)
  {
    Qassert(low_only == 0 and qinv.eigen_with_nvec == false and restart_cg == 0);
    qinv.inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
    dirac_pc = qinv.dirac_pc_odd;
  }

  if(low_only == 0 or low_only == -1)
  {
    //TIMER_FLOPS("get_staggered_prop_group_inv");
    //qinv.setup_CG();
    quda::ColorSpinorParam param = qinv.cs_gpu;
    param.is_composite  = false;
    param.is_component  = false;
    ///param.composite_dim = 1;
    param.setPrecision(qinv.inv_param.cuda_prec, qinv.inv_param.cuda_prec, true);
    //param.setPrecision(QUDA_DOUBLE_PRECISION, qinv.inv_param.cuda_prec, true);
    param.create = QUDA_REFERENCE_FIELD_CREATE;
    //
    //qlat::vector<qlat::ComplexT<double >* > Abuf;Abuf.resize(nsrc*2);
    //qlat::vector<qlat::ComplexT<float  >* > Bbuf;Bbuf.resize(nsrc*2);
    Int buf_prec = 0;
    if(qinv.inv_param.cuda_prec == QUDA_DOUBLE_PRECISION){buf_prec = 0;}
    if(qinv.inv_param.cuda_prec == QUDA_SINGLE_PRECISION){buf_prec = 1;}
    if(buf_prec == -1){abort_r("QUDA buffer prec not supported yet!");}
    //
    //quda::ColorSpinorParam param_src = qinv.cs_gpu;
    //param_src.is_composite  = false;
    //param_src.is_component  = false;
    //if(get_data_type_is_double<Ty >()){
    //  param_src.setPrecision(QUDA_DOUBLE_PRECISION, QUDA_DOUBLE_PRECISION, true);
    //}else{
    //  param_src.setPrecision(QUDA_SINGLE_PRECISION, QUDA_SINGLE_PRECISION, true);
    //}
    //param_src.create = QUDA_REFERENCE_FIELD_CREATE;
    //
    //using internal solver which always accept double input
    //Qassert(qinv.inv_param.cuda_prec == QUDA_DOUBLE_PRECISION);
    //Qassert(get_data_type_is_double<Tk >());//buf must be double !
    /////TODO need to fix the precision issue here
    //
    //prop have to be double if buf_prec is double
    if(buf_prec == 0){
      Qassert( get_data_type_is_double<Ty >());
      //Qassert(prec_type == 0 or prec_type == 2 or prec_type == 10)
    }
    //if(!get_data_type_is_double<Ty >()){
    //  Qassert(!( prec_type == 0 or prec_type == 2 or prec_type == 10 ) );
    //}
    //if(buf_prec == 1){
    //  Qassert(!get_data_type_is_double<Ty >());
    //  Qassert(prec_type == 1 or prec_type == 11)
    //}

    // double precion goal with only single output ......
    if(Is_data_double<Ty >() == 0){
      Qassert(!( prec_type == 0 or prec_type == 2 or prec_type == 10 ) );
    }

    std::vector<quda::ColorSpinorField* > Qvec;Qvec.resize(nsrc * 2);
    //quda::ColorSpinorField* Qsrc_tmp = NULL;
    ///0--nsrc for Qvec --> src; nsrc -- 2 nsrc --> res
    for(Int vi=0;vi<nsrc;vi++){
      param.v = (void*) prop[vi];
      Qvec[vi] = new quda::ColorSpinorField(param);
      qlat_cf_to_quda_cf(*Qvec[vi], src[vi], geo, map);
      //qlat_cf_to_quda_cf(buf[vi], src[vi], DIM, geo, map);
      //src and prop may be the same
      //Qsrc_tmp = new quda::ColorSpinorField(param); 
      //qlat_cf_to_quda_cf()
      //delete Qsrc_tmp;Qsrc_tmp=NULL;
      //if(buf_prec==0)cpy_GPU((qlat::ComplexT<double >*) Qvec[vi]->data(), buf[vi], Nd);
      //if(buf_prec==1)cpy_GPU((qlat::ComplexT<float  >*) Qvec[vi]->data(), buf[vi], Nd);

      param.v = (void*) buf[vi];
      Qvec[nsrc + vi] = new quda::ColorSpinorField(param);
      quda::blas::zero(*Qvec[nsrc + vi]);
    }

    //if(qinv.gtmp_invG0 == NULL){
    //  param.create = QUDA_ZERO_FIELD_CREATE;
    //  qinv.gtmp_invG0 = new quda::ColorSpinorField(param);
    //}

    qlat::vector<void* > Psrc;
    qlat::vector<void* > Pres;
    qlat::vector<void* > Pbuf;///even of (*Qvec[vi]) which is not needed later
    //qlat::ComplexT<double> rD;
    //qlat::ComplexT<float > rF;

    if(low_only == 0)
    {
      Psrc.resize(nsrc);
      Pres.resize(nsrc);
      Pbuf.resize(nsrc);
      for(Int vi=0;vi<nsrc;vi++)
      {
        quda::ColorSpinorField srcP;
        quda::ColorSpinorField sol1;
        dirac_pc->prepare(    srcP, sol1, (*Qvec[nsrc + vi]), (*Qvec[vi]), QUDA_MAT_SOLUTION);
        Psrc[vi] = Qvec[nsrc + vi]->Odd().data();
        Pres[vi] = Qvec[nsrc + vi]->Even().data();
        Pbuf[vi] = Qvec[       vi]->Even().data();
      }
      qinv.inv_param.cpu_prec = Qvec[0]->Precision();

      if(restart_cg == 0)
      {
        if(inv_even_even){
          qinv.deflate_Ty(Pres, Psrc, mass, buf_prec, 0);
        }
        //else{
        //  for(Int vi=0;vi<nsrc;vi++){
        //    quda::blas::zero(Qvec[nsrc + vi]->Odd());
        //  }
        //}
        qinv.inv_param.solution_type = QUDA_MATPC_SOLUTION;

        qlat::vector<void* > srcI;srcI.resize(nsrc);
        qlat::vector<void* > resI;resI.resize(nsrc);
        for(Int vi=0;vi<nsrc;vi++)
        {
          if(inv_even_even){
            srcI[vi] = Qvec[nsrc + vi]->Odd().data();
            resI[vi] = Qvec[nsrc + vi]->Even().data();
          }else{
            srcI[vi] = Qvec[nsrc + vi]->Even().data();
            resI[vi] = Qvec[nsrc + vi]->Odd().data();
          }
        }

        {
          qinv.callMultiSrcQuda(resI, srcI);
        }

        for(Int vi=0;vi<nsrc;vi++)
        {
          dirac_pc->reconstruct(    (*Qvec[nsrc + vi]), (*Qvec[vi]), QUDA_MAT_SOLUTION);
          //if(inv_even_even){
          //  qinv.dirac_pc->reconstruct(    (*Qvec[nsrc + vi]), (*Qvec[vi]), QUDA_MAT_SOLUTION);
          //}else{
          //  qinv.dirac_pc_odd->reconstruct((*Qvec[nsrc + vi]), (*Qvec[vi]), QUDA_MAT_SOLUTION);
          //}
        }
        qinv.inv_param.solution_type = QUDA_MAT_SOLUTION;
      }else{
        //const Int solve_mode = 0;
        ////qacc_barrier(dummy);
        ////deflate on CPU with Ebuf???
        //const Int mode_def = 1;
        qinv.inv_param.solution_type = QUDA_MATPC_SOLUTION;

        const Int iter_group =  200;
        //const Int iter_group = qinv.inv_param.maxiter;
        const Int maxN = (qinv.inv_param.maxiter + iter_group - 1) / iter_group;
        qinv.inv_param.maxiter = iter_group;

        quda::ColorSpinorParam cs_tmp(*qinv.gsrcH);
        cs_tmp.is_composite  = false;
        cs_tmp.is_component  = false;
        cs_tmp.composite_dim = 1;

        std::vector<quda::ColorSpinorField > gsumV;gsumV.resize(nsrc);
        std::vector<quda::ColorSpinorField > gresV;gresV.resize(nsrc);
        std::vector<quda::ColorSpinorField > gsrcV;gsrcV.resize(nsrc);
        std::vector<quda::ColorSpinorField > gbufV;gbufV.resize(nsrc);
        quda::ColorSpinorField tmp0(cs_tmp);
        quda::ColorSpinorField tmp1(cs_tmp);
        std::vector<double > residual_vecs;residual_vecs.resize(nsrc);
        std::vector<double > norm_vecs;norm_vecs.resize(nsrc);
        //std::vector<double > norm_vecs_buf;norm_vecs_buf.resize(nsrc);
        for(Int vi=0;vi<nsrc;vi++){
          residual_vecs[vi] = 1.0;
          gresV[vi] = quda::ColorSpinorField(cs_tmp); 
          gsumV[vi] = quda::ColorSpinorField(cs_tmp); 
          gsrcV[vi] = quda::ColorSpinorField(cs_tmp); 
          gbufV[vi] = quda::ColorSpinorField(cs_tmp); 

          quda::blas::zero(gsumV[vi]);
          gresV[vi] = Qvec[nsrc + vi]->Even();
          gsrcV[vi] = Qvec[nsrc + vi]->Odd();

          Psrc[vi]  = gsrcV[vi].data(); // redirect bufs
          Pres[vi]  = gresV[vi].data();
          Pbuf[vi]  = gbufV[vi].data();
          //norm_vecs_buf[vi] = std::sqrt( quda::blas::norm2(gsrcV[vi]) );
          //norm_vecs[vi] = std::sqrt( quda::blas::norm2(gsrcV[vi]) );
          norm_vecs[vi] = std::sqrt( quda::blas::norm2(gsrcV[vi]) ) ;
          //qmessage("vec %d, norm %.8e \n", vi, norm_vecs[vi]);
        }

        quda::Dirac& xpay = *dirac_pc;
        double mass2 = 4 * xpay.Mass() * xpay.Mass();
        quda::Complex One = quda::Complex(1.0, 0.0);
        Int iter_total    = 0;
        Int count_diverge = 0;

        for(Int it =0;it<maxN;it++)
        {
          //qinv.deflate_Ty(Pres, Psrc, mass, buf_prec, 0);
          //for(Int vi=0;vi<nsrc;vi++){
          //  quda::blas::ax( 1.0 / norm_vecs_buf[vi], gsrcV[vi]);
          //}
          Int converge = 0;
          for(Int vi=0;vi<nsrc;vi++){
            if(residual_vecs[vi] < err){converge += 1; }
          }
          if(converge == nsrc){continue ;}

          if(count_diverge / nsrc >= 2){
            qmessage("Too many inversion diverge %d, nsrc %d!", count_diverge, nsrc);
            break;
          }
          for(Int vi=0;vi<nsrc;vi++){
            qinv.inv_param.tol = err / residual_vecs[vi];
          }

          qinv.deflate_Ty(Pres, Psrc, mass, buf_prec, 0);

          //qlat::vector<void* > srcI;srcI.resize(nsrc);
          //qlat::vector<void* > resI;resI.resize(nsrc);
          //for(Int vi=0;vi<nsrc;vi++)
          //{
          //  srcI[vi] = Qvec[nsrc + vi]->Odd().data();
          //  resI[vi] = Qvec[nsrc + vi]->Even().data();
          //}

          {
            TIMER("invertQuda group");
            qinv.inv_param.num_src = nsrc;
            qinv.callMultiSrcQuda(Pres, Psrc);
            //invertMultiSrcQuda(Pres.data(), Psrc.data(), &qinv.inv_param);
            //invertQuda(Qvec[nsrc + vi]->Even().data(), Qvec[nsrc + vi]->Odd().data(), &qinv.inv_param);
          }

          //for(Int vi=0;vi<nsrc;vi++){
          //  ////TIMER("invertQuda group");
          //  invertQuda(gresV[vi].data(), gsrcV[vi].data(), &qinv.inv_param);
          //}
          iter_total += iter_group;

          for(Int vi=0;vi<nsrc;vi++){
            quda::blas::caxpy( 1.0*One , gresV[vi], gsumV[vi]);
            //quda::blas::caxpy( 1.0*One * norm_vecs[vi], gresV[vi], gsumV[vi]);

            //invertQuda((void*)Qvec[nsrc + vi]->Even().data(), (void*)Qvec[nsrc + vi]->Odd().data(), &qinv.inv_param);
 
            gsrcV[vi] = Qvec[nsrc + vi]->Odd();
            xpay.Dslash(tmp0, gsumV[vi], QUDA_ODD_PARITY);
            xpay.DslashXpay(tmp1, tmp0, QUDA_EVEN_PARITY, gsumV[vi], mass2);
            quda::blas::caxpy(-1.0*One, tmp1, gsrcV[vi]);
            //norm_vecs[vi] = std::sqrt( quda::blas::norm2(gsrcV[vi]) );
            //double na = norm_vecs[vi];

            double na = std::sqrt( quda::blas::norm2(gsrcV[vi]) ) / norm_vecs[vi];
            double Tmp= std::sqrt( quda::blas::norm2(gsrcV[vi]) );
            //qmessage("check res %.8e, self %.8e, res norm %.8e srcnorm %.8e \n", qinv.inv_param.true_res[0], na, Tmp,  norm_vecs[vi]);
            if(na > residual_vecs[vi] * 1.3){
              qmessage("residual diverge current %.8e previous %.8e! \n", na, residual_vecs[vi]);
              count_diverge += 1;
            }
            residual_vecs[vi] = na;
          }
        }

        for(Int vi=0;vi<nsrc;vi++)
        {
          gsrcV[vi] = Qvec[nsrc + vi]->Odd();
          xpay.Dslash(tmp0, gsumV[vi], QUDA_ODD_PARITY);
          xpay.DslashXpay(tmp1, tmp0, QUDA_EVEN_PARITY, gsumV[vi], mass2);
          //tmp0 = gsrcV[vi];
          quda::blas::caxpy(-1.0*One, tmp1, gsrcV[vi]);
          double na = std::sqrt( quda::blas::norm2(gsrcV[vi]) ) / norm_vecs[vi] ;
          qmessage("===check res self %.8e, niter %d \n", na, iter_total);
        }

        //Qvec[nsrc + vi]->Odd()  = (*qinv.gsrcH).Component(0);

        for(Int vi=0;vi<nsrc;vi++){
          Qvec[nsrc + vi]->Even() = gsumV[vi];
          dirac_pc->reconstruct((*Qvec[nsrc + vi]), (*Qvec[vi]), QUDA_MAT_SOLUTION);
          //delete res;delete src;
        }
        qinv.inv_param.solution_type = QUDA_MAT_SOLUTION;
      }

    }

    if(low_only == -1)
    {
      Psrc.resize(nsrc);
      Pres.resize(nsrc);
      Pbuf.resize(nsrc);

      //const Int solve_mode = 0;
      for(Int vi=0;vi<nsrc;vi++)
      {
        quda::ColorSpinorField srcP;
        quda::ColorSpinorField sol1;
        dirac_pc->prepare(srcP, sol1, (*Qvec[nsrc + vi]), (*Qvec[vi]), QUDA_MAT_SOLUTION);
        Psrc[vi] = Qvec[nsrc + vi]->Odd().data();
        Pres[vi] = Qvec[nsrc + vi]->Even().data();
        Pbuf[vi] = Qvec[       vi]->Even().data();
        Qvec[       vi]->Even() =  Qvec[nsrc + vi]->Odd();
        quda::blas::zero(Qvec[nsrc + vi]->Even());
      }

      //remove vectors from src
      qinv.deflate_Ty(Psrc, Psrc, mass, buf_prec, 2);
      //{
      //  if(buf_prec == 0){eigen.deflate(Psrc, Psrc, rD, mass, 2);}
      //  if(buf_prec == 1){eigen.deflate(Psrc, Psrc, rF, mass, 2);}
      //  if(eigenD != NULL){if(buf_prec == 0){eigenD->deflateA(Psrc, Psrc, rD, mass,  2);}}
      //  if(eigenD != NULL){if(buf_prec == 1){eigenD->deflateA(Psrc, Psrc, rF, mass,  2);}}
      //}

      qinv.inv_param.solution_type = QUDA_MATPC_SOLUTION;
      for(Int vi=0;vi<nsrc;vi++)
      {
        // the COPY could use single precison as input
        qinv.invertQuda_COPY_single(Qvec[nsrc + vi]->Even(), Qvec[nsrc + vi]->Odd());
      }

      //add low to results
      qinv.deflate_Ty(Pres, Pbuf, mass, buf_prec, 0, 0);
      //{
      //  if(buf_prec == 0){eigen.deflateA(Pres, Pbuf, rD, mass, 0);}
      //  if(buf_prec == 1){eigen.deflateA(Pres, Pbuf, rF, mass, 0);}
      //  if(eigenD != NULL){if(buf_prec == 0){eigenD->deflateA(Pres, Pbuf, rD, mass,  0);}}
      //  if(eigenD != NULL){if(buf_prec == 1){eigenD->deflateA(Pres, Pbuf, rF, mass,  0);}}
      //}

      for(Int vi=0;vi<nsrc;vi++){
        dirac_pc->reconstruct((*Qvec[nsrc + vi]), (*Qvec[vi]), QUDA_MAT_SOLUTION);
      }
      qinv.inv_param.solution_type = QUDA_MAT_SOLUTION;
    }

    for(Int vi=0;vi<nsrc;vi++){
      quda_cf_to_qlat_cf(prop[vi], *Qvec[nsrc + vi], geo, map);
      //if(buf_prec==0){
      //  quda_cf_to_qlat_cf(prop[vi], (qlat::ComplexT<double >*) (*Qvec[nsrc + vi]).data(), DIM, geo, map);
      //}
      //if(buf_prec==1){
      //  quda_cf_to_qlat_cf(prop[vi], (qlat::ComplexT<float  >*) (*Qvec[nsrc + vi]).data(), DIM, geo, map);
      //}
    }

    for(Int vi=0;vi<2*nsrc;vi++){
      delete Qvec[vi];
    }
    //gettimeofday(&tm1, NULL);double time0 = tm1.tv_sec - tm0.tv_sec;time0 += (tm1.tv_usec - tm0.tv_usec)/1000000.0;
    qinv.inv_param.secs += 1e-25;
    //timer.flops += qinv.inv_param.gflops;
    if((qinv.quda_verbos >= 0 or qinv.quda_verbos == -2) and qinv.inv_param.iter > 0)
    {
      qmessage("Done: %8d iter / %.6f secs = %.3f Gflops %.3f Gflops/GPU, Cost %.3f Gflops \n",
            qinv.inv_param.iter, qinv.inv_param.secs, qinv.inv_param.gflops / qinv.inv_param.secs,
            qinv.inv_param.gflops / (qinv.inv_param.secs * qlat::get_num_node()), qinv.inv_param.gflops);
    }

  }

  if(low_only == 1)
  {
    qinv.get_low_prop(prop, src, buf);
  }

  if(!inv_even_even)
  {
    qinv.inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
  }

}

template<typename Ty>
void get_staggered_prop_group(quda_inverter& qinv, std::vector<colorFT >& src, std::vector<colorFT >& prop,
    const double mass, const double err, const Int niter, Int low_only = 0, const Int prec_type = 0)
{
  qlat::vector<Ty* > s0;
  qlat::vector<Ty* > r0;
  const Int nsrc = src.size();
  Qassert(nsrc == Long(prop.size()));
  s0.resize(nsrc);
  r0.resize(nsrc);
  for(Int i=0;i<nsrc;i++)
  {
    Qassert(src[i].initialized and prop[i].initialized);
    s0[i] = (Ty*) qlat::get_data(src[i] ).data();
    r0[i] = (Ty*) qlat::get_data(prop[i]).data();
  }
  get_staggered_prop_group(qinv, s0, r0, mass, err, niter, low_only, prec_type);
}

// avoid using it, only SRH
template<typename Ty>
void get_staggered_prop(quda_inverter& qinv, Ty* src, Ty* prop,
    const double mass, const double err, const Int niter, Int low_only = 0, const Int prec_type = 0)
{
  qinv.setup_mat_mass(mass);
  qlat::vector<Ty* > sP;
  qlat::vector<Ty* > rP;
  sP.resize(1);
  rP.resize(1);
  sP[0] = src ;
  rP[0] = prop;

  get_staggered_prop_group(qinv, sP, rP, mass, err, niter, low_only, prec_type);

  //TIMER("QUDA inversions");
  ////if(qinv.fermion_type == 0)qinv.setup_inv_mass(mass);
  ////const Int Dim = 3;
  ////qlat_cf_to_quda_cf((qlat::ComplexD*) qinv.csrc->data(), src, geo, Dim);

  ////{
  ////const Long V = qinv.geoB().local_volume();
  ////qlat::vector_gpu<Ty > tmp;tmp.resize(V*3);
  ////tmp.copy_from(src, V*3);
  ////Ty norm = tmp.norm();
  ////qmessage("==src norm %.8e \n", norm.real());
  ////}

  ////Ty norm = norm_FieldM(prop);
  ////qmessage("normp %.3e %.3e \n", norm.real(), norm.imag());
  //if(low_only == 0){
  //  //qinv.do_inv(qinv.cres->data(), qinv.csrc->data(), mass, err, niter, prec_type);
  //  qinv.do_inv(prop, src, mass, err, niter, prec_type);
  //}

  //if(low_only == 1){
  //  qlat_cf_to_quda_cf((*qinv.gsrc), src, qinv.geoB(), qinv.map_index);
  //  qinv.prepare_low_prop();
  //  quda_cf_to_qlat_cf(prop, (*qinv.gres), qinv.geoB(), qinv.map_index);
  //}

  //if(low_only == 2){
  //  if(qinv.gadd == NULL){
  //    quda::ColorSpinorParam param(*qinv.gres);
  //    qinv.gadd  = quda::ColorSpinorField::Create(param);
  //  }
  //  qlat_cf_to_quda_cf((*qinv.gtmp0), src, qinv.geoB(), qinv.map_index);

  //  qinv.do_inv(prop, src, mass, err, niter, prec_type);
  //  //qinv.do_inv(qinv.cres->data(), qinv.csrc->data(), mass, err, niter, prec_type);
  //  (*qinv.gadd) = (*qinv.gres);
  //  (*qinv.gsrc) = (*qinv.gtmp0);
  //  qinv.prepare_low_prop();

  //  quda::Complex n_unit(+1.0, 0.0);
  //  quda::Complex p_unit(-1.0, 0.0);
  //  quda::blas::caxpby(n_unit, *qinv.gadd  , p_unit, *qinv.gres);
  //  //(*qinv.cres) = (*qinv.gres);
  //  quda_cf_to_qlat_cf(prop, (*qinv.gres), qinv.geoB(), qinv.map_index);
  //}
}

template<typename Ty>
void get_staggered_prop(quda_inverter& qinv, qlat::FieldM<Ty, 3>& src, qlat::FieldM<Ty, 3>& prop
    , const double mass, const double err, const Int niter, Int low_only = 0, const Int prec_type = 0)
{
  Ty* srcP = (Ty*) qlat::get_data(src).data();

  const Geometry& geo = src.geo();
  if(!prop.initialized){prop.init(geo);} ////allocate mem for prop
  Ty* propP = (Ty*) qlat::get_data(prop).data();

  //get_staggered_prop(qinv, srcP, propP, mass, err, niter, low_only, prec_type);

  qlat::vector<Ty* > sP;
  qlat::vector<Ty* > rP;
  sP.resize(1);
  rP.resize(1);
  for(Int iv=0;iv<1;iv++){
    sP[iv] = srcP;
    rP[iv] = propP;
  }
  get_staggered_prop_group(qinv, sP, rP, mass, err, niter, low_only, prec_type);

  //   norm = norm_FieldM(prop);
  //qmessage("normp %.3e %.3e \n", norm.real(), norm.imag());
}

template<typename Ty>
void get_staggered_prop(quda_inverter& qinv, std::vector<qlat::FieldM<Ty, 3> >& srcL, std::vector<qlat::FieldM<Ty, 3> >& propL
    , const double mass, const double err, const Int niter, Int low_only = 0, const Int prec_type = 0)
{
  const Int Ninv = srcL.size();
  if(Ninv <= 0){propL.resize(0);return ;}

  qlat::vector<Ty* > sP;
  qlat::vector<Ty* > rP;
  sP.resize(Ninv);
  rP.resize(Ninv);
  for(Int iv=0;iv<Ninv;iv++){
    sP[iv] = (Ty*) qlat::get_data(srcL[iv] ).data();
    rP[iv] = (Ty*) qlat::get_data(propL[iv]).data();
  }
  get_staggered_prop_group(qinv, sP, rP, mass, err, niter, low_only, prec_type);

  //for(Int iv=0;iv<Ninv;iv++)
  //{
  //  get_staggered_prop(qinv,srcL[iv], propL[iv], mass, err, niter, low_only, prec_type);
  //}
}

template<typename Ty>
void get_staggered_prop(quda_inverter& qinv, qlat::vector_gpu<Ty >& src, qlat::vector_gpu<Ty >& prop,
    const double mass, const double err, const Int niter, Int low_only = 0, const Int prec_type = 0)
{
  const size_t Vl = qinv.V * 3;
  Qassert(src.size() % Vl == 0);
  prop.resize(src.size());
  const Int Ninv = src.size() / Vl;
  if(Ninv <= 0){prop.resize(0);return ;}

  //Qassert(Ninv <= qinv.num_src);
  //int tmp_inv = qinv.num_src_inv;
  //qinv.num_src_inv = Ninv;

  Ty* srcP = src.data();
  Ty* resP = prop.data();
  move_index mv_civ;
  mv_civ.move_civ_out(srcP, resP, 1, Vl, 3, 1, true);

  qlat::vector<Ty* > sP;
  qlat::vector<Ty* > rP;
  sP.resize(Ninv);
  rP.resize(Ninv);
  for(Int iv=0;iv<Ninv;iv++){
    sP[iv] = srcP[iv*Vl];
    rP[iv] = resP[iv*Vl];
  }

  get_staggered_prop_group(qinv, sP, rP, mass, err, niter, low_only, prec_type);
  //   get_staggered_prop(qinv, resP, resP,  mass, err, niter, low_only, prec_type);
  mv_civ.move_civ_in( resP, resP, 1, 3, Vl, 1, true);

}

////even dslash flops
inline double get_stagger_even_flops(const Geometry& geo)
{
  Int nv[4];
  for(Int i=0;i<4;i++){
    //nv[i] = geo.node_site[i] * geo.geon.size_node[i];
    nv[i] = geo.node_site[i];
  }
  const double fac = 670.00;
  double flops     = nv[0]  * nv[1] * nv[2] * nv[3] * fac;
  return flops;
}

/*
  get even-even solutions
  still need qinv buffers
  be carefull to check the inversion with single precision input
  not used ?
*/
template<typename Ty>
void get_staggered_multishift(quda_inverter& qinv,
  std::vector<qlat::FieldM<Ty, 3> >& res,
  std::vector<qlat::FieldM<Ty, 3> >& src,
  const std::vector<double >& masses2,
  const double err, const Int niter)
{
  TIMER("get_staggered_multishift");
  Qassert(src.size() != 0);

  //const Int DIM = 3;
  const Geometry& geo = qinv.geoB();
  qlat::vector<Long >& map = qinv.map_index;

  const Int Nsrc = src.size();
  const Int Nres = res.size();
  const Int multishift = masses2.size();

  Qassert(Nres  == Nsrc * multishift);
  for(unsigned int i=0;i<src.size();i++){Qassert(src[i].initialized);}
  for(unsigned int i=0;i<res.size();i++){Qassert(res[i].initialized);}
  const Geometry& geo_src = src[0].geo();
  Qassert(geo_src.node_site == geo.node_site and geo_src.geon.size_node == geo.geon.size_node);

  ///qinv.setup_mat_mass(0.0);
  //qinv.setup_CG();
  ////reload QUDA links
  //qinv.setup_inv_param_prec(prec_type);

  //QudaInvertParam  inv_param_copy = qinv.inv_param;

  //std::vector<quda::ColorSpinorField> out_multishift_gpu(multishift);

  ////const Int Cdim = ref.CompositeDim();
  quda::ColorSpinorField& ref = *qinv.gsrc;
  const size_t Lsize = ref.TotalBytes();

  //int buf_prec = -1;
  //if(qinv.inv_param.cuda_prec == QUDA_DOUBLE_PRECISION){buf_prec = 0;}
  //if(qinv.inv_param.cuda_prec == QUDA_SINGLE_PRECISION){buf_prec = 1;}
  //if(buf_prec == -1){abort_r("QUDA buffer prec not supported yet!");}

  if(get_data_type_is_double<Ty >()){
    Qassert(qinv.inv_param.cuda_prec == QUDA_DOUBLE_PRECISION);
  }else{
    Qassert(qinv.inv_param.cuda_prec == QUDA_SINGLE_PRECISION);
  }

  std::vector<quda::ColorSpinorField>& out_multishift_gpu = qinv.out_multishift_gpu;
  std::vector<quda::ColorSpinorField>& src_multishift_gpu = qinv.src_multishift_gpu;

  bool init = false;
  if(out_multishift_gpu.size() != multishift){init = true;}
  if(src_multishift_gpu.size() != 1){init = true;}
  if(init == false){
    for(unsigned int i=0;i<out_multishift_gpu.size();i++){
      if(out_multishift_gpu[i].TotalBytes() != Lsize){init = true;}
    }
    for(unsigned int i=0;i<src_multishift_gpu.size();i++){
      if(src_multishift_gpu[i].TotalBytes() != Lsize){init = true;}
    }
  }
  if(init == true){
    quda::ColorSpinorParam cs_param = quda::ColorSpinorParam(ref);
    cs_param.is_composite  = false;cs_param.is_component  = false;
    cs_param.location = QUDA_CUDA_FIELD_LOCATION;
    out_multishift_gpu.resize(multishift);
    src_multishift_gpu.resize(1);
    for(unsigned int i=0;i<out_multishift_gpu.size();i++){
      out_multishift_gpu[i] = quda::ColorSpinorField(cs_param);
    }
    for(unsigned int i=0;i<src_multishift_gpu.size();i++){
      src_multishift_gpu[i] = quda::ColorSpinorField(cs_param);
    }
  }

  QudaInvertParam inv_param = qinv.inv_param;
  inv_param.secs = 0;
  inv_param.gflops = 0;
  inv_param.iter = 0;

  inv_param.tol = err; 
  inv_param.maxiter = niter;

  ////check whether multishift is needed
  if(multishift == 1){
    inv_param.mass = std::sqrt(masses2[0] / 4.0);
  }else{
    inv_param.num_offset = multishift;
    inv_param.compute_true_res = 0;
    //inv_param.mass = masses[0]; 
    for (Int i = 0; i < multishift; i++) {
      // Set masses and offsets
      //inv_param.offset[i] = 4 * masses[i] * masses[i] - 4 * masses[0] * masses[0];
      inv_param.offset[i] = masses2[i];
      inv_param.tol_offset[i]    = inv_param.tol;
      inv_param.tol_hq_offset[i] = 1e-5; ///donot check heavy quark, tol_hq
    }
  }

  inv_param.solution_type = QUDA_MATPC_SOLUTION;
  inv_param.solve_type    = QUDA_DIRECT_PC_SOLVE;
  inv_param.input_location  = QUDA_CUDA_FIELD_LOCATION;
  inv_param.output_location = QUDA_CUDA_FIELD_LOCATION;
  inv_param.use_init_guess = QUDA_USE_INIT_GUESS_NO;

  std::vector<void *> _hp_multi_x(multishift);
  for(Int di=0;di<multishift;di++)
  {
    _hp_multi_x[di] = out_multishift_gpu[di].Even().data();
  }

  for(Int srci=0;srci<Nsrc;srci++)
  {
    qlat_cf_to_quda_cf(src_multishift_gpu[0], (Ty*) qlat::get_data(src[srci]).data(), geo, map);
    //{
    //double na = quda::blas::norm2(src_multishift_gpu[0].Even());
    //qmessage("src %5d, %+.8e, offset %.8e, err %.3e, iter %8d \n", srci, na, inv_param.offset[0], inv_param.tol, inv_param.maxiter);
    //}

    if(multishift == 1){
      invertQuda(_hp_multi_x[0],               src_multishift_gpu[0].Even().data(), &inv_param);
    }
    else{
      invertMultiShiftQuda(_hp_multi_x.data(), src_multishift_gpu[0].Even().data(), &inv_param);
    }

    for(Int mi = 0; mi < multishift; mi++)
    {
      quda_cf_to_qlat_cf((Ty*) qlat::get_data(res[mi*Nsrc + srci]).data(), out_multishift_gpu[mi], geo, map);
    }
  }

  ////restore default inv_param
  //inv_param = inv_param_copy;
  if(qinv.quda_verbos >= 0 or qinv.quda_verbos == -2)
  {
    qmessage("Done multishift: %8d iter / %.6f secs = %.3f Gflops %.3f Gflops/GPU, Cost %.3f Gflops \n",
          inv_param.iter, inv_param.secs, inv_param.gflops / inv_param.secs,
          inv_param.gflops / (inv_param.secs * qlat::get_num_node()), inv_param.gflops);
  }
  qinv.inv_param.secs   = inv_param.secs;
  qinv.inv_param.gflops = inv_param.gflops;
  qinv.inv_param.iter   = inv_param.iter;
}

////get even-even solutions
// minimum qinv gpu memeory usage 
// GPU memory input and output
template<typename Ty>
void get_staggered_multishift_even(quda_inverter& qinv,
  std::vector<Ty* > res, std::vector<Ty* > src,
  const std::vector<double >& masses2,
  const double err, const Int niter)
{
  TIMER_FLOPS("get_staggered_multishift_even");
  Qassert(src.size() != 0);
  if(qinv.quda_verbos >= 0){print_mem_info();}

  //const Int DIM = 3;

  const Int Nsrc = src.size();
  const Int Nres = res.size();
  const Int multishift = masses2.size();

  Qassert(Nres  == Nsrc * multishift);
  //for(unsigned int i=0;i<src.size();i++){Qassert(src[i].initialized);}
  //for(unsigned int i=0;i<res.size();i++){Qassert(res[i].initialized);}
  //const Geometry& geo_src = src[0].geo();
  //Qassert(geo_src.node_site == geo.node_site and geo_src.geon.size_node == geo.geon.size_node);

  ///qinv.setup_mat_mass(0.0);
  //qinv.setup_CG();
  ////reload QUDA links
  //qinv.setup_inv_param_prec(prec_type);

  //QudaInvertParam  inv_param_copy = qinv.inv_param;

  //std::vector<quda::ColorSpinorField> out_multishift_gpu(multishift);

  //quda::ColorSpinorField& ref = *qinv.gsrc;
  //const Int Cdim = ref.CompositeDim();
  //const size_t Lsize = ref.TotalBytes() / Cdim;

  //int buf_prec = -1;
  //if(qinv.inv_param.cuda_prec == QUDA_DOUBLE_PRECISION){buf_prec = 0;}
  //if(qinv.inv_param.cuda_prec == QUDA_SINGLE_PRECISION){buf_prec = 1;}
  //if(buf_prec == -1){abort_r("QUDA buffer prec not supported yet!");}
  //Qassert(qinv.inv_param.cuda_prec == QUDA_DOUBLE_PRECISION or qinv.inv_param.cuda_prec == QUDA_SINGLE_PRECISION);

  if(get_data_type_is_double<Ty >()){
    Qassert(qinv.inv_param.cuda_prec == QUDA_DOUBLE_PRECISION);
  }else{
    Qassert(qinv.inv_param.cuda_prec == QUDA_SINGLE_PRECISION);
    Qassert(false);// default inverter could not copy single source correcttly
  }

  QudaInvertParam  inv_param_bak = qinv.inv_param;
  QudaInvertParam& inv_param     = qinv.inv_param;
  inv_param.secs = 0;
  inv_param.gflops = 0;
  inv_param.iter = 0;
  inv_param.use_init_guess = QUDA_USE_INIT_GUESS_NO;

  inv_param.tol = err; 
  inv_param.maxiter = niter;

  ////check whether multishift is needed
  if(multishift == 1){
    inv_param.mass = std::sqrt(masses2[0] / 4.0);
  }else{
    inv_param.num_offset = multishift;

    //temporary hack for refinement
    Int qlat_multi_refine = 1;
    std::string val = qlat::get_env(std::string("qlat_multi_refine"));
    if(val != ""){qlat_multi_refine = stringtonum(val);}

    inv_param.compute_true_res = qlat_multi_refine;
    //inv_param.mass = masses[0]; 
    double max_tol = 1e-13;
    //inv_param.cuda_prec_sloppy              = QUDA_DOUBLE_PRECISION;
    if(inv_param.cuda_prec_sloppy == QUDA_SINGLE_PRECISION and inv_param.cuda_prec == QUDA_DOUBLE_PRECISION){
      max_tol = 1e-8;
    }
    if(inv_param.cuda_prec_sloppy == QUDA_SINGLE_PRECISION and inv_param.cuda_prec == QUDA_SINGLE_PRECISION){
      max_tol = 1e-6;
    }

    double gap = 0.0;
    const double tol  = inv_param.tol;
    const double Etol = std::log10(tol);

    //treak the inv precision for double dslash
    if(tol > max_tol and tol > 0){
      double Emax = std::log10(max_tol);
      gap = (Etol - Emax)/ (multishift - 1.0);
    }

    //(inv_param.tol - max_tol) / std::sqrt(multishift - 1.0);
    //if(gap < 0 or inv_param.tol < 1e-11){gap = 0.0;}
    //if(qlat_multi_refine == 0){gap = 0.0;}//to zero if no refinement

    for (Int i = 0; i < multishift; i++) {
      // Set masses and offsets
      //inv_param.offset[i] = 4 * masses[i] * masses[i] - 4 * masses[0] * masses[0];
      inv_param.offset[i] = masses2[i];
      inv_param.tol_offset[i]    = std::pow(10, Etol - gap * i);

      inv_param.tol_hq_offset[i] = 5e-8; ///donot check heavy quark, tol_hq, not used currently
    }
  }

  inv_param.solution_type = QUDA_MATPC_SOLUTION;
  inv_param.solve_type    = QUDA_DIRECT_PC_SOLVE;
  inv_param.input_location  = QUDA_CUDA_FIELD_LOCATION;
  inv_param.output_location = QUDA_CUDA_FIELD_LOCATION;
  inv_param.use_init_guess = QUDA_USE_INIT_GUESS_NO;

  std::vector<void *> _hp_multi_x(multishift);

  if(multishift == 1)
  {
    qlat::vector<void *> srcP;srcP.resize(Nsrc);
    qlat::vector<void *> resP;resP.resize(Nsrc);
    for(Int srci=0;srci<Nsrc;srci++)
    {
      srcP[srci] = src[srci];
      resP[srci] = res[srci];
    }
    qinv.callMultiSrcQuda(resP, srcP);
  }
  else{
    for(Int srci=0;srci<Nsrc;srci++)
    {
      for(Int mi=0;mi<multishift;mi++)
      {
        _hp_multi_x[mi] = res[mi*Nsrc + srci];
      }
      invertMultiShiftQuda(_hp_multi_x.data(), src[srci], &inv_param);
    }
  }

  ////restore default inv_param
  if(qinv.quda_verbos >= 0 or qinv.quda_verbos == -2)
  {
    qmessage("Done multishift: %8d iter / %.6f secs = %.3f Gflops %.3f Gflops/GPU, Cost %.3f Gflops \n",
          inv_param.iter, inv_param.secs, inv_param.gflops / inv_param.secs,
          inv_param.gflops / (inv_param.secs * qlat::get_num_node()), inv_param.gflops);
  }
  qinv.inv_param =  inv_param_bak;
  qinv.inv_param.secs   = inv_param.secs;
  qinv.inv_param.gflops = inv_param.gflops;
  qinv.inv_param.iter   = inv_param.iter;
  timer.flops += inv_param.gflops * 1024 * 1024 * 1024 / qlat::get_num_node();
}

// Int gf_gauge_dir 3=Coulomb, 4=Landau
inline void quda_inverter::gauge_fix(qlat::ComplexD* quda_gf, Int gf_gauge_dir, Int gf_maxiter, double gf_tolerance, Int fix_type)
{
  TIMER("quda gauge fix");
  /////load gauge to quda GPU default position
  const double gf_ovr_relaxation_boost = 1.5;
  const Int gf_reunit_interval = 10;
  const bool gf_theta_condition = false;
  const Int gf_verbosity_interval = 100;
  const double gf_fft_alpha = 0.8;
  const bool gf_fft_autotune = false;

  QudaGaugeParam  copy_param= gauge_param;
  copy_param.reconstruct                   = QUDA_RECONSTRUCT_NO; 
  copy_param.reconstruct_sloppy            = QUDA_RECONSTRUCT_NO;
  copy_param.t_boundary = QUDA_PERIODIC_T;
  //copy_param.reconstruct                   = QUDA_RECONSTRUCT_12;
  //copy_param.reconstruct_sloppy            = QUDA_RECONSTRUCT_12;

  quda::GaugeFieldParam cpu_param(copy_param, quda_gf);
  cpu_param.location = QUDA_CPU_FIELD_LOCATION;
  if(cpu_param.order <= 4) cpu_param.ghostExchange = QUDA_GHOST_EXCHANGE_NO;
  quda::GaugeField *in = quda::GaugeField::Create(cpu_param);

  quda::GaugeFieldParam gpu_param(copy_param, quda_gf);
  // switch the parameters for creating the mirror precise cuda gauge field
  gpu_param.location = QUDA_CUDA_FIELD_LOCATION;
  gpu_param.create = QUDA_NULL_FIELD_CREATE;
  gpu_param.reconstruct = copy_param.reconstruct;
  gpu_param.setPrecision(copy_param.cuda_prec, true);
  //if(gpu_param.order <= 4) gpu_param.ghostExchange = QUDA_GHOST_EXCHANGE_NO;
  //gpu_param.ghostExchange = QUDA_GHOST_EXCHANGE_NO;
  //gpu_param.ghostExchange = QUDA_GHOST_EXCHANGE_EXTENDED;
  //gpu_param.pad = copy_param.ga_pad;

  //gpu_param.ghostExchange = QUDA_GHOST_EXCHANGE_PAD;
  gpu_param.pad = copy_param.ga_pad;

  quda::GaugeField* U = NULL;

  // need extended field
  if (quda::comm_partitioned()) {
    gpu_param.ghostExchange = QUDA_GHOST_EXCHANGE_PAD;
    //gpu_param.ghostExchange = QUDA_GHOST_EXCHANGE_NO;
    //quda::GaugeField *tmp = new quda::GaugeField(gpu_param);
    //tmp->copy(*in);

    //quda::GaugeFieldParam gpu_param(copy_param, quda_gf);
    //U = new quda::GaugeField(gpu_param);
    //copyExtendedGauge(*U, *tmp, QUDA_CUDA_FIELD_LOCATION);
    //delete tmp;
    quda::GaugeField* tmp = new quda::GaugeField(gpu_param);
    tmp->copy(*in);
    //*tmp = *in;

    quda::lat_dim_t R = {0, 0, 0, 0};
    for (Int d = 0; d < 4; d++)
      if (quda::comm_dim_partitioned(d)) R[d] = 2;
    static quda::TimeProfile GaugeFix("GaugeFix");
    U = quda::createExtendedGauge(*tmp, R, GaugeFix);
    delete tmp;
  } else {
    gpu_param.ghostExchange = QUDA_GHOST_EXCHANGE_NO;
    U = new quda::GaugeField(gpu_param);
    U->copy(*in);
    //*U = *in;
  }

  // Reunitarization
  // Int num_failures_h = 0;
  // quda::unitarizeLinks(*U, &num_failures_h);
  // qudaDeviceSynchronize();
  // if (num_failures_h > 0) errorQuda("Error in the unitarization (%d errors)", num_failures_h);

  if(fix_type == 0){
    quda::gaugeFixingOVR(*U, gf_gauge_dir, gf_maxiter, gf_verbosity_interval, gf_ovr_relaxation_boost, gf_tolerance, gf_reunit_interval, gf_theta_condition);
  }

  if(fix_type == 1){
    quda::gaugeFixingFFT(*U, gf_gauge_dir, gf_maxiter, gf_verbosity_interval, gf_fft_alpha, gf_fft_autotune, gf_tolerance, gf_theta_condition);
  }

  {
    if (quda::comm_partitioned()) {
      gpu_param.ghostExchange = QUDA_GHOST_EXCHANGE_PAD;
    }else{
      gpu_param.ghostExchange = QUDA_GHOST_EXCHANGE_NO;
    }
    quda::GaugeField *tmp = new quda::GaugeField(gpu_param);

    tmp->copy(*U);//extend copy on GPU
    in->copy(*tmp);

    delete tmp;
    ///in->copy(*U);
  }
  //*in = *U;
  delete in;
  delete  U;
}


}  // namespace qlat



#endif
