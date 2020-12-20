#pragma once

#include <qlat/qlat.h>

#include <quda.h>
#include <invert_quda.h>

#include <cstdlib>

namespace qlat
{  //

static int mpi_rank_from_coords_x(const int* coords, void* fdata)
{
  int* dims = reinterpret_cast<int*>(fdata);
  //
  int rank;
  rank = coords[3];
  for (int i = 2; i >= 0; i--) {
    rank = dims[i] * rank + coords[i];
  }
  return rank;
}
//
static int mpi_rank_from_coords_t(const int* coords, void* fdata)
{
  int* dims = reinterpret_cast<int*>(fdata);
  //
  int rank;
  rank = coords[0];
  for (int i = 1; i <= 3; i++) {
    rank = dims[i] * rank + coords[i];
  }
  return rank;
}
//
inline void quda_begin(int mpi_layout[4], bool t = false)
{
  using namespace quda;
  // The following sets the MPI comm stuff.
  if (t) {
    initCommsGridQuda(4, mpi_layout, mpi_rank_from_coords_t,
                      reinterpret_cast<void*>(mpi_layout));
  } else {
    initCommsGridQuda(4, mpi_layout, mpi_rank_from_coords_x,
                      reinterpret_cast<void*>(mpi_layout));
  }
  // comm_set_gridsize(mpi_layout);
  initQuda(-1000);
  printf(
      "initialized on quda rank #%03d (%03d,%03d,%03d,%03d), qlat rank #%03d "
      "(%03d,%03d,%03d,%03d).\n",
      comm_rank(), comm_coord(0), comm_coord(1), comm_coord(2), comm_coord(3),
      get_id_node(), get_coor_node()[0], get_coor_node()[1], get_coor_node()[2],
      get_coor_node()[3]);
  // Make sure there is no mismatch
  qassert(comm_rank() == get_id_node());
  for (int d = 0; d < 4; d++) {
    qassert(comm_coord(d) == get_coor_node()[d]);
  }
}

inline void quda_end()
{
  using namespace quda;
  endQuda();
}

template <class T>
void quda_convert_gauge(std::vector<T>& qgf, const GaugeField& gf)
{
  TIMER("quda_convert_gauge(qgf,gf)");
  const Geometry& geo = gf.geo();
  ColorMatrix* quda_pt = reinterpret_cast<ColorMatrix*>(qgf.data());
  qassert(geo.multiplicity == 4);
  long V = geo.local_volume();
  long Vh = V / 2;
  for (int qlat_idx = 0; qlat_idx < V; qlat_idx++) {
    Coordinate xl = geo.coordinate_from_index(qlat_idx);
    const Vector<ColorMatrix> ms = gf.get_elems_const(xl);
    int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    for (int mu = 0; mu < 4; mu++) {
      int quda_idx = (qlat_idx / 2 + eo * Vh) * 4 + mu;
      quda_pt[quda_idx] = ms[mu];
    }
  }
}

template <class T>
void quda_convert_fermion(FermionField5d& ff, const std::vector<T>& qff)
{
  TIMER("quda_convert_fermion(ff,qff)");
  const Geometry& geo = ff.geo();
  const WilsonVector* quda_pt =
      reinterpret_cast<const WilsonVector*>(qff.data());
  int Ls = geo.multiplicity;
  qassert(Ls > 0);
  long V = geo.local_volume();
  long Vh = V / 2;
//
#pragma omp parallel for
  for (long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
    const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
    int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    Vector<WilsonVector> wvs = ff.get_elems(xl);
    for (int s = 0; s < Ls; s++) {
      int quda_idx = eo * Vh * Ls + s * Vh + qlat_idx_4d / 2;
      wvs[s] = quda_pt[quda_idx];
    }
  }
}

template <class T>
void quda_convert_fermion(std::vector<T>& qff, const FermionField5d& ff)
{
  TIMER("quda_convert_fermion(qff,ff)");
  const Geometry& geo = ff.geo();
  WilsonVector* quda_pt = reinterpret_cast<WilsonVector*>(qff.data());
  int Ls = geo.multiplicity;
  qassert(Ls > 0);
  long V = geo.local_volume();
  long Vh = V / 2;
//
#pragma omp parallel for
  for (long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
    const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
    int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    const Vector<WilsonVector> wvs = ff.get_elems_const(xl);
    for (int s = 0; s < Ls; s++) {
      int quda_idx = eo * Vh * Ls + s * Vh + qlat_idx_4d / 2;
      quda_pt[quda_idx] = wvs[s];
    }
  }
}

QudaPrecision get_quda_precision(int byte)
{
  switch (byte) {
    case 8:
      return QUDA_DOUBLE_PRECISION;
      break;
    case 4:
      return QUDA_SINGLE_PRECISION;
      break;
    case 2:
      return QUDA_HALF_PRECISION;
      break;
    default:
      qassert(false);
      return QUDA_INVALID_PRECISION;
  }
}

void set_deflation_param(QudaEigParam& df_param)
{
  df_param.import_vectors = QUDA_BOOLEAN_NO;
  df_param.run_verify = QUDA_BOOLEAN_NO;

  df_param.nk = df_param.invert_param->nev;
  df_param.np =
      df_param.invert_param->nev * df_param.invert_param->deflation_grid;
  df_param.extlib_type = QUDA_EIGEN_EXTLIB;

  df_param.cuda_prec_ritz = df_param.invert_param->cuda_prec;
  df_param.location = QUDA_CUDA_FIELD_LOCATION;
  df_param.mem_type_ritz = QUDA_MEMORY_DEVICE;

  // set file i/o parameters
  strcpy(df_param.vec_infile, "/ccs/home/jiquntu/meas-qlat/vec_outfile.txt");
  strcpy(df_param.vec_outfile, "/ccs/home/jiquntu/meas-qlat/vec_outfile.txt");
}

void set_gauge_param(QudaGaugeParam& gauge_param, const Geometry& geo,
                     const InverterParams& ip)
{
  //
  for (int mu = 0; mu < 4; mu++) {
    gauge_param.X[mu] = geo.node_site[mu];
  }
  //
  // ... OK. I don't know what this means
  gauge_param.type = QUDA_WILSON_LINKS;
  //
  // Slowest changing to fastest changing: even-odd, mu, x_cb_4d, row, column,
  // complex See the code later in this file to see the conversion between
  // Grid inde and Quda index.
  gauge_param.gauge_order = QUDA_MILC_GAUGE_ORDER;
  //
  // The precision used here should be the same as those set in the inv_param,
  // i.e. gauge_param.cuda_prec = inv_param.cuda_prec
  // gauge_param.cuda_prec_sloppy = inv_param.cuda_prec_sloppy
  gauge_param.cpu_prec = get_quda_precision(ip.higher_precision);
  gauge_param.cuda_prec = get_quda_precision(ip.higher_precision);
  gauge_param.reconstruct = QUDA_RECONSTRUCT_NO;
  gauge_param.cuda_prec_sloppy = get_quda_precision(ip.lower_precision);
  gauge_param.cuda_prec_precondition = get_quda_precision(ip.lower_precision);
  gauge_param.cuda_prec_refinement_sloppy =
      get_quda_precision(ip.lower_precision);
  gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_NO;
  //
  gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;
  //
  gauge_param.anisotropy = 1.0;
  gauge_param.t_boundary = QUDA_PERIODIC_T;
  //
  int x_face_size = gauge_param.X[1] * gauge_param.X[2] * gauge_param.X[3] / 2;
  int y_face_size = gauge_param.X[0] * gauge_param.X[2] * gauge_param.X[3] / 2;
  int z_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[3] / 2;
  int t_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[2] / 2;
  int pad_size = std::max(x_face_size, y_face_size);
  pad_size = std::max(pad_size, z_face_size);
  pad_size = std::max(pad_size, t_face_size);
  gauge_param.ga_pad = pad_size;
}

void set_inv_param(QudaInvertParam& inv_param, const FermionAction& fa,
                   const InverterParams& ip)
{
  inv_param.Ls = fa.ls;
  inv_param.dslash_type = QUDA_MOBIUS_DWF_DSLASH;
  inv_param.mass = fa.mass;
  // Note that Quda uses -M5 as M5 ...
  inv_param.m5 = -fa.m5;
  if (fa.is_using_zmobius) {
    static_assert(sizeof(__complex__ double) == sizeof(Complex));
    memcpy(inv_param.b_5, fa.bs.data(), fa.ls * sizeof(Complex));
    memcpy(inv_param.c_5, fa.cs.data(), fa.ls * sizeof(Complex));
  } else {
    for (int s = 0; s < fa.ls; s++) {
      inv_param.b_5[s] = 0.5 * fa.mobius_scale + 0.5;
      inv_param.c_5[s] = 0.5 * fa.mobius_scale - 0.5;
    }
  }
  // kappa is irrelevant for Mobius/DWF but you have to set it.
  inv_param.kappa = 1. / (2. * (1. + 3. / 1. + fa.mass));
  inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;
  inv_param.solver_normalization = QUDA_DEFAULT_NORMALIZATION;
  //
  // Whether or not content of your input void* pointer will be modified
  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
  //
  // I don't know what these are but you have to set them.
  inv_param.use_sloppy_partial_accumulator = 0;
  inv_param.solution_accumulator_pipeline = 1;
  //
  // This is for the reliable update. Just set it to some large number.
  inv_param.max_res_increase = 20000;
  inv_param.max_res_increase_total = 200000;
  inv_param.maxiter_precondition = 6;
  //
  inv_param.mq1 = fa.mass;
  inv_param.mq2 = fa.mass;
  inv_param.mq3 = 0.01;
  inv_param.eofa_shift = +0.0;
  inv_param.eofa_pm = 1;
  //
  // The solver tolerance, i.e. |MdagM * x - b| < tol * |b|
  inv_param.tol = ip.stop_rsd;
  inv_param.tol_restart = 0.0005;
  //
  // The maximum number of iterations.
  inv_param.maxiter = ip.max_num_iter;
  //
  // This is for Quda's sophisticated reliable update. 0.1 should be good.
  inv_param.reliable_delta = 0.1;
  inv_param.use_alternative_reliable = true;
  //
  // NORMOP_PC means preconditioned normal operator MdagM
  inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
  //
  // QUDA_MATPC_EVEN_EVEN means we solve on even sites and use symmetric
  // preconditioning The other options are: QUDA_MATPC_ODD_ODD,
  // QUDA_MATPC_EVEN_EVEN_ASYMMETRIC,
  // QUDA_MATPC_ODD_ODD_ASYMMETRIC,
  //
  // There might be a performance difference.
  inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
  //
  // Eventually we want the unpreconditioned solution.
  inv_param.solution_type = QUDA_MAT_SOLUTION;
  //
  inv_param.dagger = QUDA_DAG_NO;
  //
  // The precision used to correct the inner solver.
  inv_param.cpu_prec = get_quda_precision(ip.higher_precision);
  ;
  inv_param.cuda_prec = get_quda_precision(ip.higher_precision);
  ;
  // The sloppy(inner) solver precision
  inv_param.cuda_prec_sloppy = get_quda_precision(ip.lower_precision);
  inv_param.cuda_prec_precondition = get_quda_precision(ip.lower_precision);
  inv_param.cuda_prec_refinement_sloppy =
      get_quda_precision(ip.lower_precision);
  inv_param.cuda_prec_ritz = get_quda_precision(ip.lower_precision);
  //
  inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
  inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
  //
  // I don't know what these are but you have to set them.
  inv_param.sp_pad = 0;
  inv_param.cl_pad = 0;
  //
  // Both CPS and Grid use this gamma matrix representation
  inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  //
  // Slowest changing to fastest changing: even-odd, Ls, x_cb_4d, spin, color,
  // complex See the code later in this file to see the conversion between
  // Grid inde and Quda index.
  inv_param.dirac_order = QUDA_DIRAC_ORDER;
  //
  // QUDA_DEBUG_VERBOSE is too nasty.
  inv_param.verbosity = QUDA_VERBOSE;
  //
  // It seems the initial value of this is undefined so it's better to set it
  // here. If set to QUDA_USE_INIT_GUESS_NO Quda will zero the input solution
  // pointer before the solve.
  inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;
  //
  if (ip.solver_type == 0) {
    inv_param.inv_type = QUDA_CG_INVERTER;
  } else if (ip.solver_type == 1) {
    inv_param.inv_type = QUDA_INC_EIGCG_INVERTER;
    inv_param.nev = 16;
    inv_param.deflation_grid = 8;
    inv_param.rhs_idx = 0;
    inv_param.eigcg_max_restarts = 3;
    inv_param.max_restart_num = 3;
  } else if (ip.solver_type == 2) {
    inv_param.inv_type = QUDA_MSPCG_INVERTER;
  } else {
    qassert(false);
  }
}

struct InverterDomainWallQuda : InverterDomainWall {
  // Now setup all the QUDA parameters
  QudaGaugeParam gauge_param;
  QudaInvertParam inv_param;
  // static void* df_preconditioner;
  // static QudaEigParam df_param;
  //
  std::vector<double> qgf;
  //
  bool qlat_check = true;
  //
  InverterDomainWallQuda() : qgf(0) { init(); }
  ~InverterDomainWallQuda() { init(); }
  //
  void init()
  {
    free();
    InverterDomainWall::init();
  }
  //
  void load_gauge()
  {  // Load gauge field to Quda.
    // initialize the std::vectors that holds the gauge field.
    TIMER("InvDWQuda::load_gauge()");
    size_t qgf_size = geo.local_volume() * 4 * 18;
    qgf.resize(qgf_size);
    quda_convert_gauge(qgf, this->gf);
    loadGaugeQuda((void*)qgf.data(), &gauge_param);
    double plaq[3];
    plaqQuda(plaq);
    printfQuda(
        "Computed plaquette is %16.12e (spatial = %16.12e, temporal = "
        "%16.12e)\n",
        plaq[0], plaq[1], plaq[2]);
  }
  //
  void setup()
  {
    TIMER("InvDWQuda::setup()");
    using namespace quda;
    free();
    //
    gauge_param = newQudaGaugeParam();
    inv_param = newQudaInvertParam();
    // Now setup all the QUDA parameters
    set_gauge_param(gauge_param, geo, ip);
    //
    load_gauge();
    //
    set_inv_param(inv_param, fa, ip);
    //
    printQudaInvertParam(&inv_param);
  }
  //
  void setup(const GaugeField& gf_, const FermionAction& fa_)
  {
    InverterDomainWall::setup(gf_, fa_);
    setup();
  }
  void setup(const GaugeField& gf_, const FermionAction& fa_,
             const InverterParams& ip_)
  {
    InverterDomainWall::setup(gf_, fa_, ip_);
    setup();
  }
  void setup(const GaugeField& gf_, const FermionAction& fa_,
             const LowModes& lm_)
  {
    InverterDomainWall::setup(gf_, fa_, lm_);
    setup();
  }
  //
  void free() {}
  //
  void allocate_eigcg()
  {
    TIMER("InvDWQuda::allocate_eigcg()");
    QudaEigParam* df_param = get_df_param();
    df_param->invert_param = &inv_param;
    set_deflation_param(*df_param);
    inv_param.deflation_op = newDeflationQuda(df_param);
  }
  //
  void deallocate_eigcg()
  {
    TIMER("InvDWQuda::deallocate_eigcg()");
    if (inv_param.deflation_op) {
      destroyDeflationQuda(inv_param.deflation_op);
    } else {
      Printf("inv_param.deflation_op is null.");
    }
  }

 private:
  QudaEigParam* get_df_param()
  {
    // Ugly trick to have a static object in a header only library.
    // This is potentially not thread safe but who cares ...
    static QudaEigParam df_param_ = newQudaEigParam();
    return &df_param_;
  }
};

inline void setup_inverter(InverterDomainWallQuda& inv) { inv.setup(); }

inline void setup_inverter(InverterDomainWallQuda& inv, const GaugeField& gf,
                           const FermionAction& fa)
{
  inv.setup(gf, fa);
}

inline void setup_inverter(InverterDomainWallQuda& inv, const GaugeField& gf,
                           const FermionAction& fa, const InverterParams& ip)
{
  inv.setup(gf, fa, ip);
}

inline void setup_inverter(InverterDomainWallQuda& inv, const GaugeField& gf,
                           const FermionAction& fa, const LowModes& lm)
{
  qassert(false);  // Quda inverter does NOT support deflation, yet.
  // inv.setup(gf, fa, lm);
}

inline void invert(FermionField5d& sol, const FermionField5d& src,
                   const InverterDomainWallQuda& inv)
{
  // initialize the std::vectors that hold source and solution vectors.
  size_t qff_size = inv.geo().local_volume() * inv.fa.ls * 24;
  std::vector<double> qff_src(qff_size);
  std::vector<double> qff_sol(qff_size);
  // inv_param_dup.deflation_op = inv.df_preconditioner;
  // Quda does not have D_minus built in.
  FermionField5d dm_in, check;
  if (inv.fa.is_multiplying_dminus) {
    multiply_d_minus(dm_in, src, inv);
  } else {
    dm_in = src;
  }
  // inverse_with_cg(sol, src, inv, cg_with_herm_sym_2);
  Printf("Input  5d vector norm2 = %16.12e.\n", qnorm(dm_in));
  quda_convert_fermion(qff_src, dm_in);
  // Make a copy of the inv_param
  QudaInvertParam inv_param_dup = inv.inv_param;
  // Perform the actual inversion
  invertQuda(qff_sol.data(), qff_src.data(), &inv_param_dup);
  // Printf("inv_param rhs index    = %03d.\n",
  // InverterDomainWallQuda::inv_param.rhs_idx);
  quda_convert_fermion(sol, qff_sol);
  // The difference between Quda and CPS/Grid
  sol *= 1. / ((0.5 * inv.fa.mobius_scale + 0.5) * (4. - inv.fa.m5) + 1.);
  Printf("Output 5d vector norm2 = %16.12e.\n", qnorm(sol));
  if (inv.qlat_check) {
    check.init(geo_resize(src.geo()));
    multiply_m_full(check, sol, inv);
    check -= dm_in;
    Printf("Check  5d vector norm2 = %16.12e.\n", qnorm(check));
  }
}

inline void invert(FermionField4d& sol, const FermionField4d& src,
                   const InverterDomainWallQuda& inv)
{
  TIMER_VERBOSE("invert(sol4d,src4d,inv_quda)");
  invert_dwf(sol, src, inv);
}

}  // namespace qlat
