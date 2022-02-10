#ifndef UTILS_QUDA_INVERTER_H
#define UTILS_QUDA_INVERTER_H

#pragma once

#include <qlat/qlat.h>

#include <quda.h>
#include <invert_quda.h>

#include <cstdlib>
#include "float_type.h"
#include "quda_para.h"

static TimeProfile profileEigensolve("eigensolveQuda");

namespace qlat
{  //

struct quda_inverter {
  long V;
  int X[4];
  QudaGaugeParam  gauge_param;
  QudaInvertParam inv_param;

  QudaEigParam    eig_param = newQudaEigParam();
  int nvec;
  ////void **host_evecs = NULL;
  ////double _Complex *host_evals = NULL;

  std::vector<std::vector<quda::Complex > > evecs;
  std::vector<quda::Complex > evals;

  //void *eig_preconditioner  = NULL;
  EigenSolver *eig_solve;
  DiracMatrix* mat;
  DiracMatrix* mat_Mdag;
  DiracMatrix* mat_MMdag;
  DiracMatrix* mat_MdagM;

  Dirac *dirac;
  Dirac *dSloppy;
  Dirac *dPre;
  Dirac *dEig;
  std::vector<ColorSpinorField *> kSpace;
  std::vector<ColorSpinorField *> evecs_;


  quda::ColorSpinorParam cs_cpu;
  quda::ColorSpinorParam cs_gpu;

  quda::ColorSpinorField *csrc, *cres;
  quda::ColorSpinorField *gsrc, *gres;

  quda::ColorSpinorField *ctmp0, *ctmp1, *ctmp2;
  quda::ColorSpinorField *gtmp0, *gtmp1, *gtmp2;

  ////0 for wilson, clover, 1 for stagger
  int fermion_type;
  int use_eigen_pc;
  int check_residue;
  int apply_stag_phase;

  qlat::vector<qlat::Complex > quda_clover    ;
  qlat::vector<qlat::Complex > quda_clover_inv;


  /////load gauge and set default parameters
  quda_inverter(const Geometry& geo, qlat::vector<qlat::Complex >& quda_gf, int apply_stag=0)
  {
    /////set up gauge parameters
    V = geo.local_volume();
    qassert(quda_gf.size() == V *4*3*3);
    apply_stag_phase = apply_stag;

    ////===Start of gauge_param
    gauge_param = newQudaGaugeParam();
    inv_param   = newQudaInvertParam();


    for (int mu = 0; mu < 4; mu++) {gauge_param.X[mu] = geo.node_site[mu];X[mu] = geo.node_site[mu];}


    gauge_param.type = QUDA_WILSON_LINKS; //// or QUDA_SU3_LINKS
    // Slowest changing to fastest changing: even-odd, mu, x_cb_4d, row, column,
    // complex See the code later in this file to see the conversion between
    // Grid inde and Quda index.
    gauge_param.gauge_order = QUDA_MILC_GAUGE_ORDER;

    gauge_param.t_boundary = QUDA_PERIODIC_T; ////fermion_t_boundary

    ////fine tune the precisions if needed
    gauge_param.cpu_prec               = QUDA_DOUBLE_PRECISION;
    gauge_param.cuda_prec              = QUDA_DOUBLE_PRECISION;

    gauge_param.cuda_prec_sloppy       = QUDA_SINGLE_PRECISION;
    gauge_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
    gauge_param.cuda_prec_eigensolver  = QUDA_DOUBLE_PRECISION;

    
    gauge_param.reconstruct                   = QUDA_RECONSTRUCT_NO;  ////gauge_param.reconstruct = link_recon;
    gauge_param.reconstruct_sloppy            = QUDA_RECONSTRUCT_NO;
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

    gauge_param.anisotropy = 1.0; /////anisotropy

    int pad_size = 0;
    // For multi-GPU, ga_pad must be large enough to store a time-slice
    int x_face_size = gauge_param.X[1] * gauge_param.X[2] * gauge_param.X[3] / 2;
    int y_face_size = gauge_param.X[0] * gauge_param.X[2] * gauge_param.X[3] / 2;
    int z_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[3] / 2;
    int t_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[2] / 2;
    pad_size = std::max({x_face_size, y_face_size, z_face_size, t_face_size});
    gauge_param.ga_pad = pad_size;
    gauge_param.struct_size = sizeof(gauge_param);

    if(apply_stag_phase == 1){
      gauge_param.scale = 1.0;
      gauge_param.staggered_phase_type = QUDA_STAGGERED_PHASE_MILC;
      applyGaugeFieldScaling_long((qlat::Complex*) quda_gf.data(), V/2, &gauge_param);
    }
    /////load gauge to quda GPU default position
    loadGaugeQuda((void *) quda_gf.data(), &gauge_param);

    // Compute plaquette as a sanity check from interal GPU gauge
    double plaq[3];
    plaqQuda(plaq);
    print0("Computed plaquette is %.8e (spatial = %.8e, temporal = %.8e)\n", plaq[0], plaq[1], plaq[2]);

    ////===END of gauge_param

    csrc = NULL; cres = NULL;
    gsrc = NULL; gres = NULL;

    ctmp0 = NULL; ctmp1 = NULL; ctmp2 = NULL;
    gtmp0 = NULL; gtmp1 = NULL; gtmp2 = NULL;

    eig_solve = NULL;
    mat = NULL;
    mat_Mdag  = NULL;
    mat_MMdag = NULL;
    mat_MdagM = NULL;

    dirac = NULL;
    dSloppy = NULL;
    dPre    = NULL;
    dEig    = NULL;

    nvec = 0;
    use_eigen_pc = 0;
    check_residue = 0;

  }

  void free_mem(){
    if(csrc != NULL){delete csrc;csrc=NULL;}if(cres != NULL){delete cres;cres=NULL;}
    if(gsrc != NULL){delete gsrc;gsrc=NULL;}if(gres != NULL){delete gres;gres=NULL;}
    if(ctmp0 != NULL){delete ctmp0;ctmp0=NULL;}
    if(ctmp1 != NULL){delete ctmp1;ctmp1=NULL;}
    if(ctmp2 != NULL){delete ctmp2;ctmp2=NULL;}
    if(gtmp0 != NULL){delete gtmp0;gtmp0=NULL;}
    if(gtmp1 != NULL){delete gtmp1;gtmp1=NULL;}
    if(gtmp2 != NULL){delete gtmp2;gtmp2=NULL;}

    evecs.resize(0);evals.resize(0);

    for (int i = 0; i < nvec; i++) delete evecs_[i];
    for (int i = 0; i < nvec; i++) delete kSpace[i];
    evecs_.resize(0);kSpace.resize(0);

    if(dirac != NULL){delete dirac;dirac=NULL;}
    if(dSloppy != NULL){delete dSloppy;dSloppy=NULL;}
    if(dPre != NULL){delete dPre;dPre=NULL;}
    if(dEig != NULL){delete dEig;dEig=NULL;}

    if(eig_solve != NULL){delete eig_solve;eig_solve=NULL;}
    if(mat       != NULL){delete mat      ;  mat       = NULL;}
    if(mat_Mdag  != NULL){delete mat_Mdag ;  mat_Mdag  = NULL;}
    if(mat_MMdag != NULL){delete mat_MMdag;  mat_MMdag = NULL;}
    if(mat_MdagM != NULL){delete mat_MdagM;  mat_MdagM = NULL;}

    nvec = 0;

  }

  void setup_clover(double kappa, double clover_csw, double err = 1e-15, int niter = 10000)
  {

    fermion_type = 0;
    qassert(apply_stag_phase == 0);

    /////===Start of Inv parameters
    inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH;
    /////double kappa      = kappa;
    double anisotropy = gauge_param.anisotropy;
    inv_param.kappa = kappa;
    inv_param.mass = 0.5 / kappa - (1.0 + 3.0 / anisotropy);

    printfQuda("Kappa = %.8f Mass = %.8f\n", inv_param.kappa, inv_param.mass);

    // Use 3D or 4D laplace
    //===inv_param.laplace3D = laplace3D;
    inv_param.Ls = 1;

    inv_param.cpu_prec                      = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec                     = QUDA_DOUBLE_PRECISION;

    inv_param.cuda_prec_sloppy              = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_refinement_sloppy   = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_precondition   = QUDA_DOUBLE_PRECISION;
    inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
    inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
    inv_param.dirac_order = QUDA_DIRAC_ORDER;

    ////related to eigensystem
    inv_param.cuda_prec_eigensolver =  QUDA_DOUBLE_PRECISION;


    inv_param.clover_cpu_prec               = QUDA_DOUBLE_PRECISION;
    inv_param.clover_cuda_prec              = QUDA_DOUBLE_PRECISION;
    inv_param.clover_cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
    inv_param.clover_cuda_prec_eigensolver  = QUDA_DOUBLE_PRECISION;

    inv_param.clover_cuda_prec_sloppy       = QUDA_SINGLE_PRECISION;
    inv_param.clover_cuda_prec_refinement_sloppy = QUDA_SINGLE_PRECISION;

    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
    // Use kappa * csw or supplied clover_coeff
    /////double clover_csw = clover_csw;
    bool compute_clover_trlog = false;
    //bool compute_clover_trlog = true;
    inv_param.clover_csw = clover_csw;
    inv_param.clover_coeff = clover_csw * inv_param.kappa;
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
    inv_param.solver_normalization = QUDA_DEFAULT_NORMALIZATION;
    //inv_param.solver_normalization = QUDA_KAPPA_NORMALIZATION;

    ///unknown
    int gcrNkrylov = 10;
    QudaCABasis ca_basis = QUDA_POWER_BASIS;
    double ca_lambda_min = 0.0;
    double ca_lambda_max = -1.0;

    inv_param.pipeline = 1;
    inv_param.Nsteps = 2;
    inv_param.gcrNkrylov = gcrNkrylov;
    inv_param.ca_basis = ca_basis;
    inv_param.ca_lambda_min = ca_lambda_min;
    inv_param.ca_lambda_max = ca_lambda_max;

    ////ouble err = 1e-10;
    ///int niter  = 100000;

    inv_param.tol = err;
    inv_param.tol_restart = 0.0005;
    //if (tol_hq == 0 && tol == 0) {
    //  errorQuda("qudaInvert: requesting zero residual\n");
    //  exit(1);
    //}

    // require both L2 relative and heavy quark residual to determine convergence
    inv_param.residual_type = static_cast<QudaResidualType_s>(0);
    inv_param.residual_type = static_cast<QudaResidualType_s>(inv_param.residual_type | QUDA_L2_RELATIVE_RESIDUAL);

    inv_param.tol_hq = err; // specify a tolerance for the residual for heavy quark residual

    //// Offsets used only by multi-shift solver
    //// These should be set in the application code. We set the them here by way of
    //// example
    //inv_param.num_offset = multishift;
    //for (int i = 0; i < inv_param.num_offset; i++) inv_param.offset[i] = 0.06 + i * i * 0.1;
    //// these can be set individually
    //for (int i = 0; i < inv_param.num_offset; i++) {
    //  inv_param.tol_offset[i] = inv_param.tol;
    //  inv_param.tol_hq_offset[i] = inv_param.tol_hq;
    //}
    inv_param.maxiter = niter;

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
    //inv_param.verbosity_precondition = mg_verbosity[0];
    //inv_param.omega = 1.0;

    inv_param.input_location  = QUDA_CPU_FIELD_LOCATION;
    inv_param.output_location = QUDA_CPU_FIELD_LOCATION;

    // QUDA_DEBUG_VERBOSE is too nasty.
    inv_param.verbosity   = QUDA_VERBOSE;

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

    ////int setup_clover = 1;
    ///int do_inv = 1;
    /////if(setup_clover == 1)
    {
    //////===operator define
    //void *clover = nullptr;
    //void *clover_inv = nullptr;

    int clover_site_size           = 72; // real numbers per block-diagonal clover matrix
    ////size_t host_clover_data_type_size = (cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
    int host_clover_data_type_size = sizeof(qlat::Complex)/2;
    ////size_t host_spinor_data_type_size = (cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
    int host_spinor_data_type_size = sizeof(qlat::Complex)/2;

    quda_clover.resize(    V * clover_site_size * host_clover_data_type_size);
    quda_clover_inv.resize(V * clover_site_size * host_spinor_data_type_size);

    //constructHostCloverField((void*) quda_clover.data(), (void*) quda_clover_inv.data(), inv_param);
    /////===host

    bool compute_clover = true;
    //if(in.clover_csw == 0){compute_clover = false;}
    //double norm = 0.00; // clover components are random numbers in the range (-norm, norm)
    //double diag = 1.0;  // constant added to the diagonal
    ////if (!compute_clover) constructQudaCloverField((void*) quda_clover.data(), norm, diag, inv_param.clover_cpu_prec, V);
    inv_param.compute_clover = compute_clover;
    if (compute_clover) inv_param.return_clover = 1;
    inv_param.compute_clover_inverse = 1;
    inv_param.return_clover_inverse  = 1;
    //if(in.clover_csw == 0){
    //inv_param.compute_clover_inverse = 0;
    //inv_param.return_clover_inverse  = 0;
    //}
    /////===host
    ///// Load the clover terms to the device
    loadCloverQuda((void*) quda_clover.data(), (void*)  quda_clover_inv.data(), &inv_param);
    //////===operator define
    }

    bool pc_solve = (inv_param.solve_type == QUDA_DIRECT_PC_SOLVE) ||
      (inv_param.solve_type == QUDA_NORMOP_PC_SOLVE) || (inv_param.solve_type == QUDA_NORMERR_PC_SOLVE);

    /////print0("requested precision %d\n",inv_param.cuda_prec);
    createDiracWithEig(dirac, dSloppy, dPre, dEig, inv_param, pc_solve);
    //createDirac(dirac, dSloppy, dPre, inv_param, pc_solve);

    mat        = new DiracM(*dirac);
    mat_MdagM  = new DiracMdagM(dirac);
    mat_MMdag  = new DiracMMdag(*dirac);
    mat_Mdag   = new DiracMdag(*dirac);

    constructWilsonTestSpinorParam(&cs_cpu, &inv_param, &gauge_param);

    cs_gpu = quda::ColorSpinorParam(cs_cpu);cs_gpu.location = QUDA_CUDA_FIELD_LOCATION;
    cs_gpu.create = QUDA_ZERO_FIELD_CREATE;
    cs_gpu.setPrecision(inv_param.cuda_prec_eigensolver, inv_param.cuda_prec_eigensolver, true);
    if(cs_gpu.nSpin != 1) cs_gpu.gammaBasis = QUDA_UKQCD_GAMMA_BASIS;

    csrc  = quda::ColorSpinorField::Create(cs_cpu);
    cres  = quda::ColorSpinorField::Create(cs_cpu);
    ctmp0 = quda::ColorSpinorField::Create(cs_cpu);
    ctmp1 = quda::ColorSpinorField::Create(cs_cpu);
    ctmp2 = quda::ColorSpinorField::Create(cs_cpu);

    gsrc  = quda::ColorSpinorField::Create(cs_gpu);
    gres  = quda::ColorSpinorField::Create(cs_gpu);
    gtmp0 = quda::ColorSpinorField::Create(cs_gpu);
    gtmp1 = quda::ColorSpinorField::Create(cs_gpu);
    gtmp2 = quda::ColorSpinorField::Create(cs_gpu);


  }

  void setup_stagger(double fermion_mass, double err = 1e-15, int niter = 10000)
  {
    fermion_type = 1;

    /////===Start of Inv parameters
    inv_param.dslash_type = QUDA_STAGGERED_DSLASH;
    /////double kappa      = kappa;
    //double anisotropy = gauge_param.anisotropy;
    //inv_param.kappa = kappa;
    //inv_param.mass = 0.5 / kappa - (1.0 + 3.0 / anisotropy);
    //printfQuda("Kappa = %.8f Mass = %.8f\n", inv_param.kappa, inv_param.mass);

    inv_param.mass  = fermion_mass;

    // Use 3D or 4D laplace
    //===inv_param.laplace3D = laplace3D;
    inv_param.Ls = 1;

    inv_param.cpu_prec                      = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec                     = QUDA_DOUBLE_PRECISION;

    inv_param.cuda_prec_sloppy              = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_refinement_sloppy   = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_precondition   = QUDA_DOUBLE_PRECISION;
    inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
    inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
    inv_param.dirac_order = QUDA_DIRAC_ORDER;

    ////related to eigensystem
    inv_param.cuda_prec_eigensolver =  QUDA_DOUBLE_PRECISION;

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
    //if(in.solver_type == 0){inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;}
    //if(in.solver_type == 1){inv_param.solve_type = QUDA_NORMERR_PC_SOLVE;}
    //if(in.solver_type == 2){inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;}


    // There might be a performance difference.
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;

    inv_param.dagger = QUDA_DAG_NO;

    //inv_param.mass_normalization   = QUDA_KAPPA_NORMALIZATION;
    inv_param.mass_normalization   = QUDA_MASS_NORMALIZATION;
    inv_param.solver_normalization = QUDA_DEFAULT_NORMALIZATION;
    //inv_param.solver_normalization = QUDA_KAPPA_NORMALIZATION;

    ///unknown
    int gcrNkrylov = 10;
    QudaCABasis ca_basis = QUDA_POWER_BASIS;
    double ca_lambda_min = 0.0;
    double ca_lambda_max = -1.0;

    inv_param.pipeline = 1;
    inv_param.Nsteps = 2;
    inv_param.gcrNkrylov = gcrNkrylov;
    inv_param.ca_basis = ca_basis;
    inv_param.ca_lambda_min = ca_lambda_min;
    inv_param.ca_lambda_max = ca_lambda_max;

    ////ouble err = 1e-10;
    ///int niter  = 100000;

    inv_param.tol = err;
    inv_param.tol_restart = 0.0005;
    //if (tol_hq == 0 && tol == 0) {
    //  errorQuda("qudaInvert: requesting zero residual\n");
    //  exit(1);
    //}

    // require both L2 relative and heavy quark residual to determine convergence
    inv_param.residual_type = static_cast<QudaResidualType_s>(0);
    inv_param.residual_type = static_cast<QudaResidualType_s>(inv_param.residual_type | QUDA_L2_RELATIVE_RESIDUAL);

    inv_param.tol_hq = err; // specify a tolerance for the residual for heavy quark residual

    //// Offsets used only by multi-shift solver
    //// These should be set in the application code. We set the them here by way of
    //// example
    //inv_param.num_offset = multishift;
    //for (int i = 0; i < inv_param.num_offset; i++) inv_param.offset[i] = 0.06 + i * i * 0.1;
    //// these can be set individually
    //for (int i = 0; i < inv_param.num_offset; i++) {
    //  inv_param.tol_offset[i] = inv_param.tol;
    //  inv_param.tol_hq_offset[i] = inv_param.tol_hq;
    //}
    inv_param.maxiter = niter;

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
    //inv_param.verbosity_precondition = mg_verbosity[0];
    //inv_param.omega = 1.0;

    inv_param.input_location  = QUDA_CPU_FIELD_LOCATION;
    inv_param.output_location = QUDA_CPU_FIELD_LOCATION;

    // QUDA_DEBUG_VERBOSE is too nasty.
    inv_param.verbosity   = QUDA_VERBOSE;

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

    bool pc_solve = (inv_param.solve_type == QUDA_DIRECT_PC_SOLVE) ||
      (inv_param.solve_type == QUDA_NORMOP_PC_SOLVE) || (inv_param.solve_type == QUDA_NORMERR_PC_SOLVE);

    ///////print0("requested precision %d\n",inv_param.cuda_prec);
    createDiracWithEig(dirac, dSloppy, dPre, dEig, inv_param, pc_solve);
    ////createDirac(dirac, dSloppy, dPre, inv_param, pc_solve);

    mat        = new DiracM(*dirac);
    mat_MdagM  = new DiracMdagM(dirac);
    mat_MMdag  = new DiracMMdag(*dirac);
    mat_Mdag   = new DiracMdag(*dirac);

    ////=====START construct Staggered color spin parameters
    cs_cpu.nColor = 3;
    cs_cpu.nSpin = 1;
    cs_cpu.nDim = 5;
    for (int d = 0; d < 4; d++) cs_cpu.x[d] = gauge_param.X[d];
    bool pc = isPCSolution(inv_param.solution_type);
    if (pc) cs_cpu.x[0] /= 2;
    cs_cpu.x[4] = 1;
    cs_cpu.pc_type = QUDA_4D_PC;
    cs_cpu.siteSubset = pc ? QUDA_PARITY_SITE_SUBSET : QUDA_FULL_SITE_SUBSET;
    // Lattice vector data properties
    cs_cpu.setPrecision(inv_param.cpu_prec);
    cs_cpu.pad = 0;
    cs_cpu.siteOrder = QUDA_EVEN_ODD_SITE_ORDER;
    cs_cpu.fieldOrder = QUDA_SPACE_SPIN_COLOR_FIELD_ORDER;
    cs_cpu.gammaBasis = inv_param.gamma_basis;
    cs_cpu.create = QUDA_ZERO_FIELD_CREATE;
    cs_cpu.location = QUDA_CPU_FIELD_LOCATION;
    ////=====END construct Staggered color spin parameters

    ////constructWilsonTestSpinorParam(&cs_cpu, &inv_param, &gauge_param);

    cs_gpu = quda::ColorSpinorParam(cs_cpu);cs_gpu.location = QUDA_CUDA_FIELD_LOCATION;
    cs_gpu.create = QUDA_ZERO_FIELD_CREATE;
    cs_gpu.setPrecision(inv_param.cuda_prec_eigensolver, inv_param.cuda_prec_eigensolver, true);

    csrc  = quda::ColorSpinorField::Create(cs_cpu);
    cres  = quda::ColorSpinorField::Create(cs_cpu);
    ctmp0 = quda::ColorSpinorField::Create(cs_cpu);
    ctmp1 = quda::ColorSpinorField::Create(cs_cpu);
    ctmp2 = quda::ColorSpinorField::Create(cs_cpu);

    gsrc  = quda::ColorSpinorField::Create(cs_gpu);
    gres  = quda::ColorSpinorField::Create(cs_gpu);
    gtmp0 = quda::ColorSpinorField::Create(cs_gpu);
    gtmp1 = quda::ColorSpinorField::Create(cs_gpu);
    gtmp2 = quda::ColorSpinorField::Create(cs_gpu);

  }

  void setup_eigen(int num_eigensys, double err = 1e-15)
  {
    /////may need settings
    nvec = num_eigensys;
    if(nvec == 0){return ;}
    use_eigen_pc = 1;
    //quda::setVerbosity(QUDA_VERBOSE);
    setVerbosity(QUDA_VERBOSE);

    int eig_n_kr = 3*nvec;
    double eig_tol    = err;
    double eig_qr_tol = err*0.1;
    int eig_batched_rotate = 0; // If unchanged, will be set to maximum
    int eig_check_interval = 10;
    int eig_max_restarts = 1000;
    bool eig_use_eigen_qr = true;
    bool eig_use_poly_acc = true;
    int eig_poly_deg = 100;
    double eig_amin = 0.1;
    double eig_amax = 0.0; // If zero is passed to the solver, an estimate will be computed

    ////preconditionor for the inverter
    inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;

    qassert(nvec > 0);
    printf("number of eigen vectors %d .\n", nvec);
    eig_param.invert_param = &inv_param;

    ////===basice eig parameters
    eig_param.eig_type = QUDA_EIG_TR_LANCZOS; ////or QUDA_EIG_IR_ARNOLDI
    eig_param.spectrum = QUDA_SPECTRUM_SR_EIG;
    if ((eig_param.eig_type == QUDA_EIG_TR_LANCZOS || eig_param.eig_type == QUDA_EIG_BLK_TR_LANCZOS)
        && !(eig_param.spectrum == QUDA_SPECTRUM_LR_EIG || eig_param.spectrum == QUDA_SPECTRUM_SR_EIG)) {
      errorQuda("Only real spectrum type (LR or SR) can be passed to Lanczos type solver.");
    }

    // The solver will exit when n_conv extremal eigenpairs have converged
    eig_param.n_conv       = nvec;
    // Inverters will deflate only this number of vectors.
    eig_param.n_ev_deflate = nvec;

    //eig_param.block_size
    //  = (eig_param.eig_type == QUDA_EIG_TR_LANCZOS || eig_param.eig_type == QUDA_EIG_IR_ARNOLDI) ? 1 : eig_block_size;
    eig_param.block_size = 1;

    eig_param.n_ev = nvec;
    eig_param.n_kr = eig_n_kr;
    eig_param.tol  = eig_tol;
    eig_param.qr_tol = eig_qr_tol;
    eig_param.batched_rotate = eig_batched_rotate;

    eig_param.require_convergence = QUDA_BOOLEAN_TRUE;
    eig_param.check_interval = eig_check_interval;
    eig_param.max_restarts = eig_max_restarts;


    //setEigParam(eig_param);
    //eig_param.use_norm_op = QUDA_BOOLEAN_TRUE;
    //eig_param.use_dagger  = QUDA_BOOLEAN_TRUE;

    ////** What type of Dirac operator we are using **/
    ////** If !(use_norm_op) && !(use_dagger) use M. **/
    ////** If use_dagger, use Mdag **/
    ////** If use_norm_op, use MdagM **/
    ////** If use_norm_op && use_dagger use MMdag. **/

    eig_param.use_norm_op = QUDA_BOOLEAN_TRUE;
    eig_param.use_dagger  = QUDA_BOOLEAN_FALSE;

    eig_param.compute_gamma5  = QUDA_BOOLEAN_FALSE;
    eig_param.compute_svd     = QUDA_BOOLEAN_FALSE;
    eig_param.arpack_check    = QUDA_BOOLEAN_FALSE;

    if (eig_param.compute_svd) {
      eig_param.use_dagger = QUDA_BOOLEAN_FALSE;
      eig_param.use_norm_op = QUDA_BOOLEAN_TRUE;
    }

    eig_param.use_eigen_qr = eig_use_eigen_qr ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
    eig_param.use_poly_acc = eig_use_poly_acc ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
    eig_param.poly_deg = eig_poly_deg;
    eig_param.a_min = eig_amin;
    eig_param.a_max = eig_amax;

    eig_param.arpack_check = QUDA_BOOLEAN_FALSE;

    strcpy(eig_param.vec_infile, "");
    strcpy(eig_param.vec_outfile, "");
    //eig_param.save_prec = eig_save_prec;
    bool eig_io_parity_inflate = false;
    eig_param.io_parity_inflate = eig_io_parity_inflate ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;

    eig_param.struct_size = sizeof(eig_param);
    /////===END of eig_param

    ///evecs = (void **)safe_malloc(nvec * sizeof(void *));
    long volume_size = V/2;
    if(inv_param.solution_type == QUDA_MAT_SOLUTION or inv_param.solution_type == QUDA_MATDAG_MAT_SOLUTION)
    {volume_size = V;}

    evecs.resize(nvec);evals.resize(nvec, 0.0);
    int spinor_site_size = 12;
    if(fermion_type == 1){spinor_site_size = 3;}

    for (int i = 0; i < nvec; i++) {evecs[i].resize(volume_size * spinor_site_size);}

    ColorSpinorParam cpuParam(&evecs[0][0], inv_param, &X[0], inv_param.solution_type, inv_param.input_location);
    for (int i = 0; i < num_eigensys; i++) {
      cpuParam.v = &evecs[i][0];
      evecs_.push_back(ColorSpinorField::Create(cpuParam));
    }

    for (int i = 0; i < nvec; i++) { kSpace.push_back(ColorSpinorField::Create(cs_gpu)); }

    // If you attempt to compute part of the imaginary spectrum of a symmetric matrix,
    // the solver will fail.
    if ((eig_param.spectrum == QUDA_SPECTRUM_LI_EIG || eig_param.spectrum == QUDA_SPECTRUM_SI_EIG)
        && ((eig_param.use_norm_op || (inv_param.dslash_type == QUDA_LAPLACE_DSLASH))
            || ((inv_param.dslash_type == QUDA_STAGGERED_DSLASH || inv_param.dslash_type == QUDA_ASQTAD_DSLASH)
                && inv_param.solve_type == QUDA_DIRECT_PC_SOLVE))) {
      errorQuda("Cannot compute imaginary spectra with a hermitian operator");
    }
    // Gamma5 pre-multiplication is only supported for the M type operator
    if (eig_param.compute_gamma5) {
      if (eig_param.use_norm_op || eig_param.use_dagger) {
        errorQuda("gamma5 premultiplication is only supported for M type operators: dag = %s, normop = %s",
                  eig_param.use_dagger ? "true" : "false", eig_param.use_norm_op ? "true" : "false");
      }
    }

    ///gtmp1 = quda::ColorSpinorField::Create(cs_gpu);
    ///d.prepare(in_b, out_b, *gtmp0, *gtmp1, inv_param.solution_type);

    eig_solve = EigenSolver::create(&eig_param, *mat_MMdag, profileEigensolve);

    ////eig_solve->printEigensolverSetup();
    (*eig_solve)(kSpace, evals);
    //eig_solve->computeEvals(*mat, kSpace, evals, eig_param.n_conv);
    ////printf("===Deflation size n_ev_deflate %d, n_conv %d \n", eig_param->n_ev_deflate, eig_param->n_conv );

    //eig_solve->orthoCheck(kSpace, eig_param.n_ev_deflate);
    ////eig_solve->deflate(*tmp_out, *tmp_in, kSpace, evals, true);

    ////copy vectors to host
    for (int i = 0; i < eig_param.n_conv; i++) *evecs_[i] = *kSpace[i];

  }

  void do_inv(void* res, void* src)
  {
    ColorSpinorParam cpuParam(src, inv_param, &X[0], inv_param.solution_type, cs_cpu.location);

    cpuParam.v = src;
    ColorSpinorField* cpu_src = ColorSpinorField::Create(cpuParam);
    cpuParam.v = res;
    ColorSpinorField* cpu_res = ColorSpinorField::Create(cpuParam);

    if(use_eigen_pc == 1){

      *gsrc = *cpu_src;
      *gres = *cpu_res;

      blas::ax(0.0, *gres);
      eig_solve->deflate(*gres, *gsrc, kSpace, evals, true);
      (*mat_Mdag)(*gtmp1, *gres, *gtmp0, *gtmp2);

      ////normalization of operators
      if(fermion_type == 0){blas::ax((2.0*inv_param.kappa), *gtmp1);}

      *cpu_res = *gtmp1;

    }
    //blas::ax(0.0, *cpu_res);

    invertQuda(res, src, &inv_param);

    if(check_residue == 1)
    {
      /////===check residue
      //quda::Complex zero( 0.0, 0.0);
      quda::Complex n_unit(-1.0, 0.0);
      //quda::Complex p_unit(+1.0, 0.0);

      *gsrc = *cpu_src;
      *gres = *cpu_res;

      //blas::caxpby(zero, *tmp_in , zero, *temp);
      //blas::ax(0.0, *gmtp1);
      //(*mat)(*gtmp1, *gres, *gtmp0, *gtmp2);
      //blas::ax(1.0/(2.0*inv_param.kappa), *gtmp1);

      //blas::caxpby(p_unit, *gsrc , n_unit, *gtmp1);
      //Complex residual = sqrt(blas::norm2(*gtmp1));
      //printf("===solution residual %.8e %.8e \n", residual.real(), residual.imag());

      blas::ax(0.0, *gtmp1);
      (*mat)(*gtmp1, *gres, *gtmp0, *gtmp2);

      if(fermion_type == 0){blas::ax(1.0/(2.0*inv_param.kappa), *gtmp1);}

      quda::Complex evals  = blas::cDotProduct(*gtmp1, *gsrc) / sqrt(blas::norm2(*gsrc));
      quda::Complex factor = sqrt(blas::norm2(*gsrc)) / sqrt(blas::norm2(*gtmp1));
      blas::caxpby(evals, *gsrc , n_unit, *gtmp1);
      quda::Complex residual = sqrt(blas::norm2(*gtmp1));
      printf("===solution residual %.8e %.8e, factor %.8e %.8e \n", residual.real(), residual.imag(), factor.real(), factor.imag());
      /////===check residue
    }

    delete cpu_src;delete cpu_res;

  }

  ~quda_inverter()
  {
    V = 0;
    free_mem();
  }
};


}  // namespace qlat

#endif
