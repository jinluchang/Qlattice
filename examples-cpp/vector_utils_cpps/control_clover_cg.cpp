#include <sys/sysinfo.h>
#include <unistd.h>

//#include <qutils/vector.h>
#include "quda_para.h"
#include "general_funs.h"
#include <qlat/qcd-smear.h>
//#include "utils_clover_inverter.h"

int main(int argc, char* argv[])
{
  using namespace qlat;
  inputpara in; 
  begin_Lat(&argc, &argv, in);

  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;
  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site); 
  fflush_MPI();

  int Nsrc =  12;
  int bfac =   2;
  int prec_type = 0;
  double err = 1e-15;
  int niter  = 100000;
  int quda_verbos = 0;

  in.find_para(std::string("bfac" ), bfac);
  in.find_para(std::string("prec_type" ), prec_type);
  in.find_para(std::string("err" ), err);
  in.find_para(std::string("niter" ), niter);
  in.find_para(std::string("quda_verbos" ), quda_verbos);

  //std::vector<int > sp;sp.resize(4);
  //for(int i=0;i<4;i++){sp[i] = 0;}

  int mpi_layout[4]={0,0,0,0};
  qlat::GeometryNode geon = qlat::get_geometry_node();for(int i=0;i<4;i++){mpi_layout[i] = geon.size_node[i];}

  qlat::GaugeField gf;gf.init(geo);
  qlat::GaugeField gf1;gf1.init(geo);
  random_link(gf1, in.seed);
  qlat::gf_ape_smear(gf, gf1, 0.125, 2);

  quda_begin(mpi_layout);

  long V = geo.local_volume();
  qlat::vector<qlat::Complex > quda_gf;quda_gf.resize(V * 4 * 3*3);

  quda_convert_gauge(quda_gf, gf);

  //////quda containers
  QudaGaugeParam gauge_param = newQudaGaugeParam();

  //////setWilsonGaugeParam(gauge_param);

  ////===Start of gauge_param
  for (int mu = 0; mu < 4; mu++) {gauge_param.X[mu] = geo.node_site[mu];}

  gauge_param.type = QUDA_WILSON_LINKS; //// or QUDA_SU3_LINKS
  // Slowest changing to fastest changing: even-odd, mu, x_cb_4d, row, column,
  // complex See the code later in this file to see the conversion between
  // Grid inde and Quda index.
  gauge_param.gauge_order = QUDA_MILC_GAUGE_ORDER;
  gauge_param.t_boundary = QUDA_PERIODIC_T; ////fermion_t_boundary

  /////===Start of Inv parameters
  QudaInvertParam inv_param = newQudaInvertParam();

  ///////setGaugeParam(gauge_param);
  
  if(prec_type == 0){
    gauge_param.cuda_prec              = QUDA_DOUBLE_PRECISION;
    gauge_param.cuda_prec_sloppy       = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec                          = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec_sloppy                   = QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec_refinement_sloppy        = QUDA_DOUBLE_PRECISION;

    inv_param.clover_cuda_prec                   = QUDA_DOUBLE_PRECISION;
    inv_param.clover_cuda_prec_sloppy            = QUDA_DOUBLE_PRECISION;
    inv_param.clover_cuda_prec_refinement_sloppy = QUDA_DOUBLE_PRECISION;
  }

  if(prec_type == 1){
    gauge_param.cuda_prec              = QUDA_SINGLE_PRECISION;
    gauge_param.cuda_prec_sloppy       = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec                          = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_sloppy                   = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_refinement_sloppy        = QUDA_SINGLE_PRECISION;


    inv_param.clover_cuda_prec                   = QUDA_SINGLE_PRECISION;
    inv_param.clover_cuda_prec_sloppy            = QUDA_SINGLE_PRECISION;
    inv_param.clover_cuda_prec_refinement_sloppy = QUDA_SINGLE_PRECISION;
  }

  if(prec_type == 2){
    gauge_param.cuda_prec                        = QUDA_HALF_PRECISION;
    gauge_param.cuda_prec_sloppy                 = QUDA_HALF_PRECISION;
    inv_param.cuda_prec                          = QUDA_HALF_PRECISION;
    inv_param.cuda_prec_sloppy                   = QUDA_HALF_PRECISION;
    inv_param.cuda_prec_refinement_sloppy        = QUDA_HALF_PRECISION;

    inv_param.clover_cuda_prec                   = QUDA_HALF_PRECISION;
    inv_param.clover_cuda_prec_sloppy            = QUDA_HALF_PRECISION;
    inv_param.clover_cuda_prec_refinement_sloppy = QUDA_HALF_PRECISION;
  }

  gauge_param.cpu_prec               = QUDA_DOUBLE_PRECISION;
  gauge_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
  gauge_param.cuda_prec_eigensolver  = QUDA_DOUBLE_PRECISION;

  gauge_param.reconstruct                   = QUDA_RECONSTRUCT_NO;  ////gauge_param.reconstruct = link_recon;
  gauge_param.reconstruct_sloppy            = QUDA_RECONSTRUCT_NO;
  gauge_param.reconstruct_precondition      = QUDA_RECONSTRUCT_NO;
  gauge_param.reconstruct_eigensolver       = QUDA_RECONSTRUCT_NO;
  gauge_param.reconstruct_refinement_sloppy = QUDA_RECONSTRUCT_NO;  /////link_recon
  
  gauge_param.anisotropy = 1.0;
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

  //// set parameters for the reference Dslash, and prepare fields to be loaded
  //setDims(gauge_param.X);

  loadGaugeQuda((void *) quda_gf.data(), &gauge_param);

  // Compute plaquette as a sanity check from interal GPU gauge
  //double plaq[3];
  //plaqQuda(plaq);
  //qmessage("Computed plaquette is %.8e (spatial = %.8e, temporal = %.8e)\n", plaq[0], plaq[1], plaq[2]);

  ////===END of gauge_param

  //setInvertParam(inv_param);
  // Set dslash type
  //if (dslash_type == QUDA_CLOVER_WILSON_DSLASH || dslash_type == QUDA_TWISTED_CLOVER_DSLASH) {

  inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH;
  double kappa      = in.kappa;
  double anisotropy = gauge_param.anisotropy;
  inv_param.kappa = kappa;
  inv_param.mass = 0.5 / kappa - (1.0 + 3.0 / anisotropy);

  qmessage("Kappa = %.8f Mass = %.8f\n", inv_param.kappa, inv_param.mass);

  // Use 3D or 4D laplace
  //===inv_param.laplace3D = laplace3D;
  inv_param.Ls = 1;

  inv_param.clover_cpu_prec               = QUDA_DOUBLE_PRECISION;

  inv_param.clover_cuda_prec_precondition = QUDA_SINGLE_PRECISION;
  inv_param.clover_cuda_prec_eigensolver  = QUDA_DOUBLE_PRECISION;
  inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
  // Use kappa * csw or supplied clover_coeff
  double clover_csw = in.clover_csw;
  bool compute_clover_trlog = false;
  //bool compute_clover_trlog = true;
  inv_param.clover_csw = clover_csw;
  inv_param.clover_coeff = clover_csw * inv_param.kappa;

  /////===unknow para
  inv_param.compute_clover_trlog = compute_clover_trlog ? 1 : 0;

  // General parameter setup
  inv_param.inv_type = QUDA_CG_INVERTER;

  // Eventually we want the unpreconditioned solution.
  inv_param.solution_type = QUDA_MAT_SOLUTION;

  // NORMOP_PC means preconditioned normal operator MdagM
  inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;

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
  //inv_param.cuda_prec_precondition = cuda_prec_precondition;
  //inv_param.cuda_prec_eigensolver = cuda_prec_eigensolver;
  //inv_param.omega = 1.0;

  //inv_param.cuda_prec                     = QUDA_DOUBLE_PRECISION;
  //inv_param.cuda_prec_sloppy              = QUDA_SINGLE_PRECISION;
  //inv_param.cuda_prec_refinement_sloppy   = QUDA_SINGLE_PRECISION;
  inv_param.cpu_prec                      = QUDA_DOUBLE_PRECISION;
  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
  inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  inv_param.dirac_order = QUDA_DIRAC_ORDER;

  inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
  inv_param.output_location = QUDA_CPU_FIELD_LOCATION;

  //inv_param.sp_pad = 0;
  //inv_param.cl_pad = 0;

  // QUDA_DEBUG_VERBOSE is too nasty.

  if(quda_verbos == 2){
    setVerbosity(QUDA_DEBUG_VERBOSE);
    inv_param.verbosity   = QUDA_VERBOSE;
  }

  if(quda_verbos == 1){
    inv_param.verbosity = QUDA_VERBOSE;
    setVerbosity(QUDA_VERBOSE);
  }

  if(quda_verbos == 0){
    inv_param.verbosity = QUDA_SUMMARIZE;
  }


  if(quda_verbos == -1){
    inv_param.verbosity = QUDA_SILENT;
    setVerbosity(QUDA_SILENT);
  }

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

  int setup_clover = 1;
  int do_inv = 1;

  if(setup_clover == 1){
  //////===operator define
  //void *clover = nullptr;
  //void *clover_inv = nullptr;

  int clover_site_size           = 72; // real numbers per block-diagonal clover matrix
  ////size_t host_clover_data_type_size = (cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
  int host_clover_data_type_size = sizeof(qlat::Complex)/2;
  ////size_t host_spinor_data_type_size = (cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
  int host_spinor_data_type_size = sizeof(qlat::Complex)/2;
  qlat::vector<char > quda_clover    ;quda_clover.resize(    V * clover_site_size * host_clover_data_type_size);
  qlat::vector<char > quda_clover_inv;quda_clover_inv.resize(V * clover_site_size * host_spinor_data_type_size);

  //constructHostCloverField((void*) quda_clover.data(), (void*) quda_clover_inv.data(), inv_param);
  /////===host

  bool compute_clover = true;
  //if(in.clover_csw == 0){compute_clover = false;}
  //double norm = 0.00; // clover components are random numbers in the range (-norm, norm)
  //double diag = 1.0;  // constant added to the diagonal
  //if (!compute_clover) constructQudaCloverField((void*) quda_clover.data(), norm, diag, inv_param.clover_cpu_prec, V);
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


  //int Nsrc = 12;
  //std::vector<qlat::FermionField4d > qlat_ff;qlat_ff.resize(Nsrc);
  //for(int i=0;i<qlat_ff.size();i++){qlat_ff[i].init(geo);}
  ////===Invertion part
  if(do_inv == 1)
  {
    std::vector<quda::ColorSpinorField *> quda_in(Nsrc);
    std::vector<quda::ColorSpinorField *> quda_out(Nsrc);
    quda::ColorSpinorParam cs_param;
    constructWilsonTestSpinorParam(&cs_param, &inv_param, &gauge_param);
    for (int i = 0; i < Nsrc; i++) {
      // Populate the host spinor with random numbers.
      quda_in[i] = quda::ColorSpinorField::Create(cs_param);
      quda_out[i] = quda::ColorSpinorField::Create(cs_param);
    }

    for(int bi=0;bi<bfac;bi++)
    {
      //// Vector construct START
      ////-----------------------------------------------------------------------------------

      inv_param.secs = 0;
      inv_param.gflops = 0;
      inv_param.iter = 0;

      //quda::ColorSpinorField *check;
      //check = quda::ColorSpinorField::Create(cs_param);

      for (int i = 0; i < Nsrc; i++) {
        // Populate the host spinor with random numbers.
        //quda_in[i]->Source(QUDA_RANDOM_SOURCE);

        //////set point src at zero
        qlat::Complex* src = (qlat::Complex*) (quda_in[i]->data());
        //long Vh = V / 2;
        random_Ty(src, V*12, 0);
        //#pragma omp parallel for
        //for (long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
        //  const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
        //  const Coordinate xg = geo.coordinate_g_from_l(xl);
        //  int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
        //  int quda_idx = eo * Vh + qlat_idx_4d / 2;

        //  for(int dc=0;dc<12;dc++){
        //  if(xg[0] == sp[0] and xg[1] == sp[1] and xg[2] == sp[2] and xg[3] == sp[3] and dc == i){
        //    res[quda_idx*12 + dc] = qlat::Complex(1.0, 0.0);
        //  }
        //  else{
        //    res[quda_idx*12 + dc] = qlat::Complex(0.0, 0.0);
        //  }
        //  }
        //}
        //////set point src at zero

      }

      std::vector<double> time(Nsrc);
      std::vector<double> gflops(Nsrc);
      std::vector<int> iter(Nsrc);

      for (int i = 0; i < Nsrc; i++) {
        invertQuda(quda_out[i]->data(), quda_in[i]->data(), &inv_param);
        time[i] = inv_param.secs;
        gflops[i] = inv_param.gflops / inv_param.secs;
        iter[i] = inv_param.iter;
        //printfQuda("Done: %i iter / %g secs = %g Gflops\n\n", inv_param.iter, inv_param.secs,
        //           inv_param.gflops / inv_param.secs);
        //quda_ff_to_Ffield4d(qlat_ff[i], (qlat::Complex*) quda_out[i]->data());
      }

      qmessage("Done g%03d: %8d iter / %.6f secs = %.3f Gflops, Cost %.3f Gflops \n", bi,
            inv_param.iter, inv_param.secs, inv_param.gflops / inv_param.secs, inv_param.gflops);

    }
    for (auto p : quda_in) { delete p; }
    for (auto p : quda_out) { delete p; }
  }

  quda_end();

  fflush_MPI();
  qlat::Timer::display();
  qlat::end();
  return 0;
}

