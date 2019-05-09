#pragma once

#include <qlat/qlat.h>

// Quda
#include <quda.h>
#include <invert_quda.h>

#include <cstdlib>

QLAT_START_NAMESPACE

static int mpi_rank_from_coords(const int* coords, void* fdata)
{
  int* dims = reinterpret_cast<int*>(fdata);

  int rank = coords[3];
  for (int i = 2; i >= 0; i--) {
    rank = dims[i] * rank + coords[i];
  }
  return rank;
}

static void comm_set_gridsize(int* grid)
{
  initCommsGridQuda(4, grid, mpi_rank_from_coords, reinterpret_cast<void*>(grid));
}

namespace qlat_quda {

  inline void quda_begin(int mpi_layout[4])
  {  
    using namespace quda;
    // The following sets the MPI comm stuff.
    comm_set_gridsize(mpi_layout);
    initQuda(-1000);
    printfQuda("Quda initialized on rank #%03d.\n", comm_rank());
  }
  
  inline void quda_end()
  {
    using namespace quda;
    endQuda();
  }
  
  template<class T>
  void quda_convert_gauge(std::vector<T>& qgf, const GaugeField& gf)
  {
    TIMER_VERBOSE("quda_convert_gauge(qgf,gf)");
    const Geometry& geo = gf.geo;
    ColorMatrix* quda_pt = reinterpret_cast<ColorMatrix*>(qgf.data());
    qassert(geo.multiplicity == 4);
    long V = geo.local_volume();
    long Vh = V/2;
    
    for(int qlat_idx = 0; qlat_idx < V; qlat_idx++){
      Coordinate xl = geo.coordinate_from_index(qlat_idx );
      const Vector<ColorMatrix> ms = gf.get_elems_const(xl);
      int eo = (xl[0]+xl[1]+xl[2]+xl[3])%2;
      for(int mu = 0; mu < 4; mu++){
        int quda_idx = (qlat_idx/2 + eo*Vh)*4 + mu;
        quda_pt[quda_idx] = ms[mu];
      }
    }
  }
  
  template<class T>
  void quda_convert_fermion(FermionField5d& ff, const std::vector<T>& qff)
  {
    TIMER_VERBOSE("quda_convert_fermion(ff,qff)");
    const Geometry& geo = ff.geo;
    const WilsonVector* quda_pt = reinterpret_cast<const WilsonVector*>(qff.data());
    int Ls = geo.multiplicity;
    qassert(Ls > 0);
    long V = geo.local_volume();
    long Vh = V/2;
    
    #pragma omp parallel for
    for(long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++){
      const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
      int eo = (xl[0]+xl[1]+xl[2]+xl[3])%2;
      Vector<WilsonVector> wvs = ff.get_elems(xl);
      for (int s = 0; s < Ls; s++) {
        int quda_idx = eo*Vh*Ls + s*Vh + qlat_idx_4d/2; 
        wvs[s] = quda_pt[quda_idx];
      }
    }
  }
  
  template<class T>
  void quda_convert_fermion(std::vector<T>& qff, const FermionField5d& ff)
  {
    TIMER_VERBOSE("quda_convert_fermion(qff,ff)");
    const Geometry& geo = ff.geo;
    WilsonVector* quda_pt = reinterpret_cast<WilsonVector*>(qff.data());
    int Ls = geo.multiplicity;
    qassert(Ls > 0);
    long V = geo.local_volume();
    long Vh = V/2;
    
    #pragma omp parallel for
    for(long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++){
      const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
      int eo = (xl[0]+xl[1]+xl[2]+xl[3])%2;
      const Vector<WilsonVector> wvs = ff.get_elems_const(xl);
      for (int s = 0; s < Ls; s++) {
        int quda_idx = eo*Vh*Ls + s*Vh + qlat_idx_4d/2; 
        quda_pt[quda_idx] = wvs[s];
      }
    }
  }
  
  struct InverterDomainWallQuda : InverterDomainWall {
    // Now setup all the QUDA parameters
    QudaGaugeParam gauge_param;
    // newQudaGaugeParam();
    QudaInvertParam inv_param; 
    // newQudaInvertParam();
   
    std::vector<double> qff_src;
    std::vector<double> qff_sol;
  
    std::vector<double> qgf;
  
    InverterDomainWallQuda(): qff_src(0), qff_sol(0), qgf(0) { init(); }
    ~InverterDomainWallQuda() { init(); }
    //
    void init()
    {
      free();
      InverterDomainWall::init();
    }
    //
    void setup()
    {
      TIMER_VERBOSE("InvDWQuda::setup");
      using namespace quda;
      free();
  
      // Now setup all the QUDA parameters
      gauge_param = newQudaGaugeParam();
      inv_param = newQudaInvertParam();
      
      for(int mu = 0; mu < 4; mu++){
        gauge_param.X[mu]     = geo.node_site[mu];
      }
      
      // ... OK. I don't know what this means
      gauge_param.type        = QUDA_WILSON_LINKS;
      
      // Slowest changing to fastest changing: even-odd, mu, x_cb_4d, row, column, complex 
      // See the code later in this file to see the conversion between Grid inde and Quda index.
      gauge_param.gauge_order = QUDA_MILC_GAUGE_ORDER;
  
      // The precision used here should be the same as those set in the inv_param, i.e.
      // gauge_param.cuda_prec = inv_param.cuda_prec
      // gauge_param.cuda_prec_sloppy = inv_param.cuda_prec_sloppy
      gauge_param.cpu_prec    = QUDA_DOUBLE_PRECISION;
      gauge_param.cuda_prec   = QUDA_DOUBLE_PRECISION;
      gauge_param.reconstruct = QUDA_RECONSTRUCT_NO;
      gauge_param.cuda_prec_sloppy    
                              = QUDA_HALF_PRECISION;
      gauge_param.reconstruct_sloppy  
                              = QUDA_RECONSTRUCT_NO;
      
      gauge_param.gauge_fix   = QUDA_GAUGE_FIXED_NO;
  
      gauge_param.anisotropy  = 1.0;
      gauge_param.t_boundary  = QUDA_PERIODIC_T;
  
      int x_face_size = gauge_param.X[1] * gauge_param.X[2] * gauge_param.X[3] / 2;
      int y_face_size = gauge_param.X[0] * gauge_param.X[2] * gauge_param.X[3] / 2;
      int z_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[3] / 2;
      int t_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[2] / 2;
      int pad_size = std::max(x_face_size, y_face_size);
          pad_size = std::max(pad_size, z_face_size);
          pad_size = std::max(pad_size, t_face_size);
      gauge_param.ga_pad      = pad_size;
     
      // initialize the std::vectors that holds the gauge field.
      size_t qgf_size = geo.local_volume() * 4 * 18;
      qgf.resize(qgf_size);
      quda_convert_gauge(qgf, this->gf);
      loadGaugeQuda((void*)qgf.data(), &gauge_param);
      double plaq[3];
      plaqQuda(plaq);
      printfQuda("Computed plaquette is %16.12e (spatial = %16.12e, temporal = %16.12e)\n", plaq[0], plaq[1], plaq[2]);
  
      inv_param.Ls            = fa.ls;
      inv_param.dslash_type   = QUDA_MOBIUS_DWF_EOFA_DSLASH;
      inv_param.mass          = fa.mass;
      // Note that Quda uses -M5 as M5 ...
      inv_param.m5            = -fa.m5;
      if(fa.is_using_zmobius){
        // TODO: Error!
      }else{
        for(int s = 0; s < fa.ls; s++){
          inv_param.b_5[s]    = 0.5 * fa.mobius_scale + 0.5;
          inv_param.c_5[s]    = 0.5 * fa.mobius_scale - 0.5;
        }
      }
      // kappa is irrelevant for Mobius/DWF but you have to set it.
      inv_param.kappa         = 1./(2.*(1.+3./1.+fa.mass));
      inv_param.mass_normalization    
                              = QUDA_KAPPA_NORMALIZATION; 
      inv_param.solver_normalization  
                              = QUDA_DEFAULT_NORMALIZATION;
      
      // Whether or not content of your input void* pointer will be modified
      inv_param.preserve_source       
                              = QUDA_PRESERVE_SOURCE_YES;
  
      // I don't know what these are but you have to set them.
      inv_param.use_sloppy_partial_accumulator 
                              = 0;
      inv_param.solution_accumulator_pipeline  
                              = 1;
      
      // This is for the reliable update. Just set it to some large number.
      inv_param.max_res_increase 
                              = 20000;
  
      inv_param.mq1           = fa.mass;
      inv_param.mq2           = fa.mass;
      inv_param.mq3           = 0.01;
      inv_param.eofa_shift    = -0.12345678;
      inv_param.eofa_pm       = 1;
   
      // The solver tolerance, i.e. |MdagM * x - b| < tol * |b|
      inv_param.tol           = 1e-12;
      inv_param.tol_restart   = 1e-3;
      
      // The maximum number of iterations.
      inv_param.maxiter       = 50000;
  
      // This is for Quda's sophisticated reliable update. 0.1 should be good.
      inv_param.reliable_delta
                              = 0.1;
  
      // NORMOP_PC means preconditioned normal operator MdagM
      inv_param.solve_type    = QUDA_NORMOP_PC_SOLVE;
      
      // QUDA_MATPC_EVEN_EVEN means we solve on even sites and use symmetric preconditioning
      // The other options are:
      // QUDA_MATPC_ODD_ODD,
      // QUDA_MATPC_EVEN_EVEN_ASYMMETRIC,
      // QUDA_MATPC_ODD_ODD_ASYMMETRIC,
      //
      // There might be a performance difference.
      inv_param.matpc_type    = QUDA_MATPC_EVEN_EVEN;
  
      // Eventually we want the unpreconditioned solution.
      inv_param.solution_type = QUDA_MAT_SOLUTION;
      
      // MSPCG does NOT support EOFA, yet.
      inv_param.inv_type      = QUDA_CG_INVERTER;
  
      inv_param.dagger        = QUDA_DAG_NO;
  
      // The precision used to correct the inner solver.
      inv_param.cpu_prec      = QUDA_DOUBLE_PRECISION;;
      inv_param.cuda_prec     = QUDA_DOUBLE_PRECISION;;
      // The sloppy(inner) solver precision 
      inv_param.cuda_prec_sloppy 
                              = QUDA_HALF_PRECISION;
  
      inv_param.input_location  
                              = QUDA_CPU_FIELD_LOCATION;
      inv_param.output_location 
                              = QUDA_CPU_FIELD_LOCATION;
      
      // I don't know what these are but you have to set them.
      inv_param.sp_pad = 0;
      inv_param.cl_pad = 0;
  
      // Both CPS and Grid use this gamma matrix representation
      inv_param.gamma_basis   = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
      
      // Slowest changing to fastest changing: even-odd, Ls, x_cb_4d, spin, color, complex
      // See the code later in this file to see the conversion between Grid inde and Quda index.
      inv_param.dirac_order   = QUDA_DIRAC_ORDER;
  
      // QUDA_DEBUG_VERBOSE is too nasty.
      inv_param.verbosity     = QUDA_VERBOSE;
  
      // It seems the initial value of this is undefined so it's better to set it here. 
      // If set to QUDA_USE_INIT_GUESS_NO Quda will zero the input solution pointer before the solve.
      inv_param.use_init_guess 
                              = QUDA_USE_INIT_GUESS_YES;
    
      // initialize the std::vectors that hold source and solution vectors.
      size_t qff_size = geo.local_volume() * fa.ls * 24;
      qff_src.resize(qff_size);
      qff_sol.resize(qff_size);
    }
    void setup(const GaugeField& gf_, const FermionAction& fa_)
    {
      InverterDomainWall::setup(gf_, fa_);
      setup();
    }
    void setup(const GaugeField& gf_, const FermionAction& fa_,
               const LowModes& lm_)
    {
      InverterDomainWall::setup(gf_, fa_, lm_);
      setup();
    }
    //
    void free() { }
  };
  
  inline void setup_inverter(InverterDomainWallQuda& inv) { inv.setup(); }
  
  inline void setup_inverter(InverterDomainWallQuda& inv, const GaugeField& gf,
                             const FermionAction& fa)
  {
    inv.setup(gf, fa);
  }
  
  inline void setup_inverter(InverterDomainWallQuda& inv, const GaugeField& gf,
                             const FermionAction& fa, const LowModes& lm)
  {
    inv.setup(gf, fa, lm);
  }
  
  inline void invert(FermionField5d& sol, const FermionField5d& src,
                      const InverterDomainWallQuda& inv)
  {
    // inverse_with_cg(sol, src, inv, cg_with_herm_sym_2);
    // TODO
  }
  
  inline void invert(FermionField4d& sol, const FermionField4d& src,
                      const InverterDomainWallQuda& inv)
  {
    invert_dwf(sol, src, inv);
  }

} // namespace qlat_quda

QLAT_END_NAMESPACE
