#pragma once

// can only be used with luchang's version of CPS

#include <qlat/qlat.h>

#ifdef NO_CPS

namespace qlat
{  //

typedef InverterDomainWall InverterDomainWallCPS;

inline void cps_begin(int* argc, char** argv[], const Coordinate& total_site)
{
  begin(argc, argv);
}

void cps_end() { end(); }

}  // namespace qlat

#else

#define QLAT_CPS

#include <gf/tools.h>

#include <gf/gauge_field.h>

#include <gf/fermion_field.h>

#include <gf/qed.h>

#include <gf/inverter.h>

#include <gf/rng_state.h>

#include <alg/alg_fix_gauge.h>

#include <alg/alg_rnd_gauge.h>

extern MPI_Comm QMP_COMM_WORLD;

#include "cps-utils.h"

#include "cps-pio.h"

#include "cps-lanc.h"

namespace qlat
{  //

void set_do_arg(cps::DoArg& do_arg, const Coordinate& total_site)
{
  using namespace cps;
  do_arg.x_sites = total_site[0];
  do_arg.y_sites = total_site[1];
  do_arg.z_sites = total_site[2];
  do_arg.t_sites = total_site[3];
  do_arg.s_sites = 2;
  do_arg.dwf_height = 1.0;
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_INPUT;
  do_arg.start_seed_value = 123121;
  do_arg.x_nodes = 0;
  do_arg.y_nodes = 0;
  do_arg.z_nodes = 0;
  do_arg.t_nodes = 0;
  do_arg.s_nodes = 0;
  do_arg.x_node_sites = 0;
  do_arg.y_node_sites = 0;
  do_arg.z_node_sites = 0;
  do_arg.t_node_sites = 0;
  do_arg.s_node_sites = 0;
  do_arg.gfix_chkb = 1;
}

void cps_begin(int* argc, char** argv[], const Coordinate& total_site)
{
  cps::Start(argc, argv);
  cps::DoArg do_arg;
  set_do_arg(do_arg, total_site);
  cps::GJP.Initialize(do_arg);
  cps::LRG.Initialize();
  Coordinate size_node(cps::SizeX(), cps::SizeY(), cps::SizeZ(), cps::SizeT());
  begin(cps::UniqueID(), size_node);
  sync_node();
  cps::dataWriteParNumber() = 256;
  cps::dataReadParNumber() = 256;
}

void cps_end() { cps::End(); }

inline cps::Geometry geo_convert(const Geometry& geo)
{
  Coordinate total_site = geo.total_site();
  cps::Geometry cgeo;
  cgeo.init(total_site.data(), geo.multiplicity);
  return cgeo;
}

inline Geometry geo_convert(const cps::Geometry& cgeo)
{
  qlat::Coordinate total_site(cgeo.totalSite(0), cgeo.totalSite(1),
                              cgeo.totalSite(2), cgeo.totalSite(3));
  qlat::Geometry geo;
  geo.init(total_site, cgeo.multiplicity);
  return geo;
}

template <class M, class N>
inline void value_convert(M& x, const N& y)
{
  qassert(sizeof(M) == sizeof(N));
  std::memcpy(&x, &y, sizeof(M));
}

template <class M, class N>
inline void field_convert(cps::GridComm<M>& gc, const Field<N>& f)
{
  TIMER_VERBOSE("field_convert");
  const Geometry& geo = f.geo();
  cps::Geometry cgeo = geo_convert(geo);
  gc.init(cgeo);
  qassert(check_matching_geo_mult(f.geo(), geo_convert(gc.getGeometry())));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<N> v = f.get_elems_const(xl);
    M* gcv = gc.getElem(xl.data(), 0);
    for (int m = 0; m < geo.multiplicity; ++m) {
      value_convert(gcv[m], v[m]);
    }
  }
}

template <class M, class N>
inline void field_convert(Field<N>& f, const cps::GridComm<M>& gc)
{
  TIMER_VERBOSE("field_convert");
  const cps::Geometry cgeo = gc.getGeometry();
  Geometry geo = geo_convert(cgeo);
  f.init(geo);
  qassert(check_matching_geo_mult(f.geo(), geo));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const M* gcv = gc.getElem(xl.data(), 0);
    Vector<N> v = f.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      value_convert(v[m], gcv[m]);
    }
  }
}

inline void gf_fix_gauge_landau(GaugeField& gf,
                                const double stop_cond = 1.0e-12,
                                const double max_iter_num = 500000)
{
  TIMER_VERBOSE("gf_fix_gauge_landau")
  cps::GaugeField cgf;
  field_convert(cgf, gf);
  cgf.exportLattice();
  cps::FixGaugeArg fix_gauge_arg;
  fix_gauge_arg.fix_gauge_kind = cps::FIX_GAUGE_LANDAU;
  fix_gauge_arg.stop_cond = stop_cond;
  fix_gauge_arg.max_iter_num = max_iter_num;
  cps::CommonArg common_arg;
  cps::Lattice& lat =
      cps::LatticeFactory::Create(cps::F_CLASS_NONE, cps::G_CLASS_NONE);
  cps::AlgFixGauge fg(lat, &common_arg, &fix_gauge_arg);
  fg.run();
  cps::AlgRotateGauge rg(lat, &common_arg);
  rg.run();
  fg.free();
  cps::LatticeFactory::Destroy();
  cgf.importLattice();
  field_convert(gf, cgf);
}

inline void gt_gf_fix_gauge_coulomb(GaugeTransform& gt, const GaugeField& gf,
                                    const double stop_cond = 1.0e-12,
                                    const double max_iter_num = 500000)
{
  TIMER_VERBOSE("gt_gf_fix_gauge_coulomb")
  cps::GaugeField cgf;
  field_convert(cgf, gf);
  cgf.exportLattice();
  cps::FixGaugeArg fix_gauge_arg;
  fix_gauge_arg.fix_gauge_kind = cps::FIX_GAUGE_COULOMB_T;
  fix_gauge_arg.stop_cond = stop_cond;
  fix_gauge_arg.max_iter_num = max_iter_num;
  fix_gauge_arg.hyperplane_start = 0;
  fix_gauge_arg.hyperplane_step = 1;
  fix_gauge_arg.hyperplane_num = gf.geo().total_site()[3];
  cps::CommonArg common_arg;
  cps::Lattice& lat =
      cps::LatticeFactory::Create(cps::F_CLASS_NONE, cps::G_CLASS_NONE);
  cps::AlgFixGauge fg(lat, &common_arg, &fix_gauge_arg);
  fg.run();
  gt.init(gf.geo());
#pragma omp parallel for
  for (long index = 0; index < gt.geo().local_volume(); ++index) {
    const Coordinate xl = gt.geo().coordinate_from_index(index);
    const int k =
        ((xl[2] * gt.geo().node_site[1]) + xl[1]) * gt.geo().node_site[0] + xl[0];
    value_convert(gt.get_elem(xl), lat.fix_gauge_ptr[xl[3]][k]);
  }
  fg.free();
  cps::LatticeFactory::Destroy();
}

inline void gf_fix_gauge_coulomb(GaugeField& gf,
                                 const double stop_cond = 1.0e-12,
                                 const double max_iter_num = 500000)
{
  TIMER_VERBOSE("gf_fix_gauge_coulomb")
  GaugeTransform gt;
  gt_gf_fix_gauge_coulomb(gt, gf, stop_cond, max_iter_num);
  gf_apply_gauge_transformation(gf, gf, gt);
}

inline cps::FermionActionDomainWall fa_convert(const FermionAction& fa)
{
  cps::FermionActionDomainWall cfa(fa.mass, fa.ls, fa.m5, true, fa.mobius_scale,
                                   fa.is_multiplying_dminus);
  if (fa.is_using_zmobius) {
    cfa.is_using_zmobius_action = true;
    cfa.CGdiagonalMee = fa.cg_diagonal_mee;
    cfa.bs = fa.bs;
    cfa.cs = fa.cs;
  }
  return cfa;
}

inline cps::GaugeActionImprRect ga_convert(const GaugeAction& ga)
{
  cps::GaugeActionImprRect cga(ga.beta, ga.c1);
  return cga;
}

inline cps::LancArg lanc_arg_convert(const LancArg& la, const FermionAction& fa)
{
  qassert(la.initialized);
  cps::LancArg cla;
  cla.mass = fa.mass;
  cla.N_use = la.n_use;
  cla.N_get = la.n_get;
  cla.N_true_get = la.n_true_get;
  cla.ch_beta = la.ch_beta;
  cla.ch_alpha = la.ch_alpha;
  cla.ch_ord = la.ch_ord;
  cla.stop_rsd = 1e-8;
  cla.qr_rsd = 1e-14;
  cla.EigenOper = cps::DDAGD;
  cla.precon = true;
  cla.ch_sh = false;
  cla.ch_mu = 0.0;
  cla.lock = false;
  cla.maxits = 10000;
  cla.fname = "lanc.dat";
  return cla;
}

struct LowModesCPS {
  bool initialized;
  Geometry geo;
  cps::LanczosDefault lanc;
  //
  void init()
  {
    initialized = false;
    lanc.free_evecs();
  }
  void init(const Geometry& geo_, const FermionAction& fa)
  {
    if (initialized) {
      return;
    }
    initialized = true;
    geo = geo_;
    lanc.init(geo_convert(geo), fa_convert(fa));
  }
  //
  LowModesCPS() { init(); }
};

inline void run_lanc(LowModesCPS& lm, const GaugeField& gf,
                     const FermionAction& fa, const LancArg& la)
{
  TIMER_VERBOSE("run_lanc");
  if (!la.initialized) {
    displayln_info("WARNING: " + std::string(fname) +
                   ": la is not initialized.");
    return;
  }
  qassert(!lm.initialized);
  cps::GaugeField cgf;
  field_convert(cgf, gf);
  cps::LancArg cla = lanc_arg_convert(la, fa);
  lm.init(geo_reform(gf.geo(), fa.ls, 0), fa);
  cps::LanczosRun<float>::run(lm.lanc, cgf, fa_convert(fa), cla);
}

inline long read_low_modes_compressed(LowModesCPS& lm, const std::string& path)
{
  if (!does_file_exist_sync_node(path + "/metadata.txt")) {
    displayln_info(
        ssprintf("load_from_compressed_eigen_vectors: '%s' do not exist.",
                 path.c_str()));
    return 0;
  }
  TIMER_VERBOSE_FLOPS("read_low_modes_compressed");
  cps::LanczosDefault& lanc = lm.lanc;
  lanc.free_evecs();
  vector<double> vals;
  std::vector<HalfVector> hvs;
  {
    CompressedEigenSystemInfo cesi;
    CompressedEigenSystemBases cesb;
    CompressedEigenSystemCoefs cesc;
    if (0 == load_compressed_eigen_vectors(vals, cesi, cesb, cesc, path)) {
      return 0;
    }
    std::vector<BlockedHalfVector> bhvs;
    decompress_eigen_system(bhvs, cesb, cesc);
    convert_half_vectors(hvs, bhvs);
  }
  const long nvec = hvs.size();
  const long vec_size = get_data(hvs[0]).data_size();
  timer.flops += nvec * vec_size * get_num_node();
  qassert(nvec == vals.size());
  qassert(vec_size == lanc.vec_size);
  lanc.resize(nvec);
  for (int k = 0; k < nvec; k++) {
    lanc.getVal(k) = vals[k];
  }
  {
    TIMER_VERBOSE_FLOPS("load_from_compressed_eigen_vectors-alloc");
    timer.flops += nvec * vec_size;
    for (int k = 0; k < nvec; k++) {
      lanc.alloc(k);
      convert_half_vector_bfm_format(
          Vector<ComplexF>((ComplexF*)lanc.getVec(k),
                           vec_size / sizeof(ComplexF)),
          hvs[k]);
      hvs[k].init();
    }
  }
  displayln_info(fname + ssprintf(": Reading %ld vectors, done.", nvec));
  return nvec * vec_size;
}

inline long read_low_modes(LowModesCPS& lm, const std::string& path)
// lm must be initialized
{
  TIMER_VERBOSE("read_low_modes");
  qassert(lm.initialized);
  long size = read_low_modes_compressed(lm, path);
  if (size == 0) {
    size = cps::lanczosReadParNode(lm.lanc, path);
  }
  return size;
}

inline long load_low_modes(LowModesCPS& lm, const std::string& path)
{
  TIMER_VERBOSE("load_low_modes");
  if (not lm.initialized) {
    return 0;
  }
  const long total_bytes = read_low_modes(lm, path);
  if (0 != total_bytes) {
    lm.initialized = true;
  }
  return total_bytes;
}

inline long write_low_modes(const LowModesCPS& lm, const std::string& path)
{
  TIMER_VERBOSE("write_low_modes");
  qassert(lm.initialized);
  return cps::lanczosWriteParNode(lm.lanc, path);
}

inline void load_or_compute_low_modes(LowModesCPS& lm, const std::string& path,
                                      const GaugeField& gf,
                                      const FermionAction& fa,
                                      const LancArg& la)
{
  TIMER_VERBOSE("load_or_compute_low_modes");
  qassert(!lm.initialized);
  if (path == "") {
    run_lanc(lm, gf, fa, la);
  } else {
    lm.init(geo_reform(gf.geo(), fa.ls, 0), fa);
    if (0 == read_low_modes(lm, path)) {
      lm.init();
      run_lanc(lm, gf, fa, la);
      if (lm.initialized) {
        write_low_modes(lm, path);
      }
    }
  }
}

struct InverterDomainWallCPS {
  box_acc<Geometry> geo;
  FermionAction fa;
  GaugeField gf;
  //
  cps::InverterDomainWallDefault inverter;
  ConstHandle<LowModesCPS> lm;
  //
  InverterDomainWallCPS() { init(); }
  //
  void init()
  {
    stop_rsd() = 1.0e-8;
    max_num_iter() = 50000;
    max_mixed_precision_cycle() = 100;
  }
  //
  void setup() { inverter.reinit(); }
  void setup(const GaugeField& gf_, const FermionAction& fa_)
  {
    geo = geo_reform(gf_.geo(), fa_.ls);
    fa = fa_;
    gf = gf_;
    cps::FermionActionDomainWall cfa = fa_convert(fa);
    cps::GaugeField cgf;
    field_convert(cgf, gf);
    inverter.geo.initialized = false;
    inverter.init(cgf, cfa);
  }
  //
  double& stop_rsd() { return inverter.stop_rsd; }
  //
  int& max_num_iter() { return inverter.max_num_iter; }
  //
  int& max_mixed_precision_cycle()
  {
    return inverter.max_mixed_precision_cycle;
  }
};

inline void setup_inverter(InverterDomainWallCPS& inverter)
{
  inverter.setup();
  inverter.lm.init();
}

inline void setup_inverter(InverterDomainWallCPS& inverter,
                           const GaugeField& gf, const FermionAction& fa)
{
  inverter.setup(gf, fa);
  inverter.lm.init();
}

inline void setup_inverter(InverterDomainWallCPS& inverter,
                           const GaugeField& gf, const FermionAction& fa,
                           const LowModesCPS& lm)
{
  setup_inverter(inverter, gf, fa);
  inverter.lm.init(lm);
}

inline long invert(FermionField5d& sol, const FermionField5d& src,
                   const InverterDomainWallCPS& inverter)
// sol do not need to be initialized
{
  TIMER_VERBOSE("invert(5d,5d,IDWCPS)");
  const Geometry& geo = src.geo();
  sol.init(geo);
  cps::FermionField5d csol, csrc;
  field_convert(csrc, src);
  field_convert(csol, sol);
  InverterDomainWallCPS& inv = *((InverterDomainWallCPS*)(&inverter));
  if (inverter.lm.null() or (not inverter.lm().initialized)) {
    inv.inverter.inv(csol, csrc);
  } else {
    inv.inverter.inv(csol, csrc, NULL, &(inverter.lm().lanc));
  }
  field_convert(sol, csol);
  return 0;
}

inline long invert(FermionField4d& sol, const FermionField4d& src,
                   const InverterDomainWallCPS& inverter)
{
  return invert_dwf(sol, src, inverter);
}

}  // namespace qlat

#include <gf/eff_overlap.h>

#endif
