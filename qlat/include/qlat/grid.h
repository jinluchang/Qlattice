#pragma once


#ifdef NO_GRID

#include <qlat/qlat.h>

namespace qlat
{  //

typedef InverterDomainWall InverterDomainWallGrid;

inline void grid_begin(
    int* argc, char** argv[],
    const std::vector<Coordinate>& size_node_list = std::vector<Coordinate>())
{
  begin(argc, argv, size_node_list);
}

void grid_end() { end(); }

}  // namespace qlat

#else

#define QLAT_GRID

#include <Grid/Grid.h>
#include <qlat/qlat.h>

#include <cstdlib>

namespace qlat
{  //

inline Grid::Coordinate grid_convert(const Coordinate& x)
{
  std::vector<int> ret(4);
  for (int mu = 0; mu < 4; ++mu) {
    ret[mu] = x[mu];
  }
  return Grid::Coordinate(ret);
}

inline Grid::Coordinate grid_convert(const Coordinate& x, const int m)
{
  std::vector<int> ret(5);
  for (int mu = 0; mu < 4; ++mu) {
    ret[mu + 1] = x[mu];
  }
  ret[0] = m;
  return Grid::Coordinate(ret);
}

inline Coordinate grid_convert(const Grid::Coordinate& x)
{
  if (x.size() == 4) {
    return Coordinate(x[0], x[1], x[2], x[3]);
  } else if (x.size() == 5) {
    return Coordinate(x[1], x[2], x[3], x[4]);
  } else {
    qassert(false);
  }
}

inline int id_node_from_grid(const Grid::GridCartesian* UGrid)
{
  using namespace Grid;
  const Grid::Coordinate& mpi_layout = UGrid->_processors;
  const Grid::Coordinate& mpi_corr = UGrid->_processor_coor;
  const Coordinate size_node = grid_convert(mpi_layout);
  const Coordinate coor_node = grid_convert(mpi_corr);
  return index_from_coordinate(coor_node, size_node);
}

inline void grid_begin(
    int* argc, char** argv[],
    const std::vector<Coordinate>& size_node_list = std::vector<Coordinate>())
{
  using namespace Grid;
  system("rm /dev/shm/Grid*");
  Grid_init(argc, argv);
  const int num_node = init_mpi(argc, argv);
  Coordinate size_node;
  for (int i = 0; i < (int)size_node_list.size(); ++i) {
    size_node = size_node_list[i];
    if (num_node == product(size_node)) {
      break;
    }
  }
  if (num_node != product(size_node)) {
    size_node = plan_size_node(num_node);
  }
  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      grid_convert(size_node * 4), GridDefaultSimd(Nd, vComplexF::Nsimd()),
      grid_convert(size_node));
  UGrid->show_decomposition();
  displayln_info(
      ssprintf("GridThreads::GridThreads() = %d", GridThread::GetThreads()));
  const int id_node = id_node_from_grid(UGrid);
  delete UGrid;
  begin(id_node, size_node);
}

void grid_end()
{
  Grid::Grid_finalize();
  system("rm /dev/shm/Grid*");
  end();
}

void grid_convert(Grid::LatticeGaugeFieldF& ggf, const GaugeField& gf)
{
  TIMER_VERBOSE("grid_convert(ggf,gf)");
  using namespace Grid;
  const Geometry& geo = gf.geo();
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Grid::Coordinate coor = grid_convert(xl);
    const Vector<ColorMatrix> ms = gf.get_elems_const(xl);
    LorentzColourMatrixF gms;
    array<ComplexF, sizeof(LorentzColourMatrixF) / sizeof(ComplexF)>& fs =
        (array<ComplexF, sizeof(LorentzColourMatrixF) / sizeof(ComplexF)>&)
            gms;
    array<ComplexT, sizeof(LorentzColourMatrixF) / sizeof(ComplexF)>& ds =
        *((array<ComplexT, sizeof(LorentzColourMatrixF) /
                                    sizeof(ComplexF)>*)ms.data());
    qassert(sizeof(LorentzColourMatrixF) ==
            ms.data_size() / sizeof(ComplexT) * sizeof(ComplexF));
    qassert(fs.size() * sizeof(ComplexT) == ms.data_size());
    qassert(fs.size() == ds.size());
    for (int i = 0; i < fs.size(); ++i) {
      fs[i] = ds[i];
    }
    pokeLocalSite(gms, ggf, coor);
  }
}

void grid_convert(FermionField5d& ff, const Grid::LatticeFermionF& gff)
// ff need to be initialized with correct geo
{
  TIMER_VERBOSE("grid_convert(ff,gff)");
  using namespace Grid;
  const Geometry& geo = ff.geo();
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Grid::Coordinate coor = grid_convert(xl, 0);
    Vector<WilsonVector> wvs = ff.get_elems(xl);
    array<ComplexF, sizeof(WilsonVector) / sizeof(ComplexT)> fs;
    for (int m = 0; m < geo.multiplicity; ++m) {
      coor[0] = m;
      peekLocalSite(fs, gff, coor);
      array<ComplexT, sizeof(WilsonVector) / sizeof(ComplexT)>& ds =
          (array<ComplexT, sizeof(WilsonVector) / sizeof(ComplexT)>&)
              wvs[m];
      for (int k = 0; k < sizeof(WilsonVector) / sizeof(ComplexT); ++k) {
        ds[k] = fs[k];
      }
    }
  }
}

void grid_convert(Grid::LatticeFermionF& gff, const FermionField5d& ff)
{
  TIMER_VERBOSE("grid_convert(gff,ff)");
  using namespace Grid;
  const Geometry& geo = ff.geo();
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Grid::Coordinate coor = grid_convert(xl, 0);
    const Vector<WilsonVector> wvs = ff.get_elems_const(xl);
    array<ComplexF, sizeof(WilsonVector) / sizeof(ComplexT)> fs;
    for (int m = 0; m < geo.multiplicity; ++m) {
      coor[0] = m;
      const array<ComplexT, sizeof(WilsonVector) / sizeof(ComplexT)>& ds =
          (const array<ComplexT, sizeof(WilsonVector) / sizeof(ComplexT)>&)
              wvs[m];
      for (int k = 0; k < sizeof(WilsonVector) / sizeof(ComplexT); ++k) {
        fs[k] = ds[k];
      }
      pokeLocalSite(fs, gff, coor);
    }
  }
}

struct InverterDomainWallGrid : InverterDomainWall {
  Grid::GridCartesian* UGrid = NULL;
  Grid::GridRedBlackCartesian* UrbGrid = NULL;
  Grid::GridCartesian* FGrid = NULL;
  Grid::GridRedBlackCartesian* FrbGrid = NULL;
  Grid::LatticeGaugeFieldF* Umu = NULL;
  Grid::ZMobiusFermionF* Ddwf = NULL;
  Grid::SchurDiagTwoOperator<Grid::ZMobiusFermionF, Grid::LatticeFermionF>*
      HermOp = NULL;
  //
  InverterDomainWallGrid() { init(); }
  ~InverterDomainWallGrid() { init(); }
  //
  void init()
  {
    free();
    InverterDomainWall::init();
  }
  //
  void setup()
  {
    TIMER_VERBOSE("InvGrid::setup");
    using namespace Grid;
    free();
    const Coordinate total_site = geo.total_site();
    const Coordinate size_node = geo.geon.size_node;
    UGrid = SpaceTimeGrid::makeFourDimGrid(
        grid_convert(total_site), GridDefaultSimd(Nd, vComplexF::Nsimd()),
        grid_convert(size_node));
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
    FGrid = SpaceTimeGrid::makeFiveDimGrid(fa.ls, UGrid);
    FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(fa.ls, UGrid);
    qassert(geo.geon.id_node == id_node_from_grid(UGrid));
    Umu = new LatticeGaugeFieldF(UGrid);
    grid_convert(*Umu, gf);
    std::vector<Complex> omega(fa.ls, 0.0);
    for (int i = 0; i < fa.ls; ++i) {
      omega[i] = 1.0 / (fa.bs[i] + fa.cs[i]);
    }
    Ddwf = new ZMobiusFermionF(*Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid,
                               fa.mass, fa.m5, omega, 1.0, 0.0);
    HermOp =
        new Grid::SchurDiagTwoOperator<ZMobiusFermionF, LatticeFermionF>(*Ddwf);
  }
  void setup(const GaugeField& gf_, const FermionAction& fa_)
  {
    InverterDomainWall::setup(gf_, fa_);
    setup();
  }
  void setup(const GaugeField& gf_, const FermionAction& fa_, LowModes& lm_)
  {
    InverterDomainWall::setup(gf_, fa_, lm_);
    setup();
  }
  //
  void free()
  {
    TIMER_VERBOSE("InvGrid::free");
    if (HermOp != NULL) {
      delete HermOp;
      HermOp = NULL;
    }
    if (Ddwf != NULL) {
      delete Ddwf;
      Ddwf = NULL;
    }
    if (Umu != NULL) {
      delete Umu;
      Umu = NULL;
    }
    if (FrbGrid != NULL) {
      delete FrbGrid;
      FrbGrid = NULL;
    }
    if (FGrid != NULL) {
      delete FGrid;
      FGrid = NULL;
    }
    if (UrbGrid != NULL) {
      delete UrbGrid;
      UrbGrid = NULL;
    }
    if (UGrid != NULL) {
      delete UGrid;
      UGrid = NULL;
    }
  }
};

inline void setup_inverter(InverterDomainWallGrid& inv) { inv.setup(); }

inline void setup_inverter(InverterDomainWallGrid& inv, const GaugeField& gf,
                           const FermionAction& fa)
{
  inv.setup(gf, fa);
}

inline void setup_inverter(InverterDomainWallGrid& inv, const GaugeField& gf,
                           const FermionAction& fa, LowModes& lm)
{
  inv.setup(gf, fa, lm);
}

inline void multiply_m_grid(FermionField5d& out, const FermionField5d& in,
                            const InverterDomainWallGrid& inv)
{
  TIMER("multiply_m_grid(5d,5d,InvGrid)");
  out.init(geo_resize(in.geo()));
  using namespace Grid;
  GridCartesian* FGrid = inv.FGrid;
  LatticeFermionF gin(FGrid), gout(FGrid);
  grid_convert(gin, in);
  (*inv.Ddwf).M(gin, gout);
  grid_convert(out, gout);
}

inline void invert_grid_no_dminus(FermionField5d& sol,
                                  const FermionField5d& src,
                                  const InverterDomainWallGrid& inv)
// sol do not need to be initialized
{
  TIMER_VERBOSE("invert_grid_no_dminus(5d,5d,InvGrid)");
  sol.init(geo_resize(src.geo()));
  using namespace Grid;
  GridCartesian* FGrid = inv.FGrid;
  LatticeFermionF gsrc(FGrid), gsol(FGrid);
  grid_convert(gsrc, src);
  grid_convert(gsol, sol);
  if (is_checking_invert()) {
    FermionField5d ff;
    ff.init(geo_resize(src.geo()));
    grid_convert(ff, gsrc);
    displayln_info(fname + ssprintf(": src qnorm = %24.17E", qnorm(ff)));
    grid_convert(ff, gsol);
    displayln_info(fname + ssprintf(": sol qnorm = %24.17E", qnorm(ff)));
  }
  ConjugateGradient<LatticeFermionF> CG(inv.ip.stop_rsd, inv.ip.max_num_iter,
                                        false);
  SchurRedBlackDiagMooeeSolve<LatticeFermionF> SchurSolver(CG);
  SchurSolver(*inv.Ddwf, gsrc, gsol);
  grid_convert(sol, gsol);
  if (is_checking_invert()) {
    FermionField5d src1;
    src1.init(geo_resize(src.geo()));
    multiply_m(src1, sol, inv);
    src1 -= src;
    displayln_info(fname + ssprintf(": diff qnorm = %24.17E", qnorm(src1)));
    displayln_info(fname + ssprintf(": sol after qnorm = %24.17E", qnorm(sol)));
  }
}

inline long cg_with_herm_sym_2(FermionField5d& sol, const FermionField5d& src,
                               const InverterDomainWallGrid& inv,
                               const double stop_rsd = 1e-8,
                               const long max_num_iter = 50000)
{
  TIMER_VERBOSE_FLOPS("cg_with_herm_sym_2(5d,5d,InvGrid)");
  displayln_info(fname + ssprintf(": stop_rsd=%8.2E ; max_num_iter=%5ld",
                                  stop_rsd, max_num_iter));
  sol.init(geo_resize(src.geo()));
  qassert(sol.geo().eo == 1);
  qassert(src.geo().eo == 1);
  using namespace Grid;
  GridRedBlackCartesian* FrbGrid = inv.FrbGrid;
  LatticeFermionF gsrc(FrbGrid), gsol(FrbGrid);
  gsrc.Checkerboard() = Odd;
  gsol.Checkerboard() = Odd;
  grid_convert(gsrc, src);
  grid_convert(gsol, sol);
  ConjugateGradient<LatticeFermionF> CG(stop_rsd, max_num_iter, false);
  CG(*inv.HermOp, gsrc, gsol);
  grid_convert(sol, gsol);
  const long iter = CG.IterationsToComplete;
  timer.flops += 5500 * iter * inv.fa.ls * inv.geo().local_volume();
  return iter;
}

inline long invert(FermionField5d& sol, const FermionField5d& src,
                   const InverterDomainWallGrid& inv)
{
  return invert_with_cg(sol, src, inv, cg_with_herm_sym_2);
}

inline long invert(FermionField4d& sol, const FermionField4d& src,
                   const InverterDomainWallGrid& inv)
{
  return invert_dwf(sol, src, inv);
}

}  // namespace qlat

#endif
