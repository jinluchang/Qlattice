#pragma once

#define QLAT_GRID

#include <qlat/qlat.h>

#include <Grid/Grid.h>

#include <cstdlib>

QLAT_START_NAMESPACE

inline std::vector<int> grid_convert(const Coordinate& x)
{
  std::vector<int> ret(4);
  for (int mu = 0; mu < 4; ++mu) {
    ret[mu] = x[mu];
  }
  return ret;
}

inline std::vector<int> grid_convert(const Coordinate& x, const int m)
{
  std::vector<int> ret(5);
  for (int mu = 0; mu < 4; ++mu) {
    ret[mu + 1] = x[mu];
  }
  ret[0] = m;
  return ret;
}

inline Coordinate grid_convert(const std::vector<int>& x)
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
  using namespace Grid::QCD;
  const std::vector<int>& mpi_layout = UGrid->_processors;
  const std::vector<int>& mpi_corr = UGrid->_processor_coor;
  const Coordinate size_node = grid_convert(mpi_layout);
  const Coordinate coor_node = grid_convert(mpi_corr);
  return index_from_coordinate(coor_node, size_node);
}

inline void grid_begin(int* argc, char** argv[], const Coordinate& size_node_ = Coordinate())
{
  using namespace Grid;
  using namespace Grid::QCD;
  system("rm /dev/shm/Grid*");
  Grid_init(argc, argv);
  const int num_node = init_mpi(argc, argv);
  const Coordinate size_node = size_node_ == Coordinate() ? plan_size_node(num_node) : size_node_;
  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      grid_convert(size_node * 4),
      GridDefaultSimd(Nd, vComplex::Nsimd()),
      grid_convert(size_node));
  UGrid->show_decomposition();
  displayln_info(ssprintf("GridThreads::GridThreads() = %d", GridThread::GetThreads()));
  const int id_node = id_node_from_grid(UGrid);
  delete UGrid;
  begin(id_node, size_node);
}

void grid_end()
{
  Grid::Grid_finalize();
  system("rm /dev/shm/Grid*");
}

void grid_convert(Grid::QCD::LatticeGaugeField& ggf, const GaugeField& gf)
{
  TIMER_VERBOSE("grid_convert(ggf,gf)");
  using namespace Grid;
  using namespace Grid::QCD;
  const Geometry& geo = gf.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    std::vector<int> coor = grid_convert(xl);
    const Vector<ColorMatrix> ms = gf.get_elems_const(xl);
    LorentzColourMatrix gms;
    qassert(sizeof(LorentzColourMatrix) == ms.data_size());
    memcpy(&gms, ms.data(), ms.data_size());
    pokeLocalSite(gms, ggf, coor);
  }
}

void grid_convert(FermionField5d& ff, const Grid::QCD::LatticeFermion& gff)
  // ff need to be initialized with correct geo
{
  TIMER_VERBOSE("grid_convert(ff,gff)");
  using namespace Grid;
  using namespace Grid::QCD;
  const Geometry& geo = ff.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    std::vector<int> coor = grid_convert(xl, 0);
    Vector<WilsonVector> wvs = ff.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      coor[0] = m;
      peekLocalSite(wvs[m], gff, coor);
    }
  }
}

void grid_convert(Grid::QCD::LatticeFermion& gff, const FermionField5d& ff)
{
  TIMER_VERBOSE("grid_convert(gff,ff)");
  using namespace Grid;
  using namespace Grid::QCD;
  const Geometry& geo = ff.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    std::vector<int> coor = grid_convert(xl, 0);
    const Vector<WilsonVector> wvs = ff.get_elems_const(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      coor[0] = m;
      pokeLocalSite(wvs[m], gff, coor);
    }
  }
}

struct InverterDomainWallGrid
{
  Geometry geo;
  FermionAction fa;
  GaugeField gf;
  InverterParams ip;
  //
  Grid::GridCartesian* UGrid = NULL;
  Grid::GridRedBlackCartesian* UrbGrid = NULL;
  Grid::GridCartesian* FGrid = NULL;
  Grid::GridRedBlackCartesian* FrbGrid = NULL;
  Grid::QCD::LatticeGaugeField* Umu = NULL;
  Grid::QCD::DomainWallFermionR* Ddwf = NULL;
  Grid::ConjugateGradient<Grid::QCD::LatticeFermion>* CG = NULL;
  Grid::QCD::SchurRedBlackDiagMooeeSolve<Grid::QCD::LatticeFermion>* SchurSolver = NULL;
  //
  InverterDomainWallGrid()
  {
    init();
  }
  ~InverterDomainWallGrid()
  {
    init();
  }
  //
  void init()
  {
    ip.init();
    free();
  }
  //
  void setup()
  {
    TIMER_VERBOSE("InvGrid::setup");
    using namespace Grid;
    using namespace Grid::QCD;
    free();
    const Coordinate total_site = geo.total_site();
    const Coordinate size_node = geo.geon.size_node;
    UGrid = SpaceTimeGrid::makeFourDimGrid(grid_convert(total_site),
        GridDefaultSimd(Nd,vComplex::Nsimd()), grid_convert(size_node));
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
    FGrid = SpaceTimeGrid::makeFiveDimGrid(fa.ls, UGrid);
    FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(fa.ls, UGrid);
    qassert(geo.geon.id_node == id_node_from_grid(UGrid));
    Umu = new LatticeGaugeField(UGrid);
    grid_convert(*Umu, gf);
    Ddwf = new DomainWallFermionR(*Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, fa.mass, fa.m5);
    CG = new ConjugateGradient<LatticeFermion>(ip.stop_rsd, ip.max_num_iter);
    SchurSolver = new SchurRedBlackDiagMooeeSolve<LatticeFermion>(*CG);
  }
  void setup(const GaugeField& gf_, const FermionAction& fa_)
  {
    geo = geo_reform(gf_.geo);
    gf.init();
    gf.init(geo);
    gf = gf_;
    fa = fa_;
    setup();
  }
  void setup(const GaugeField& gf_, const FermionAction& fa_, const LowModes& lm_)
  {
    setup(gf_, fa_);
    // TODO
  }
  //
  void free()
  {
    if (UGrid != NULL) {
      delete UGrid;
      UGrid = NULL;
    }
    if (UrbGrid != NULL) {
      delete UrbGrid;
      UrbGrid = NULL;
    }
    if (FGrid != NULL) {
      delete FGrid;
      FGrid = NULL;
    }
    if (FrbGrid != NULL) {
      delete FrbGrid;
      FrbGrid = NULL;
    }
    if (Umu != NULL) {
      delete Umu;
      Umu = NULL;
    }
    if (Ddwf != NULL) {
      delete Ddwf;
      Ddwf = NULL;
    }
    if (CG != NULL) {
      delete CG;
      CG = NULL;
    }
    if (SchurSolver != NULL) {
      delete SchurSolver;
      SchurSolver = NULL;
    }
  }
  //
  double& stop_rsd()
  {
    return ip.stop_rsd;
  }
  //
  long& max_num_iter()
  {
    return ip.max_num_iter;
  }
  //
  long& max_mixed_precision_cycle()
  {
    return ip.max_mixed_precision_cycle;
  }
};

inline void setup_inverter(InverterDomainWallGrid& inv)
{
  inv.setup();
}

inline void setup_inverter(InverterDomainWallGrid& inv, const GaugeField& gf, const FermionAction& fa)
{
  inv.setup(gf, fa);
}

inline void setup_inverter(InverterDomainWallGrid& inv, const GaugeField& gf, const FermionAction& fa, const LowModes& lm)
{
  inv.setup(gf, fa, lm);
}

inline void inverse(FermionField5d& sol, const FermionField5d& src, const InverterDomainWallGrid& inv)
  // sol do not need to be initialized
{
  TIMER_VERBOSE("inverse(5d,5d,InvGrid)");
  sol.init(geo_resize(src.geo));
  using namespace Grid;
  using namespace Grid::QCD;
  GridCartesian* FGrid = inv.FGrid;
  LatticeFermion gsrc(FGrid), gsol(FGrid);
  grid_convert(gsrc, src);
  grid_convert(gsol, sol);
  if (true) {
    FermionField5d ff;
    ff.init(geo_resize(src.geo));
    grid_convert(ff, gsrc);
    displayln_info(fname + ssprintf(": src norm = %24.17E", norm(ff)));
    grid_convert(ff, gsol);
    displayln_info(fname + ssprintf(": sol norm = %24.17E", norm(ff)));
  }
  (*inv.SchurSolver)(*inv.Ddwf, gsrc, gsol);
  grid_convert(sol, gsol);
}

inline void inverse(FermionField4d& sol, const FermionField4d& src, const InverterDomainWallGrid& inv)
{
  inverse_dwf(sol, src, inv);
}

inline void multiply_m(FermionField5d& out, const FermionField5d& in, const InverterDomainWallGrid& inv)
{
  TIMER("multiply_m(5d,5d,InvGrid)");
  out.init(geo_resize(in.geo));
  using namespace Grid;
  using namespace Grid::QCD;
  GridCartesian* FGrid = inv.FGrid;
  LatticeFermion gin(FGrid), gout(FGrid);
  grid_convert(gin, in);
  (*inv.Ddwf).M(gin, gout);
  grid_convert(out, gout);
}

QLAT_END_NAMESPACE

