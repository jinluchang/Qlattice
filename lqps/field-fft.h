#pragma once

#include <lqps/config.h>
#include <lqps/utils.h>
#include <lqps/mpi.h>
#include <lqps/geometry.h>
#include <lqps/field.h>

#include <fftw3.h>
#include <omp.h>

#include <vector>

LQPS_START_NAMESPACE

struct fftComplexFieldPlan
{
  Geometry geo;    // geo.isOnlyLocal == true
  int mc;          // geo.multiplicity * sizeof(M) / sizeof(Complex)
  Coordinate dirs; // 0 is no transform, 1 is forward transform, -1 is backward transform
  //
  virtual const char* cname()
  {
    return "fftComplexFieldPlan";
  }
  //
  bool isMatch(const Geometry& geo_, const int mc_, const Coordinate& dirs_)
  {
    return geo_ == geo && mc_ == mc && dirs_ == dirs;
  }
  //
  static bool check(const Geometry& geo_, const int mc_, const Coordinate& dirs_)
  {
    assert(0 < geo_.multiplicity);
    bool b = true;
    b = b && geo_.isOnlyLocal();
    b = b && mc_ % geo_.multiplicity == 0;
    for (int mu = 0; mu < 4; mu++) {
      b = b && (dirs_[mu] == -1 || dirs_[mu] == 0 || dirs_[mu] == 1);
    }
    return b;
  }
  //
  static fftComplexFieldPlan& getPlan(const Geometry& geo_, const int mc_, const Coordinate dirs_)
  {
    TIMER("fftComplexFieldPlan::getPlan");
    assert(check(geo_, mc_, dirs_));
    static std::vector<fftComplexFieldPlan> planV(100);
    static int next_plan_index = 0;
    for (int i = 0; i < planV.size(); i++) {
      if (planV[i].isMatch(geo_, mc_, dirs_)) {
        return planV[i];
      }
    }
    DisplayInfo("fftComplexFieldPlan", "getPlan", "start to make a new fft plan with id = %d\n", next_plan_index);
    fftComplexFieldPlan& plan = planV[next_plan_index];
    next_plan_index++;
    next_plan_index %= planV.size();
    plan.end();
    plan.init(geo_, mc_, dirs_);
    return plan;
  }
  //
  fftComplexFieldPlan()
  {
  }
  //
  ~fftComplexFieldPlan()
  {
    end();
  }
  //
  void end()
  {
    if (geo.initialized) {
      DisplayInfo(cname(), "end", "free a plan\n");
      fftw_destroy_plan(fftplan);
      geo.initialized = false;
    }
  }
  //
  void init(const Geometry& geo_, const int mc_, const Coordinate dirs_)
  {
    TIMER_VERBOSE("fftComplexFieldPlan::init");
    assert(check(geo_, mc_, dirs_));
    geo = geo_;
    mc = mc_;
    dirs = dirs_;
    // FIXME currently can only transform in one direction
    int dir = 0;
    bool isForward = true;
    for (int mu = 0; mu < 4; mu++) {
      if (0 != dirs[mu]) {
        dir = mu;
        isForward = dirs[mu] == 1;
        break;
      }
    }
    const int sizec = geo.totalSite(dir);
    const int nc = geo.localVolume() / geo.nodeSite[dir] * mc;
    const int chunk = (nc-1) / geo.geon.sizeNode[dir]+1;
    const int nc_start = std::min(nc, geo.geon.coorNode[dir] * chunk);
    const int nc_stop = std::min(nc, nc_start + chunk);
    const int nc_size = nc_stop - nc_start;
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
    DisplayInfo("fftComplexFieldPlan", "init", "malloc %d\n", nc_size * sizec * sizeof(Complex));
    Complex* fftdatac = (Complex*)fftw_malloc(nc_size * sizec * sizeof(Complex));
    const int rank = 1;
    const int n[1] = { sizec };
    const int howmany = nc_size;
    const int dist = 1;
    const int stride = nc_size;
    fftplan = fftw_plan_many_dft(rank, n, howmany,
        (fftw_complex*)fftdatac, n, stride, dist,
        (fftw_complex*)fftdatac, n, stride, dist,
        isForward ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_free(fftdatac);
    DisplayInfo("fftComplexFieldPlan", "init", "free %d\n", nc_size * sizec * sizeof(Complex));
  }
  //
  fftw_plan fftplan;
};

template<class M>
void fftComplexFieldDirs(Field<M>& field, const Coordinate& dirs)
{
  TIMER("fftComplexFieldDirs");
  Geometry geo = field.geo; geo.resize(0);
  const int mc = geo.multiplicity * sizeof(M) / sizeof(Complex);
  fftComplexFieldPlan& plan = fftComplexFieldPlan::getPlan(geo, mc, dirs);
  fftw_plan& fftplan = plan.fftplan;
  // FIXME currently can only transform in one direction
  int dir = 0;
  bool isForward = true;
  for (int mu = 0; mu < 4; mu++) {
    if (0 != dirs[mu]) {
      dir = mu;
      isForward = dirs[mu] == 1;
      break;
    }
  }
  const int sizec = geo.totalSite(dir);
  const int nc = geo.localVolume() / geo.nodeSite[dir] * mc;
  const int chunk = (nc-1)/geo.geon.sizeNode[dir]+1;
  const int nc_start = std::min(nc, geo.geon.coorNode[dir] * chunk);
  const int nc_stop = std::min(nc, nc_start + chunk);
  const int nc_size = nc_stop - nc_start;
  Complex* fftdatac = (Complex*)fftw_malloc(nc_size * sizec * sizeof(Complex));
  Field<M> fields; fields.init(geo);
  Field<M> fieldr; fieldr.init(geo);
  Geometry geos; geos.init(geo);
  const int fieldsize = getDataSize(fields) / sizeof(double);
  fields = field;
  for (int i = 0; i < geos.geon.sizeNode[dir]; i++) {
#pragma omp parallel for
    for (long index = 0; index < geos.localVolume(); index++) {
      Coordinate xl; geos.coordinateFromIndex(xl, index);
      Coordinate xg; geos.coordinateGfL(xg, xl);
      int nc_index = 0;
      int nc_offset = mc;
      for (int mu = 0; mu < 4; mu++) {
        if (dir != mu) {
          nc_index += xl[mu] * nc_offset;
          nc_offset *= geos.nodeSite[mu];
        }
      }
      Complex* pc = (Complex*)fields.getElems(xl).data();
      for (int nci = std::max(0, nc_start-nc_index); nci < std::min(mc, nc_stop-nc_index); nci++) {
        fftdatac[nci+nc_index-nc_start + nc_size*xg[dir]] = pc[nci];
      }
    }
    if (i == geos.geon.sizeNode[dir] - 1) {
      break;
    }
    {
      TIMER_FLOPS("fftComplexFieldDirs_getData");
      timer.flops += getDataSize(fields);
      getDataPlusMu(getData(fieldr), getData(fields), dir);
    }
    swap(fields, fieldr);
    geos.geon.coorNode[dir] = mod(geos.geon.coorNode[dir] + 1, geos.geon.sizeNode[dir]);
  }
  {
    TIMER("fftComplexFieldDirs_fftw");
    fftw_execute_dft(fftplan, (fftw_complex*)fftdatac, (fftw_complex*)fftdatac);
  }
  geos.geon.coorNode[dir] = mod(geo.geon.coorNode[dir] + 1, geos.geon.sizeNode[dir]);
  for (int i = 0; i < geos.geon.sizeNode[dir]; i++) {
#pragma omp parallel for
    for (long index = 0; index < geos.localVolume(); index++) {
      Coordinate xl; geos.coordinateFromIndex(xl, index);
      Coordinate xg; geos.coordinateGfL(xg, xl);
      int nc_index = 0;
      int nc_offset = mc;
      for (int mu = 0; mu < 4; mu++) {
        if (dir != mu) {
          nc_index += xl[mu] * nc_offset;
          nc_offset *= geos.nodeSite[mu];
        }
      }
      Complex* pc = (Complex*)fields.getElems(xl).data();
      for (int nci = std::max(0, nc_start-nc_index); nci < std::min(mc, nc_stop-nc_index); nci++) {
        pc[nci] = fftdatac[nci+nc_index-nc_start + nc_size*xg[dir]];
      }
    }
    if (i == geos.geon.sizeNode[dir] - 1) {
      break;
    }
    {
      TIMER_FLOPS("fftComplexFieldDirs_getData");
      timer.flops += getDataSize(fields);
      getDataPlusMu(getData(fieldr), getData(fields), dir);
    }
    swap(fields, fieldr);
    geos.geon.coorNode[dir] = mod(geos.geon.coorNode[dir] + 1, geos.geon.sizeNode[dir]);
  }
  field = fields;
  fftw_free(fftdatac);
}

template<class M>
void fftComplexField(Field<M>& field, const bool isForward = true)
{
  TIMER("fftComplexField");
  // forward compute
  // field(k) <- \sum_{x} exp( - ii * 2 pi * k * x ) field(x)
  // backwards compute
  // field(x) <- \sum_{k} exp( + ii * 2 pi * k * x ) field(k)
  Coordinate dirs;
  for (int dir = 0; dir < 4; ++dir) {
    setZero(dirs);
    dirs[dir] = isForward ? 1 : -1;
    fftComplexFieldDirs(field, dirs);
  }
}

LQPS_END_NAMESPACE
