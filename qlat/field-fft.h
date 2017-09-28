#pragma once

#include <qlat/config.h>
#include <qlat/utils.h>
#include <qlat/utils-coordinate.h>
#include <qlat/mpi.h>
#include <qlat/geometry.h>
#include <qlat/field.h>

#include <fftw3.h>
#include <omp.h>

#include <vector>

QLAT_START_NAMESPACE

struct fft_complex_field_plan
{
  Geometry geo;    // geo.is_only_local == true
  int mc;          // geo.multiplicity * sizeof(M) / sizeof(Complex)
  Coordinate dirs; // 0 is no transform, 1 is forward transform, -1 is backward transform
  //
  virtual const std::string& cname()
  {
    static const std::string s = "fft_complex_field_plan";
    return s;
  }
  //
  bool is_match(const Geometry& geo_, const int mc_, const Coordinate& dirs_)
  {
    return geo_ == geo && mc_ == mc && dirs_ == dirs;
  }
  //
  static bool check(const Geometry& geo_, const int mc_, const Coordinate& dirs_)
  {
    qassert(0 < geo_.multiplicity);
    bool b = true;
    b = b && geo_.is_only_local();
    b = b && mc_ % geo_.multiplicity == 0;
    for (int mu = 0; mu < 4; mu++) {
      b = b && (dirs_[mu] == -1 || dirs_[mu] == 0 || dirs_[mu] == 1);
    }
    return b;
  }
  //
  static fft_complex_field_plan& get_plan(const Geometry& geo_, const int mc_, const Coordinate dirs_)
  {
    TIMER("fft_complex_field_plan::get_plan");
    qassert(check(geo_, mc_, dirs_));
    static std::vector<fft_complex_field_plan> planV(100);
    static int next_plan_index = 0;
    for (int i = 0; i < (int)planV.size(); i++) {
      if (planV[i].is_match(geo_, mc_, dirs_)) {
        return planV[i];
      }
    }
    DisplayInfo("", fname, "start to make a new fft plan with id = %d\n", next_plan_index);
    fft_complex_field_plan& plan = planV[next_plan_index];
    next_plan_index++;
    next_plan_index %= planV.size();
    plan.end();
    plan.init(geo_, mc_, dirs_);
    return plan;
  }
  //
  fft_complex_field_plan()
  {
  }
  //
  ~fft_complex_field_plan()
  {
    end();
  }
  //
  void end()
  {
    if (geo.initialized) {
      displayln_info(cname() + "::end(): free a plan.");
      fftw_destroy_plan(fftplan);
      geo.initialized = false;
    }
  }
  //
  void init(const Geometry& geo_, const int mc_, const Coordinate dirs_)
  {
    TIMER_VERBOSE("fft_complex_field_plan::init");
    qassert(check(geo_, mc_, dirs_));
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
    const int sizec = geo.total_site()[dir];
    const int nc = geo.local_volume() / geo.node_site[dir] * mc;
    const int chunk = (nc-1) / geo.geon.size_node[dir]+1;
    const int nc_start = std::min(nc, geo.geon.coor_node[dir] * chunk);
    const int nc_stop = std::min(nc, nc_start + chunk);
    const int nc_size = nc_stop - nc_start;
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
    DisplayInfo("fft_complex_field_plan", "init", "malloc %d\n", nc_size * sizec * sizeof(Complex));
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
    DisplayInfo("fft_complex_field_plan", "init", "free %d\n", nc_size * sizec * sizeof(Complex));
  }
  //
  fftw_plan fftplan;
};

template<class M>
void fft_complex_field_dirs(Field<M>& field, const Coordinate& dirs)
{
  TIMER("fft_complex_field_dirs");
  Geometry geo = field.geo; geo.resize(0);
  const int mc = geo.multiplicity * sizeof(M) / sizeof(Complex);
  fft_complex_field_plan& plan = fft_complex_field_plan::get_plan(geo, mc, dirs);
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
  const int sizec = geo.total_site()[dir];
  const int nc = geo.local_volume() / geo.node_site[dir] * mc;
  const int chunk = (nc-1)/geo.geon.size_node[dir]+1;
  const int nc_start = std::min(nc, geo.geon.coor_node[dir] * chunk);
  const int nc_stop = std::min(nc, nc_start + chunk);
  const int nc_size = nc_stop - nc_start;
  Complex* fftdatac = (Complex*)fftw_malloc(nc_size * sizec * sizeof(Complex));
  Field<M> fields; fields.init(geo);
  Field<M> fieldr; fieldr.init(geo);
  Geometry geos = geo;
  // FIXME my compiler says unused variable 'fieldsize'
  // const int fieldsize = get_data_size(fields) / sizeof(double);
  fields = field;
  for (int i = 0; i < geos.geon.size_node[dir]; i++) {
#pragma omp parallel for
    for (long index = 0; index < geos.local_volume(); index++) {
      Coordinate xl = geos.coordinate_from_index(index);
      Coordinate xg = geos.coordinate_g_from_l(xl);
      int nc_index = 0;
      int nc_offset = mc;
      for (int mu = 0; mu < 4; mu++) {
        if (dir != mu) {
          nc_index += xl[mu] * nc_offset;
          nc_offset *= geos.node_site[mu];
        }
      }
      Complex* pc = (Complex*)fields.get_elems(xl).data();
      for (int nci = std::max(0, nc_start-nc_index); nci < std::min(mc, nc_stop-nc_index); nci++) {
        fftdatac[nci+nc_index-nc_start + nc_size*xg[dir]] = pc[nci];
      }
    }
    if (i == geos.geon.size_node[dir] - 1) {
      break;
    }
    {
      TIMER_FLOPS("fft_complex_field_dirs-get_data");
      timer.flops += get_data_size(fields);
      get_data_plus_mu(get_data(fieldr), get_data(fields), dir);
    }
    std::swap(fields, fieldr);
    geos.geon.coor_node[dir] = mod(geos.geon.coor_node[dir] + 1, geos.geon.size_node[dir]);
  }
  {
    TIMER("fft_complex_field_dirs-fftw");
    fftw_execute_dft(fftplan, (fftw_complex*)fftdatac, (fftw_complex*)fftdatac);
  }
  geos.geon.coor_node[dir] = mod(geo.geon.coor_node[dir] + 1, geos.geon.size_node[dir]);
  for (int i = 0; i < geos.geon.size_node[dir]; i++) {
#pragma omp parallel for
    for (long index = 0; index < geos.local_volume(); index++) {
      Coordinate xl = geos.coordinate_from_index(index);
      Coordinate xg = geos.coordinate_g_from_l(xl);
      int nc_index = 0;
      int nc_offset = mc;
      for (int mu = 0; mu < 4; mu++) {
        if (dir != mu) {
          nc_index += xl[mu] * nc_offset;
          nc_offset *= geos.node_site[mu];
        }
      }
      Complex* pc = (Complex*)fields.get_elems(xl).data();
      for (int nci = std::max(0, nc_start-nc_index); nci < std::min(mc, nc_stop-nc_index); nci++) {
        pc[nci] = fftdatac[nci+nc_index-nc_start + nc_size*xg[dir]];
      }
    }
    if (i == geos.geon.size_node[dir] - 1) {
      break;
    }
    {
      TIMER_FLOPS("fft_complex_field_dirs-get_data");
      timer.flops += get_data_size(fields);
      get_data_plus_mu(get_data(fieldr), get_data(fields), dir);
    }
    std::swap(fields, fieldr);
    geos.geon.coor_node[dir] = mod(geos.geon.coor_node[dir] + 1, geos.geon.size_node[dir]);
  }
  field = fields;
  fftw_free(fftdatac);
}

template<class M>
void fft_complex_field(Field<M>& field, const bool isForward = true)
{
  TIMER("fft_complex_field");
  // forward compute
  // field(k) <- \sum_{x} exp( - ii * 2 pi * k * x ) field(x)
  // backwards compute
  // field(x) <- \sum_{k} exp( + ii * 2 pi * k * x ) field(k)
  Coordinate dirs;
  for (int dir = 0; dir < 4; ++dir) {
    set_zero(dirs);
    dirs[dir] = isForward ? 1 : -1;
    fft_complex_field_dirs(field, dirs);
  }
}

QLAT_END_NAMESPACE
