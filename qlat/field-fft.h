#pragma once

#include <qlat/config.h>
#include <qlat/field.h>
#include <qlat/geometry.h>
#include <qlat/mpi.h>
#include <qlat/utils-coordinate.h>
#include <qlat/utils.h>

#include <fftw3.h>

#include <vector>

QLAT_START_NAMESPACE

struct fft_complex_field_plan {
  Geometry geo;     // geo.is_only_local == true
  int mc;           // geo.multiplicity * sizeof(M) / sizeof(Complex)
  int dir;          // direction of the fft
  bool is_forward;  // is forward fft (forward \sum_x f[x] exp(-i k x))
  //
  virtual const std::string& cname()
  {
    static const std::string s = "fft_complex_field_plan";
    return s;
  }
  //
  bool is_match(const Geometry& geo_, const int mc_, const int dir_, const bool is_forward_)
  {
    return geo_ == geo && mc_ == mc && dir_ == dir && is_forward_ == is_forward;
  }
  //
  static bool check(const Geometry& geo_, const int mc_,
                    const int dir_, const bool is_forward_)
  {
    qassert(0 < geo_.multiplicity);
    bool b = true;
    b = b && geo_.is_only_local();
    b = b && mc_ % geo_.multiplicity == 0;
    b = b && 0 <= dir_ and dir_ < 4;
    b = b && (is_forward_ == true or is_forward_ == false);
    return b;
  }
  //
  static fft_complex_field_plan& get_plan(const Geometry& geo_, const int mc_,
                                          const int dir_,
                                          const bool is_forward_)
  {
    TIMER("fft_complex_field_plan::get_plan");
    qassert(check(geo_, mc_, dir_, is_forward_));
    static std::vector<fft_complex_field_plan> planV(100);
    static int next_plan_index = 0;
    for (int i = 0; i < (int)planV.size(); i++) {
      if (planV[i].is_match(geo_, mc_, dir_, is_forward_)) {
        return planV[i];
      }
    }
    DisplayInfo("", fname, "start to make a new fft plan with id = %d\n",
                next_plan_index);
    fft_complex_field_plan& plan = planV[next_plan_index];
    next_plan_index++;
    next_plan_index %= planV.size();
    plan.end();
    plan.init(geo_, mc_, dir_, is_forward_);
    return plan;
  }
  //
  fft_complex_field_plan() {}
  //
  ~fft_complex_field_plan() { end(); }
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
  void init(const Geometry& geo_, const int mc_, const int dir_, const bool is_forward_)
  {
    TIMER_VERBOSE("fft_complex_field_plan::init");
    qassert(check(geo_, mc_, dir_, is_forward_));
    geo = geo_;
    mc = mc_;
    dir = dir_;
    is_forward = is_forward_;
    const int sizec = geo.total_site()[dir];
    const int nc = geo.local_volume() / geo.node_site[dir] * mc;
    const int chunk = (nc - 1) / geo.geon.size_node[dir] + 1;
    const int nc_start = std::min(nc, geo.geon.coor_node[dir] * chunk);
    const int nc_stop = std::min(nc, nc_start + chunk);
    const int nc_size = nc_stop - nc_start;
    // fftw_init_threads();
    // fftw_plan_with_nthreads(omp_get_max_threads());
    DisplayInfo("fft_complex_field_plan", "init", "malloc %d\n",
                nc_size * sizec * sizeof(Complex));
    Complex* fftdatac =
        (Complex*)fftw_malloc(nc_size * sizec * sizeof(Complex));
    const int rank = 1;
    const int n[1] = {sizec};
    const int howmany = nc_size;
    const int dist = 1;
    const int stride = nc_size;
    fftplan = fftw_plan_many_dft(
        rank, n, howmany, (fftw_complex*)fftdatac, n, stride, dist,
        (fftw_complex*)fftdatac, n, stride, dist,
        is_forward ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_free(fftdatac);
    DisplayInfo("fft_complex_field_plan", "init", "free %d\n",
                nc_size * sizec * sizeof(Complex));
  }
  //
  fftw_plan fftplan;
};

template <class M>
void fft_complex_field_dir(Field<M>& field, const int dir, const bool is_forward)
{
  TIMER("fft_complex_field_dir");
  Geometry geo = field.geo;
  geo.resize(0);
  const int mc = geo.multiplicity * sizeof(M) / sizeof(Complex);
  fft_complex_field_plan& plan =
      fft_complex_field_plan::get_plan(geo, mc, dir, is_forward);
  fftw_plan& fftplan = plan.fftplan;
  const int sizec = geo.total_site()[dir];
  const int nc = geo.local_volume() / geo.node_site[dir] * mc;
  const int chunk = (nc - 1) / geo.geon.size_node[dir] + 1;
  const int nc_start = std::min(nc, geo.geon.coor_node[dir] * chunk);
  const int nc_stop = std::min(nc, nc_start + chunk);
  const int nc_size = nc_stop - nc_start;
  Complex* fftdatac = (Complex*)fftw_malloc(nc_size * sizec * sizeof(Complex));
  Field<M> fields;
  fields.init(geo);
  Field<M> fieldr;
  fieldr.init(geo);
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
      for (int nci = std::max(0, nc_start - nc_index);
           nci < std::min(mc, nc_stop - nc_index); nci++) {
        fftdatac[nci + nc_index - nc_start + nc_size * xg[dir]] = pc[nci];
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
    qswap(fields, fieldr);
    geos.geon.coor_node[dir] =
        mod(geos.geon.coor_node[dir] + 1, geos.geon.size_node[dir]);
  }
  {
    TIMER("fft_complex_field_dirs-fftw");
    fftw_execute_dft(fftplan, (fftw_complex*)fftdatac, (fftw_complex*)fftdatac);
  }
  geos.geon.coor_node[dir] =
      mod(geo.geon.coor_node[dir] + 1, geos.geon.size_node[dir]);
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
      for (int nci = std::max(0, nc_start - nc_index);
           nci < std::min(mc, nc_stop - nc_index); nci++) {
        pc[nci] = fftdatac[nci + nc_index - nc_start + nc_size * xg[dir]];
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
    qswap(fields, fieldr);
    geos.geon.coor_node[dir] =
        mod(geos.geon.coor_node[dir] + 1, geos.geon.size_node[dir]);
  }
  field = fields;
  fftw_free(fftdatac);
}

template <class M>
void fft_complex_field(Field<M>& field, const bool is_forward = true)
{
  TIMER_FLOPS("fft_complex_field");
  timer.flops += get_data(field).data_size() * get_num_node();
  // forward compute
  // field(k) <- \sum_{x} exp( - ii * 2 pi * k * x ) field(x)
  // backwards compute
  // field(x) <- \sum_{k} exp( + ii * 2 pi * k * x ) field(k)
  for (int dir = 0; dir < 4; ++dir) {
    fft_complex_field_dir(field, dir, is_forward);
  }
}

template <class M>
void fft_complex_field_spatial(Field<M>& field, const bool is_forward = true)
{
  TIMER_FLOPS("fft_complex_field_spatial");
  timer.flops += get_data(field).data_size() * get_num_node();
  // forward compute
  // field(k) <- \sum_{x} exp( - ii * 2 pi * k * x ) field(x)
  // backwards compute
  // field(x) <- \sum_{k} exp( + ii * 2 pi * k * x ) field(k)
  for (int dir = 0; dir < 3; ++dir) {
    fft_complex_field_dirs(field, dir, is_forward);
  }
}

QLAT_END_NAMESPACE
