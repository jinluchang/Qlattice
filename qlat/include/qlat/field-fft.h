#pragma once

#include <qlat/setup.h>
#include <qlat/field.h>
#include <qlat/geometry.h>
#include <qlat/mpi.h>
#include <qlat/utils-coordinate.h>

#include <qlat/field-shuffle.h>

#include <fftw3.h>

#include <vector>

namespace qlat
{  //

inline bool check_fft_plan_key(const Geometry& geo, const int mc, const int dir,
                               const bool is_forward)
{
  qassert(0 < geo.multiplicity);
  bool b = true;
  b = b && geo.eo == 0;
  b = b && geo.is_only_local;
  b = b && mc % geo.multiplicity == 0;
  b = b && 0 <= dir and dir < 4;
  b = b && (is_forward == true or is_forward == false);
  return b;
}

struct API FftComplexFieldPlan {
  Geometry geo;     // geo.is_only_local == true
  int mc;           // geo.multiplicity * sizeof(M) / sizeof(Complex)
  int dir;          // direction of the fft
  bool is_forward;  // is forward fft (forward \sum_x f[x] exp(-i k x))
  //
  fftw_plan fftplan;
  ShufflePlan sp;
  //
  FftComplexFieldPlan() {}
  //
  ~FftComplexFieldPlan() { end(); }
  //
  void end()
  {
    if (geo.initialized) {
      displayln_info("FftComplexFieldPlan::end(): free a plan.");
      fftw_destroy_plan(fftplan);
      geo.initialized = false;
    }
  }
  //
  void init(const Geometry& geo_, const int mc_, const int dir_,
            const bool is_forward_)
  {
    TIMER_VERBOSE("FftComplexFieldPlan::init");
    qassert(check_fft_plan_key(geo_, mc_, dir_, is_forward_));
    geo = geo_;
    mc = mc_;
    dir = dir_;
    is_forward = is_forward_;
    const int sizec = geo.total_site()[dir];
    const long nc = geo.local_volume() / geo.node_site[dir] * mc;
    const long chunk = ((nc / mc - 1) / geo.geon.size_node[dir] + 1) * mc;
    const long nc_start = std::min(nc, geo.geon.coor_node[dir] * chunk);
    const long nc_stop = std::min(nc, nc_start + chunk);
    const long nc_size = nc_stop - nc_start;
    // fftw_init_threads();
    // fftw_plan_with_nthreads(omp_get_max_threads());
    displayln_info(ssprintf("FftComplexFieldPlan::init: malloc %ld",
                            nc_size * sizec * sizeof(Complex)));
    Complex* fftdatac =
        (Complex*)fftw_malloc(nc_size * sizec * sizeof(Complex));
    const int rank = 1;
    const int n[1] = {sizec};
    const long howmany = nc_size;
    const long dist = 1;
    const long stride = nc_size;
    fftplan = fftw_plan_many_dft(
        rank, n, howmany, (fftw_complex*)fftdatac, n, stride, dist,
        (fftw_complex*)fftdatac, n, stride, dist,
        is_forward ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_free(fftdatac);
    displayln_info(ssprintf("FftComplexFieldPlan::init: free %ld",
                            nc_size * sizec * sizeof(Complex)));
    sp = make_shuffle_plan_fft(geo.total_site(), dir);
  }
};

API inline Cache<std::string, FftComplexFieldPlan>& get_fft_plan_cache()
{
  static Cache<std::string, FftComplexFieldPlan> cache("FftComplexFieldPlanCache", 32);
  return cache;
}

inline FftComplexFieldPlan& get_fft_plan(const Geometry& geo, const int mc,
                                         const int dir, const bool is_forward)
{
  TIMER("get_fft_plan");
  qassert(check_fft_plan_key(geo, mc, dir, is_forward));
  Cache<std::string, FftComplexFieldPlan>& cache = get_fft_plan_cache();
  const std::string key =
      ssprintf("%s %s %d %d %d %d %d", show(geo.node_site).c_str(),
               show(geo.geon.size_node).c_str(), geo.geon.id_node,
               geo.multiplicity, mc, dir, is_forward ? 1 : 0);
  if (cache.has(key)) {
    return cache[key];
  }
  FftComplexFieldPlan& plan = cache[key];
  plan.init(geo, mc, dir, is_forward);
  return plan;
}

template <class M>
void fft_complex_field_dir(Field<M>& field1, const Field<M>& field, const int dir,
                           const bool is_forward)
{
  TIMER("fft_complex_field_dir");
  Geometry geo = field.geo();
  geo.resize(0);
  const int mc = geo.multiplicity * sizeof(M) / sizeof(Complex);
  FftComplexFieldPlan& plan = get_fft_plan(geo, mc, dir, is_forward);
  fftw_plan& fftplan = plan.fftplan;
  const int sizec = geo.total_site()[dir];
  const long nc = geo.local_volume() / geo.node_site[dir] * mc;
  const long chunk = ((nc / mc - 1) / geo.geon.size_node[dir] + 1) * mc;
  const long nc_start = std::min(nc, geo.geon.coor_node[dir] * chunk);
  const long nc_stop = std::min(nc, nc_start + chunk);
  const long nc_size = nc_stop - nc_start;
  qassert(nc_size >= 0);
  const ShufflePlan& sp = plan.sp;
  std::vector<Field<M> > fft_fields;
  shuffle_field(fft_fields, field, sp);
  field1.init();
  Complex* fftdatac = (Complex*)fftw_malloc(nc_size * sizec * sizeof(Complex));
#pragma omp parallel for
  for (int i = 0; i < (int)fft_fields.size(); ++i) {
    if (not(get_data_size(fft_fields[i]) == nc_size * (int)sizeof(Complex))) {
      displayln(fname +
                ssprintf(": get_data_size=%d ; nc_size*sizeof(Complex)=%d",
                         get_data_size(fft_fields[i]),
                         nc_size * (int)sizeof(Complex)));
      qassert(get_data_size(fft_fields[i]) == nc_size * (int)sizeof(Complex));
    }
    std::memcpy((void*)&fftdatac[nc_size * i],
                (void*)get_data(fft_fields[i]).data(),
                get_data_size(fft_fields[i]));
  }
  {
    TIMER("fft_complex_field_dir-fftw");
    fftw_execute_dft(fftplan, (fftw_complex*)fftdatac, (fftw_complex*)fftdatac);
  }
#pragma omp parallel for
  for (int i = 0; i < (int)fft_fields.size(); ++i) {
    std::memcpy((void*)get_data(fft_fields[i]).data(),
                (void*)&fftdatac[nc_size * i], get_data_size(fft_fields[i]));
  }
  fftw_free(fftdatac);
  field1.init(geo);
  shuffle_field_back(field1, fft_fields, sp);
}

template <class M>
void fft_complex_field_dir(Field<M>& field, const int dir,
                           const bool is_forward)
{
  fft_complex_field_dir(field, field, dir, is_forward);
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
    fft_complex_field_dir(field, dir, is_forward);
  }
}

}  // namespace qlat
