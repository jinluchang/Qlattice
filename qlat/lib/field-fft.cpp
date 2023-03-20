#define QLAT_INSTANTIATE_FIELD_FFT

#include <qlat/field-fft.h>

namespace qlat
{  //

bool check_fft_plan_key(const Geometry& geo, const int mc, const int dir,
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

void FftComplexFieldPlan::init()
{
  if (geo.initialized) {
    displayln_info("FftComplexFieldPlan::end(): free a plan.");
    fftw_destroy_plan(fftplan);
    geo.init();
  }
}

void FftComplexFieldPlan::init(const Geometry& geo_, const int mc_,
                               const int dir_, const bool is_forward_)
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
  Complex* fftdatac = (Complex*)fftw_malloc(nc_size * sizec * sizeof(Complex));
  const int rank = 1;
  const int n[1] = {sizec};
  const long howmany = nc_size;
  const long dist = 1;
  const long stride = nc_size;
  fftplan = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)fftdatac, n,
                               stride, dist, (fftw_complex*)fftdatac, n, stride,
                               dist, is_forward ? FFTW_FORWARD : FFTW_BACKWARD,
                               FFTW_ESTIMATE);
  fftw_free(fftdatac);
  displayln_info(ssprintf("FftComplexFieldPlan::init: free %ld",
                          nc_size * sizec * sizeof(Complex)));
  sp = make_shuffle_plan_fft(geo.total_site(), dir);
}

FftComplexFieldPlan& get_fft_plan(const Geometry& geo, const int mc,
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

}  // namespace qlat
