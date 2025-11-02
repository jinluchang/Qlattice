#define QLAT_INSTANTIATE_FIELD_FFT

#include <qlat/field-fft.h>

namespace qlat
{  //

bool check_fft_plan_key(const Geometry& geo, const Int mc, const Int dir,
                        const bool is_forward)
{
  (void)mc;
  bool b = true;
  b = b && geo.eo == 0;
  b = b && geo.is_only_local;
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

void FftComplexFieldPlan::init(const Geometry& geo_, const Int mc_,
                               const Int dir_, const bool is_forward_)
{
  TIMER_VERBOSE("FftComplexFieldPlan::init");
  qassert(check_fft_plan_key(geo_, mc_, dir_, is_forward_));
  geo = geo_;
  mc = mc_;
  dir = dir_;
  is_forward = is_forward_;
  const Int sizec = geo.total_site()[dir];
  const Long nc = geo.local_volume() / geo.node_site[dir] * mc;
  const Long chunk = ((nc / mc - 1) / geo.geon.size_node[dir] + 1) * mc;
  const Long nc_start = std::min(nc, geo.geon.coor_node[dir] * chunk);
  const Long nc_stop = std::min(nc, nc_start + chunk);
  const Long nc_size = nc_stop - nc_start;
  // fftw_init_threads();
  // fftw_plan_with_nthreads(omp_get_max_threads());
  displayln_info(ssprintf("FftComplexFieldPlan::init: malloc %ld",
                          nc_size * sizec * sizeof(ComplexD)));
  ComplexD* fftdatac = (ComplexD*)fftw_malloc(nc_size * sizec * sizeof(ComplexD));
  const Int rank = 1;
  const Int n[1] = {sizec};
  const Long howmany = nc_size;
  const Long dist = 1;
  const Long stride = nc_size;
  unsigned fftw_plan_flag = FFTW_ESTIMATE;
  if (get_fftw_plan_flag() == "measure") {
    fftw_plan_flag = FFTW_MEASURE;
  } else {
    qassert(get_fftw_plan_flag() == "estimate");
    fftw_plan_flag = FFTW_ESTIMATE;
  }
  fftplan = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)fftdatac, n,
                               stride, dist, (fftw_complex*)fftdatac, n, stride,
                               dist, is_forward ? FFTW_FORWARD : FFTW_BACKWARD,
                               fftw_plan_flag);
  fftw_free(fftdatac);
  displayln_info(ssprintf("FftComplexFieldPlan::init: free %ld",
                          nc_size * sizec * sizeof(ComplexD)));
  sp = make_shuffle_plan_fft(geo.total_site(), dir);
}

FftComplexFieldPlan& get_fft_plan(const Geometry& geo, const Int mc,
                                  const Int dir, const bool is_forward)
{
  TIMER("get_fft_plan");
  qassert(check_fft_plan_key(geo, mc, dir, is_forward));
  Cache<std::string, FftComplexFieldPlan>& cache = get_fft_plan_cache();
  const std::string key =
      ssprintf("%s %s %d %d %d %d", show(geo.node_site).c_str(),
               show(geo.geon.size_node).c_str(), geo.geon.id_node, mc, dir,
               is_forward ? 1 : 0);
  if (cache.has(key)) {
    return cache[key];
  }
  FftComplexFieldPlan& plan = cache[key];
  plan.init(geo, mc, dir, is_forward);
  return plan;
}

}  // namespace qlat
