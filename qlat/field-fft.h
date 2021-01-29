#pragma once

#include <qlat/config.h>
#include <qlat/field.h>
#include <qlat/geometry.h>
#include <qlat/mpi.h>
#include <qlat/utils-coordinate.h>
#include <qlat/utils.h>
#include <qlat/field-shuffle.h>

#include <fftw3.h>

#include <vector>

namespace qlat
{  //

struct fft_complex_field_plan {
  Geometry geo;     // geo.is_only_local == true
  int mc;           // geo.multiplicity * sizeof(M) / sizeof(Complex)
  int dir;          // direction of the fft
  bool is_forward;  // is forward fft (forward \sum_x f[x] exp(-i k x))
  //
  fftw_plan fftplan;
  ShufflePlan sp;
  //
  bool is_match(const Geometry& geo_, const int mc_, const int dir_,
                const bool is_forward_)
  {
    return geo_ == geo && mc_ == mc && dir_ == dir && is_forward_ == is_forward;
  }
  //
  static bool check(const Geometry& geo_, const int mc_, const int dir_,
                    const bool is_forward_)
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
    displayln_info(fname +
                   ssprintf(": start to make a new fft plan with id = %d",
                            next_plan_index));
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
      displayln_info("fft_complex_field_plan::end(): free a plan.");
      fftw_destroy_plan(fftplan);
      geo.initialized = false;
    }
  }
  //
  void init(const Geometry& geo_, const int mc_, const int dir_,
            const bool is_forward_)
  {
    TIMER_VERBOSE("fft_complex_field_plan::init");
    qassert(check(geo_, mc_, dir_, is_forward_));
    geo = geo_;
    mc = mc_;
    dir = dir_;
    is_forward = is_forward_;
    const int sizec = geo.total_site()[dir];
    const int nc = geo.local_volume() / geo.node_site[dir] * mc;
    const int chunk = ((nc / mc - 1) / geo.geon.size_node[dir] + 1) * mc;
    const int nc_start = std::min(nc, geo.geon.coor_node[dir] * chunk);
    const int nc_stop = std::min(nc, nc_start + chunk);
    const int nc_size = nc_stop - nc_start;
    // fftw_init_threads();
    // fftw_plan_with_nthreads(omp_get_max_threads());
    displayln_info(ssprintf("fft_complex_field_plan::init: malloc %d",
                            nc_size * sizec * sizeof(Complex)));
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
    displayln_info(ssprintf("fft_complex_field_plan::init: free %d",
                            nc_size * sizec * sizeof(Complex)));
    sp = make_shuffle_plan_fft(geo.total_site(), dir);
  }
};

template <class M>
void fft_complex_field_dir(Field<M>& field, const int dir,
                           const bool is_forward)
{
  TIMER("fft_complex_field_dir");
  Geometry geo = field.geo();
  geo.resize(0);
  const int mc = geo.multiplicity * sizeof(M) / sizeof(Complex);
  fft_complex_field_plan& plan =
      fft_complex_field_plan::get_plan(geo, mc, dir, is_forward);
  fftw_plan& fftplan = plan.fftplan;
  const int sizec = geo.total_site()[dir];
  const int nc = geo.local_volume() / geo.node_site[dir] * mc;
  const int chunk = ((nc / mc - 1) / geo.geon.size_node[dir] + 1) * mc;
  const int nc_start = std::min(nc, geo.geon.coor_node[dir] * chunk);
  const int nc_stop = std::min(nc, nc_start + chunk);
  const int nc_size = nc_stop - nc_start;
  qassert(nc_size >= 0);
  Complex* fftdatac = (Complex*)fftw_malloc(nc_size * sizec * sizeof(Complex));
  const ShufflePlan& sp = plan.sp;
  std::vector<Field<M> > fft_fields;
  shuffle_field(fft_fields, field, sp);
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
  shuffle_field_back(field, fft_fields, sp);
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
    fft_complex_field_dir(field, dir, is_forward);
  }
}

}  // namespace qlat
