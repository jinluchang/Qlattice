#pragma once

#include <fftw3.h>
#include <qlat/core.h>
#include <qlat/field-shuffle.h>
#include <qlat/utils-coordinate.h>
#include <qlat/vector_utils/utils_FFT_GPU.h>

#include <vector>

namespace qlat
{  //

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
  ~FftComplexFieldPlan() { init(); }
  //
  void init();
  void init(const Geometry& geo_, const int mc_, const int dir_,
            const bool is_forward_);
};

API inline Cache<std::string, FftComplexFieldPlan>& get_fft_plan_cache()
{
  static Cache<std::string, FftComplexFieldPlan> cache(
      "FftComplexFieldPlanCache", 32);
  return cache;
}

FftComplexFieldPlan& get_fft_plan(const Geometry& geo, const int mc,
                                  const int dir, const bool is_forward);

template <class M>
void fft_complex_field_dir(Field<M>& field1, const Field<M>& field,
                           const int dir, const bool is_forward)
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

template <class M>
void fft_complex_fields(std::vector<Handle<Field<M> > >& vec,
                        const std::vector<int>& fft_dirs,
                        const std::vector<bool>& fft_is_forwards,
                        int mode_fft = 1)
{
  const long n_field = vec.size();
  int use_plan = 0;
  bool fft_direction = false;
  bool ft4D = false;
  if (mode_fft == 1) {
    ////check 3D
    if (fft_dirs.size() == 3 and fft_dirs[0] == 0 and fft_dirs[1] == 1 and
        fft_dirs[2] == 2) {
      if (fft_is_forwards[0] == fft_is_forwards[1] and
          fft_is_forwards[0] == fft_is_forwards[2]) {
        use_plan = 1;
        fft_direction = fft_is_forwards[0];
        ft4D = false;
      }
    }
    ////check 4D
    if (fft_dirs.size() == 4 and fft_dirs[0] == 0 and fft_dirs[1] == 1 and
        fft_dirs[2] == 2 and fft_dirs[3] == 3) {
      if (fft_is_forwards[0] == fft_is_forwards[1] and
          fft_is_forwards[0] == fft_is_forwards[2])
        if (fft_is_forwards[0] == fft_is_forwards[3]) {
          use_plan = 1;
          fft_direction = fft_is_forwards[0];
          ft4D = true;
        }
    }
  }
  if (use_plan == 0) {
    for (long i = 0; i < n_field; ++i) {
      Field<M> ft;
      for (long k = 0; k < (long)fft_dirs.size(); ++k) {
        ft = vec[i]();
        const int fft_dir = fft_dirs[k];
        const bool is_forward = fft_is_forwards[k];
        fft_complex_field_dir(vec[i](), ft, fft_dir, is_forward);
      }
    }
  } else if (use_plan == 1) {
    fft_fieldM(vec, fft_direction, ft4D);
  } else {
    qassert(false);
  }
}

// --------------------

#ifdef QLAT_INSTANTIATE_FIELD_FFT
#define QLAT_EXTERN
#else
#define QLAT_EXTERN extern
#endif

#define QLAT_EXTERN_TEMPLATE(TYPENAME)                                       \
                                                                             \
  QLAT_EXTERN template void fft_complex_field_dir<TYPENAME>(                 \
      Field<TYPENAME> & field1, const Field<TYPENAME>& field, const int dir, \
      const bool is_forward);                                                \
                                                                             \
  QLAT_EXTERN template void fft_complex_field_dir<TYPENAME>(                 \
      Field<TYPENAME> & field, const int dir, const bool is_forward);        \
                                                                             \
  QLAT_EXTERN template void fft_complex_field<TYPENAME>(                     \
      Field<TYPENAME> & field, const bool is_forward);                       \
                                                                             \
  QLAT_EXTERN template void fft_complex_field_spatial<TYPENAME>(             \
      Field<TYPENAME> & field, const bool is_forward);                       \
                                                                             \
  QLAT_EXTERN template void fft_complex_fields<TYPENAME>(                    \
      std::vector<Handle<Field<TYPENAME> > > & vec,                          \
      const std::vector<int>& fft_dirs,                                      \
      const std::vector<bool>& fft_is_forwards, int mode_fft)

QLAT_CALL_WITH_TYPES(QLAT_EXTERN_TEMPLATE);

#undef QLAT_EXTERN
#undef QLAT_EXTERN_TEMPLATE

}  // namespace qlat
