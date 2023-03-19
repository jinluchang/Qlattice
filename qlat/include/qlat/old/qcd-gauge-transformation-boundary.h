#pragma once

#include <qlat/qcd-gauge-transformation.h>

namespace qlat
{  //

struct FourInterval : std::pair<Coordinate, Coordinate> {
};

inline void make_local_deflation_plan(
    std::vector<U1GaugeTransform>& u1gts,
    std::vector<FourInterval>& global_partition, const Geometry& geo,
    const Coordinate& tw_par)
// tw_par indicates the partition number in each of the directions.
// Should be large or equal than 1. 1 means boundary condition is not changed in
// this direction. u1gts constains the U1 gauge transformation(s) that move the
// boundaries to approriate places.
{
  // total number of partition of the global lattice
  int Np = product(tw_par);
  Printf("Number of partitions = %d\n", Np);

  Coordinate global_size = geo.global_size();
  assert(global_size % tw_par == Coordinate(0, 0, 0, 0));

  assert(geo.eo == 0);            // Only work with odd sites, for now
  assert(geo.multiplicity == 1);  // Only work with odd sites, for now

  Coordinate partition_size = global_size / tw_par;

  // initialize the U1 field
  u1gts.resize(Np);
  for (int i = 0; i < Np; i++) {
    u1gts[i].init(geo);
    for (long j = 0; j < u1gts[i].field.size(); j++) {
      u1gts[i].field[j] = +1.;
    }
  }

  // initialize the global_partition.
  global_partition.resize(Np);

  for (int i = 0; i < Np; i++) {
    Printf("partition #%04d:\n", i);
    Coordinate partition_coor = qlat::coordinate_from_index(i, tw_par);
    for (int mu = 0; mu < 4; mu++) {
      global_partition[i].first[mu] = partition_coor[mu] * partition_size[mu];
      global_partition[i].second[mu] =
          (partition_coor[mu] + 1) *
          partition_size[mu];  // left(first) including, right(second)
                               // excluding. [first, second)

      // Find the border that is fartest from the source(partition)
      int target_border = (global_partition[i].second[mu] +
                           (global_size[mu] - partition_size[mu]) / 2) %
                          global_size[mu];

      Printf("direction = %03d, first = %03d, second = %03d, target = %03d\n",
             mu, global_partition[i].first[mu], global_partition[i].second[mu],
             target_border);

      // Flip all gt link in [0, target_border) to -1
      long num_flip = 0;
      // #pragma omp parallel for
      for (long index = 0; index < geo.local_volume(); index++) {
        Coordinate local_coor = geo.coordinate_from_index(index);
        Coordinate global_coor = geo.coordinate_g_from_l(local_coor);
        if ((target_border < global_size[mu] / 2 && global_coor[mu] >= 0 &&
             global_coor[mu] < target_border) ||
            (target_border >= global_size[mu] / 2 && global_coor[mu] >= 0 &&
             global_coor[mu] >= target_border)) {
          u1gts[i].field[index * 1] *= -1.;
          num_flip++;
          // printf("(%02d,%02d,%02d,%02d): YES.\n", global_coor[0],
          // global_coor[1], global_coor[2], global_coor[3]);
        } else {
          // if(i == 31 and mu == 3) printf("(%02d,%02d,%02d,%02d)[%04d]: NO.
          // \n", global_coor[0], global_coor[1], global_coor[2],
          // global_coor[3], count);
        }
      }
      glb_sum(num_flip);
      Printf("flipped %08d times\n", num_flip);
    }
    double sum_real = 0.;
    double sum_imag = 0.;
    // TODO: Test!!!
    for (long index = 0; index < geo.local_volume();
         index++) {  // We are only working with vectors so it seems we don't
                     // need to worry about communication?
      sum_real += u1gts[i].field[index].real();
      sum_imag += u1gts[i].field[index].imag();
    }
    glb_sum(sum_real);
    glb_sum(sum_imag);
    Printf("sum_real = %.8E\n", sum_real);
    Printf("sum_imag = %.8E\n", sum_imag);
  }
}

inline void apply_u1_gauge_tranform_on_bfm_vct(void* bfm_vct, size_t Ls,
                                               const U1GaugeTransform& u1gt)
{
  // Assuming single precision.
  // Since this is only a U1 tranformation we treat spin and color equally.
  // Assuming bfm_vct is using even-odd preconditioning(checkerboarding).

  assert(u1gt.geo().eo == 0);  // NO checkerboarding for qlat
  assert(u1gt.geo().multiplicity == 1);

  size_t bfm_vct_block_size =
      Ls * 12;  // 12 = 3 * 4. bfm_vct_block_size is the number single precision
                // complex number on one 4d site.
  ComplexF* bfm_vct_ComplexF = (ComplexF*)bfm_vct;
  size_t site_size_4d = u1gt.geo().local_volume();
  // Printf("[Apply U1 bfm]: Ls = %d, bfm_vct_block_size = %d, site_size_4d =
  // %d\n", Ls, bfm_vct_block_size, site_size_4d);
#pragma omp for
  for (size_t m = 0; m < site_size_4d / 2; m++) {
    Coordinate local_coor = u1gt.geo().coordinate_from_index(m);
    if (sum(local_coor) % 2 == 0)
      continue;  // TODO: temporary fix. Fix me!!! We don't want even sites.

    size_t m1 = m;
    size_t m2 = m + site_size_4d / 2;
    size_t b1, b2;
    for (size_t s = 0; s < bfm_vct_block_size; s++) {
      b1 = (m / 2 * bfm_vct_block_size + s) * 2;
      b2 = (m / 2 * bfm_vct_block_size + s) * 2 + 1;
      bfm_vct_ComplexF[b1] *= u1gt.field[m1];
      bfm_vct_ComplexF[b2] *= u1gt.field[m2];
      // Printf("(bfm_idx, cps_idx) = (%06d, %06d)\n", b1*bfm_vct_block_size+s,
      // m1); Printf("(bfm_idx, cps_idx) = (%06d, %06d)\n",
      // b2*bfm_vct_block_size+s, m2);
    }
  }
}

inline bool is_inside(const Coordinate& coor, const Coordinate& lower,
                      const Coordinate& upper)
{
  return (lower[0] <= coor[0] and coor[0] < upper[0]) and
         (lower[1] <= coor[1] and coor[1] < upper[1]) and
         (lower[2] <= coor[2] and coor[2] < upper[2]) and
         (lower[3] <= coor[3] and coor[3] < upper[3]);
}

inline void extract_par_vct_from_bfm_vct(void* par_vct, const void* bfm_vct,
                                         size_t Ls, const FourInterval& par,
                                         const Geometry& geo)
{
  // Assuming single precision.
  // Assuming bfm_vct is using even-odd preconditioning(checkerboarding).

  assert(geo.eo == 0);  // NO checkerboarding for qlat
  assert(geo.multiplicity == 1);

  ComplexF* parp = (ComplexF*)par_vct;
  ComplexF* bfmp = (ComplexF*)bfm_vct;

  size_t bfm_vct_block_size =
      Ls * 12;  // 12 = 3 * 4. bfm_vct_block_size is the number of single
                // precision complex number on one 4d site.
  size_t site_size_4d = geo.local_volume();
// #pragma omp parallel for
#pragma omp for
  for (size_t m = 0; m < site_size_4d / 2; m++) {
    size_t m1 = m;
    size_t b1;
    Coordinate local_coor1 = geo.coordinate_from_index(m1);
    Coordinate global_coor1 = geo.coordinate_g_from_l(local_coor1);

    if (sum(local_coor1) % 2 == 0)
      continue;  // TODO: temporary fix. Fix me!!! We don't want even sites.

    if (is_inside(global_coor1, par.first, par.second)) {
      for (size_t s = 0; s < bfm_vct_block_size; s++) {
        b1 = (m / 2 * bfm_vct_block_size + s) * 2;
        parp[b1] = bfmp[b1];
      }
    } else {
      for (size_t s = 0; s < bfm_vct_block_size; s++) {
        b1 = (m / 2 * bfm_vct_block_size + s) * 2;
        parp[b1] = 0.;
      }
    }

    size_t m2 = m + site_size_4d / 2;
    size_t b2;
    Coordinate local_coor2 = geo.coordinate_from_index(m2);
    Coordinate global_coor2 = geo.coordinate_g_from_l(local_coor2);
    if (is_inside(global_coor2, par.first, par.second)) {
      for (size_t s = 0; s < bfm_vct_block_size; s++) {
        b2 = (m / 2 * bfm_vct_block_size + s) * 2 + 1;
        parp[b2] = bfmp[b2];
      }
    } else {
      for (size_t s = 0; s < bfm_vct_block_size; s++) {
        b2 = (m / 2 * bfm_vct_block_size + s) * 2 + 1;
        parp[b2] = 0.;
      }
    }
  }
}

inline void scalar_multiplication_by_partition(void* out_vct,
                                               const void* bfm_vct,
                                               const std::vector<Complex>& b,
                                               int Ls, const Coordinate& tw_par,
                                               const Geometry& geo)
{
  // Assuming single precision.
  // Assuming bfm_vct is using even-odd preconditioning(checkerboarding).

  Coordinate global_size = geo.global_size();
  Coordinate partition_size = global_size / tw_par;

  ComplexF* bfmp = (ComplexF*)bfm_vct;
  ComplexF* outp = (ComplexF*)out_vct;

  size_t bfm_vct_block_size =
      Ls * 12;  // 12 = 3 * 4. bfm_vct_block_size is the number of single
                // precision complex number on one 4d site.
  size_t site_size_4d = geo.local_volume();
// #pragma omp parallel for
#pragma omp for
  for (size_t m = 0; m < site_size_4d / 2; m++) {
    size_t m1 = m;
    size_t b1;
    Coordinate local_coor1 = geo.coordinate_from_index(m1);
    if (sum(local_coor1) % 2 == 0)
      continue;  // TODO: temporary fix. Fix me!!! We don't want even sites.

    Coordinate global_coor1 = geo.coordinate_g_from_l(local_coor1);
    int p = qlat::index_from_coordinate(global_coor1 / partition_size, tw_par);

    for (size_t s = 0; s < bfm_vct_block_size; s++) {
      b1 = (m / 2 * bfm_vct_block_size + s) * 2;
      outp[b1] = bfmp[b1] * (ComplexF)b[p];
    }

    // size_t m2 = m + site_size_4d / 2;
    size_t b2;
    // Coordinate local_coor2 = geo.coordinate_from_index(m2);
    // Coordinate global_coor2 = geo.coordinate_g_from_l(local_coor2);
    p = qlat::index_from_coordinate(global_coor1 / partition_size, tw_par);
    for (size_t s = 0; s < bfm_vct_block_size; s++) {
      b2 = (m / 2 * bfm_vct_block_size + s) * 2 + 1;
      outp[b2] = bfmp[b2] * (ComplexF)b[p];
    }
  }
}

inline void fft_convolution(std::vector<Complex>& out,
                            const std::vector<Complex>& x,
                            const std::vector<Complex>& y)
{
  // global calculation: same on all nodes.
  // single thread!!!

  TIMER("fft_convolution()");

  assert(x.size() == y.size());
  out.resize(x.size());

  // static fftw_complex* x_in;
  static fftw_complex* x_out;
  // static fftw_complex* y_in;
  static fftw_complex* y_out;

  static fftw_complex* z_in;
  static fftw_complex* z_out;

  static fftw_plan p_forward;
  static fftw_plan p_backward;

  static fftw_complex* f_in;
  static fftw_complex* f_out;

  static bool initialized = false;
  static int N;

  if (not initialized) {
    N = x.size();

    // x_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    x_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    // y_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    y_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

    z_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    z_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

    f_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    f_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

    p_forward = fftw_plan_dft_1d(N, f_in, f_out, FFTW_FORWARD, FFTW_MEASURE);
    p_backward = fftw_plan_dft_1d(N, z_in, z_out, FFTW_BACKWARD, FFTW_MEASURE);

    initialized = true;
  }

  assert(N == (int)x.size());

  std::memcpy(f_in, x.data(), sizeof(fftw_complex) * N);
  fftw_execute(p_forward);
  std::memcpy(x_out, f_out, sizeof(fftw_complex) * N);

  std::memcpy(f_in, y.data(), sizeof(fftw_complex) * N);
  fftw_execute(p_forward);
  std::memcpy(y_out, f_out, sizeof(fftw_complex) * N);

  for (int i = 0; i < N; i++) {
    z_in[i][0] = x_out[i][0] * y_out[i][0] - x_out[i][1] * y_out[i][1];
    z_in[i][1] = x_out[i][1] * y_out[i][0] + x_out[i][0] * y_out[i][1];
    // Printf("x[%d] = %.8E + i %.8E, y[%d] = %.8E + i %.8E\n", i, x[i].real(),
    // x[i].imag(), i, y[i].real(), y[i].imag());
  }

  fftw_execute(p_backward);

  std::memcpy((void*)out.data(), (void*)z_out, sizeof(fftw_complex) * N);

  for (int i = 0; i < N; i++) {
    out[i] /= (double)N;
  }

  return;
}

}  // namespace qlat
