#pragma once

#include <qlat-utils/matrix.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>

namespace qlat
{  //

void gf_ape_smear(GaugeField& gf, const GaugeField& gf0, const double alpha,
                  const Long steps = 1);

void gf_spatial_ape_smear(GaugeField& gf, const GaugeField& gf0,
                          const double alpha, const Long steps = 1);

void gf_hyp_smear(GaugeField& gf, const GaugeField& gf0, const double alpha1,
                  const double alpha2, const double alpha3);

template <class T>
void prop_smear(Propagator4dT<T>& prop, const GaugeFieldT<T>& gf1,
                const double coef, const int step,
                const CoordinateD& mom = CoordinateD(),
                const bool smear_in_time_dir = false)
// gf1 is left_expanded and refreshed
// set_left_expanded_gauge_field(gf1, gf)
// prop is of normal size
{
  TIMER_FLOPS("prop_smear");
  const int n_avg = smear_in_time_dir ? 8 : 6;
  const Long v_gb = prop.geo().local_volume() * 12 * 4;
  timer.flops += v_gb * step * n_avg * (3 * (3 * 6 + 2 * 2));
  if (0 == step) {
    return;
  }
  const Geometry& geo = prop.geo();
  const Geometry geo1 =
      smear_in_time_dir
          ? geo_resize(geo, 1)
          : geo_resize(geo, Coordinate(1, 1, 1, 0), Coordinate(1, 1, 1, 0));
  const int dir_limit = smear_in_time_dir ? 4 : 3;
  array<ComplexD, 8> mom_factors_v;
  box<array<ComplexD, 8>> mom_factors(
      mom_factors_v);  // (array<ComplexD, 8>());
  for (int i = 0; i < 8; ++i) {
    const int dir = i - 4;
    const double phase = dir >= 0 ? mom[dir] : -mom[-dir - 1];
    mom_factors()[i] = qpolar(coef / n_avg, -phase);
  }
  Propagator4dT<T> prop1;
  prop1.init(geo1);
  for (int i = 0; i < step; ++i) {
    prop1 = prop;
    refresh_expanded_1(prop1);
    qacc_for(index, geo.local_volume(), {
      const Coordinate xl = prop.geo().coordinate_from_index(index);
      WilsonMatrixT<T>& wm = prop.get_elem(xl);
      wm *= 1 - coef;
      for (int dir = -dir_limit; dir < dir_limit; ++dir) {
        const Coordinate xl1 = coordinate_shifts(xl, dir);
        ColorMatrixT<T> link =
            dir >= 0
                ? gf1.get_elem(xl, dir)
                : (ColorMatrixT<T>)matrix_adjoint(gf1.get_elem(xl1, -dir - 1));
        link *= mom_factors()[dir + 4];
        wm += link * prop1.get_elem(xl1);
      }
    });
  }
}

void prop_spatial_smear_no_comm(std::vector<FermionField4d>& ff_vec,
                                const GaugeField& gf, const RealD coef,
                                const Long step,
                                const CoordinateD& mom = CoordinateD());

#ifdef QLAT_INSTANTIATE_SMEAR
#define QLAT_EXTERN
#else
#define QLAT_EXTERN extern
#endif

QLAT_EXTERN template void prop_smear<Real>(Propagator4d&, const GaugeField&,
                                           const double, const int,
                                           const CoordinateD&, const bool);

#undef QLAT_EXTERN

}  // namespace qlat
