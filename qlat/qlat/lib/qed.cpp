#include <qlat/core.h>
#include <qlat/qed.h>
#include "qlat-utils/matrix.h"
#include "qlat-utils/qacc-func.h"
#include "qlat-utils/types.h"

namespace qlat
{  //

void set_left_expanded_gauge_field(Field<ComplexD>& gf1,
                                   const Field<ComplexD>& gf)
{
  TIMER("set_left_expanded_gauge_field");
  const Coordinate expansion_left(1, 1, 1, 1);
  const Coordinate expansion_right(0, 0, 0, 0);
  const Geometry geo1 = geo_resize(gf.geo(), expansion_left, expansion_right);
  gf1.init(geo1, 4);
  Qassert(gf1.geo() == geo1);
  gf1 = gf;
  refresh_expanded_1(gf1);
}

void multiply_m_dwf_qed(Field<ComplexD>& out, const Field<ComplexD>& in,
                        const Field<ComplexD>& gf1, const RealD mass,
                        const RealD m5, const Int ls)
// set_left_expanded_gauge_field(gf1, gf);
// in.geo() should not be expanded.
// out can be the same object as in.
// mass is the masss of the fermion, for example, its value can be `0.1`.
// m5 should typically to 1.0
// ls should typically to a large even integer, such as 64.
{
  TIMER("multiply_m_dwf_qed");
  Qassert(mass > 0.0);
  Qassert(ls > 0);
  const Geometry& geo = in.geo();
  const Geometry& geo_gf1 = gf1.geo();
  Qassert(is_matching_geo(geo, geo_gf1));
  Qassert(geo_gf1.expansion_left == Coordinate(1, 1, 1, 1));
  Qassert(geo_gf1.expansion_right == Coordinate(0, 0, 0, 0));
  Qassert(gf1.multiplicity == 4);
  Qassert(in.multiplicity == 4 * ls);
  Qassert(geo.is_only_local == true);
  const Coordinate expansion_left(1, 1, 1, 1);
  const Coordinate expansion_right(1, 1, 1, 1);
  const Geometry geo1 = geo_resize(geo, expansion_left, expansion_right);
  Qassert(geo.eo == 0);
  Qassert(geo1.eo == 0);
  Qassert(geo_gf1.eo == 0);
  const box<SpinMatrixConstants>& smc = get_spin_matrix_constants();
  Field<ComplexD> in1;
  in1.init(geo1, 4 * ls);
  in1 = in;
  refresh_expanded_1(in1);
  out.init(geo, 4 * ls);
  Qassert(out.geo().is_only_local);
  set_zero(out);
  out = in1;
  qacc_for(index, geo.local_volume(), {
    const SpinMatrix& gamma_5 = smc().gamma5;
    const array<SpinMatrix, 4>& gammas = smc().cps_gammas;
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ComplexD> v = out.get_elems(index);
    {
      const Vector<ComplexD> iv = in1.get_elems_const(xl);
      mul_vec_plus(v, (ComplexD)(5.0 - m5), iv);
      const Vector<ComplexD> iv_m(iv.p + 4, 4 * (ls - 1));
      const Vector<ComplexD> iv_m0(iv.p, 4);
      const Vector<ComplexD> iv_p(iv.p, 4 * (ls - 1));
      const Vector<ComplexD> iv_p0(iv.p + 4 * (ls - 1), 4);
      Vector<ComplexD> v_m(v.p, 4 * (ls - 1));
      Vector<ComplexD> v_m0(v.p + 4 * (ls - 1), 4);
      Vector<ComplexD> v_p(v.p + 4, 4 * (ls - 1));
      Vector<ComplexD> v_p0(v.p, 4);
      mul_vec_plus(v_m, (ComplexD)(-0.5), iv_m);
      mul_vec_plus(v_m0, (ComplexD)(0.5 * mass), iv_m0);
      mat_mul_multi_vec_plus(v_m, (ComplexD)(0.5), gamma_5, iv_m);
      mat_mul_multi_vec_plus(v_m0, (ComplexD)(-0.5 * mass), gamma_5, iv_m0);
      mul_vec_plus(v_p, (ComplexD)(-0.5), iv_p);
      mul_vec_plus(v_p0, (ComplexD)(0.5 * mass), iv_p0);
      mat_mul_multi_vec_plus(v_p, (ComplexD)(-0.5), gamma_5, iv_p);
      mat_mul_multi_vec_plus(v_p0, (ComplexD)(0.5 * mass), gamma_5, iv_p0);
    }
    for (Int mu = 0; mu < 4; ++mu) {
      const Coordinate xl_p = coordinate_shifts(xl, mu);
      const Coordinate xl_m = coordinate_shifts(xl, -mu - 1);
      const ComplexD u_p = gf1.get_elem(xl, mu);
      const ComplexD u_m = 1.0 / gf1.get_elem(xl_m, mu);
      const Vector<ComplexD> iv_p = in1.get_elems_const(xl_p);
      const Vector<ComplexD> iv_m = in1.get_elems_const(xl_m);
      mat_mul_multi_vec_plus(v, (ComplexD)(0.5 * u_p), gammas[mu], iv_p);
      mat_mul_multi_vec_plus(v, (ComplexD)(-0.5 * u_m), gammas[mu], iv_m);
      mul_vec_plus(v, (ComplexD)(-0.5 * u_p), iv_p);
      mul_vec_plus(v, (ComplexD)(-0.5 * u_m), iv_m);
    }
  });
}

}  // namespace qlat
