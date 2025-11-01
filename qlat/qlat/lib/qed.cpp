#include <qlat/core.h>
#include <qlat/qed.h>

namespace qlat
{  //

void set_left_expanded_gauge_field(Field<ComplexD>& gf1, const Field<ComplexD>& gf)
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
                        const Int ls)
// set_left_expanded_gauge_field(gf1, gf);
// in.geo() should not be expanded.
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
  Field<ComplexD> in1;
  in1.init(geo1, 4 * ls);
  in1 = in;
  refresh_expanded_1(in1);
  out.init(geo, 4 * ls);
  set_zero(out);
  // Placeholder implementation: copy input to output.
  out = in1;
}

}  // namespace qlat
