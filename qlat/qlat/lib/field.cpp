#define QLAT_INSTANTIATE_FIELD

#include <qlat/field.h>

namespace qlat
{  //

void set_sqrt_field(Field<RealD>& f, const Field<RealD>& f1)
{
  TIMER("set_sqrt_field(f,f1)");
  const Geometry geo = geo_resize(f1.geo());
  const Int multiplicity = f1.multiplicity;
  f.init(geo, multiplicity);
  qacc_for(index, geo.local_volume(), {
    const Vector<RealD> f1v = f1.get_elems_const(index);
    Vector<RealD> fv = f.get_elems(index);
    for (Int m = 0; m < multiplicity; ++m) {
      fv[m] = std::sqrt(f1v[m]);
    }
  });
}

void set_mom_phase_field(FieldM<ComplexD, 1>& f, const CoordinateD& mom)
// mom is in lattice unit (1/a)
// exp(i * mom \cdot xg )
{
  TIMER("set_mom_phase_field");
  const Geometry& geo = f.geo();
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = f.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    double phase = 0;
    for (Int k = 0; k < DIMN; ++k) {
      phase += mom[k] * xg[k];
    }
    f.get_elem(xl) = qpolar(1.0, phase);
  });
}

void set_phase_field(FieldM<ComplexD, 1>& f, const CoordinateD& lmom)
// lmom is in lattice momentum unit
// exp(i * 2*pi/L * lmom \cdot xg )
{
  TIMER("set_phase_field");
  const CoordinateD mom = lmom * lattice_mom_mult(f.geo());
  set_mom_phase_field(f, mom);
}

void set_xg_field(Field<Int>& f, const Geometry& geo_)
{
  TIMER("set_xg_field(f,geo)");
  const Geometry geo = geo_resize(geo_);
  f.init(geo, DIMN);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    Vector<Int> fv = f.get_elems(xl);
    for (Int mu = 0; mu < DIMN; ++mu) {
      fv[mu] = xg[mu];
    }
  });
}

}  // namespace qlat
