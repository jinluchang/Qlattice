#pragma once

#include <qlat/contract-wall-src-prop.h>

namespace qlat
{  //

inline LatData mk_hvp_table(const Coordinate& total_site)
{
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_number("mu", 0, 2));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  return ld;
}

inline LatData contract_chvp3(const SelProp& prop1, const SelProp& prop2,
                              const Int tslice_src, const FieldSelection& fsel)
{
  TIMER_VERBOSE("contract_chvp3");
  const Geometry& geo = prop1.geo();
  const Coordinate total_site = geo.total_site();
  const array<SpinMatrix, 4>& v_ms = SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  qassert(fsel.n_elems == prop1.n_elems);
  qassert(fsel.n_elems == prop2.n_elems);
  const Int multiplicity = 3;
  SelectedField<ComplexD> chvp;
  chvp.init(fsel, multiplicity);
  qacc_for(idx, fsel.n_elems, {
    const WilsonMatrix& wm1_x_y = prop1.get_elem(idx);
    const WilsonMatrix& wm2_x_y = prop2.get_elem(idx);
    const WilsonMatrix wm2_y_x =
        gamma5 * (WilsonMatrix)matrix_adjoint(wm2_x_y) * gamma5;
    Vector<ComplexD> chvp_v = chvp.get_elems(idx);
    for (Int mu = 0; mu < 3; ++mu) {
      const WilsonMatrix wm = wm2_y_x * v_ms[mu] * wm1_x_y;
      chvp_v[mu] = matrix_trace(wm, v_ms[mu]);
    }
  });
  LatData ld = mk_hvp_table(total_site);
  for (Long idx = 0; idx < fsel.n_elems; ++idx) {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Int tsep = mod(xg[3] - tslice_src, total_site[3]);
    const Vector<ComplexD> chvp_v = chvp.get_elems_const(idx);
    Vector<ComplexD> ldv = lat_data_cget(ld, make_array<int>(tsep));
    ldv += chvp_v;
  }
  glb_sum(ld);
  ld *= 1.0 / get_fsel_prob(fsel);
  return ld;
}

}  // namespace qlat
