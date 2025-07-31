#define QLAT_INSTANTIATE_SMEAR

#include <qlat/qcd-smear.h>
#include <qlat/vector_utils/utils_smear_vecs.h>

namespace qlat
{  //

ColorMatrix gf_link_ape_smear_no_comm(const GaugeField& gf,
                                      const Coordinate& xl, const int mu,
                                      const double alpha)
{
  return color_matrix_su_projection(
      (ComplexD)(1.0 - alpha) * gf.get_elem(xl, mu) +
      (ComplexD)(alpha / 6.0) * gf_staple_no_comm(gf, xl, mu));
}

void gf_ape_smear_no_comm(GaugeField& gf, const GaugeField& gf0,
                          const double alpha)
{
  TIMER_VERBOSE("gf_ape_smear_no_comm");
  qassert(&gf != &gf0);
  const Geometry& geo = gf0.geo();
  gf.init(geo_resize(geo));
  qassert(is_matching_geo(geo, gf.geo()));
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int mu = 0; mu < DIMN; ++mu) {
      v[mu] = gf_link_ape_smear_no_comm(gf0, xl, mu, alpha);
    }
  }
}

void gf_ape_smear(GaugeField& gf, const GaugeField& gf0, const double alpha,
                  const Long steps)
{
  TIMER_VERBOSE("gf_ape_smear");
  gf = gf0;
  GaugeField gf1;
  gf1.init(geo_resize(gf0.geo(), 1));
  for (Long i = 0; i < steps; ++i) {
    gf1 = gf;
    refresh_expanded(gf1);
    gf_ape_smear_no_comm(gf, gf1, alpha);
  }
}

ColorMatrix gf_link_spatial_ape_smear_no_comm(const GaugeField& gf,
                                              const Coordinate& xl,
                                              const int mu, const double alpha)
{
  const double multi = mu == 3 ? 6.0 : 4.0;
  return color_matrix_su_projection(
      (ComplexD)(1.0 - alpha) * gf.get_elem(xl, mu) +
      (ComplexD)(alpha / multi) * gf_spatial_staple_no_comm(gf, xl, mu));
}

void gf_spatial_ape_smear_no_comm(GaugeField& gf, const GaugeField& gf0,
                                  const double alpha)
{
  TIMER_VERBOSE("gf_spatial_ape_smear_no_comm");
  qassert(&gf != &gf0);
  const Geometry& geo = gf0.geo();
  gf.init(geo_resize(geo));
  qassert(is_matching_geo(geo, gf.geo()));
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int mu = 0; mu < 3; ++mu) {
      // No need to smear the temperal link (mu == 3)
      v[mu] = gf_link_spatial_ape_smear_no_comm(gf0, xl, mu, alpha);
    }
  }
}

void gf_spatial_ape_smear(GaugeField& gf, const GaugeField& gf0,
                          const double alpha, const Long steps)
{
  TIMER_VERBOSE("gf_spatial_ape_smear");
  gf = gf0;
  const Coordinate expan_left(1, 1, 1, 0);
  const Coordinate expan_right(1, 1, 1, 0);
  GaugeField gf1;
  gf1.init(geo_resize(gf0.geo(), expan_left, expan_right));
  for (Long i = 0; i < steps; ++i) {
    gf1 = gf;
    refresh_expanded(gf1);
    gf_spatial_ape_smear_no_comm(gf, gf1, alpha);
  }
}

ColorMatrix gf_link_hyp_smear_3_no_comm(const GaugeField& gf,
                                        const Coordinate& xl, const int mu,
                                        const int nu, const int rho,
                                        const double alpha3)
{
  ColorMatrix ret;
  set_zero(ret);
  const Coordinate xl_mu = coordinate_shifts(xl, mu);
  for (int m = 0; m < DIMN; ++m) {
    if (mu != m && nu != m && rho != m) {
      ret += gf.get_elem(xl, m) * gf.get_elem(coordinate_shifts(xl, m), mu) *
             matrix_adjoint(gf.get_elem(xl_mu, m));
      ret += matrix_adjoint(gf.get_elem(coordinate_shifts(xl, -m - 1), m)) *
             gf.get_elem(coordinate_shifts(xl, -m - 1), mu) *
             gf.get_elem(coordinate_shifts(xl_mu, -m - 1), m);
    }
  }
  ret = (ComplexD)(1.0 - alpha3) * gf.get_elem(xl, mu) +
        (ComplexD)(alpha3 / 2.0) * ret;
  return color_matrix_su_projection(ret);
}

ColorMatrix gf_link_hyp_smear_2_no_comm(const GaugeField& gf,
                                        const Coordinate& xl, const int mu,
                                        const int nu, const double alpha2,
                                        const double alpha3)
{
  ColorMatrix ret;
  set_zero(ret);
  const Coordinate xl_mu = coordinate_shifts(xl, mu);
  for (int m = 0; m < DIMN; ++m) {
    if (mu != m && nu != m) {
      ret += gf_link_hyp_smear_3_no_comm(gf, xl, m, mu, nu, alpha3) *
             gf_link_hyp_smear_3_no_comm(gf, coordinate_shifts(xl, m), mu, m,
                                         nu, alpha3) *
             matrix_adjoint(
                 gf_link_hyp_smear_3_no_comm(gf, xl_mu, m, mu, nu, alpha3));
      ret += matrix_adjoint(gf_link_hyp_smear_3_no_comm(
                 gf, coordinate_shifts(xl, -m - 1), m, mu, nu, alpha3)) *
             gf_link_hyp_smear_3_no_comm(gf, coordinate_shifts(xl, -m - 1), mu,
                                         m, nu, alpha3) *
             gf_link_hyp_smear_3_no_comm(gf, coordinate_shifts(xl_mu, -m - 1),
                                         m, mu, nu, alpha3);
    }
  }
  ret = (ComplexD)(1.0 - alpha2) * gf.get_elem(xl, mu) +
        (ComplexD)(alpha2 / 4.0) * ret;
  return color_matrix_su_projection(ret);
}

ColorMatrix gf_link_hyp_smear_1_no_comm(const GaugeField& gf,
                                        const Coordinate& xl, const int mu,
                                        const double alpha1,
                                        const double alpha2,
                                        const double alpha3)
{
  ColorMatrix ret;
  set_zero(ret);
  const Coordinate xl_mu = coordinate_shifts(xl, mu);
  for (int m = 0; m < DIMN; ++m) {
    if (mu != m) {
      ret += gf_link_hyp_smear_2_no_comm(gf, xl, m, mu, alpha2, alpha3) *
             gf_link_hyp_smear_2_no_comm(gf, coordinate_shifts(xl, m), mu, m,
                                         alpha2, alpha3) *
             matrix_adjoint(
                 gf_link_hyp_smear_2_no_comm(gf, xl_mu, m, mu, alpha2, alpha3));
      ret += matrix_adjoint(gf_link_hyp_smear_2_no_comm(
                 gf, coordinate_shifts(xl, -m - 1), m, mu, alpha2, alpha3)) *
             gf_link_hyp_smear_2_no_comm(gf, coordinate_shifts(xl, -m - 1), mu,
                                         m, alpha2, alpha3) *
             gf_link_hyp_smear_2_no_comm(gf, coordinate_shifts(xl_mu, -m - 1),
                                         m, mu, alpha2, alpha3);
    }
  }
  ret = (ComplexD)(1.0 - alpha1) * gf.get_elem(xl, mu) +
        (ComplexD)(alpha1 / 6.0) * ret;
  return color_matrix_su_projection(ret);
}

ColorMatrix gf_link_hyp_smear_no_comm(const GaugeField& gf,
                                      const Coordinate& xl, const int mu,
                                      const double alpha1, const double alpha2,
                                      const double alpha3)
{
  return gf_link_hyp_smear_1_no_comm(gf, xl, mu, alpha1, alpha2, alpha3);
}

void gf_hyp_smear_no_comm(GaugeField& gf, const GaugeField& gf0,
                          const double alpha1, const double alpha2,
                          const double alpha3)
{
  TIMER_VERBOSE("gf_hyp_smear_no_comm");
  qassert(&gf != &gf0);
  const Geometry& geo = gf0.geo();
  gf.init(geo_resize(geo));
  qassert(is_matching_geo(geo, gf.geo()));
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int mu = 0; mu < DIMN; ++mu) {
      v[mu] = gf_link_hyp_smear_no_comm(gf0, xl, mu, alpha1, alpha2, alpha3);
    }
  }
}

void gf_hyp_smear(GaugeField& gf, const GaugeField& gf0, const double alpha1,
                  const double alpha2, const double alpha3)
// values in paper is 0.75 0.6 0.3
// 10.1103/PhysRevD.64.034504 Eq(4)
{
  TIMER_VERBOSE("gf_hyp_smear");
  GaugeField gf1;
  gf1.init(geo_resize(gf0.geo(), 1));
  gf1 = gf0;
  refresh_expanded(gf1);
  gf_hyp_smear_no_comm(gf, gf1, alpha1, alpha2, alpha3);
}

void prop_spatial_smear_no_comm(std::vector<FermionField4d>& ff_vec,
                                const GaugeField& gf, const RealD coef,
                                const Long step, const CoordinateD& mom)
// `gf` and each of `ff_vec` should contain entire time slices.
// No communication will be performed.
{
  TIMER_FLOPS("prop_spatial_smear_no_comm(ff_vec,gf,coef,step,mom)");
  const Geometry geo = gf.geo.get();
  const Int num_field = ff_vec.size();
  const Int n_avg = 6;
  const Long v_gb = geo.local_volume() * num_field * 4;
  timer.flops += v_gb * step * n_avg * (3 * (3 * 6 + 2 * 2));
  if (0 == step) {
    return;
  }
  const Int dir_limit = 3;
  array<ComplexD, 6> mom_factors_v;
  for (int i = 0; i < 6; ++i) {
    const Int dir = i - 3;
    const RealD phase = dir >= 0 ? mom[dir] : -mom[-dir - 1];
    mom_factors_v[i] = qpolar(coef / n_avg, -phase);
  }
  box<array<ComplexD, 6>> mom_factors(mom_factors_v,
                                      MemType::Acc);  // (array<ComplexD, 8>());
  const Int t_size = geo.total_site()[3];
  qassert(geo.geon.size_node == Coordinate(1, 1, 1, t_size));
  qassert(geo.is_only_local);
  qassert(num_field >= 0);
  if (num_field == 0) {
    return;
  }
  vector<FermionField4d> ffv_vec(num_field, MemType::Cpu);
  set_zero(ffv_vec);
  qfor(id_field, num_field, {
    qassert(ff_vec[id_field].geo.get() == geo);
    ff_vec[id_field].set_mem_type(MemType::Acc);
    ffv_vec[id_field].set_view(ff_vec[id_field]);
  });
  ffv_vec.set_mem_type(MemType::Acc);
  Field<ComplexD> ff, ff1;
  ff.set_mem_type(MemType::Acc);
  ff1.set_mem_type(MemType::Acc);
  ff.init(geo, num_field * 12);
  ff1.init(geo, num_field * 12);
  Field<ColorMatrix> gf_spatial;
  gf_spatial.set_mem_type(MemType::Acc);
  gf_spatial.init(geo, 6);
#ifdef QLAT_USE_ACC
  constexpr Int order_type = 1;
#else
  constexpr Int order_type = 0;
#endif
  const Int num_color_vec = num_field * 4;
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf_spatial.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const array<ComplexD, 6>& mfv = mom_factors();
    Vector<ComplexD> v = ff.get_elems(index);
    for (Int id_field = 0; id_field < num_field; ++id_field) {
      const WilsonVector& wv = ffv_vec[id_field].get_elem(index);
      for (Int s = 0; s < 4; ++s) {
        for (Int c = 0; c < 3; ++c) {
          if (order_type == 0) {
            v.p[c * num_color_vec + id_field * 4 + s] = wv.p[s * 3 + c];
          } else if (order_type == 1) {
            v.p[id_field * 12 + s * 3 + c] = wv.p[s * 3 + c];
          } else {
            qassert(false);
          }
        }
      }
    }
    Vector<ColorMatrix> gfv = gf_spatial.get_elems(index);
    for (Int dir = -dir_limit; dir < dir_limit; ++dir) {
      const Coordinate xl1 = coordinate_shifts(xl, dir);
      const Long index1 = geo.index_from_coordinate(xl1);
      ColorMatrix link = dir >= 0
                             ? gf.get_elem(index, dir)
                             : matrix_adjoint(gf.get_elem(index1, -dir - 1));
      link *= mfv[dir + 3];
      if (order_type == 0) {
        gfv[dir + 3] = link;
      } else if (order_type == 1) {
        gfv[dir + 3] = link;
      } else {
        qassert(false);
      }
    }
  });
  const RealD one_minus_coef = 1.0 - coef;
  for (Int iter = 0; iter < step; ++iter) {
    qswap(ff, ff1);
    qacc_for(index, geo.local_volume(), {
      const Geometry& geo = gf_spatial.geo();
      const Coordinate xl = geo.coordinate_from_index(index);
      const Vector<ColorMatrix> gfv = gf_spatial.get_elems_const(index);
      const Vector<ComplexD> v1 = ff1.get_elems_const(index);
      Vector<ComplexD> v = ff.get_elems(index);
      qassert(v.size() == v1.size());
      qassert(v.size() == num_color_vec * 3);
      qfor(i, v.size(), { v.p[i] = one_minus_coef * v1.p[i]; });
      for (Int dir = -dir_limit; dir < dir_limit; ++dir) {
        const Coordinate xl1 = coordinate_shifts(xl, dir);
        const Long index1 = geo.index_from_coordinate(xl1);
        const Vector<ComplexD> v11 = ff1.get_elems_const(index1);
        const ColorMatrix& link = gfv[dir + 3];
        alignas(16) const ComplexD* pl = link.p;
        alignas(64) const ComplexD* p11 = v11.p;
        alignas(64) ComplexD* p = v.p;
        if (order_type == 0) {
          for (Int c2 = 0; c2 < 3; ++c2) {
            alignas(64) const ComplexD* pp11 = &(p11[c2 * num_color_vec]);
            for (Int c1 = 0; c1 < 3; ++c1) {
              alignas(64) ComplexD* pp = &(p[c1 * num_color_vec]);
              const ComplexD lc = pl[c1 * 3 + c2];
              for (Int is = 0; is < num_color_vec; is += 4) {
                pp[is + 0] += lc * pp11[is + 0];
                pp[is + 1] += lc * pp11[is + 1];
                pp[is + 2] += lc * pp11[is + 2];
                pp[is + 3] += lc * pp11[is + 3];
              }
            }
          }
        } else if (order_type == 1) {
          for (Int ss = 0; ss < num_color_vec; ++ss) {
            alignas(16) const ComplexD* pp11 = &(p11[ss * 3]);
            alignas(16) ComplexD* pp = &(p[ss * 3]);
            pp[0] += pl[0] * pp11[0] + pl[0 + 1] * pp11[1] + pl[0 + 2] * pp11[2];
            pp[1] += pl[3] * pp11[0] + pl[3 + 1] * pp11[1] + pl[3 + 2] * pp11[2];
            pp[2] += pl[6] * pp11[0] + pl[6 + 1] * pp11[1] + pl[6 + 2] * pp11[2];
          }
        } else {
          qassert(false);
        }
      }
    });
  }
  qacc_for(index, geo.local_volume(), {
    const Vector<ComplexD> v = ff.get_elems_const(index);
    for (Int id_field = 0; id_field < num_field; ++id_field) {
      WilsonVector& wv = ffv_vec[id_field].get_elem(index);
      for (Int s = 0; s < 4; ++s) {
        for (Int c = 0; c < 3; ++c) {
          if (order_type == 0) {
            wv.p[s * 3 + c] = v.p[c * num_color_vec + id_field * 4 + s];
          } else if (order_type == 1) {
            wv.p[s * 3 + c] = v.p[id_field * 12 + s * 3 + c];
          } else {
            qassert(false);
          }
        }
      }
    }
  });
  // Can remove this when set_mem_type is used extensively.
  qfor(id_field, num_field,
       { ff_vec[id_field].set_mem_type(get_default_mem_type()); });
}

}  // namespace qlat
