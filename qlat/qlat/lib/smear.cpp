#include <qlat/qcd-acc.h>
#include <qlat/qcd-smear.h>
#include <qlat/vector_utils/utils_smear_vecs.h>

namespace qlat
{  //

static qacc ColorMatrix color_matrix_sub_invert(const ColorMatrix& x,
                                                const Int ind)
// get su2 submatrix of x and return the su3 matrix that
// has the inverse of this matrix in the relevant row and column
{
  const Int su2_index[][3] = {{0, 1, 2}, {0, 2, 1}, {1, 2, 0}};
  const Int i1 = su2_index[ind][0];
  const Int i2 = su2_index[ind][1];
  // const Int zero_rc = su2_index[ind][2];
  // project onto SU(2)
  double p0 = x(i1, i1).real() + x(i2, i2).real();
  double p1 = x(i1, i2).imag() + x(i2, i1).imag();
  double p2 = x(i1, i2).real() - x(i2, i1).real();
  double p3 = x(i1, i1).imag() - x(i2, i2).imag();
  const double psqr = sqrt(p0 * p0 + p1 * p1 + p2 * p2 + p3 * p3);
  ColorMatrix y;
  set_unit(y);
  if (psqr == 0.0) {
    return y;
  } else {
    double ipsqr = 1.0 / psqr;
    p0 *= ipsqr;
    p1 *= ipsqr;
    p2 *= ipsqr;
    p3 *= ipsqr;
    // fill with inverse
    y(i1, i1) = ComplexD(p0, -p3);
    y(i2, i2) = ComplexD(p0, p3);
    y(i1, i2) = ComplexD(-p2, -p1);
    y(i2, i1) = ComplexD(p2, -p1);
    return y;
  }
}

static qacc ColorMatrix color_matrix_su_projection(
    const ColorMatrix& x, const double tolerance = 1.0e-8)
{
  // usually takes ~5 hits, so just exit
  // if hits the max, as something is
  // probably very wrong.
  const Int max_iter = 10000;
  const ColorMatrix xdag = matrix_adjoint(x);
  ColorMatrix tmp = xdag;
  ColorMatrix y;
  set_unit(y);
  double old_tr = matrix_trace(xdag).real();
  for (Int i = 0; i < max_iter; i++) {
    // loop over su2 subgroups
    double diff = 0.0;
    for (Int j = 0; j < 3; j++) {
      const ColorMatrix inv = color_matrix_sub_invert(tmp, j);
      // y  .DotMEqual( inv, ycopy );
      y = inv * y;
      // tmp.DotMEqual( y, xdag );
      tmp = y * xdag;
      const double tr = matrix_trace(tmp).real();
      const double dtr = tr - old_tr;
      if (dtr > diff) {
        diff = dtr;
      }
      old_tr = tr;
    }
    // for single precision the difference seems
    // to never get below 1e-7 (not too suprising)
    if (diff < tolerance) {
      break;
    }
    qassert(i < max_iter - 1);
  }
  unitarize(y);
  return y;
}

static qacc ColorMatrix gf_link_ape_smear_no_comm(const GaugeField& gf,
                                                  const Coordinate& xl,
                                                  const Int mu,
                                                  const double alpha)
{
  return color_matrix_su_projection(
      (ComplexD)(1.0 - alpha) * gf.get_elem(xl, mu) +
      (ComplexD)(alpha / 6.0) * gf_staple_no_comm(gf, xl, mu));
}

static void gf_ape_smear_no_comm(GaugeField& gf, const GaugeField& gf0,
                                 const double alpha)
{
  TIMER_VERBOSE("gf_ape_smear_no_comm");
  Qassert(&gf != &gf0);
  const Geometry geo0 = gf0.get_geo();
  const Geometry geo = geo_resize(geo0);
  gf.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (Int mu = 0; mu < DIMN; ++mu) {
      v[mu] = gf_link_ape_smear_no_comm(gf0, xl, mu, alpha);
    }
  });
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

static qacc ColorMatrix gf_link_spatial_ape_smear_no_comm(const GaugeField& gf,
                                                          const Coordinate& xl,
                                                          const Int mu,
                                                          const double alpha)
{
  const double multi = mu == 3 ? 6.0 : 4.0;
  return color_matrix_su_projection(
      (ComplexD)(1.0 - alpha) * gf.get_elem(xl, mu) +
      (ComplexD)(alpha / multi) * gf_spatial_staple_no_comm(gf, xl, mu));
}

static void gf_spatial_ape_smear_no_comm(GaugeField& gf, const GaugeField& gf0,
                                         const double alpha)
{
  TIMER_VERBOSE("gf_spatial_ape_smear_no_comm");
  Qassert(&gf != &gf0);
  const Geometry geo0 = gf0.get_geo();
  const Geometry geo = geo_resize(geo0);
  gf.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (Int mu = 0; mu < 3; ++mu) {
      // No need to smear the temperal link (mu == 3)
      v[mu] = gf_link_spatial_ape_smear_no_comm(gf0, xl, mu, alpha);
    }
  });
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
                                        const Coordinate& xl, const Int mu,
                                        const Int nu, const Int rho,
                                        const double alpha3)
{
  ColorMatrix ret;
  set_zero(ret);
  const Coordinate xl_mu = coordinate_shifts(xl, mu);
  for (Int m = 0; m < DIMN; ++m) {
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
                                        const Coordinate& xl, const Int mu,
                                        const Int nu, const double alpha2,
                                        const double alpha3)
{
  ColorMatrix ret;
  set_zero(ret);
  const Coordinate xl_mu = coordinate_shifts(xl, mu);
  for (Int m = 0; m < DIMN; ++m) {
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
                                        const Coordinate& xl, const Int mu,
                                        const double alpha1,
                                        const double alpha2,
                                        const double alpha3)
{
  ColorMatrix ret;
  set_zero(ret);
  const Coordinate xl_mu = coordinate_shifts(xl, mu);
  for (Int m = 0; m < DIMN; ++m) {
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
                                      const Coordinate& xl, const Int mu,
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
  Qassert(&gf != &gf0);
  const Geometry& geo = gf0.geo();
  gf.init(geo_resize(geo));
  Qassert(is_matching_geo(geo, gf.geo()));
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (Int mu = 0; mu < DIMN; ++mu) {
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

template <class T>
static void prop_smear(Propagator4dT<T>& prop, const GaugeFieldT<T>& gf1,
                       const double coef, const Int step,
                       const CoordinateD& mom = CoordinateD(),
                       const bool smear_in_time_dir = false)
// gf1 is left_expanded and refreshed
// set_left_expanded_gauge_field(gf1, gf)
// prop is of normal size
{
  TIMER_FLOPS("prop_smear");
  const Int n_avg = smear_in_time_dir ? 8 : 6;
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
  const Int dir_limit = smear_in_time_dir ? 4 : 3;
  array<ComplexD, 8> mom_factors_v;
  box<array<ComplexD, 8>> mom_factors(
      mom_factors_v);  // (array<ComplexD, 8>());
  for (Int i = 0; i < 8; ++i) {
    const Int dir = i - 4;
    const double phase = dir >= 0 ? mom[dir] : -mom[-dir - 1];
    mom_factors()[i] = qpolar(coef / n_avg, -phase);
  }
  Propagator4dT<T> prop1;
  prop1.init(geo1);
  for (Int i = 0; i < step; ++i) {
    prop1 = prop;
    refresh_expanded_1(prop1);
    qacc_for(index, geo.local_volume(), {
      const Coordinate xl = prop.geo().coordinate_from_index(index);
      WilsonMatrixT<T>& wm = prop.get_elem(xl);
      wm *= 1 - coef;
      for (Int dir = -dir_limit; dir < dir_limit; ++dir) {
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

void prop_smear(Propagator4dT<RealD>& prop, const GaugeFieldT<RealD>& gf1,
                const double coef, const Int step, const CoordinateD& mom,
                const bool smear_in_time_dir)
{
  prop_smear<RealD>(prop, gf1, coef, step, mom, smear_in_time_dir);
}

void prop_smear_qlat_convension(Propagator4dT<RealD>& prop,
                                const GaugeFieldT<RealD>& gf, const double coef,
                                const Int step, const CoordinateD& mom,
                                const bool smear_in_time_dir, const Int mode)
{
  prop_smear_qlat_convension<RealD, RealD>(prop, gf, coef, step, mom,
                                           smear_in_time_dir, mode);
}

void prop_smear_qlat_convension(Propagator4dT<RealF>& prop,
                                const GaugeFieldT<RealF>& gf, const double coef,
                                const Int step, const CoordinateD& mom,
                                const bool smear_in_time_dir, const Int mode)
{
  prop_smear_qlat_convension<RealF, RealF>(prop, gf, coef, step, mom,
                                           smear_in_time_dir, mode);
}

static void prop_spatial_smear_no_comm_acc(std::vector<FermionField4d>& ff_vec,
                                           const GaugeField& gf,
                                           const RealD coef, const Long step,
                                           const CoordinateD& mom)
// `gf` and each of `ff_vec` should contain entire time slices.
// No communication will be performed.
// More suitable for ACC execution, but should also work on CPU.
{
  TIMER_FLOPS("prop_spatial_smear_no_comm_acc(ff_vec,gf,coef,step,mom)");
  const Geometry geo = gf.geo.get();
  const Int num_field = ff_vec.size();
  const Int n_avg = 6;
  const Long v_gb = geo.local_volume() * num_field * 4;
  const Long flops = v_gb * step * n_avg * (3 * (3 * 6 + 2 * 2));
  timer.flops += flops;
  if (0 == step) {
    return;
  }
  const Int dir_limit = 3;
  array<ComplexD, 6> mom_factors_v;
  for (Int i = 0; i < 6; ++i) {
    const Int dir = i - 3;
    const RealD phase = dir >= 0 ? mom[dir] : -mom[-dir - 1];
    mom_factors_v[i] = qpolar(coef / n_avg, -phase);
  }
  box<array<ComplexD, 6>> mom_factors(mom_factors_v,
                                      MemType::Acc);  // (array<ComplexD, 8>());
  const Int t_size = geo.total_site()[3];
  Qassert(geo.geon.size_node == Coordinate(1, 1, 1, t_size));
  Qassert(geo.is_only_local);
  Qassert(num_field >= 0);
  if (num_field == 0) {
    return;
  }
  vector<FermionField4d> ffv_vec(num_field, MemType::Cpu);
  set_zero(ffv_vec);
  qfor(id_field, num_field, {
    Qassert(ff_vec[id_field].geo.get() == geo);
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
  const Int num_color_vec = num_field * 4;
  const Int chunk_num_color_vec = get_env_long_default(
      "q_prop_spatial_smear_no_comm_chunk_num_color_vec", 1);
  const Int num_chunk_color_vec = num_color_vec / chunk_num_color_vec;
  Qassert(num_chunk_color_vec * chunk_num_color_vec == num_color_vec);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf_spatial.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const array<ComplexD, 6>& mfv = mom_factors();
    Vector<ComplexD> v = ff.get_elems(index);
    for (Int id_field = 0; id_field < num_field; ++id_field) {
      const WilsonVector& wv = ffv_vec[id_field].get_elem(index);
      for (Int s = 0; s < 4; ++s) {
        for (Int c = 0; c < 3; ++c) {
          v.p[id_field * 12 + s * 3 + c] = wv.p[s * 3 + c];
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
      gfv[dir + 3] = link;
    }
  });
  const RealD one_minus_coef = 1.0 - coef;
  {
    TIMER_FLOPS("prop_spatial_smear_no_comm_acc-smear");
    timer.flops += flops;
    for (Int iter = 0; iter < step; ++iter) {
      qswap(ff, ff1);
      qacc_for(idx, geo.local_volume() * num_chunk_color_vec, {
        const Long index = idx / num_chunk_color_vec;
        const Int id_chunk_color_vec = idx % num_chunk_color_vec;
        const Geometry& geo = gf_spatial.geo();
        const Coordinate xl = geo.coordinate_from_index(index);
        Vector<ComplexD> v = ff.get_elems(index);
        qassert(v.size() == num_color_vec * 3);
        {
          const Vector<ComplexD> v1 = ff1.get_elems_const(index);
          qassert(v1.size() == v.size());
          alignas(16) const ComplexD* p1 =
              &(v1.p[id_chunk_color_vec * chunk_num_color_vec * 3]);
          alignas(16) ComplexD* p =
              &(v.p[id_chunk_color_vec * chunk_num_color_vec * 3]);
          for (Int ss = 0; ss < chunk_num_color_vec; ++ss) {
            alignas(16) const ComplexD* pp1 = &(p1[ss * 3]);
            alignas(16) ComplexD* pp = &(p[ss * 3]);
            pp[0] = one_minus_coef * pp1[0];
            pp[1] = one_minus_coef * pp1[1];
            pp[2] = one_minus_coef * pp1[2];
          }
        }
        const Vector<ColorMatrix> gfv = gf_spatial.get_elems_const(index);
        for (Int dir = -dir_limit; dir < dir_limit; ++dir) {
          const Coordinate xl1 = coordinate_shifts(xl, dir);
          const Long index1 = geo.index_from_coordinate(xl1);
          const Vector<ComplexD> v11 = ff1.get_elems_const(index1);
          qassert(v11.size() == v.size());
          const ColorMatrix& link = gfv[dir + 3];
          alignas(16) const ComplexD* pl = link.p;
          alignas(16) const ComplexD* p11 =
              &(v11.p[id_chunk_color_vec * chunk_num_color_vec * 3]);
          alignas(16) ComplexD* p =
              &(v.p[id_chunk_color_vec * chunk_num_color_vec * 3]);
          for (Int ss = 0; ss < chunk_num_color_vec; ++ss) {
            alignas(16) const ComplexD* pp11 = &(p11[ss * 3]);
            alignas(16) ComplexD* pp = &(p[ss * 3]);
            pp[0] +=
                pl[0] * pp11[0] + pl[0 + 1] * pp11[1] + pl[0 + 2] * pp11[2];
            pp[1] +=
                pl[3] * pp11[0] + pl[3 + 1] * pp11[1] + pl[3 + 2] * pp11[2];
            pp[2] +=
                pl[6] * pp11[0] + pl[6 + 1] * pp11[1] + pl[6 + 2] * pp11[2];
          }
        }
      });
    }
  }
  ff1.init();
  gf_spatial.init();
  qacc_for(index, geo.local_volume(), {
    const Vector<ComplexD> v = ff.get_elems_const(index);
    for (Int id_field = 0; id_field < num_field; ++id_field) {
      WilsonVector& wv = ffv_vec[id_field].get_elem(index);
      for (Int s = 0; s < 4; ++s) {
        for (Int c = 0; c < 3; ++c) {
          wv.p[s * 3 + c] = v.p[id_field * 12 + s * 3 + c];
        }
      }
    }
  });
  // Can remove this when set_mem_type is used extensively.
  qfor(id_field, num_field,
       { ff_vec[id_field].set_mem_type(get_default_mem_type()); });
}

static void prop_spatial_smear_no_comm_cpu(std::vector<FermionField4d>& ff_vec,
                                           const GaugeField& gf,
                                           const RealD coef, const Long step,
                                           const CoordinateD& mom)
// `gf` and each of `ff_vec` should contain entire time slices.
// No communication will be performed.
// More suitable for CPU execution, but should also work on ACC.
{
  TIMER_FLOPS("prop_spatial_smear_no_comm_cpu(ff_vec,gf,coef,step,mom)");
  const Geometry geo = gf.geo.get();
  const Int num_field = ff_vec.size();
  const Int n_avg = 6;
  const Long v_gb = geo.local_volume() * num_field * 4;
  const Long flops = v_gb * step * n_avg * (3 * (3 * 6 + 2 * 2));
  timer.flops += flops;
  if (0 == step) {
    return;
  }
  const Int dir_limit = 3;
  array<ComplexD, 6> mom_factors_v;
  for (Int i = 0; i < 6; ++i) {
    const Int dir = i - 3;
    const RealD phase = dir >= 0 ? mom[dir] : -mom[-dir - 1];
    mom_factors_v[i] = qpolar(coef / n_avg, -phase);
  }
  box<array<ComplexD, 6>> mom_factors(mom_factors_v,
                                      MemType::Acc);  // (array<ComplexD, 8>());
  const Int t_size = geo.total_site()[3];
  Qassert(geo.geon.size_node == Coordinate(1, 1, 1, t_size));
  Qassert(geo.is_only_local);
  Qassert(num_field >= 0);
  if (num_field == 0) {
    return;
  }
  vector<FermionField4d> ffv_vec(num_field, MemType::Cpu);
  set_zero(ffv_vec);
  qfor(id_field, num_field, {
    Qassert(ff_vec[id_field].geo.get() == geo);
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
          v.p[c * num_color_vec + id_field * 4 + s] = wv.p[s * 3 + c];
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
      gfv[dir + 3] = link;
    }
  });
  const RealD one_minus_coef = 1.0 - coef;
  {
    TIMER_FLOPS("prop_spatial_smear_no_comm_cpu-smear");
    timer.flops += flops;
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
        }
      });
    }
  }
  ff1.init();
  gf_spatial.init();
  qacc_for(index, geo.local_volume(), {
    const Vector<ComplexD> v = ff.get_elems_const(index);
    for (Int id_field = 0; id_field < num_field; ++id_field) {
      WilsonVector& wv = ffv_vec[id_field].get_elem(index);
      for (Int s = 0; s < 4; ++s) {
        for (Int c = 0; c < 3; ++c) {
          wv.p[s * 3 + c] = v.p[c * num_color_vec + id_field * 4 + s];
        }
      }
    }
  });
  // Can remove this when set_mem_type is used extensively.
  qfor(id_field, num_field,
       { ff_vec[id_field].set_mem_type(get_default_mem_type()); });
}

void prop_spatial_smear_no_comm(std::vector<FermionField4d>& ff_vec,
                                const GaugeField& gf, const RealD coef,
                                const Long step, const CoordinateD& mom)
// `gf` and each of `ff_vec` should contain entire time slices.
// No communication will be performed.
{
#ifdef QLAT_USE_ACC
  static const std::string q_prop_spatial_smear_no_comm =
      get_env_default("q_prop_spatial_smear_no_comm", "acc");
#else
  static const std::string q_prop_spatial_smear_no_comm =
      get_env_default("q_prop_spatial_smear_no_comm", "cpu");
#endif
  if (q_prop_spatial_smear_no_comm == "acc") {
    prop_spatial_smear_no_comm_acc(ff_vec, gf, coef, step, mom);
  } else if (q_prop_spatial_smear_no_comm == "cpu") {
    prop_spatial_smear_no_comm_cpu(ff_vec, gf, coef, step, mom);
  } else {
    Qassert(false);
  }
}

void gf_reduce_half(GaugeField& hgf, const GaugeField& gf)
// xl = coordinate_shifts(hxl * 2, 0)
{
  TIMER_VERBOSE("gf_reduce_half(hgf,gf)");
  const Geometry geo = gf.get_geo();
  Geometry hgeo;
  hgeo.init(geo.geon, geo.node_site / 2);
  Qassert(geo.node_site == hgeo.node_site * 2);
  Qassert(DIMN == gf.multiplicity);
  hgf.init(hgeo, DIMN);
  qacc_for(hindex, hgeo.local_volume(), {
    const Geometry& hgeo = hgf.geo();
    const Coordinate hxl = hgeo.coordinate_from_index(hindex);
    const Coordinate xl = coordinate_shifts(hxl * 2, 0);
    for (Int m = 0; m < DIMN; ++m) {
      hgf.get_elem(hxl, m) =
          gf.get_elem(xl, m) * gf.get_elem(coordinate_shifts(xl, m), m);
    }
  });
}

}  // namespace qlat
