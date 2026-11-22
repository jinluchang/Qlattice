#include <qlat/gauge-fixing.h>

#include <algorithm>
#include <cmath>

namespace qlat
{

namespace
{

const Int exp_order = 20;

Geometry get_expanded_geometry(const Geometry& geo)
{
  return geo_resize(geo, Coordinate(1, 1, 1, 1), Coordinate(1, 1, 1, 1));
}

void gf_set_random_local_gradient(Field<ColorMatrix>& D, const GaugeField& gf,
                                  const GaugeField& gfe,
                                  const std::vector<Int>& dirs)
// gfe is the expanded copy of gf (neighbour links are available)
{
  const Int n_dirs = (Int)dirs.size();
  array<Int, 4> dirs_arr;
  dirs_arr.fill(0);
  for (Int i = 0; i < n_dirs; ++i) {
    dirs_arr[i] = dirs[i];
  }
  const Geometry& geo = gf.geo();
  D.init(geo, 1);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrix> v = gf.get_elems_const(xl);
    ColorMatrix K;
    set_zero(K);
    for (Int i = 0; i < n_dirs; ++i) {
      const Int mu = dirs_arr[i];
      Coordinate xlm = xl;
      xlm[mu] -= 1;
      const Vector<ColorMatrix> vm = gfe.get_elems_const(xlm);
      K += v[mu] + matrix_adjoint(vm[mu]);
    }
    // D(x) = ((K^dagger - K) / 2) with the trace part removed
    D.get_elem(xl) = make_tr_less_anti_herm_matrix(matrix_adjoint(K));
  });
}

void gf_update_gauge_field(GaugeField& gf, const Field<ColorMatrix>& g,
                           const Field<ColorMatrix>& ge)
// ge is the expanded copy of g (g(x+mu) is available)
{
  const Geometry& geo = gf.geo();
  for (Int mu = 0; mu < 4; ++mu) {
    qacc_for(index, geo.local_volume(), {
      const Geometry& geo = gf.geo();
      const Coordinate xl = geo.coordinate_from_index(index);
      Coordinate xlp = xl;
      xlp[mu] += 1;
      ColorMatrix& u = gf.get_elem(xl, mu);
      u = g.get_elem(xl) * u * matrix_adjoint(ge.get_elem(xlp));
    });
  }
}

}  // namespace

void gf_gauge_gradient(Field<ColorMatrix>& D, const GaugeField& gf,
                       const std::vector<Int>& dirs)
{
  TIMER("gf_gauge_gradient");
  GaugeField gfe;
  gfe.init(get_expanded_geometry(gf.geo()));
  gfe = gf;
  refresh_expanded(gfe);
  gf_set_random_local_gradient(D, gf, gfe, dirs);
}

RealD gf_gauge_functional(const GaugeField& gf, const std::vector<Int>& dirs)
{
  TIMER("gf_gauge_functional");
  const Int n_dirs = (Int)dirs.size();
  array<Int, 4> dirs_arr;
  dirs_arr.fill(0);
  for (Int i = 0; i < n_dirs; ++i) {
    dirs_arr[i] = dirs[i];
  }
  const Geometry& geo = gf.geo();
  Field<RealD> f;
  f.init(geo, 1);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrix> v = gf.get_elems_const(xl);
    RealD sum = 0.0;
    for (Int i = 0; i < n_dirs; ++i) {
      sum += matrix_trace(v[dirs_arr[i]]).real();
    }
    f.get_elem(xl) = sum;
  });
  RealD sum = field_glb_sum(f)[0];
  return sum / (geo.total_volume() * (RealD)n_dirs * 3.0);
}

RealD gf_gauge_divergence_from_gradient(const Field<ColorMatrix>& D)
{
  TIMER("gf_gauge_divergence_from_gradient");
  const Geometry& geo = D.geo();
  Field<RealD> f;
  f.init(geo, 1);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = D.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const ColorMatrix& d = D.get_elem(xl);
    RealD norm2 = 0.0;
    for (Int i = 0; i < 3; ++i) {
      for (Int j = 0; j < 3; ++j) {
        norm2 += qnorm(d(i, j));
      }
    }
    f.get_elem(xl) = norm2;
  });
  const RealD sum = field_glb_sum(f)[0];
  return std::sqrt(sum / (geo.total_volume() * 3.0));
}

RealD gf_gauge_divergence(const GaugeField& gf, const std::vector<Int>& dirs)
{
  Field<ColorMatrix> D;
  gf_gauge_gradient(D, gf, dirs);
  return gf_gauge_divergence_from_gradient(D);
}

RealD gf_gauge_unitarity_error(const GaugeField& gf)
{
  TIMER("gf_gauge_unitarity_error");
  const Geometry& geo = gf.geo();
  Field<RealD> f;
  f.init(geo, 1);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrix> v = gf.get_elems_const(xl);
    RealD err = 0.0;
    for (Int mu = 0; mu < 4; ++mu) {
      const ColorMatrix m = v[mu] * matrix_adjoint(v[mu]);
      for (Int i = 0; i < 3; ++i) {
        for (Int j = 0; j < 3; ++j) {
          const ComplexD delta =
              m(i, j) - (i == j ? ComplexD(1.0) : ComplexD(0.0));
          // std::abs(ComplexD) is a host-only std::complex overload that is not
          // device callable inside this qacc_for CUDA lambda;
          // std::sqrt(qnorm(delta)) is the device-safe equivalent (qnorm is the
          // squared magnitude under QLAT_USE_ACC).
          err = std::max(err, std::sqrt(qnorm(delta)));
        }
      }
    }
    f.get_elem(xl) = err;
  });
  return field_glb_max(f)[0];
}

GaugeFixingInfo gf_fix_gauge(GaugeField& gf, const std::vector<Int>& dirs,
                             const Int max_iter, const RealD eps,
                             const RealD alpha, const bool verbose)
{
  TIMER("gf_fix_gauge");
  Qassert(gf.initialized);
  Qassert(dirs.size() > 0);
  for (Int i = 0; i < (Int)dirs.size(); ++i) {
    Qassert(0 <= dirs[i] and dirs[i] < 4);
  }
  const Geometry& geo = gf.geo();
  const Geometry geo_expanded = get_expanded_geometry(geo);
  GaugeField gfe;
  gfe.init(geo_expanded);
  Field<ColorMatrix> D;
  D.init(geo, 1);
  Field<ColorMatrix> g;
  g.init(geo, 1);
  Field<ColorMatrix> ge;
  ge.init(geo_expanded, 1);
  GaugeFixingInfo info;
  info.converged = false;
  info.n_iter = 0;
  info.max_iter = max_iter;
  info.eps = eps;
  info.alpha = alpha;
  info.functional = 0.0;
  info.divergence = 0.0;
  for (Int n_iter = 1; n_iter <= max_iter; ++n_iter) {
    gfe = gf;
    refresh_expanded(gfe);
    gf_set_random_local_gradient(D, gf, gfe, dirs);
    const RealD divergence = gf_gauge_divergence_from_gradient(D);
    const RealD functional = gf_gauge_functional(gf, dirs);
    info.functional_list.push_back(functional);
    info.divergence_list.push_back(divergence);
    info.n_iter = n_iter;
    if (verbose) {
      displayln_info(
          ssprintf("gf_fix_gauge: sweep %d: functional = %.15e;"
                   " divergence = %.3e",
                   n_iter, functional, divergence));
    }
    if (divergence <= eps) {
      info.converged = true;
      break;
    }
    qacc_for(index, geo.local_volume(), {
      const Geometry& geo = gf.geo();
      const Coordinate xl = geo.coordinate_from_index(index);
      const ColorMatrix h = (alpha * D.get_elem(xl));
      g.get_elem(xl) = make_matrix_exp(h, exp_order);
    });
    ge = g;
    refresh_expanded(ge);
    gf_update_gauge_field(gf, g, ge);
  }
  gf_gauge_gradient(D, gf, dirs);
  info.functional = gf_gauge_functional(gf, dirs);
  info.divergence = gf_gauge_divergence_from_gradient(D);
  return info;
}

GaugeFixingInfo gf_fix_landau_gauge(GaugeField& gf, const Int max_iter,
                                    const RealD eps, const RealD alpha,
                                    const bool verbose)
{
  const std::vector<Int> dirs = {0, 1, 2, 3};
  return gf_fix_gauge(gf, dirs, max_iter, eps, alpha, verbose);
}

GaugeFixingInfo gf_fix_coulomb_gauge(GaugeField& gf, const Int max_iter,
                                     const RealD eps, const RealD alpha,
                                     const bool verbose)
{
  const std::vector<Int> dirs = {1, 2, 3};
  return gf_fix_gauge(gf, dirs, max_iter, eps, alpha, verbose);
}

}  // namespace qlat
