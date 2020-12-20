#include <qlat/qlat.h>

namespace qlat
{  //

inline void scalar_inversion(Field<Complex>& sol, const Field<Complex>& src,
                             const double mass,
                             const CoordinateD momtwist = CoordinateD())
// the mass is not necessarily the exponent of the exponential fall off
{
  TIMER("scalar_inversion");
  const Geometry geo = geo_resize(src.geo());
  const Coordinate total_site = geo.total_site();
  sol.init(geo);
  sol = src;
  fft_complex_field(sol, true);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    CoordinateD kk;
    double s2 = 0.0;
    for (int i = 0; i < DIMN; i++) {
      kg[i] = smod(kg[i], total_site[i]);
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (double)total_site[i];
      s2 += 4.0 * sqr(std::sin(kk[i] / 2.0));
    }
    const double fac = 1.0 / (s2 + sqr(mass));
    Vector<Complex> v = sol.get_elems(kl);
    for (int m = 0; m < v.size(); ++m) {
      v[m] *= fac;
    }
  }
  fft_complex_field(sol, false);
  sol *= 1.0 / geo.total_volume();
}

inline void scalar_derivative(Field<Complex>& sol, const Field<Complex>& src,
                              const CoordinateD momtwist = CoordinateD())
// v[m*4 + mu] = sv[m] * std::sin(kk[mu]);
{
  TIMER("scalar_derivative");
  const Geometry geo = geo_reform(src.geo(), src.geo().multiplicity * 4);
  sol.init(geo);
  qassert(sol.geo() == geo);
  const Coordinate total_site = geo.total_site();
  Field<Complex> src_mom;
  src_mom.init(geo_resize(src.geo()));
  src_mom = src;
  fft_complex_field(src_mom, true);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    CoordinateD kk;
    for (int i = 0; i < DIMN; i++) {
      kg[i] = smod(kg[i], total_site[i]);
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (double)total_site[i];
    }
    const Vector<Complex> sv = src_mom.get_elems_const(kl);
    Vector<Complex> v = sol.get_elems(kl);
    for (int m = 0; m < sv.size(); ++m) {
      for (int mu = 0; mu < 4; ++mu) {
        v[m * 4 + mu] = sv[m] * std::sin(kk[mu]);
      }
    }
  }
  fft_complex_field(sol, false);
  sol *= 1.0 / geo.total_volume();
}

inline void scalar_divergence(Field<Complex>& sol, const Field<Complex>& src,
                              const CoordinateD momtwist = CoordinateD())
// v[m] += sv[m*4+mu] * std::sin(kk[mu]);
{
  TIMER("scalar_derivative");
  const Geometry geo = geo_reform(src.geo(), src.geo().multiplicity / 4);
  sol.init(geo);
  qassert(sol.geo() == geo);
  const Coordinate total_site = geo.total_site();
  Field<Complex> src_mom;
  src_mom.init(geo_resize(src.geo()));
  src_mom = src;
  fft_complex_field(src_mom, true);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    CoordinateD kk;
    for (int i = 0; i < DIMN; i++) {
      kg[i] = smod(kg[i], total_site[i]);
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (double)total_site[i];
    }
    const Vector<Complex> sv = src_mom.get_elems_const(kl);
    Vector<Complex> v = sol.get_elems(kl);
    for (int m = 0; m < v.size(); ++m) {
      v[m] = 0;
      for (int mu = 0; mu < 4; ++mu) {
        v[m] += sv[m * 4 + mu] * std::sin(kk[mu]);
      }
    }
  }
  fft_complex_field(sol, false);
  sol *= 1.0 / geo.total_volume();
}

inline void set_pion_photon_photon_vertex_two_end(
    FieldM<Complex, 4 * 4>& pion, const FieldM<Complex, 4>& photon1,
    const FieldM<Complex, 4>& photon2, const double m_vector, const double f_pi)
{
  TIMER("set_pion_photon_photon_vertex_two_end");
  qassert(photon1.geo() == photon2.geo());
  const Geometry geo = geo_reform(photon1.geo(), 4 * 4, 0);
  FieldM<Complex, 4 * 4> pinv1, pinv2, p1, p2;
  scalar_derivative(p1, photon1);
  scalar_derivative(p2, photon2);
  scalar_inversion(pinv1, p1, m_vector);
  scalar_inversion(pinv2, p2, m_vector);
  pion.init(geo);
  set_zero(pion);
  const Complex fac = 1.0 / (4.0 * sqr(PI) * f_pi) * ii * sqr(m_vector) / 2.0;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<Complex> v1 = pinv1.get_elems_const(xl);
    const Vector<Complex> v2 = pinv2.get_elems_const(xl);
    const Vector<Complex> d1 = p1.get_elems_const(xl);
    const Vector<Complex> d2 = p2.get_elems_const(xl);
    Vector<Complex> pv = pion.get_elems(xl);
    for (int mu = 0; mu < 4; ++mu) {
      for (int nu = 0; nu < 4; ++nu) {
        if (nu == mu) {
          continue;
        }
        for (int rho = 0; rho < 4; ++rho) {
          if (rho == mu or rho == nu) {
            continue;
          }
          for (int sigma = 0; sigma < 4; ++sigma) {
            if (sigma == mu or sigma == nu or sigma == rho) {
              continue;
            }
            pv[mu * 4 + nu] += fac *
                               (Complex)epsilon_tensor(mu, nu, rho, sigma) *
                               (v1[mu * 4 + rho] * d2[nu * 4 + sigma] +
                                d1[mu * 4 + rho] * v2[nu * 4 + sigma]);
          }
        }
      }
    }
  }
}

inline void set_pion_photon_photon_vertex_vmd(FieldM<Complex, 4 * 4>& pion,
                                              const FieldM<Complex, 4>& photon1,
                                              const FieldM<Complex, 4>& photon2,
                                              const double m_vector,
                                              const double f_pi)
{
  TIMER("set_pion_photon_photon_vertex_vmd");
  qassert(photon1.geo() == photon2.geo());
  const Geometry geo = geo_reform(photon1.geo(), 4 * 4, 0);
  FieldM<Complex, 4 * 4> pinv1, pinv2, p1, p2;
  scalar_derivative(p1, photon1);
  scalar_derivative(p2, photon2);
  scalar_inversion(pinv1, p1, m_vector);
  scalar_inversion(pinv2, p2, m_vector);
  pion.init(geo);
  set_zero(pion);
  const Complex fac = 1.0 / (4.0 * sqr(PI) * f_pi) * ii * sqr(sqr(m_vector));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<Complex> v1 = pinv1.get_elems_const(xl);
    const Vector<Complex> v2 = pinv2.get_elems_const(xl);
    Vector<Complex> pv = pion.get_elems(xl);
    for (int mu = 0; mu < 4; ++mu) {
      for (int nu = 0; nu < 4; ++nu) {
        if (nu == mu) {
          continue;
        }
        for (int rho = 0; rho < 4; ++rho) {
          if (rho == mu or rho == nu) {
            continue;
          }
          for (int sigma = 0; sigma < 4; ++sigma) {
            if (sigma == mu or sigma == nu or sigma == rho) {
              continue;
            }
            pv[mu * 4 + nu] += fac *
                               (Complex)epsilon_tensor(mu, nu, rho, sigma) *
                               v1[mu * 4 + rho] * v2[nu * 4 + sigma];
          }
        }
      }
    }
  }
}

inline void set_pion_photon_photon_vertex_lmd(FieldM<Complex, 4 * 4>& pion,
                                              const FieldM<Complex, 4>& photon1,
                                              const FieldM<Complex, 4>& photon2,
                                              const double m_vector,
                                              const double f_pi)
{
  TIMER_VERBOSE("set_pion_photon_photon_vertex_lmd");
  FieldM<Complex, 4 * 4> pion_tmp;
  set_pion_photon_photon_vertex_vmd(pion, photon1, photon2, m_vector, f_pi);
  set_pion_photon_photon_vertex_two_end(pion_tmp, photon1, photon2, m_vector,
                                        f_pi);
  const double coef = 8.0 / 3.0 * sqr(PI * f_pi / m_vector);
  pion *= 1.0 - coef;
  pion_tmp *= coef;
  pion += pion_tmp;
}

inline void set_pion_photon_photon_vertex(FieldM<Complex, 4 * 4>& pion,
                                          const FieldM<Complex, 4>& photon1,
                                          const FieldM<Complex, 4>& photon2,
                                          const double m_vector,
                                          const double f_pi)
// vertex function
{
  // SADJUST ME
  // set_pion_photon_photon_vertex_vmd(pion, photon1, photon2, m_vector, f_pi);
  // set_pion_photon_photon_vertex_two_end(pion, photon1, photon2, m_vector,
  // f_pi);
  set_pion_photon_photon_vertex_lmd(pion, photon1, photon2, m_vector, f_pi);
}

inline void set_photon_pion_photon_vertex_two_end(
    FieldM<Complex, 4>& photon1, const FieldM<Complex, 1>& pion,
    const FieldM<Complex, 4>& photon2, const double m_vector, const double f_pi)
{
  TIMER("set_photon_pion_photon_vertex_two_end");
  const Geometry geo = geo_reform(photon2.geo(), 4, 0);
  qassert(is_matching_geo(pion.geo(), photon2.geo()));
  photon1.init(geo);
  qassert(photon1.geo() == geo);
  set_zero(photon1);
  FieldM<Complex, 4 * 4> pinv1, pinv2, p1, p2;
  scalar_derivative(p2, photon2);
  scalar_inversion(pinv2, p2, m_vector);
  p1.init(geo);
  pinv1.init(geo);
  const Complex fac = 1.0 / (4.0 * sqr(PI) * f_pi) * ii * sqr(m_vector) / 2.0;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<Complex> d1 = p1.get_elems(xl);
    const Vector<Complex> d2 = p2.get_elems_const(xl);
    Vector<Complex> v1 = pinv1.get_elems(xl);
    const Vector<Complex> v2 = pinv2.get_elems_const(xl);
    const Complex& v = pion.get_elem(xl);
    for (int mu = 0; mu < 4; ++mu) {
      for (int nu = 0; nu < 4; ++nu) {
        if (nu == mu) {
          continue;
        }
        for (int rho = 0; rho < 4; ++rho) {
          if (rho == mu or rho == nu) {
            continue;
          }
          for (int sigma = 0; sigma < 4; ++sigma) {
            if (sigma == mu or sigma == nu or sigma == rho) {
              continue;
            }
            const Complex c =
                -fac * (Complex)epsilon_tensor(mu, nu, rho, sigma) * v;
            d1[mu * 4 + rho] += c * v2[nu * 4 + sigma];
            v1[mu * 4 + rho] += c * d2[nu * 4 + sigma];
          }
        }
      }
    }
  }
  scalar_inversion(pinv1, pinv1, m_vector);
  p1 += pinv1;
  scalar_divergence(photon1, p1);
}

inline void set_photon_pion_photon_vertex_vmd(FieldM<Complex, 4>& photon1,
                                              const FieldM<Complex, 1>& pion,
                                              const FieldM<Complex, 4>& photon2,
                                              const double m_vector,
                                              const double f_pi)
{
  TIMER("set_photon_pion_photon_vertex_vmd");
  const Geometry geo = geo_reform(photon2.geo(), 4, 0);
  qassert(is_matching_geo(pion.geo(), photon2.geo()));
  photon1.init(geo);
  qassert(photon1.geo() == geo);
  set_zero(photon1);
  FieldM<Complex, 4 * 4> pinv1, pinv2, p1, p2;
  scalar_derivative(p2, photon2);
  scalar_inversion(pinv2, p2, m_vector);
  p1.init(geo);
  pinv1.init(geo);
  const Complex fac = 1.0 / (4.0 * sqr(PI) * f_pi) * ii * sqr(sqr(m_vector));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<Complex> v1 = pinv1.get_elems(xl);
    const Vector<Complex> v2 = pinv2.get_elems_const(xl);
    const Complex& v = pion.get_elem(xl);
    for (int mu = 0; mu < 4; ++mu) {
      for (int nu = 0; nu < 4; ++nu) {
        if (nu == mu) {
          continue;
        }
        for (int rho = 0; rho < 4; ++rho) {
          if (rho == mu or rho == nu) {
            continue;
          }
          for (int sigma = 0; sigma < 4; ++sigma) {
            if (sigma == mu or sigma == nu or sigma == rho) {
              continue;
            }
            const Complex c =
                -fac * (Complex)epsilon_tensor(mu, nu, rho, sigma) * v;
            v1[mu * 4 + rho] += c * v2[nu * 4 + sigma];
          }
        }
      }
    }
  }
  scalar_inversion(p1, pinv1, m_vector);
  scalar_divergence(photon1, p1);
}

inline void set_photon_pion_photon_vertex_lmd(FieldM<Complex, 4>& photon1,
                                              const FieldM<Complex, 1>& pion,
                                              const FieldM<Complex, 4>& photon2,
                                              const double m_vector,
                                              const double f_pi)
{
  TIMER_VERBOSE("set_photon_pion_photon_vertex_lmd");
  FieldM<Complex, 4> photon1_tmp;
  set_photon_pion_photon_vertex_vmd(photon1, pion, photon2, m_vector, f_pi);
  set_photon_pion_photon_vertex_two_end(photon1_tmp, pion, photon2, m_vector,
                                        f_pi);
  const double coef = 8.0 / 3.0 * sqr(PI * f_pi / m_vector);
  photon1 *= 1.0 - coef;
  photon1_tmp *= coef;
  photon1 += photon1_tmp;
}

inline void set_photon_pion_photon_vertex(FieldM<Complex, 4>& photon1,
                                          const FieldM<Complex, 1>& pion,
                                          const FieldM<Complex, 4>& photon2,
                                          const double m_vector,
                                          const double f_pi)
// vertex function
{
  // SADJUST ME
  // set_photon_pion_photon_vertex_two_end(photon1, pion, photon2, m_vector,
  // f_pi); set_photon_pion_photon_vertex_vmd(photon1, pion, photon2, m_vector,
  // f_pi);
  set_photon_pion_photon_vertex_lmd(photon1, pion, photon2, m_vector, f_pi);
}

inline void set_pion_gg_decay(FieldM<Complex, 4 * 4>& decay,
                              const Coordinate& total_site, const int t_sep,
                              const double m_pi, const double f_pi,
                              const double m_vector)
{
  TIMER_VERBOSE("set_pion_gg_decay");
  Geometry geo;
  geo.init(total_site, 1);
  //
  decay.init();
  decay.init(geo);
  //
  FieldM<Complex, 1> pion;
  pion.init(geo);
  set_zero(pion);
  const int t_src = mod(-t_sep, total_site[3]);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    if (xg[3] == t_src) {
      pion.get_elem(xl) = 1.0;
    }
  }
  scalar_inversion(pion, pion, m_pi);
  //
  pion *= exp(m_pi * t_sep);
  //
  for (int nu = 0; nu < 4; ++nu) {
    FieldM<Complex, 4> photon1, photon2;
    photon1.init(geo);
    photon2.init(geo);
    set_zero(photon1);
    set_zero(photon2);
    const Coordinate xl_0 = geo.coordinate_l_from_g(Coordinate());
    if (geo.is_local(xl_0)) {
      photon2.get_elem(xl_0, nu) = 1.0;
    }
    set_photon_pion_photon_vertex(photon1, pion, photon2, m_vector, f_pi);
#pragma omp parallel for
    for (long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      Vector<Complex> v = decay.get_elems(xl);
      const Vector<Complex> p1v = photon1.get_elems(xl);
      for (int mu = 0; mu < 4; ++mu) {
        v[mu * 4 + nu] = p1v[mu];
      }
    }
  }
}

inline std::vector<Complex> get_pion_corr(const Coordinate& total_site,
                                          const double m_pi)
{
  TIMER_VERBOSE("get_pion_corr");
  Geometry geo;
  geo.init(total_site, 1);
  //
  FieldM<Complex, 1> pion;
  pion.init(geo);
  set_zero(pion);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    if (xg[3] == 0) {
      pion.get_elem(xl) = 1.0;
    }
  }
  scalar_inversion(pion, pion, m_pi);
  //
  std::vector<Complex> corr(total_site[3], 0.0);
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    corr[xg[3]] += pion.get_elem(xl);
  }
  glb_sum(get_data(corr));
  return corr;
}

inline void save_pion_gg(const std::string& path, const Coordinate& total_site,
                         const int t_sep, const double ainv_gev,
                         const double m_pi_gev, const double f_pi_gev,
                         const double m_vector_gev)
{
  if (does_file_exist_sync_node(path + "/decay.field")) {
    return;
  }
  TIMER_VERBOSE("save_pion_gg");
  qmkdir_info(path);
  {
    std::ostringstream out;
    out << ssprintf("total_site = %s\n", show(total_site).c_str());
    out << ssprintf("t_sep = %d\n", t_sep);
    out << ssprintf("ainv_gev = %24.17E\n", ainv_gev);
    out << ssprintf("m_pi_gev = %24.17E\n", m_pi_gev);
    out << ssprintf("f_pi_gev = %24.17E\n", f_pi_gev);
    out << ssprintf("m_vector_gev = %24.17E\n", m_vector_gev);
    qtouch_info(path + "/info.txt", out.str());
  }
  {
    const std::vector<Complex> corr =
        get_pion_corr(total_site, m_pi_gev / ainv_gev);
    std::ostringstream out;
    out << "# t C[t].real() C[t].imag()\n";
    for (int t = 0; t < total_site[3]; ++t) {
      out << ssprintf("%4d %24.17E %24.17E\n", t, corr[t].real(),
                      corr[t].imag());
    }
    qtouch_info(path + "/pion-corr.txt", out.str());
  }
  {
    FieldM<Complex, 4 * 4> decay;
    set_pion_gg_decay(decay, total_site, t_sep, m_pi_gev / ainv_gev,
                      f_pi_gev / ainv_gev, m_vector_gev / ainv_gev);
    write_field_double(decay, path + "/decay.field");
  }
}

inline void compute_all()
{
  TIMER_VERBOSE("compute_all");
  qmkdir_info("results");
  save_pion_gg("results/physical-24nt96-1.0", Coordinate(24, 24, 24, 96), 24, 1.0, 0.1349766, 0.092, 0.77);
  // save_pion_gg("results/physical-32nt128-1.0", Coordinate(32, 32, 32, 128), 32, 1.0, 0.1349766, 0.092, 0.77);
  // save_pion_gg("results/physical-32nt128-1.3333", Coordinate(32, 32, 32, 128), 32, 1.3333, 0.1349766, 0.092, 0.77);
  // save_pion_gg("results/physical-48nt192-1.0", Coordinate(48, 48, 48, 192), 48, 1.0, 0.1349766, 0.092, 0.77);
  // save_pion_gg("results/physical-48nt192-2.0", Coordinate(48, 48, 48, 192), 48, 2.0, 0.1349766, 0.092, 0.77);
  // save_pion_gg("results/heavy-24nt96-1.0", Coordinate(24, 24, 24, 96), 24, 1.0, 0.340, 0.105, 0.83);
  // save_pion_gg("results/heavy-32nt128-1.0", Coordinate(32, 32, 32, 128), 32, 1.0, 0.340, 0.105, 0.83);
  // save_pion_gg("results/heavy-32nt128-1.3333", Coordinate(32, 32, 32, 128), 32, 1.3333, 0.340, 0.105, 0.83);
  // save_pion_gg("results/heavy-48nt192-1.0", Coordinate(48, 48, 48, 192), 48, 1.0, 0.340, 0.105, 0.83);
  // save_pion_gg("results/heavy-48nt192-2.0", Coordinate(48, 48, 48, 192), 48, 2.0, 0.340, 0.105, 0.83);
}


}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  begin(&argc, &argv);
  compute_all();
  Timer::display();
  end();
  return 0;
}
