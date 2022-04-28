#include <qlat/qlat.h>

namespace qlat
{
const double alpha_qed = 1.0 / 137.035999074;  // e_charge^2 / (4 * PI)

struct LFuncIntegrand {
  double mass;
  double ts, x;
  //
  double operator()(const double p) const
  {
    if (p == 0.0) {
      return 1.0;
    } else {
      const double ep = std::sqrt(sqr(mass) + sqr(p));
      if (x == 0.0) {
        return p / (p + ep - mass) * std::exp(-p * ts);
      } else {
        return std::sin(p * x) / ((p + ep - mass) * x) * std::exp(-p * ts);
      }
    }
  }
};

double l_func(const double ts, const double x_sq, const double mass)
// actually equals $ e_charge^2 \times L $
{
  const double coef = alpha_qed / PI;
  LFuncIntegrand f;
  f.mass = mass;
  f.ts = ts;
  f.x = sqrt(x_sq);
  return coef * adaptive_simpson(f, 0.0, 1.0 / 0.0, 1e-10);
}

void set_h_field(FieldM<Complex, 1>& h, const Coordinate& total_site,
                 const double mass)
{
  TIMER_VERBOSE("set_h_field");
  Geometry geo;
  geo.init(total_site, 1);
  h.init(geo);
  set_zero(h);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Coordinate krel = smod(xg, total_site);
    const double ts = xg[3];
    const double p_vec_sq = sqr(2.0 * PI / (double)total_site[0] * krel[0]) +
                            sqr(2.0 * PI / (double)total_site[1] * krel[1]) +
                            sqr(2.0 * PI / (double)total_site[2] * krel[2]);
    const double ep = std::sqrt(p_vec_sq + sqr(mass));
    h.get_elem(xl) = (mass + ep) / (2.0 * ep) * std::exp(-(ep - mass) * ts);
  }
  h *= (double)total_site[3] / (double)geo.total_volume();
  fft_complex_field_spatial(h, false);
}

const std::vector<std::vector<double> >& get_l_func_table(
    const double mass, const int ts_limit, const long x_vec_sq_limit)
{
  static Cache<std::string, std::vector<std::vector<double> > > cache(
      "l_func_tables", 16);
  std::vector<std::vector<double> >& l_func_table =
      cache[ssprintf("%24.17E %d %ld", mass, ts_limit, x_vec_sq_limit)];
  if (l_func_table.size() == 0) {
    TIMER_VERBOSE("acc_h_field-l_func_table");
    const int num_node = get_num_node();
    const int id_node = get_id_node();
    l_func_table.resize(ts_limit);
    for (int ts = 0; ts < ts_limit; ++ts) {
      l_func_table[ts].resize(x_vec_sq_limit, 0.0);
#pragma omp parallel for
      for (long x_vec_sq = 0; x_vec_sq < x_vec_sq_limit; ++x_vec_sq) {
        if (x_vec_sq % num_node == id_node) {
          l_func_table[ts][x_vec_sq] = l_func(ts, x_vec_sq, mass);
        } else {
          l_func_table[ts][x_vec_sq] = 0.0;
        }
      }
      glb_sum(get_data(l_func_table[ts]));
    }
  }
  return l_func_table;
}

void acc_h_field(const FieldM<Complex, 1>& h, const double mass,
                 const int ts_limit, const long x_vec_sq_limit)
{
  TIMER_VERBOSE("acc_h_field");
  const Geometry& geo = h.geo();
  const Coordinate total_site = geo.total_site();
  qassert(ts_limit >= total_site[3]);
  qassert(x_vec_sq_limit >= sqr(total_site[0] / 2) + sqr(total_site[1] / 2) +
                                sqr(total_site[2] / 2) + 1);
  const std::vector<std::vector<double> >& l_func_table = get_l_func_table(mass, ts_limit, x_vec_sq_limit);
  std::vector<double> int_short(total_site[3], 0.0);
  std::vector<double> bound_long(total_site[3], 0.0);
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Coordinate xrel = smod(xg, total_site);
    const int ts = xg[3];
    const long x_vec_sq = sqr(xrel[0]) + sqr(xrel[1]) + sqr(xrel[2]);
    const double x_sq = x_vec_sq + sqr(ts);
    const double photon_prop =
        alpha_qed / PI * (x_sq == 0.0 ? 2.76963 : 1.0 / x_sq);
    const Complex& h_val = h.get_elem(xl);
    qassert(x_vec_sq < x_vec_sq_limit);
    int_short[ts] += photon_prop * h_val.real();
    bound_long[ts] += l_func_table[ts][x_vec_sq] * h_val.real();
  }
  glb_sum(get_data(int_short));
  glb_sum(get_data(bound_long));
  std::vector<double> acc_short(total_site[3], 0.0);
  std::vector<double> total(total_site[3], 0.0);
  for (int t = 0; t < total_site[3]; ++t) {
    if (t == 0) {
      acc_short[t] = 0.5 * int_short[t];
      total[t] = 0.5 * int_short[t];
    } else {
      acc_short[t] =
          acc_short[t - 1] + 0.5 * int_short[t - 1] + 0.5 * int_short[t];
      total[t] = acc_short[t] + bound_long[t];
    }
    displayln_info(ssprintf("RESULT: %10.5f %3d %3d %3d %5d %24.17E %24.17E", mass,
                            total_site[0], total_site[1], total_site[2], t,
                            acc_short[t], total[t]));
  }
}

void compute(const int l_size, const double mass, const int ts_limit, const long x_vec_sq_limit)
{
  {
    TIMER_VERBOSE("compute");
    const Coordinate total_site(l_size, l_size, l_size, ts_limit);
    FieldM<Complex, 1> h;
    set_h_field(h, total_site, mass);
    acc_h_field(h, mass, ts_limit, x_vec_sq_limit);
  }
  Timer::display();
}

}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  begin(&argc, &argv);
  switch_monitor_file_info("log.txt");
  displayln_info(show(l_func(1.0, 1.0, 0.14)));
  displayln_info(show(l_func(2.0, 1.0, 0.14)));
  displayln_info(show(l_func(2.0, 2.0, 0.14)));
  const int ts_limit = 200;
  const long x_vec_sq_limit = 1 + 3 * sqr(48);
  for (int i = 96; i >= 2; i -= 2) {
    compute(i, 0.14, ts_limit, x_vec_sq_limit);
    compute(i, 0.07, ts_limit, x_vec_sq_limit);
    compute(i, 0.035, ts_limit, x_vec_sq_limit);
  }
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
