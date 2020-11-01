#include "qlat-setup.h"

namespace qlat
{  //


inline void reflect_and_revert_mu_nu(FieldM<Complex, 8 * 8>& f_munu)
{
  TIMER_VERBOSE("reflect_and_revert_mu_nu");
  reflect_field(f_munu);
  FieldM<Complex, 8 * 8> f_munu_tmp;
  f_munu_tmp = f_munu;
  const Geometry& geo = f_munu.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Vector<Complex> fv0 = f_munu_tmp.get_elems(index);
    Vector<Complex> fv = f_munu.get_elems(index);
    for (int mu = 0; mu < 8; ++mu) {
      for (int nu = 0; nu < 8; ++nu) {
        fv[mu * 8 + nu] = fv0[nu * 8 + mu];
      }
    }
  }
}

inline void test()
{
  TIMER_VERBOSE("test");
  FieldM<Complex, 8 * 8> meson_vv, meson_vv_old, meson_vv_diff;
  std::vector<int> trajs;
  trajs.push_back(1900);
  trajs.push_back(2260);
  trajs.push_back(2270);
  trajs.push_back(2280);
  trajs.push_back(2290);
  trajs.push_back(2300);
  trajs.push_back(2310);
  for (long i = 0; i < (long)trajs.size(); ++i) {
    const int traj = trajs[i];
    FieldM<Complex, 8 * 8> meson_vv_tmp;
    const std::string path = ssprintf(
        "analysis/field-meson-vv/24D/results=%d/decay-0-0-0.field", traj);
    read_field_double_from_float(meson_vv_tmp, path);
    meson_vv_tmp *= 1.0 / (double)trajs.size();
    meson_vv += meson_vv_tmp;
  }
  std::vector<int> trajs_old;
  trajs_old.push_back(1010);
  trajs_old.push_back(1030);
  for (long i = 0; i < (long)trajs_old.size(); ++i) {
    const int traj = trajs_old[i];
    FieldM<Complex, 8 * 8> meson_vv_tmp;
    const std::string path = ssprintf(
        "/sdcc/u/jluchang/qcdqedta/luchang/all-analysis-data/field-pion-gg/"
        "24D-0.00107/results=%d/decay_type_1.field",
        traj);
    read_field_double(meson_vv_tmp, path);
    meson_vv_tmp *= 0.5 / (double)trajs_old.size();
    reflect_and_revert_mu_nu(meson_vv_tmp);
    meson_vv_old += meson_vv_tmp;
  }
  for (long i = 0; i < (long)trajs_old.size(); ++i) {
    const int traj = trajs_old[i];
    FieldM<Complex, 8 * 8> meson_vv_tmp;
    const std::string path = ssprintf(
        "/sdcc/u/jluchang/qcdqedta/luchang/all-analysis-data/field-pion-gg/"
        "24D-0.00107/results=%d/decay_type_2.field",
        traj);
    read_field_double(meson_vv_tmp, path);
    meson_vv_tmp *= 0.5 / (double)trajs_old.size();
    meson_vv_old += meson_vv_tmp;
  }
  const Geometry& geo = meson_vv.geo;
  displayln_info(show(geo));
  const double pion_mass = 0.139;
  meson_vv_old *= std::exp(pion_mass * 2);
  displayln_info(ssprintf("meson_vv qnorm = %24.17E", qnorm(meson_vv)));
  displayln_info(ssprintf("meson_vv_old qnorm = %24.17E", qnorm(meson_vv_old)));
  meson_vv_diff = meson_vv;
  meson_vv_diff -= meson_vv_old;
  displayln_info(
      ssprintf("meson_vv_diff qnorm = %24.17E", qnorm(meson_vv_diff)));
  displayln_info(
      ssprintf("meson_vv_diff_ratio sqrt(qnorm) = %24.17E",
               std::sqrt(qnorm(meson_vv_diff) / qnorm(meson_vv_old))));
}

}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 1));
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 4));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 2, 16));
  size_node_list.push_back(Coordinate(1, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 4, 16));
  size_node_list.push_back(Coordinate(2, 4, 4, 16));
  size_node_list.push_back(Coordinate(4, 4, 4, 16));
  begin(&argc, &argv, size_node_list);
  setup();
  //
  test();
  Timer::display();
  //
  end();
  return 0;
}
