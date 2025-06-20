#include "qlat-setup.h"

namespace qlat
{  //

inline void reflect_and_revert_mu_nu(FieldM<ComplexD, 8 * 8>& f_munu)
{
  TIMER_VERBOSE("reflect_and_revert_mu_nu");
  reflect_field(f_munu);
  field_permute_mu_nu(f_munu);
}

inline void set_pfdist(FieldM<ComplexD, 1>& pfdist,
                       const std::vector<int>& traj_list)
{
  TIMER_VERBOSE("set_pfdist");
  pfdist.init();
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 1> pfdist_tmp;
    const std::string path = ssprintf(
        "analysis/field-psel-fsel-distribution/24D/results=%d/avg-0.field",
        traj);
    read_field_double(pfdist_tmp, path);
    pfdist_tmp *= 1.0 / (double)traj_list.size();
    pfdist += pfdist_tmp;
  }
}

inline void set_meson_vv(FieldM<ComplexD, 8 * 8>& meson_vv,
                         const std::vector<int>& traj_list, const std::string& tag)
{
  TIMER_VERBOSE("set_meson_vv");
  meson_vv.init();
  FieldM<ComplexD, 1> pfdist;
  set_pfdist(pfdist, traj_list);
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> meson_vv_tmp;
    const std::string path =
        ssprintf("analysis/field-meson-vv/24D/results=%d/%s-0-0-0.field", traj,
                 tag.c_str());
    read_field_double_from_float(meson_vv_tmp, path);
    meson_vv_tmp *= 1.0 / (double)traj_list.size();
    meson_vv += meson_vv_tmp;
  }
  rescale_field_with_psel_fsel_distribution(meson_vv, pfdist);
}

inline void set_meson_vv_meson(FieldM<ComplexD, 8 * 8>& meson_vv_meson,
                               const std::vector<int>& traj_list,
                               const std::string& tag)
{
  TIMER_VERBOSE("set_meson_vv_meson");
  meson_vv_meson.init();
  FieldM<ComplexD, 1> pfdist;
  set_pfdist(pfdist, traj_list);
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> tmp;
    const std::string path = ssprintf(
        "analysis/field-meson-vv-meson/24D/results=%d/%s-0-0-0-0.field", traj,
        tag.c_str());
    read_field_double_from_float(tmp, path);
    tmp *= 1.0 / (double)traj_list.size();
    meson_vv_meson += tmp;
  }
  rescale_field_with_psel_fsel_distribution(meson_vv_meson, pfdist);
}

inline void set_meson_vv_old(FieldM<ComplexD, 8 * 8>& meson_vv,
                             const std::vector<int>& traj_list,
                             const std::string& tag)
{
  TIMER_VERBOSE("set_meson_vv_old");
  meson_vv.init();
  FieldM<ComplexD, 1> pfdist;
  set_pfdist(pfdist, traj_list);
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> meson_vv_tmp;
    const std::string path = ssprintf(
        "/sdcc/u/jluchang/qcdqedta/luchang/all-analysis-data/field-pion-gg/"
        "24D-0.00107/results=%d/%s_type_1.field",
        traj, tag.c_str());
    read_field_double(meson_vv_tmp, path);
    meson_vv_tmp *= 0.5 / (double)traj_list.size();
    reflect_and_revert_mu_nu(meson_vv_tmp);
    meson_vv += meson_vv_tmp;
  }
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> meson_vv_tmp;
    const std::string path = ssprintf(
        "/sdcc/u/jluchang/qcdqedta/luchang/all-analysis-data/field-pion-gg/"
        "24D-0.00107/results=%d/%s_type_2.field",
        traj, tag.c_str());
    read_field_double(meson_vv_tmp, path);
    meson_vv_tmp *= 0.5 / (double)traj_list.size();
    meson_vv += meson_vv_tmp;
  }
  rescale_field_with_psel_fsel_distribution(meson_vv, pfdist);
}

inline void set_meson_vv_meson_old(FieldM<ComplexD, 8 * 8>& meson_vv_meson,
                                   const std::vector<int>& traj_list,
                                   const std::string& tag)
{
  TIMER_VERBOSE("set_meson_vv_meson_old");
  meson_vv_meson.init();
  FieldM<ComplexD, 1> pfdist;
  set_pfdist(pfdist, traj_list);
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> tmp;
    const std::string path = ssprintf(
        "/sdcc/u/jluchang/qcdqedta/luchang/all-analysis-data/field-ppgg/"
        "24D-0.00107/results=%d/%s_type_1.field",
        traj, tag.c_str());
    read_field_double(tmp, path);
    tmp *= 0.5 / (double)traj_list.size();
    reflect_and_revert_mu_nu(tmp);
    meson_vv_meson += tmp;
  }
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> tmp;
    const std::string path = ssprintf(
        "/sdcc/u/jluchang/qcdqedta/luchang/all-analysis-data/field-ppgg/"
        "24D-0.00107/results=%d/%s_type_2.field",
        traj, tag.c_str());
    read_field_double(tmp, path);
    tmp *= 0.5 / (double)traj_list.size();
    meson_vv_meson += tmp;
  }
  rescale_field_with_psel_fsel_distribution(meson_vv_meson, pfdist);
}

template <class M>
void set_field_range(Field<M>& f, const Long dis_sq_range)
{
  TIMER_VERBOSE("set_field_range");
  const Geometry& geo = f.geo;
  const Coordinate total_site = geo.total_site();
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Coordinate xgrel = smod(xg, total_site);
    const Long dis_sq = sqr(xgrel);
    if (dis_sq > dis_sq_range) {
      set_zero(f.get_elems(xl));
    }
  }
}

inline void test_meson_vv()
{
  TIMER_VERBOSE("test_meson_vv");
  const std::string tag = "decay";
  // const std::string tag = "fission";
  FieldM<ComplexD, 8 * 8> meson_vv, meson_vv_old, meson_vv_diff;
  std::vector<int> traj_list;
  traj_list.push_back(1010);
  traj_list.push_back(1030);
  // traj_list.push_back(1900);
  // traj_list.push_back(2260);
  // traj_list.push_back(2270);
  // traj_list.push_back(2280);
  // traj_list.push_back(2290);
  // traj_list.push_back(2300);
  // traj_list.push_back(2310);
  // traj_list.push_back(2320);
  // traj_list.push_back(2330);
  // traj_list.push_back(2340);
  // traj_list.push_back(2350);
  set_meson_vv(meson_vv, traj_list, tag);
  std::vector<int> traj_list_old;
  traj_list_old.push_back(1010);
  traj_list_old.push_back(1030);
  set_meson_vv_old(meson_vv_old, traj_list_old, tag);
  const Long dis_sq_range = sqr(3);
  set_field_range(meson_vv, dis_sq_range);
  set_field_range(meson_vv_old, dis_sq_range);
  const Geometry& geo = meson_vv.geo;
  displayln_info(show(geo));
  const double pion_mass = 0.139;
  meson_vv_old *= std::exp(pion_mass * 2);
  displayln_info(ssprintf("meson_vv qnorm = %24.17E", qnorm(meson_vv)));
  displayln_info(ssprintf("meson_vv proj to meson_vv = %24.17E",
                          qnorm_double(meson_vv, meson_vv) / qnorm(meson_vv)));
  displayln_info(ssprintf("meson_vv_old qnorm = %24.17E", qnorm(meson_vv_old)));
  displayln_info(
      ssprintf("meson_vv_old proj to meson_vv = %24.17E",
               qnorm_double(meson_vv, meson_vv_old) / qnorm(meson_vv)));
  meson_vv_diff = meson_vv;
  meson_vv_diff -= meson_vv_old;
  displayln_info(
      ssprintf("meson_vv_diff qnorm = %24.17E", qnorm(meson_vv_diff)));
  displayln_info(
      ssprintf("meson_vv_diff_ratio sqrt(qnorm) = %24.17E",
               std::sqrt(qnorm(meson_vv_diff) / qnorm(meson_vv_old))));
}

inline void test_meson_vv_meson()
{
  TIMER_VERBOSE("test_meson_vv_meson");
  // const std::string tag = "forward";
  const std::string tag = "backward";
  FieldM<ComplexD, 8 * 8> meson_vv_meson, meson_vv_meson_old,
      meson_vv_meson_diff;
  std::vector<int> traj_list;
  traj_list.push_back(1010);
  // traj_list.push_back(1030);
  // traj_list.push_back(1900);
  // traj_list.push_back(2260);
  // traj_list.push_back(2270);
  // traj_list.push_back(2280);
  // traj_list.push_back(2290);
  // traj_list.push_back(2300);
  // traj_list.push_back(2310);
  // traj_list.push_back(2320);
  // traj_list.push_back(2330);
  // traj_list.push_back(2340);
  // traj_list.push_back(2350);
  set_meson_vv_meson(meson_vv_meson, traj_list, tag);
  std::vector<int> traj_list_old;
  traj_list_old.push_back(1010);
  // traj_list_old.push_back(1030);
  set_meson_vv_meson_old(meson_vv_meson_old, traj_list_old, tag);
  const Long dis_sq_range = sqr(5);
  set_field_range(meson_vv_meson, dis_sq_range);
  set_field_range(meson_vv_meson_old, dis_sq_range);
  const Geometry& geo = meson_vv_meson.geo;
  displayln_info(show(geo));
  displayln_info(
      ssprintf("meson_vv_meson qnorm = %24.17E", qnorm(meson_vv_meson)));
  displayln_info(ssprintf(
      "meson_vv_meson proj to meson_vv = %24.17E",
      qnorm_double(meson_vv_meson, meson_vv_meson) / qnorm(meson_vv_meson)));
  displayln_info(ssprintf("meson_vv_meson_old qnorm = %24.17E",
                          qnorm(meson_vv_meson_old)));
  displayln_info(ssprintf("meson_vv_meson_old proj to meson_vv_meson = %24.17E",
                          qnorm_double(meson_vv_meson, meson_vv_meson_old) /
                              qnorm(meson_vv_meson)));
  meson_vv_meson_diff = meson_vv_meson;
  meson_vv_meson_diff -= meson_vv_meson_old;
  displayln_info(ssprintf("meson_vv_meson_diff qnorm = %24.17E",
                          qnorm(meson_vv_meson_diff)));
  displayln_info(ssprintf(
      "meson_vv_meson_diff_ratio sqrt(qnorm) = %24.17E",
      std::sqrt(qnorm(meson_vv_meson_diff) / qnorm(meson_vv_meson_old))));
}

inline void set_meson_vv_v2(FieldM<ComplexD, 8 * 8>& meson_vv,
                            const std::vector<int>& traj_list,
                            const std::string& tag)
{
  TIMER_VERBOSE("set_meson_vv_v2");
  meson_vv.init();
  FieldM<ComplexD, 1> pfdist;
  set_pfdist(pfdist, traj_list);
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> tmp;
    const std::string path =
        ssprintf("analysis-v2/field-meson-vv/24D/results=%d/%s-0.field",
                 traj, tag.c_str());
    read_field_double_from_float(tmp, path);
    tmp *= 1.0 / (double)traj_list.size();
    meson_vv += tmp;
  }
  rescale_field_with_psel_fsel_distribution(meson_vv, pfdist);
}

inline void ref_avg(FieldM<ComplexD, 8 * 8>& f)
{
  TIMER_VERBOSE("ref_avg");
  FieldM<ComplexD, 8 * 8> tmp;
  tmp = f;
  reflect_field(tmp);
  field_permute_mu_nu(tmp);
  field_conjugate_mu_nu(tmp);
  field_complex_conjugate(tmp);
  f += tmp;
  f *= 0.5;
}

inline void set_meson_vv_v3(FieldM<ComplexD, 8 * 8>& meson_vv,
                            const std::vector<int>& traj_list,
                            const std::string& tag1, const std::string& tag2)
{
  TIMER_VERBOSE("set_meson_vv_v3");
  meson_vv.init();
  FieldM<ComplexD, 1> pfdist;
  set_pfdist(pfdist, traj_list);
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> tmp;
    const std::string path =
        ssprintf("analysis-v3/field-meson-vv/24D/results=%d/%s-0.field",
                 traj, tag1.c_str());
    read_field_double_from_float(tmp, path);
    meson_vv += tmp;
  }
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> tmp;
    const std::string path =
        ssprintf("analysis-v3/field-meson-vv/24D/results=%d/%s-0.field", traj,
                 tag2.c_str());
    read_field_double_from_float(tmp, path);
    reflect_field(tmp);
    field_permute_mu_nu(tmp);
    field_conjugate_mu_nu(tmp);
    field_complex_conjugate(tmp);
    meson_vv += tmp;
  }
  meson_vv *= 0.5 / (double)traj_list.size();
  rescale_field_with_psel_fsel_distribution(meson_vv, pfdist);
}

inline void set_meson_vv_v4(FieldM<ComplexD, 8 * 8>& meson_vv,
                            const std::vector<int>& traj_list,
                            const std::string& tag)
{
  TIMER_VERBOSE("set_meson_vv_v4");
  meson_vv.init();
  FieldM<ComplexD, 1> pfdist;
  set_pfdist(pfdist, traj_list);
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> tmp;
    const std::string path =
        ssprintf("analysis/field-meson-vv/24D/results=%d/%s-0.field",
                 traj, tag.c_str());
    read_field_double_from_float(tmp, path);
    meson_vv += tmp;
  }
  meson_vv *= 1.0 / (double)traj_list.size();
  rescale_field_with_psel_fsel_distribution(meson_vv, pfdist);
}

inline void test_meson_vv_v4()
{
  TIMER_VERBOSE("test_meson_vv_v4");
  // const std::string tag = "decay";
  // const std::string tag = "fission";
  FieldM<ComplexD, 8 * 8> meson_vv, meson_vv_old, meson_vv_diff;
  std::vector<int> traj_list;
  for (Long traj = 1000; traj < 3000; traj += 10) {
    if (does_file_exist_sync_node(ssprintf(
            "analysis/field-meson-vv/24D/results=%d/decay-0-1-0.field",
            traj))) {
      traj_list.push_back(traj);
    }
    if (traj_list.size() >= 1) {
      break;
    }
  }
  set_meson_vv_v3(meson_vv_old, traj_list, "decay-0-1", "decay-1-0");
  set_meson_vv_v3(meson_vv, traj_list, "fission-0-1", "fission-1-0");
  reflect_field(meson_vv);
  meson_vv_old += meson_vv;
  meson_vv_old *= 0.5;
  set_meson_vv_v4(meson_vv, traj_list, "decay-0-1");
  const Long dis_sq_range = sqr(3);
  set_field_range(meson_vv, dis_sq_range);
  set_field_range(meson_vv_old, dis_sq_range);
  const Geometry& geo = meson_vv.geo;
  displayln_info(show(geo));
  displayln_info(ssprintf("meson_vv qnorm = %24.17E", qnorm(meson_vv)));
  displayln_info(ssprintf("meson_vv proj to meson_vv = %24.17E",
                          qnorm_double(meson_vv, meson_vv) / qnorm(meson_vv)));
  displayln_info(ssprintf("meson_vv_old qnorm = %24.17E", qnorm(meson_vv_old)));
  displayln_info(
      ssprintf("meson_vv_old proj to meson_vv = %24.17E",
               qnorm_double(meson_vv, meson_vv_old) / qnorm(meson_vv)));
  meson_vv_diff = meson_vv;
  meson_vv_diff -= meson_vv_old;
  displayln_info(
      ssprintf("meson_vv_diff qnorm = %24.17E", qnorm(meson_vv_diff)));
  displayln_info(
      ssprintf("meson_vv_diff_ratio sqrt(qnorm) = %24.17E",
               std::sqrt(qnorm(meson_vv_diff) / qnorm(meson_vv_old))));
}

inline void set_meson_vv_meson_v2(FieldM<ComplexD, 8 * 8>& meson_vv_meson,
                                  const std::vector<int>& traj_list,
                                  const std::string& tag)
{
  TIMER_VERBOSE("set_meson_vv_meson_v2");
  meson_vv_meson.init();
  FieldM<ComplexD, 1> pfdist;
  set_pfdist(pfdist, traj_list);
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> tmp;
    const std::string path = ssprintf(
        "analysis-v2/field-meson-vv-meson/24D/results=%d/%s-0-0.field",
        traj, tag.c_str());
    read_field_double_from_float(tmp, path);
    meson_vv_meson += tmp;
  }
  meson_vv_meson *= 1.0 / (double)traj_list.size();
  rescale_field_with_psel_fsel_distribution(meson_vv_meson, pfdist);
}

inline void set_meson_vv_meson_v2avg(FieldM<ComplexD, 8 * 8>& meson_vv_meson,
                                  const std::vector<int>& traj_list,
                                  const std::string& tag)
{
  TIMER_VERBOSE("set_meson_vv_meson_v2avg");
  meson_vv_meson.init();
  FieldM<ComplexD, 1> pfdist;
  set_pfdist(pfdist, traj_list);
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> tmp;
    const std::string path = ssprintf(
        "analysis-v2/field-meson-vv-meson-avg/24D/results=%d/%s-0-0.field",
        traj, tag.c_str());
    read_field_double_from_float(tmp, path);
    meson_vv_meson += tmp;
  }
  meson_vv_meson *= 1.0 / (double)traj_list.size();
  rescale_field_with_psel_fsel_distribution(meson_vv_meson, pfdist);
}

inline void set_meson_vv_meson_v3(FieldM<ComplexD, 8 * 8>& meson_vv_meson,
                                  const std::vector<int>& traj_list,
                                  const std::string& tag)
{
  TIMER_VERBOSE("set_meson_vv_meson_v3");
  meson_vv_meson.init();
  FieldM<ComplexD, 1> pfdist;
  set_pfdist(pfdist, traj_list);
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> tmp;
    const std::string path = ssprintf(
        "analysis-v3/field-meson-vv-meson/24D/results=%d/%s-0-0.field",
        traj, tag.c_str());
    read_field_double_from_float(tmp, path);
    meson_vv_meson += tmp;
  }
  meson_vv_meson *= 1.0 / (double)traj_list.size();
  rescale_field_with_psel_fsel_distribution(meson_vv_meson, pfdist);
}

inline void set_meson_vv_meson_v4(FieldM<ComplexD, 8 * 8>& meson_vv_meson,
                                  const std::vector<int>& traj_list,
                                  const std::string& tag)
{
  TIMER_VERBOSE("set_meson_vv_meson_v4");
  meson_vv_meson.init();
  FieldM<ComplexD, 1> pfdist;
  set_pfdist(pfdist, traj_list);
  for (Long i = 0; i < (Long)traj_list.size(); ++i) {
    const int traj = traj_list[i];
    FieldM<ComplexD, 8 * 8> tmp;
    const std::string path = ssprintf(
        "analysis/field-meson-vv-meson/24D/results=%d/%s-0-0.field",
        traj, tag.c_str());
    read_field_double_from_float(tmp, path);
    meson_vv_meson += tmp;
  }
  meson_vv_meson *= 1.0 / (double)traj_list.size();
  rescale_field_with_psel_fsel_distribution(meson_vv_meson, pfdist);
}

inline void test_meson_vv_meson_v4()
{
  TIMER_VERBOSE("test_meson_vv_meson_v4");
  // const std::string tag = "forward";
  // const std::string tag = "backward";
  FieldM<ComplexD, 8 * 8> meson_vv_meson, meson_vv_meson_old,
      meson_vv_meson_diff;
  std::vector<int> traj_list;
  for (Long traj = 1000; traj < 3000; traj += 10) {
    if (does_file_exist_sync_node(
            ssprintf("analysis/field-meson-vv-meson/24D/results=%d/"
                     "forward-0-0-0-0.field",
                     traj))) {
      traj_list.push_back(traj);
    }
    if (traj_list.size() >= 1) {
      break;
    }
  }
  // set_meson_vv_meson_v2(meson_vv_meson, traj_list, "forward-0-1");
  // set_meson_vv_meson_v3(meson_vv_meson, traj_list, "forward-0-1");
  // set_meson_vv_meson_v3(meson_vv_meson_old, traj_list, "backward-1-0");
  // reflect_field(meson_vv_meson_old);
  // field_permute_mu_nu(meson_vv_meson_old);
  // field_conjugate_mu_nu(meson_vv_meson_old);
  // field_complex_conjugate(meson_vv_meson_old);
  // meson_vv_meson += meson_vv_meson_old;
  // meson_vv_meson *= 0.5;
  //
  set_meson_vv_meson_v2avg(meson_vv_meson_old, traj_list, "forward-0-1");
  set_meson_vv_meson_v2avg(meson_vv_meson, traj_list, "backward-0-1");
  reflect_field(meson_vv_meson);
  meson_vv_meson_old += meson_vv_meson;
  meson_vv_meson_old *= 0.5;
  set_meson_vv_meson_v4(meson_vv_meson, traj_list, "forward-0-1");
  const Long dis_sq_range = sqr(5);
  set_field_range(meson_vv_meson, dis_sq_range);
  set_field_range(meson_vv_meson_old, dis_sq_range);
  const Geometry& geo = meson_vv_meson.geo;
  displayln_info(show(geo));
  displayln_info(
      ssprintf("meson_vv_meson qnorm = %24.17E", qnorm(meson_vv_meson)));
  displayln_info(ssprintf(
      "meson_vv_meson proj to meson_vv = %24.17E",
      qnorm_double(meson_vv_meson, meson_vv_meson) / qnorm(meson_vv_meson)));
  displayln_info(ssprintf("meson_vv_meson_old qnorm = %24.17E",
                          qnorm(meson_vv_meson_old)));
  displayln_info(ssprintf("meson_vv_meson_old proj to meson_vv_meson = %24.17E",
                          qnorm_double(meson_vv_meson, meson_vv_meson_old) /
                              qnorm(meson_vv_meson)));
  meson_vv_meson_diff = meson_vv_meson;
  meson_vv_meson_diff -= meson_vv_meson_old;
  displayln_info(ssprintf("meson_vv_meson_diff qnorm = %24.17E",
                          qnorm(meson_vv_meson_diff)));
  displayln_info(ssprintf(
      "meson_vv_meson_diff_ratio sqrt(qnorm) = %24.17E",
      std::sqrt(qnorm(meson_vv_meson_diff) / qnorm(meson_vv_meson_old))));
}

inline void test()
{
  TIMER_VERBOSE("test");
  // test_meson_vv();
  // test_meson_vv_meson();
  test_meson_vv_v4();
  test_meson_vv_meson_v4();
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
