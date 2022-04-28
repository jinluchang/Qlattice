#include <qlat/qlat.h>

#include <complex>
#include <iostream>
#include <vector>

using namespace qlat;
using namespace std;

inline std::vector<Coordinate> get_new_size_node_list()
{
  std::vector<Coordinate> new_size_nodes;
  new_size_nodes.push_back(Coordinate(1, 1, 1, 1));
  new_size_nodes.push_back(Coordinate(2, 2, 2, 2));
  // new_size_nodes.push_back(Coordinate(1,2,1,32));
  new_size_nodes.push_back(Coordinate(1, 1, 1, 8));
  new_size_nodes.push_back(Coordinate(1, 2, 4, 2));
  new_size_nodes.push_back(Coordinate(2, 2, 1, 1));
  // new_size_nodes.push_back(Coordinate(16,1,1,1));
  // new_size_nodes.push_back(Coordinate(4,8,1,1));
  new_size_nodes.push_back(Coordinate(4, 4, 4, 4));
  // new_size_nodes.push_back(Coordinate(16,16,1,1));
  // new_size_nodes.push_back(Coordinate(4,2,16,32));
  // new_size_nodes.push_back(Coordinate(4,16,16,32));
  // new_size_nodes.push_back(Coordinate(16,2,1,32));
  // new_size_nodes.push_back(Coordinate(16,16,8,32));
  return new_size_nodes;
}

inline std::vector<Coordinate> get_shift_list()
{
  std::vector<Coordinate> shift_list;
  shift_list.push_back(Coordinate(12, 34, 45, 12));
  shift_list.push_back(Coordinate(12, -34, -45, 12));
  shift_list.push_back(Coordinate(0, 0, 0, 0));
  shift_list.push_back(Coordinate(1, 1, 1, 1));
  return shift_list;
}

inline void test_shuffle()
{
  TIMER("test_shuffle");
  // Coordinate total_site(4, 4, 4, 8);
  Coordinate total_site(8, 8, 8, 8);
  // Coordinate total_site(16, 16, 16, 32);
  // Coordinate total_site(32, 32, 32, 64);
  // Coordinate total_site(48, 48, 48, 96);
  // Coordinate total_site(64, 64, 64, 64);
  Geometry geo;
  geo.init(total_site, 1);
  GaugeField gf;
  gf.init(geo);
  set_unit(gf);
  set_g_rand_color_matrix_field(gf, RngState("rgf-0.2"), 0.2);
  const crc32_t crc_0 = field_crc32(gf);
  displayln_info(ssprintf("crc32 = %08X", crc_0));
  const std::vector<Coordinate> new_size_nodes = get_new_size_node_list();
  for (size_t i = 0; i < new_size_nodes.size(); ++i) {
    const Coordinate& new_size_node = new_size_nodes[i];
    std::vector<Field<ColorMatrix> > gfs;
    shuffle_field(gfs, gf, new_size_node);
    set_unit(gf);
    shuffle_field_back(gf, gfs, new_size_node);
    const crc32_t crc_1 = field_crc32(gf);
    displayln_info(ssprintf(
        "crc32 = %08X with shuffle_field and shuffle_field_back", crc_1));
    qassert(crc_1 == crc_0);
  }
  const std::vector<Coordinate> shift_list = get_shift_list();
  for (size_t i = 0; i < shift_list.size(); ++i) {
    const Coordinate& shift = shift_list[i];
    GaugeField gf_shift;
    field_shift(gf_shift, gf, shift);
    const crc32_t crc_2 = field_crc32(gf_shift);
    displayln_info(
        ssprintf("crc32 = %08X after shift %s", crc_2, show(shift).c_str()));
    field_shift(gf_shift, gf_shift, -shift);
    const crc32_t crc_3 = field_crc32(gf_shift);
    displayln_info(
        ssprintf("crc32 = %08X after shift %s", crc_3, show(-shift).c_str()));
    qassert(crc_3 == crc_0);
  }
  {
    GaugeField gf_fft;
    gf_fft = gf;
    fft_complex_field(gf_fft, true);
    const crc32_t crc_4 = field_crc32(gf_fft);
    displayln_info(ssprintf("crc32 = %08X after fft forward", crc_4));
    displayln_info(ssprintf("qnorm = %24.17E after fft forward", qnorm(gf_fft)));
    fft_complex_field(gf_fft, false);
    const crc32_t crc_5 = field_crc32(gf_fft);
    displayln_info(ssprintf("crc32 = %08X after fft backward", crc_5));
    displayln_info(ssprintf("qnorm = %24.17E after fft backward", qnorm(gf_fft)));
    gf_fft *= 1.0 / (double)gf_fft.geo().total_volume();
    const crc32_t crc_6 = field_crc32(gf_fft);
    displayln_info(ssprintf("crc32 = %08X after rescale", crc_6));
    displayln_info(ssprintf("qnorm = %24.17E after rescale", qnorm(gf_fft)));
    gf_fft -= gf;
    const crc32_t crc_7 = field_crc32(gf_fft);
    displayln_info(ssprintf("crc32 = %08X after subtract", crc_7));
    displayln_info(ssprintf("qnorm = %24.17E after subtract", qnorm(gf_fft)));
    qassert(qnorm(gf_fft) < 1.0e-12);
  }
  {
    reflect_field(gf);
    const crc32_t crc_8 = field_crc32(gf);
    displayln_info(ssprintf("crc32 = %08X after reflect", crc_8));
    reflect_field(gf);
    const crc32_t crc_9 = field_crc32(gf);
    displayln_info(ssprintf("crc32 = %08X after reflect twice", crc_9));
    qassert(crc_9 == crc_0);
  }
}

inline void test_io()
{
  TIMER("test_io");
  qmkdir_sync_node("huge-data");
  RngState rs(get_global_rng_state(), fname);
  if (get_id_node() == 0) {
    crc32_check();
  }
  dist_write_par_limit() = 8;
  dist_read_par_limit() = 8;
  // Coordinate total_site(4, 4, 4, 8);
  Coordinate total_site(8, 8, 8, 8);
  // Coordinate total_site(16, 16, 16, 32);
  // Coordinate total_site(32, 32, 32, 64);
  // Coordinate total_site(48, 48, 48, 96);
  // Coordinate total_site(64, 64, 64, 64);
  Geometry geo;
  geo.init(total_site, 1);
  GaugeField gf;
  gf.init(geo);
  set_unit(gf);
  dist_write_field(gf, "huge-data/gauge_field_unit");
  set_unit(gf);
  set_g_rand_color_matrix_field(gf, RngState(rs, "rgf-0.1"), 0.1);
  dist_write_field(gf, "huge-data/gauge_field_rgf-0.1");
  set_unit(gf);
  dist_write_field(gf, "huge-data/gauge_field_unit_2");
  dist_read_field(gf, "huge-data/gauge_field_rgf-0.1");
  dist_write_field(gf, "huge-data/gauge_field_rgf-0.1_2");
  gf_show_info(gf);
  dist_write_field_float_from_double(gf, "huge-data/gauge_field_rgf-0.1_f");
  dist_read_field_double_from_float(gf, "huge-data/gauge_field_rgf-0.1_f");
  dist_write_field_float_from_double(gf, "huge-data/gauge_field_rgf-0.1_f_2");
  gf_show_info(gf);
  set_g_rand_color_matrix_field(gf, RngState(rs, "rgf-0.1"), 0.1);
  gf_show_info(gf);
  get_shuffle_plan_cache().limit = 16;
  displayln_info(ssprintf("dist_crc32 = %08X", field_dist_crc32(gf)));
  displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
  const crc32_t crc = field_dist_crc32(gf);
  std::vector<Coordinate> new_size_nodes = get_new_size_node_list();
  for (size_t i = 0; i < new_size_nodes.size(); ++i) {
    const Coordinate& new_size_node = new_size_nodes[i];
    std::vector<Field<ColorMatrix> > gfs;
    shuffle_field(gfs, gf, new_size_node);
    for (size_t k = 0; k < gfs.size(); ++k) {
      Field<ColorMatrix>& gfk = gfs[k];
      set_unit(gfk);
      set_g_rand_color_matrix_field(gfk, RngState(rs, "rgf-0.1"), 0.1);
    }
    set_unit(gf);
    shuffle_field_back(gf, gfs, new_size_node);
    displayln_info(ssprintf("dist_crc32 = %08X", field_dist_crc32(gf)));
    displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
    qassert(crc == field_dist_crc32(gf));
  }
  for (size_t i = 0; i < new_size_nodes.size(); ++i) {
    const Coordinate& new_size_node = new_size_nodes[i];
    displayln_info(ssprintf("dist_crc32 = %08X", field_dist_crc32(gf)));
    displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
    dist_write_field(
        gf, new_size_node,
        ssprintf("huge-data/gauge_field ; ") + show(new_size_node));
    set_unit(gf);
    dist_read_field(gf,
                    ssprintf("huge-data/gauge_field ; ") + show(new_size_node));
    displayln_info(ssprintf("dist_crc32 = %08X", field_dist_crc32(gf)));
    displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
    qassert(crc == field_dist_crc32(gf));
  }
  for (size_t i = 0; i < new_size_nodes.size(); ++i) {
    const Coordinate& new_size_node = new_size_nodes[i];
    displayln_info(ssprintf("dist_crc32 = %08X", field_dist_crc32(gf)));
    displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
    dist_write_field(
        gf, new_size_node,
        ssprintf("huge-data/gauge_field ; ") + show(new_size_node) + "_2");
    set_unit(gf);
    dist_read_field(
        gf, ssprintf("huge-data/gauge_field ; ") + show(new_size_node) + "_2");
    displayln_info(ssprintf("dist_crc32 = %08X", field_dist_crc32(gf)));
    displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
    qassert(crc == field_dist_crc32(gf));
  }
  qassert(crc == field_dist_crc32(gf));
  Field<float> fgf;
  convert_field_float_from_double(fgf, gf);
  GaugeField dfgf;
  convert_field_double_from_float(dfgf, fgf);
  const crc32_t fcrc = field_dist_crc32(dfgf);
  for (size_t i = 0; i < new_size_nodes.size(); ++i) {
    const Coordinate& new_size_node = new_size_nodes[i];
    displayln_info(ssprintf("dist_crc32 = %08X", field_dist_crc32(gf)));
    displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
    dist_write_field_float_from_double(
        gf, new_size_node,
        ssprintf("huge-data/gauge_field ; ") + show(new_size_node) + "_f");
    set_unit(gf);
    dist_read_field_double_from_float(
        gf, ssprintf("huge-data/gauge_field ; ") + show(new_size_node) + "_f");
    displayln_info(ssprintf("dist_crc32 = %08X", field_dist_crc32(gf)));
    displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
    qassert(fcrc == field_dist_crc32(gf));
  }
  set_g_rand_color_matrix_field(gf, RngState(rs, "rgf-0.1"), 0.1);
  GaugeField ugf;
  ugf = gf;
  unitarize(ugf);
  const crc32_t ucrc = field_dist_crc32(ugf);
  write_field(gf, ssprintf("huge-data/rgf-0.1.field.conf.%010d", 0));
  set_unit(gf);
  read_field(gf, ssprintf("huge-data/rgf-0.1.field.conf.%010d", 0));
  displayln_info(ssprintf("dist_crc32 = %08X", field_dist_crc32(gf)));
  displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
  qassert(crc == field_dist_crc32(gf));
  set_unit(gf);
  read_field(gf, ssprintf("huge-data/rgf-0.1.field.conf.%010d", 0));
  displayln_info(ssprintf("dist_crc32 = %08X", field_dist_crc32(gf)));
  displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
  qassert(crc == field_dist_crc32(gf));
  save_gauge_field(gf, ssprintf("huge-data/rgf-0.1.gf.conf.%010d", 0));
  set_unit(gf);
  load_gauge_field(gf, ssprintf("huge-data/rgf-0.1.gf.conf.%010d", 0));
  displayln_info(ssprintf("dist_crc32 = %08X", field_dist_crc32(gf)));
  displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
  qassert(ucrc == field_dist_crc32(gf));
  set_unit(gf);
  load_gauge_field(gf, ssprintf("huge-data/rgf-0.1.gf.conf.%010d", 0));
  displayln_info(ssprintf("dist_crc32 = %08X", field_dist_crc32(gf)));
  displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
  qassert(ucrc == field_dist_crc32(gf));
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  test_io();
  test_shuffle();
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
