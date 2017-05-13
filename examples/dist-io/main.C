#include <qlat/qlat.h>

#include <iostream>
#include <complex>
#include <vector>

using namespace qlat;
using namespace std;

void test_io()
{
  TIMER("test_io");
  qmkdir_sync_node("huge-data");
  RngState rs(get_global_rng_state(), fname);
  crc32_check();
  dist_write_par_limit() = 10;
  dist_read_par_limit() = 10;
  // Coordinate total_site(48, 48, 48, 96);
  Coordinate total_site(16, 16, 16, 32);
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
  get_shuffle_plan_cache().limit = 16;
  displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
  crc32_t crc = field_crc32(gf);
  std::vector<Coordinate> new_size_nodes;
  new_size_nodes.push_back(Coordinate(1,2,1,32));
  new_size_nodes.push_back(Coordinate(1,1,1,8));
  new_size_nodes.push_back(Coordinate(1,2,4,2));
  new_size_nodes.push_back(Coordinate(2,2,1,1));
  new_size_nodes.push_back(Coordinate(16,2,1,32));
  new_size_nodes.push_back(Coordinate(16,1,1,1));
  new_size_nodes.push_back(Coordinate(4,8,1,1));
  new_size_nodes.push_back(Coordinate(4,2,16,32));
  new_size_nodes.push_back(Coordinate(4,16,16,32));
  new_size_nodes.push_back(Coordinate(4,4,4,4));
  for (size_t i = 0; i < new_size_nodes.size(); ++i) {
    const Coordinate new_size_node = new_size_nodes[i];
    displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
    dist_write_field(gf, new_size_node, ssprintf("huge-data/gauge_field ; ") + show(new_size_node));
    set_unit(gf);
    dist_read_field(gf, new_size_node, ssprintf("huge-data/gauge_field ; ") + show(new_size_node));
    displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
    qassert(crc == field_crc32(gf));
  }
  for (size_t i = 0; i < new_size_nodes.size(); ++i) {
    const Coordinate new_size_node = new_size_nodes[i];
    displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
    dist_write_field(gf, new_size_node, ssprintf("huge-data/gauge_field ; ") + show(new_size_node) + "_2");
    set_unit(gf);
    dist_read_field(gf, new_size_node, ssprintf("huge-data/gauge_field ; ") + show(new_size_node) + "_2");
    displayln_info(ssprintf("crc32 = %08X", field_crc32(gf)));
    qassert(crc == field_crc32(gf));
  }
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  test_io();
  Timer::display();
  end();
  return 0;
}
