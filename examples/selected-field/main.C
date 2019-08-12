#include <qlat/qlat.h>

using namespace qlat;

inline void test()
{
  TIMER_VERBOSE("test");
  const Coordinate total_site(4,4,4,8);
  Geometry geo;
  geo.init(total_site, 1);
  RngState rs("selected-field");
  qmkdir_info("huge-data");
  // init f
  Field<Complex> f;
  f.init(geo, 2);
  // init fsel
  const long n_per_tslice = 16;
  FieldSelection fsel;
  set_field_selection(fsel, total_site, n_per_tslice, rs.split("free-4nt8").split(0));
  // test of partial f
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init"));
  displayln_info(ssprintf(": %06X <- f-init", field_crc32(f)));
  only_keep_selected_points(f, fsel);
  displayln_info(ssprintf("%06X <- only_keep_selected_points", field_crc32(f)));
  // test save and load fsel
  write_field_selection(fsel, "huge-data/fsel.field");
  fsel.init();
  read_field_selection(fsel, "huge-data/fsel.field", n_per_tslice);
  // test of partial f
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init"));
  displayln_info(ssprintf(": %06X <- f-init", field_crc32(f)));
  only_keep_selected_points(f, fsel);
  displayln_info(ssprintf("%06X <- only_keep_selected_points", field_crc32(f)));
  // test of reconstructed f
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init"));
  displayln_info(ssprintf(": %06X <- f-init", field_crc32(f)));
  SelectedField<Complex> sf;
  set_selected_field(sf, f, fsel);
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init-1"));
  displayln_info(ssprintf(": %06X <- f-init-1", field_crc32(f)));
  set_field_selected(f, sf, fsel);
  displayln_info(ssprintf("%06X <- set_field_selected", field_crc32(f)));
  // test of reconstructed f with slow version
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init"));
  displayln_info(ssprintf(": %06X <- f-init", field_crc32(f)));
  sf.init();
  set_selected_field_slow(sf, f, fsel);
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init-2"));
  displayln_info(ssprintf(": %06X <- f-init-2", field_crc32(f)));
  set_field_selected_slow(f, sf, fsel);
  displayln_info(ssprintf("%06X <- set_field_selected_slow", field_crc32(f)));
  // write and read field only selected
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init"));
  displayln_info(ssprintf(": %06X <- f-init", field_crc32(f)));
  write_selected_field(f, "huge-data/free-4nt8-init.sf", fsel);
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init-3"));
  displayln_info(ssprintf(": %06X <- f-init-3", field_crc32(f)));
  read_selected_field(f, "huge-data/free-4nt8-init.sf", fsel, Coordinate(1,1,1,8));
  displayln_info(ssprintf("%06X <- write and read back", field_crc32(f)));
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  test();
  Timer::display();
  end();
  return 0;
}
