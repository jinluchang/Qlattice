#include <qlat/qlat.h>

using namespace qlat;

inline void demo()
{
  TIMER_VERBOSE("demo");
  const Coordinate total_site(4,4,4,8);
  Geometry geo;
  geo.init(total_site, 1);
  qmkdir_info("huge-data");
  //
  displayln_info(fname + ssprintf(": create splittable random number generator 'rs'"));
  RngState rs = RngState("selected-field");
  //
  displayln_info(fname + ssprintf(": init field 'f'"));
  Field<Complex> f;
  f.init(geo, 2);
  set_zero(f);
  //
  displayln_info(fname + ssprintf(": compute crc32=%06X.", field_crc32(f)));
  //
  displayln_info(fname + ssprintf(": fill with random numbers"));
  set_u_rand_double(f, rs.split("f-init"));
  //
  displayln_info(fname + ssprintf(": compute crc32=%06X.", field_crc32(f)));
  //
  displayln_info(fname + ssprintf(": save to disk"));
  write_field(f, "huge-data/f.field");
  //
  displayln_info(fname + ssprintf(": clear field 'f' in memory"));
  set_zero(f);
  //
  displayln_info(fname + ssprintf(": compute crc32=%06X.", field_crc32(f)));
  //
  displayln_info(fname + ssprintf(": read from disk"));
  read_field(f, "huge-data/f.field");
  //
  displayln_info(fname + ssprintf(": compute crc32=%06X.", field_crc32(f)));
  //
  const long n_per_tslice = 16;
  displayln_info(fname + ssprintf(": random select points on each time slice"));
  displayln_info(fname + ssprintf(": select %d points per time slice", n_per_tslice));
  displayln_info(fname + ssprintf(": selection is stored in 'fsel' of type 'FieldSelection'"));
  FieldSelection fsel;
  set_field_selection(fsel, total_site, n_per_tslice, rs.split("free-4nt8").split(0));
  //
  displayln_info(fname + ssprintf(": save selection 'fsel' to disk as a field"));
  write_field_selection(fsel, "huge-data/fsel.field");
  //
  fsel.init();
  displayln_info(fname + ssprintf(": read selection 'fsel' from disk as a field"));
  read_field_selection(fsel, "huge-data/fsel.field", n_per_tslice);
  //
  displayln_info(fname + ssprintf(": possible to only keep selected points non-zero for field 'f'."));
  displayln_info(fname + ssprintf(": (DOES NOT CHANGE THE SIZE IN MOMEORY)"));
  only_keep_selected_points(f, fsel);
  //
  displayln_info(fname + ssprintf(": compute crc32=%06X.", field_crc32(f)));
  //
  displayln_info(fname + ssprintf(": declare 'sf' or type 'SelectedField<Complex>'"));
  SelectedField<Complex> sf;
  //
  displayln_info(fname + ssprintf(": 'sf' only store the selected points' data of 'f'"));
  set_selected_field(sf, f, fsel);
  //
  displayln_info(fname + ssprintf(": clear field 'f' in memory"));
  set_zero(f);
  //
  displayln_info(fname + ssprintf(": compute crc32=%06X.", field_crc32(f)));
  //
  displayln_info(fname + ssprintf(": restore the selected points data to 'f' from 'sf'"));
  set_field_selected(f, sf, fsel);
  //
  displayln_info(fname + ssprintf(": compute crc32=%06X.", field_crc32(f)));
  //
  displayln_info(fname + ssprintf(": only save the selected points data of 'f' to disk"));
  displayln_info(fname + ssprintf(": (NOTE THAT 'f', NOT 'sf', IS USED HERE)"));
  write_selected_field(f, "huge-data/f.sfield", fsel);
  //
  displayln_info(fname + ssprintf(": clear field 'f' in memory"));
  set_zero(f);
  //
  displayln_info(fname + ssprintf(": compute crc32=%06X.", field_crc32(f)));
  //
  displayln_info(fname + ssprintf(": only read the selected points data of 'f' from disk"));
  displayln_info(fname + ssprintf(": (NOTE THAT 'f', NOT 'sf', IS USED HERE)"));
  read_selected_field(f, "huge-data/f.sfield", fsel);
  //
  displayln_info(fname + ssprintf(": compute crc32=%06X.", field_crc32(f)));
  //
  displayln_info(fname + ssprintf(": compute the sum of all the selected points in 'f'"));
  Complex sum = 0.0;
  for (long index = 0; index < f.geo.local_volume(); ++index) {
    const long idx = fsel.f_local_idx.get_elems_const(index)[0];
    if (idx < 0) {
      continue;
    }
    const Coordinate xl = f.geo.coordinate_from_index(index);
    const Vector<Complex> fv = f.get_elems_const(xl);
    // f.geto.multiplicity = 2 is the number of elements per site
    for (int m = 0; m < f.geo.multiplicity; ++m) {
      sum += fv[m];
    }
  }
  displayln_info(fname + ssprintf(": sum is %s", show(sum).c_str()));
}

inline void test(const std::string& tag, const long n_per_tslice)
{
  TIMER_VERBOSE("test");
  const Coordinate total_site(4,4,4,8);
  Geometry geo;
  geo.init(total_site, 1);
  RngState rs("selected-field");
  qmkdir_info("huge-data");
  qmkdir_info("huge-data/" + tag);
  // init f
  Field<Complex> f;
  f.init(geo, 2);
  // init fsel
  FieldSelection fsel;
  set_field_selection(fsel, total_site, n_per_tslice, rs.split("free-4nt8").split(0));
  // test of partial f
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init"));
  displayln_info(ssprintf(": %06X <- f-init", field_crc32(f)));
  only_keep_selected_points(f, fsel);
  displayln_info(ssprintf("%06X <- only_keep_selected_points", field_crc32(f)));
  // test save and load fsel
  write_field_selection(fsel, "huge-data/" + tag + "/fsel.field");
  fsel.init();
  read_field_selection(fsel, "huge-data/" + tag + "/fsel.field", n_per_tslice);
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
  write_field(f, "huge-data/" + tag + "/free-4nt8-init.field");
  displayln_info(ssprintf(": %06X <- f-init", field_crc32(f)));
  write_selected_field(f, "huge-data/" + tag + "/free-4nt8-init.sfield", fsel);
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init-3"));
  displayln_info(ssprintf(": %06X <- f-init-3", field_crc32(f)));
  read_selected_field(f, "huge-data/" + tag + "/free-4nt8-init.sfield", fsel, Coordinate(1,1,1,8));
  displayln_info(ssprintf("%06X <- write and read back", field_crc32(f)));
}

inline void test_grid(const std::string& tag, const long n_per_tslice)
{
  TIMER_VERBOSE("test_grid");
  const Coordinate total_site(4, 4, 4, 8);
  Geometry geo;
  geo.init(total_site, 1);
  RngState rs("selected-field");
  qmkdir_info("huge-data");
  qmkdir_info("huge-data/" + tag);
  // init f
  Field<Complex> f;
  f.init(geo, 2);
  // init fsel
  FieldSelection fsel;
  set_grid_field_selection(fsel, total_site, n_per_tslice,
                           rs.split("free-4nt8").split(0));
  // test of partial f
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init"));
  displayln_info(ssprintf(": %06X <- f-init", field_crc32(f)));
  only_keep_selected_points(f, fsel);
  displayln_info(ssprintf("%06X <- only_keep_selected_points", field_crc32(f)));
  // test save and load fsel
  write_field_selection(fsel, "huge-data/" + tag + "/fsel.field");
  fsel.init();
  read_field_selection(fsel, "huge-data/" + tag + "/fsel.field", n_per_tslice);
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
  write_field(f, "huge-data/" + tag + "/free-4nt8-init.field");
  displayln_info(ssprintf(": %06X <- f-init", field_crc32(f)));
  write_selected_field(f, "huge-data/" + tag + "/free-4nt8-init.sfield", fsel);
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init-3"));
  displayln_info(ssprintf(": %06X <- f-init-3", field_crc32(f)));
  read_selected_field(f, "huge-data/" + tag + "/free-4nt8-init.sfield", fsel,
                      Coordinate(1, 1, 1, 8));
  displayln_info(ssprintf("%06X <- write and read back", field_crc32(f)));
}

inline void test_selected_points(const std::string& tag, const long n_points)
{
  const Coordinate total_site(4, 4, 4, 8);
  Geometry geo;
  geo.init(total_site, 1);
  const RngState rs = RngState("test_selected_points").split(tag);
  PointSelection psel;
  {
    const RngState rs_psel = rs.split("psel");
    for (int i = 0; i < n_points; ++i) {
      RngState rsi = rs_psel.split(i);
      const Coordinate xg = mod(Coordinate(rand_gen(rsi), rand_gen(rsi),
                                           rand_gen(rsi), rand_gen(rsi)),
                                total_site);
      psel.push_back(xg);
    }
  }
  qmkdir_info("results");
  save_point_selection_info(psel, "results/point-selection.txt");
  const PointSelection psel_load =
      load_point_selection_info("results/point-selection.txt");
  qassert(psel == psel_load);
  //
  Field<Complex> f;
  f.init(geo, 2);
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init"));
  displayln_info(ssprintf(": %06X <- f-init", field_crc32(f)));
  //
  Field<Complex> f1;
  f1 = f;
  only_keep_selected_points(f1, psel);
  const crc32_t crc1 = field_crc32(f1);
  displayln_info(ssprintf(": %06X <- f1 only selected", crc1));
  //
  SelectedPoints<Complex> sp;
  set_selected_points(sp, f, psel);
  save_selected_points_complex(sp, "results/f.lat");
  SelectedPoints<Complex> sp2;
  load_selected_points_complex(sp2, "results/f.lat");
  Field<Complex> f2;
  set_field_selected(f2, sp2, f.geo, psel);
  const crc32_t crc2 = field_crc32(f2);
  displayln_info(ssprintf(": %06X <- f2 set and set", crc2));
  //
  qassert(crc1 == crc2);
  //
  const Coordinate new_size_node(1, 1, 2, 4);
  const ShuffledBitSet sbs = mk_shuffled_bitset(total_site, psel, new_size_node);
  const std::string path = "results/fields";
  {
    ShuffledFieldsWriter sfw(path, new_size_node);
    write(sfw, "f.psel", f, sbs);
  }
  Field<Complex> f3;
  {
    ShuffledFieldsReader sfr(path);
    read(sfr, "f.psel", f3);
  }
  const crc32_t crc3 = field_crc32(f3);
  displayln_info(ssprintf(": %06X <- f3 write and read", crc3));
  //
  qassert(crc1 == crc3);
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  demo();
  test("16", 16);
  test("all", -1);
  test_grid("16", 16);
  test_selected_points("16", 1024);
  Timer::display();
  end();
  return 0;
}
