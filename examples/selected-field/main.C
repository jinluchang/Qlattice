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
  write_selected_field(sf, "huge-data/f.sfield", fsel);
  //
  displayln_info(fname + ssprintf(": clear field 'f' in memory"));
  set_zero(f);
  //
  displayln_info(fname + ssprintf(": compute crc32=%06X.", field_crc32(f)));
  //
  displayln_info(fname + ssprintf(": only read the selected points data of 'f' from disk"));
  read_selected_field(sf, "huge-data/f.sfield", fsel);
  set_field_selected(f, sf, fsel);
  //
  displayln_info(fname + ssprintf(": compute crc32=%06X.", field_crc32(f)));
  //
  displayln_info(fname + ssprintf(": compute the sum of all the selected points in 'f'"));
  {
    Complex sum = 0.0;
    for (long index = 0; index < f.geo().local_volume(); ++index) {
      const long idx = fsel.f_local_idx.get_elems_const(index)[0];
      if (idx < 0) {
        continue;
      }
      const Coordinate xl = f.geo().coordinate_from_index(index);
      const Vector<Complex> fv = f.get_elems_const(xl);
      // f.geto.multiplicity = 2 is the number of elements per site
      for (int m = 0; m < f.geo().multiplicity; ++m) {
        sum += fv[m];
      }
    }
    displayln_info(fname + ssprintf(": v1 sum is %s", show(sum).c_str()));
  }
  {
    Complex sum = 0.0;
    for (long index = 0; index < f.geo().local_volume(); ++index) {
      const long idx = fsel.f_local_idx.get_elems_const(index)[0];
      if (idx < 0) {
        continue;
      }
      const Vector<Complex> fv = sf.get_elems_const(idx);
      // f.geto.multiplicity = 2 is the number of elements per site
      for (int m = 0; m < f.geo().multiplicity; ++m) {
        sum += fv[m];
      }
    }
    displayln_info(fname + ssprintf(": v2 sum is %s", show(sum).c_str()));
  }
  {
    Complex sum = 0.0;
    for (long idx = 0; idx < sf.n_elems; ++idx) {
      // const long index = fsel.indices[idx];
      const Vector<Complex> fv = sf.get_elems_const(idx);
      // f.geto.multiplicity = 2 is the number of elements per site
      for (int m = 0; m < f.geo().multiplicity; ++m) {
        sum += fv[m];
      }
    }
    displayln_info(fname + ssprintf(": v3 sum is %s", show(sum).c_str()));
  }
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
  set_zero(f);
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
  // write and read field only selected different format
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init"));
  const crc32_t crc_0 = field_crc32(f);
  displayln_info(ssprintf("%06X full", crc_0));
  //
  set_u_rand_double(f, rs.split("f-init"));
  write_field_64(f, "huge-data/" + tag + "/f-init.64.field");
  write_field_double(f, "huge-data/" + tag + "/f-init.double.field");
  write_field_float_from_double(f, "huge-data/" + tag + "/f-init.float.field");
  //
  only_keep_selected_points(f, fsel);
  const crc32_t crc_1 = field_crc32(f);
  displayln_info(ssprintf("%06X selected", crc_1));
  //
  set_u_rand_double(f, rs.split("f-init"));
  write_selected_field_64(f, "huge-data/" + tag + "/sf-init.64.field", fsel);
  write_selected_field_double(f, "huge-data/" + tag + "/sf-init.double.field", fsel);
  write_selected_field_float_from_double(f, "huge-data/" + tag + "/sf-init.float.field", fsel);
  //
  set_zero(f);
  read_field_64(f, "huge-data/" + tag + "/f-init.64.field");
  const crc32_t crc_0_0 = field_crc32(f);
  displayln_info(ssprintf("%06X full read 64", crc_0_0));
  qassert(crc_0_0 == crc_0);
  set_zero(f);
  read_field_double(f, "huge-data/" + tag + "/f-init.double.field");
  const crc32_t crc_0_1 = field_crc32(f);
  displayln_info(ssprintf("%06X full read double", crc_0_1));
  qassert(crc_0_1 == crc_0);
  set_zero(f);
  read_field_double_from_float(f, "huge-data/" + tag + "/f-init.float.field");
  const crc32_t crc_0_2 = field_crc32(f);
  displayln_info(ssprintf("%06X full read float", crc_0_2));
  //
  set_zero(f);
  read_selected_field_64(f, "huge-data/" + tag + "/sf-init.64.field", fsel);
  const crc32_t crc_1_0 = field_crc32(f);
  displayln_info(ssprintf("%06X selected read 64", crc_1_0));
  qassert(crc_1_0 == crc_1);
  set_zero(f);
  read_selected_field_double(f, "huge-data/" + tag + "/sf-init.double.field", fsel);
  const crc32_t crc_1_1 = field_crc32(f);
  displayln_info(ssprintf("%06X selected read double", crc_1_1));
  qassert(crc_1_1 == crc_1);
  set_zero(f);
  read_selected_field_double_from_float(f, "huge-data/" + tag + "/sf-init.float.field", fsel);
  const crc32_t crc_1_2 = field_crc32(f);
  displayln_info(ssprintf("%06X selected read float", crc_1_2));
  check_all_files_crc32_info("huge-data/" + tag);
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
  set_zero(f);
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
  check_all_files_crc32_info("huge-data/" + tag);
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
      Coordinate xgr;
      for (int m = 0; m < 4; ++m) {
        xgr[m] = rand_gen(rsi);
      }
      const Coordinate xg = mod(xgr, total_site);
      psel.push_back(xg);
    }
  }
  qmkdir_info("huge-data");
  qmkdir_info("huge-data/" + tag);
  save_point_selection_info(psel, "huge-data/" + tag + "/point-selection.txt");
  const PointSelection psel_load =
      load_point_selection_info("huge-data/" + tag + "/point-selection.txt");
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
  save_selected_points_complex(sp, "huge-data/" + tag + "/f.lat");
  SelectedPoints<Complex> sp2;
  load_selected_points_complex(sp2, "huge-data/" + tag + "/f.lat");
  Field<Complex> f2;
  set_field_selected(f2, sp2, f.geo(), psel);
  const crc32_t crc2 = field_crc32(f2);
  displayln_info(ssprintf(": %06X <- f2 set and set", crc2));
  //
  qassert(crc1 == crc2);
  //
  const Coordinate new_size_node(1, 1, 2, 4);
  const ShuffledBitSet sbs = mk_shuffled_bitset(total_site, psel, new_size_node);
  const std::string path = "huge-data/" + tag + "/fields";
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
  check_all_files_crc32_info("huge-data/" + tag);
}

inline void test_shift(const std::string& tag, const long n_per_tslice, const long n_points)
{
  TIMER_VERBOSE("test_shift");
  const Coordinate total_site(4,4,4,8);
  Geometry geo;
  geo.init(total_site, 1);
  RngState rs = RngState("test_shift");
  Field<Complex> f;
  f.init(geo, 2);
  set_zero(f);
  set_u_rand_double(f, rs.split("f-init"));
  displayln_info(fname + ssprintf(": f crc32 = %08X", field_crc32(f)));
  SelectedField<Complex> sf;
  const PointSelection psel =
      mk_random_point_selection(total_site, n_points, rs.split("psel"));
  FieldSelection fsel;
  set_field_selection(fsel, total_site, n_per_tslice, rs.split("fsel"), psel);
  set_selected_field(sf, f, fsel);
  SelectedField<Complex> sf0;
  sf0 = sf;
  Field<Complex> f0;
  set_field_selected(f0, sf, fsel);
  displayln_info(fname + ssprintf(": f0 (with sparse) crc32 = %08X", field_crc32(f0)));
  std::vector<Coordinate> shift_list;
  shift_list.push_back(Coordinate());
  shift_list.push_back(Coordinate(12, 32, 22, -123));
  shift_list.push_back(Coordinate(0, 0, 0, -123));
  shift_list.push_back(Coordinate(0, 2, 2, 0));
  for (int i = 0; i < (int)shift_list.size(); ++i) {
    const Coordinate& shift = shift_list[i];
    f.init();
    f = f0;
    field_shift(f, f, shift);
    const crc32_t crc0 = field_crc32(f);
    displayln_info(fname +
                   ssprintf(": f (with sparse and shift) crc32 = %08X ; shift=%s", crc0, show(shift).c_str()));
    const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, shift);
    f.init();
    set_zero(f);
    sf.init();
    sf = sf0;
    field_shift(sf, sf, ssp);
    set_field_selected(f, sf, ssp.fsel);
    const crc32_t crc1 = field_crc32(f);
    displayln_info(fname +
                   ssprintf(": f (from shifted sparse field) crc32 = %08X ; shift=%s", crc1, show(shift).c_str()));
    qassert(crc1 == crc0);
  }
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  demo();
  test("16", 16);
  test("all", -1);
  test_grid("grid-16", 16);
  test_selected_points("points-1024", 1024);
  test_shift("shift-16-1024", 2, 8);
  test_shift("shift-0-1024", 0, 8);
  test_shift("shift-16-0", 2, 0);
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
