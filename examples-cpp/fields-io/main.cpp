#include <qlat/qlat.h>

using namespace qlat;

inline void test_read(const std::string& path, const std::string& fn)
{
  TIMER_VERBOSE("test_read");
  Field<Complex> f;
  read_field_double_from_float(f, path, fn);
  const crc32_t crc = field_crc32(f);
  displayln_info(fname +
                 ssprintf(": compute crc32=%08X for fn='%s' from path '%s'.",
                          crc, fn.c_str(), path.c_str()));
}

inline void test_read_sbs(const std::string& path, const std::string& fn, const ShuffledBitSet& sbs)
{
  TIMER_VERBOSE("test_read_sbs");
  Field<Complex> f;
  read_field_double_from_float(f, path, fn);
  const crc32_t crc = field_crc32(f);
  displayln_info(fname +
                 ssprintf(": compute crc32=%08X for fn='%s' from path '%s'.",
                          crc, fn.c_str(), path.c_str()));
  SelectedField<Complex> sf;
  read_field_double_from_float(sf, path, fn, sbs);
  f.init();
  set_field_selected(f, sf, sbs.fsel);
  const crc32_t crc_1 = field_crc32(f);
  displayln_info(fname +
                 ssprintf(": compute crc32=%08X for fn='%s' from path '%s' with sbs.",
                          crc_1, fn.c_str(), path.c_str()));
  qassert(crc_1 == crc);
}

inline void demo(const std::string& tag, const Coordinate& total_site,
                 const long n_random_points, const long n_per_tslice,
                 const Coordinate& new_size_node)
{
  TIMER_VERBOSE("demo");
  Geometry geo;
  geo.init(total_site, 1);
  qmkdir_info("huge-data");
  qmkdir_info("huge-data/" + tag);
  //
  displayln_info(fname + ssprintf(": create splittable random number generator 'rs'"));
  const RngState rs = RngState("selected-field");
  //
  displayln_info(fname + ssprintf(": random select %d points on entire lattice",
                                  n_random_points));
  const PointSelection psel =
      mk_random_point_selection(total_site, n_random_points, rs.split("psel"));
  //
  displayln_info(fname + ssprintf(": random select points on each time slice"));
  displayln_info(fname + ssprintf(": select %d points per time slice", n_per_tslice));
  displayln_info(fname + ssprintf(": selection is stored in 'fsel' of type 'FieldSelection'"));
  FieldSelection fsel;
  mk_field_selection(fsel.f_rank, total_site, n_per_tslice, rs.split("demo"));
  add_field_selection(fsel.f_rank, psel);
  update_field_selection(fsel);
  update_field_selection(fsel, n_per_tslice);
  //
  displayln_info(fname + ssprintf(": save selection 'fsel' to disk as a field"));
  write_field_selection(fsel, "huge-data/" + tag + "/fsel.field");
  //
  fsel.init();
  displayln_info(fname + ssprintf(": read selection 'fsel' from disk as a field"));
  read_field_selection(fsel, "huge-data/" + tag + "/fsel.field", n_per_tslice);
  //
  displayln_info(fname + ssprintf(": init field 'f'"));
  Field<Complex> f, sf, rf;
  f.init(geo, 2);
  set_zero(f);
  //
  const crc32_t crc_1 = field_crc32(f);
  displayln_info(fname + ssprintf(": compute crc32=%08X.", crc_1));
  //
  displayln_info(fname + ssprintf(": fill with random numbers"));
  set_u_rand_double(f, rs.split("f-init"));
  //
  const crc32_t crc_2 = field_crc32(f);
  displayln_info(fname + ssprintf(": compute crc32=%08X.", crc_2));
  //
  displayln_info(
      fname +
      ssprintf(
          ": possible to only keep selected points non-zero for field 'f'."));
  displayln_info(fname + ssprintf(": (DOES NOT CHANGE THE SIZE IN MOMEORY)"));
  sf = f;
  only_keep_selected_points(sf, fsel);
  //
  const crc32_t crc_3 = field_crc32(sf);
  displayln_info(fname + ssprintf(": compute crc32=%08X.", crc_3));
  //
  {
    ShuffledFieldsWriter sfw("huge-data/" + tag + "/demo.lfs", new_size_node);
    //
    const ShuffledBitSet sbs = mk_shuffled_bitset(fsel, new_size_node);
    //
    displayln_info(fname + ssprintf(": save to disk"));
    write(sfw, "f.field", f);
    write(sfw, "f.sfield", f, sbs);
    write_float_from_double(sfw, "f.float.field", f);
    write_float_from_double(sfw, "f.float.sfield", f, sbs);
    //
    //
    displayln_info(fname + ssprintf(": save sf to disk"));
    write(sfw, "sf.field", sf);
    write(sfw, "sf.sfield", sf, sbs);
    write_float_from_double(sfw, "sf.float.field", sf);
    write_float_from_double(sfw, "sf.float.sfield", sf, sbs);
  }
  {
    ShuffledFieldsWriter sfw("huge-data/" + tag + "/demo.lfs", new_size_node, true);
    //
    const ShuffledBitSet sbs = mk_shuffled_bitset(fsel, new_size_node);
    //
    displayln_info(fname + ssprintf(": save to disk append"));
    write(sfw, "fa.field", f);
    write(sfw, "fa.sfield", f, sbs);
    write_float_from_double(sfw, "fa.float.field", f);
    write_float_from_double(sfw, "fa.float.sfield", f, sbs);
    //
    //
    displayln_info(fname + ssprintf(": save sf to disk append"));
    write(sfw, "sfa.field", sf);
    write(sfw, "sfa.sfield", sf, sbs);
    write_float_from_double(sfw, "sfa.float.field", sf);
    write_float_from_double(sfw, "sfa.float.sfield", sf, sbs);
  }
  //
  {
    const std::string path = "huge-data/" + tag + "/demo.lfs";
    const std::vector<std::string> fns = list_fields(path);
    for (long i = 0; i < (long)fns.size(); ++i) {
      displayln_info(fname + ssprintf(": %5d : '%s' from '%s'.", i, fns[i].c_str(), path.c_str()));
    }
  }
  //
  {
    const ShuffledBitSet sbs = mk_shuffled_bitset(fsel, new_size_node);
    test_read("huge-data/" + tag + "/demo.lfs", "f.float.field");
    test_read_sbs("huge-data/" + tag + "/demo.lfs", "f.float.sfield", sbs);
    test_read("huge-data/" + tag + "/demo.lfs", "sf.float.field");
    test_read_sbs("huge-data/" + tag + "/demo.lfs", "sf.float.sfield", sbs);
    test_read("huge-data/" + tag + "/demo.lfs", "fa.float.field");
    test_read_sbs("huge-data/" + tag + "/demo.lfs", "fa.float.sfield", sbs);
    test_read("huge-data/" + tag + "/demo.lfs", "sfa.float.field");
    test_read_sbs("huge-data/" + tag + "/demo.lfs", "sfa.float.sfield", sbs);
  }
  //
  {
    ShuffledFieldsReader sfr("huge-data/" + tag + "/demo.lfs");
    for (int i = 0; i < 4; ++i) {
      displayln_info(fname + ssprintf(": read from disk 'f.field'"));
      //
      rf.init();
      read(sfr, "f.field", rf);
      const crc32_t crc_4 = field_crc32(rf);
      displayln_info(fname +
                     ssprintf(": compute crc32=%08X.", crc_4));
      qassert(crc_4 == crc_2)
      //
      displayln_info(fname + ssprintf(": read from disk 'sf.field'"));
      //
      rf.init();
      read(sfr, "sf.field", rf);
      const crc32_t crc_5 = field_crc32(rf);
      displayln_info(fname +
                     ssprintf(": compute crc32=%08X.", crc_5));
      qassert(crc_5 == crc_3)
      //
      rf.init();
      read(sfr, "f.sfield", rf);
      const crc32_t crc_6 = field_crc32(rf);
      displayln_info(fname +
                     ssprintf(": compute crc32=%08X.", crc_6));
      qassert(crc_6 == crc_3)
      //
      rf.init();
      read(sfr, "sf.sfield", rf);
      const crc32_t crc_7 = field_crc32(rf);
      displayln_info(fname +
                     ssprintf(": compute crc32=%08X.", crc_7));
      qassert(crc_7 == crc_3)
      //
      rf.init();
      read(sfr, "fa.field", rf);
      const crc32_t crc_8 = field_crc32(rf);
      displayln_info(fname +
                     ssprintf(": compute crc32=%08X.", crc_4));
      qassert(crc_8 == crc_2)
      //
      rf.init();
      read(sfr, "sfa.field", rf);
      const crc32_t crc_9 = field_crc32(rf);
      displayln_info(fname +
                     ssprintf(": compute crc32=%08X.", crc_5));
      qassert(crc_9 == crc_3)
      //
      rf.init();
      read(sfr, "fa.sfield", rf);
      const crc32_t crc_10 = field_crc32(rf);
      displayln_info(fname +
                     ssprintf(": compute crc32=%08X.", crc_6));
      qassert(crc_10 == crc_3)
      //
      rf.init();
      read(sfr, "sfa.sfield", rf);
      const crc32_t crc_11 = field_crc32(rf);
      displayln_info(fname +
                     ssprintf(": compute crc32=%08X.", crc_7));
      qassert(crc_11 == crc_3)
    }
  }
  check_all_files_crc32("huge-data/" + tag);
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  demo("t0", Coordinate(6, 6, 6, 8), 0, 16, Coordinate(2, 2, 2, 8));
  demo("t1", Coordinate(6, 6, 6, 8), 4, 2, Coordinate(2, 2, 2, 8));
  demo("t2", Coordinate(6, 6, 6, 8), 8, 0, Coordinate(2, 2, 2, 8));
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
