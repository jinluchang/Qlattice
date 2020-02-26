#include <qlat/qlat.h>

using namespace qlat;

inline void demo()
{
  TIMER_VERBOSE("demo");
  const Coordinate total_site(6,6,6,8);
  Geometry geo;
  geo.init(total_site, 1);
  qmkdir_info("huge-data");
  //
  displayln_info(fname + ssprintf(": create splittable random number generator 'rs'"));
  const RngState rs = RngState("selected-field");
  //
  const long n_per_tslice = 16;
  displayln_info(fname + ssprintf(": random select points on each time slice"));
  displayln_info(fname + ssprintf(": select %d points per time slice", n_per_tslice));
  displayln_info(fname + ssprintf(": selection is stored in 'fsel' of type 'FieldSelection'"));
  FieldSelection fsel;
  set_field_selection(fsel, total_site, n_per_tslice, rs.split("demo"));
  //
  displayln_info(fname + ssprintf(": save selection 'fsel' to disk as a field"));
  write_field_selection(fsel, "huge-data/fsel.field");
  //
  fsel.init();
  displayln_info(fname + ssprintf(": read selection 'fsel' from disk as a field"));
  read_field_selection(fsel, "huge-data/fsel.field", n_per_tslice);
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
  const Coordinate new_size_node(2, 2, 2, 8);
  //
  {
    ShuffledFieldsWriter sfw("huge-data/demo.lfs", new_size_node);
    //
    const ShuffledBitSet sbs = mk_shuffled_bitset(fsel, new_size_node);
    //
    displayln_info(fname + ssprintf(": save to disk"));
    write(sfw, "f.field", f);
    write(sfw, "f.sfield", f, sbs);
    //
    //
    displayln_info(fname + ssprintf(": save sf to disk"));
    write(sfw, "sf.field", sf);
    write(sfw, "sf.sfield", sf, sbs);
  }
  {
    ShuffledFieldsWriter sfw("huge-data/demo.lfs", new_size_node, true);
    //
    const ShuffledBitSet sbs = mk_shuffled_bitset(fsel, new_size_node);
    //
    displayln_info(fname + ssprintf(": save to disk"));
    write(sfw, "fa.field", f);
    write(sfw, "fa.sfield", f, sbs);
    //
    //
    displayln_info(fname + ssprintf(": save sf to disk"));
    write(sfw, "sfa.field", sf);
    write(sfw, "sfa.sfield", sf, sbs);
  }
  //
  {
    ShuffledFieldsReader sfr("huge-data/demo.lfs");
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
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  demo();
  Timer::display();
  end();
  return 0;
}
