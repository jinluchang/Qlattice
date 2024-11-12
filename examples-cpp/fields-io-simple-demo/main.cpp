#include <qlat/qlat.h>

using namespace qlat;

inline void demo()
{
  const std::string path = "huge-data/fields";
  const Coordinate new_size_node = Coordinate(2, 2, 2, 4);
  const bool is_append = false;
  const bool is_removing_old = true;
  get_shuffled_fields_writer(path, new_size_node, is_append, is_removing_old);
  const Coordinate total_site = Coordinate(4, 4, 4, 8);
  const Geometry geo(total_site);
  const RngState rs = RngState("seed");
  //
  Field<Complex> f1;
  f1.init(geo, 2);
  set_zero(f1);
  set_u_rand(f1, rs.split("f1"));
  const crc32_t crc_f1 = field_crc32(f1);
  Field<ComplexF> f2;
  f2.init(geo, 3);
  set_zero(f2);
  set_u_rand(f2, rs.split("f2"));
  const crc32_t crc_f2 = field_crc32(f2);
  //
  write_field(f1, path, "f1");
  write_field(f2, path, "f2");
  close_shuffled_fields_writer(path);
  //
  Field<Complex> f1r;
  Field<ComplexF> f2r;
  //
  read_field(f1r, path, "f1");
  read_field(f2r, path, "f2");
  close_shuffled_fields_reader(path);
  //
  const crc32_t crc_f1r = field_crc32(f1r);
  qassert(crc_f1 == crc_f1r);
  const crc32_t crc_f2r = field_crc32(f2r);
  qassert(crc_f2 == crc_f2r);
  //
  clear_shuffled_fields_writer_cache();
  clear_shuffled_fields_reader_cache();
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  demo();
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
