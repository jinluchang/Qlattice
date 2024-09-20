#include <qlat/qlat.h>

using namespace qlat;

inline void demo()
{
  const std::string path = "huge-data/fields";
  const Coordinate new_size_node = Coordinate(2, 2, 2, 4);
  const bool is_append = false;
  get_shuffled_fields_writer(path, new_size_node, is_append);
  const Coordinate total_site = Coordinate(4, 4, 4, 8);
  const Geometry geo(total_site);
  Field<Complex> f1;
  f1.init(geo, 2);
  write_field(f1, path, "f1");
  Field<ComplexF> f2;
  f2.init(geo, 3);
  write_field(f2, path, "f2");
  clear_shuffled_fields_writer_cache();
  clear_shuffled_fields_reader_cache();
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
