#include <qlat/qlat.h>

int main(int argc, char* argv[])
{
  using namespace qlat;
  begin(&argc, &argv);
  displayln_info("hello world.");
  displayln(ssprintf("%d / %d", get_id_node(), get_num_node()));
  Geometry geo;
  geo.init(Coordinate(8,8,8,8), 1);
  FieldM<Complex,1> field;
  field.init(geo);
  set_zero(field);
  end();
  return 0;
}
