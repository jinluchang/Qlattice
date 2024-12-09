#include "qlat/grid.h"

namespace qlat
{  //

inline void test()
{
  const Coordinate total_site(4, 4, 4, 8);
  const Geometry geo(total_site);
  GaugeField gf;
  gf.init(geo);
  const RngState rs = RngState().split("load_configuration");
  set_g_rand_color_matrix_field(gf, rs, 1.0);
  gf_show_info(gf);
  Grid::GridCartesian* UGrid = Grid::SpaceTimeGrid::makeFourDimGrid(
      grid_convert(total_site),
      Grid::GridDefaultSimd(Grid::Nd, Grid::vComplexD::Nsimd()),
      grid_convert(geo.geon.size_node));
  qassert(geo.geon.id_node == id_node_from_grid(UGrid));
  Grid::LatticeGaugeField ggf(UGrid);
  grid_convert(ggf, gf);
  delete UGrid;
}

}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 1));
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 2, 2));
  size_node_list.push_back(Coordinate(1, 2, 2, 2));
  size_node_list.push_back(Coordinate(2, 2, 2, 2));
  size_node_list.push_back(Coordinate(2, 2, 2, 4));
  size_node_list.push_back(Coordinate(2, 2, 4, 4));
  size_node_list.push_back(Coordinate(2, 4, 4, 4));
  size_node_list.push_back(Coordinate(4, 4, 4, 4));
  size_node_list.push_back(Coordinate(4, 4, 4, 8));
  grid_begin(&argc, &argv, size_node_list);
  test();
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  grid_end();
  return 0;
}
