#include <qlat/qcd.h>
#include <sys/sysinfo.h>



int main(int argc, char* argv[])
{
  using namespace qlat;

  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1,  1));
  size_node_list.push_back(Coordinate(1, 1, 1,  2));
  size_node_list.push_back(Coordinate(1, 1, 1,  4));
  size_node_list.push_back(Coordinate(1, 1, 1,  8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 1, 32));
  size_node_list.push_back(Coordinate(1, 1, 1, 64));
  size_node_list.push_back(Coordinate(4, 4, 8, 16));
  size_node_list.push_back(Coordinate(4, 8, 8, 16));

  begin(&argc, &argv, size_node_list);

  int nx,ny,nz,nt;
  nx = 24;
  ny = 24;
  nz = 24;
  nt = 48;

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 

  qlat::vector_acc<Complex > src;src.resize(nt);
  qacc_for(isp , nt, {
    src[isp] += 1;
  });
  for(int isp=0;isp<src.size();isp++)
  {
    qlat::displayln_info(qlat::ssprintf("t %5d %.8e \n", isp, src[isp].real()));
  }




  qlat::Timer::display();
  qlat::end();
  return 0;
}

