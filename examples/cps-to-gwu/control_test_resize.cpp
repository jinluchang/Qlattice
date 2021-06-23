#include <qlat/qlat.h>
#include <sys/sysinfo.h>

void resize_EigenM(std::vector< qlat::vector<qlat::Complex > >& a, size_t n0, size_t n1)
{
  /////a.resize(0);
  a.resize(n0);
  for(size_t iv=0;iv<n0;iv++)
  {
    a[iv].resize(n1);
    set_zero(a[iv]);
  }
}

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
  nt = 64;


  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 

  std::vector<Propagator4d > propS;
  propS.resize(1);for(unsigned int i=0;i<propS.size();i++)propS[i].init(geo);

  ////propS.resize(0);
  propS.resize(5);for(unsigned int i=0;i<propS.size();i++)propS[i].init(geo);

  propS.resize(2);for(unsigned int i=0;i<propS.size();i++)propS[i].init(geo);

  std::vector< qlat::vector<qlat::Complex > > src;

  resize_EigenM(src, 2, 3);
  resize_EigenM(src,10,10);


  qlat::Timer::display();
  qlat::end();
  return 0;
}

