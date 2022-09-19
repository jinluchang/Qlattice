#include <qlat/qcd.h>
#include <sys/sysinfo.h>

void resize_EigenM(std::vector< qlat::vector<qlat::Complex > >& a, size_t n0, size_t n1)
{
  TIMER_VERBOSE("resize_EigenM");
  /////a.resize(0);
  a.resize(n0);
  for(size_t iv=0;iv<n0;iv++)
  {
    qlat::displayln_info(fname + qlat::ssprintf(": %ld", iv));
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

  std::vector< qlat::vector<qlat::Complex > > src;
  resize_EigenM(src, 2, 3);
  resize_EigenM(src,10,10);

  int nmass = 2;
  std::vector<Propagator4d > prop0;
  prop0.resize(nmass);for(unsigned int i=0;i<prop0.size();i++)prop0[i].init(geo);
  qacc_for(isp, long(geo.local_volume()),{
    for(int im=0;im<nmass;im++)
    {   
    qlat::WilsonMatrix& v0 =  prop0[im].get_elem(isp);
    for(int d0=0;d0<12;d0++)
    for(int d1=0;d1<12;d1++)
    {   
      v0(d0,d1) = im*7 + isp*0.5 + d0*7 + d1*8;
    }   
    }   
  }); 

  std::vector<Propagator4d > propS;propS.resize(1);for(unsigned int i=0;i<propS.size();i++)propS[i].init(geo);
  ////propS.resize(0);
  propS.resize(nmass);for(unsigned int i=0;i<propS.size();i++)propS[i].init(geo);
  for(int im=0;im<nmass;im++){propS[im] = prop0[im];}

  std::vector<Propagator4d > propD;propD.resize(nmass);for(unsigned int i=0;i<propD.size();i++)propD[i].init(geo);
  for(int im=0;im<nmass;im++){propD[im] = prop0[im];}

  double diff = 0.0;
  for(int im=0;im<nmass;im++)
  {
  qacc_for(isp, long(geo.local_volume()), {
      qlat::WilsonMatrix&  vs =  propS[im].get_elem(isp);
      qlat::WilsonMatrix&  vd =  propD[im].get_elem(isp);
      for(int d0=0;d0<12;d0++)
      for(int d1=0;d1<12;d1++)
      {
        qlat::Complex tem = vs(d0,d1) - vd(d0,d1);
        diff += tem.real()*tem.real() + tem.imag()*tem.imag();
      }
  });
  }

  qlat::displayln_info(qlat::ssprintf("Prop ori diff %.5e \n",diff));




  qlat::Timer::display();
  qlat::end();
  return 0;
}

