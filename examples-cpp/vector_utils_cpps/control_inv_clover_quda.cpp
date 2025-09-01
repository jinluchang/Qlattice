#include <sys/sysinfo.h>
#include <unistd.h>

//#include <qutils/vector.h>
#include "general_funs.h"
#include "utils_io_vec.h"
#include "utils_clover_inverter.h"
#include "utils_construction.h"

int main(int argc, char* argv[])
{
  using namespace qlat;
  inputpara in;
  begin_Lat(&argc, &argv, in);

  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;

  {
  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site); 
  fflush_MPI();

  std::vector<int > sp;sp.resize(4);
  for(int i=0;i<4;i++){sp[i] = 0;}
  if(in.paraI != "None"){
    std::vector<std::string > Li = stringtolist(in.paraI);
    Qassert(sp.size() == Li.size());
    for(int i=0;i<4;i++){sp[i] = stringtonum(Li[i]);}
  }


  int mpi_layout[4]={0,0,0,0};
  qlat::GeometryNode geon = qlat::get_geometry_node();for(int i=0;i<4;i++){mpi_layout[i] = geon.size_node[i];}

  qlat::GaugeField gf;gf.init(geo);

  char rbc_conf[500];
  sprintf(rbc_conf,in.Link_name.c_str(), in.icfg);
  qlat::load_gwu_link(rbc_conf, gf);

  {
    double splaq = gf_avg_spatial_plaq(gf);
    double gplaq = gf_avg_plaq(gf);
    qmessage("spatial plaquette %.8e , plaquette %.8e \n", splaq, gplaq);
  }

  quda_begin(mpi_layout);

  Long V = geo.local_volume();
  qlat::vector<qlat::ComplexD > quda_gf;quda_gf.resize(V * 4 * 3*3);
  quda_convert_gauge(quda_gf, gf);

  const int Nsrc = 12;
  //int Nsrc = 1;
  //std::vector<qlat::FermionField4d > qlat_ff;qlat_ff.resize(Nsrc);
  //for(int i=0;i<qlat_ff.size();i++){qlat_ff[i].init(geo);}

  //std::vector<qlat::FermionField4d > qlat_fq;qlat_fq.resize(Nsrc);
  //for(int i=0;i<qlat_fq.size();i++){qlat_fq[i].init(geo);}

  quda_clover_inverter qinv(geo, QUDA_PERIODIC_T);
  //quda_clover_inverter qinv(geo, QUDA_ANTI_PERIODIC_T);
  qinv.setup_link(quda_gf.data());
  qinv.setup_clover(in.kappa, in.clover_csw);

  Propagator4d qlat_prop;qlat_prop.init(geo);
  Propagator4d qlat_src; qlat_src.init(geo);

  {

  qlat::set_zero(qlat_src);
  qlat::ComplexD* srcP = (qlat::ComplexD*) (qlat::get_data(qlat_src).data());
  for (int is = 0; is < Nsrc; is++) {
    #pragma omp parallel for
    for (Long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
      const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
      const Coordinate xg = geo.coordinate_g_from_l(xl);

      for(int dc=0;dc<12;dc++){
      if(xg[0] == sp[0] and xg[1] == sp[1] and xg[2] == sp[2] and xg[3] == sp[3] and dc == is){
        srcP[qlat_idx_4d*12*12 + dc * 12 + dc] = qlat::ComplexD(1.0, 0.0);
      }
      }
    }
  }
  }

  get_clover_prop(qinv, qlat_src, qlat_prop, in.kappa, 1e-12, 10000);

  //qlat_prop.init(geo);

  //////qlat_prop *= (2*in.kappa);

  char namew[500];
  sprintf(namew,in.Pname.c_str(), in.icfg);
  //////save_qlat_prop(namew, qlat_prop , false);
  save_gwu_prop(namew, qlat_prop);

  print_pion(qlat_prop, qlat_prop, std::string("quda ") );

  qinv.free_mem();
  quda_end();
  }

  fflush_MPI();
  qlat::Timer::display();
  qlat::end();
  return 0;
}

