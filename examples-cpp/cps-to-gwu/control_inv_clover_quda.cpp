#include <sys/sysinfo.h>
#include <unistd.h>

#include <qutils/vector.h>
#include "general_funs.h"
#include "utils_io_vec.h"
#include "utils_construction.h"
#include "quda_para.h"
#include "utils_quda_inverter.h"

int main(int argc, char* argv[])
{
  using namespace qlat;
  inputpara in;int mode_dis = -1; 
  begin_Lat(&argc, &argv, in, mode_dis);

  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;
  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 
  fflush_MPI();

  std::vector<int > sp;sp.resize(4);
  for(int i=0;i<4;i++){sp[i] = 0;}
  if(in.paraI != "None"){
    std::vector<std::string > Li = stringtolist(in.paraI);
    qassert(sp.size() == Li.size());
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
    print0("spatial plaquette %.8e , plaquette %.8e \n", splaq, gplaq);
  }

  quda_begin(mpi_layout);

  long V = geo.local_volume();
  qlat::vector_acc<qlat::Complex > quda_gf;quda_gf.resize(V * 4 * 3*3);
  quda_convert_gauge(quda_gf, gf);

  int Nsrc = 12;
  //int Nsrc = 1;
  std::vector<qlat::FermionField4d > qlat_ff;qlat_ff.resize(Nsrc);
  for(int i=0;i<qlat_ff.size();i++){qlat_ff[i].init(geo);}

  quda_inverter qinv(geo, quda_gf);
  qinv.setup_clover(in.kappa, in.clover_csw);
  qinv.check_residue = 1;

  ////////quda containers

  //////===Invertion part
  {
  //std::vector<double> time(Nsrc);
  //std::vector<double> gflops(Nsrc);
  //std::vector<int> iter(Nsrc);

  qinv.setup_eigen(in.nvec);

  for (int i = 0; i < Nsrc; i++) {
    //////set point src at zero
    qlat::Complex* res = (qlat::Complex*) (qinv.csrc->V());
    long Vh = V / 2;
    #pragma omp parallel for
    for (long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
      const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
      int quda_idx = eo * Vh + qlat_idx_4d / 2;

      for(int dc=0;dc<12;dc++){
      if(xg[0] == sp[0] and xg[1] == sp[1] and xg[2] == sp[2] and xg[3] == sp[3] and dc == i){
        res[quda_idx*12 + dc] = qlat::Complex(1.0, 0.0);
      }
      else{
        res[quda_idx*12 + dc] = qlat::Complex(0.0, 0.0);
      }
      }
    }
    //////set point src at zero

    ///quda_out[i] = quda::ColorSpinorField::Create(cs_param);
    qinv.do_inv(qinv.cres->V(), qinv.csrc->V());

    //time[i] = quda_inverter.inv_param.secs;
    //gflops[i] = quda_inverter.inv_param.gflops / quda_inverter.inv_param.secs;
    //iter[i] = quda_inverter.inv_param.iter;
    printfQuda("Done: %i iter / %g secs = %g Gflops\n\n", qinv.inv_param.iter, qinv.inv_param.secs,
               qinv.inv_param.gflops / qinv.inv_param.secs);

    quda_ff_to_Ffield4d(qlat_ff[i], (qlat::Complex*) qinv.cres->V());

  }

  }

  Propagator4d qlat_prop;
  Fermion_to_prop4d(qlat_prop, qlat_ff);
  //qlat_prop.init(geo);

  //////qlat_prop *= (2*in.kappa);

  char namew[500];
  sprintf(namew,in.Pname.c_str(), in.icfg);
  //////save_qlat_prop(namew, qlat_prop , false);
  save_gwu_prop(namew, qlat_prop);

  print_meson(qlat_prop, qlat_prop,   "quda ");


  qinv.free_mem();

  quda_end();

  fflush_MPI();
  qlat::Timer::display();
  qlat::end();
  return 0;
}

