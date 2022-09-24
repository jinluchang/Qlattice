#include <sys/sysinfo.h>
#include <unistd.h>

#include <qutils/vector.h>
#include "general_funs.h"
#include "utils_io_vec.h"
#include "utils_construction.h"
#include "utils_stagger_contractions.h"
#include "quda_para.h"
#include "utils_quda_inverter.h"
#include "utils_gaugefield.h"

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

  Coordinate sp = string_to_Coordinate(in.paraI);

  int mpi_layout[4]={0,0,0,0};
  qlat::GeometryNode geon = qlat::get_geometry_node();for(int i=0;i<4;i++){mpi_layout[i] = geon.size_node[i];}

  qlat::GaugeField gf;gf.init(geo);

  if(in.Link_name != "NONE"){
    char rbc_conf[500];
    sprintf(rbc_conf,in.Link_name.c_str(), in.icfg);
    qlat::load_gwu_link(rbc_conf, gf);
  }else{
    random_link(gf, in.seed);
    qlat::gf_ape_smear(gf, gf, 0.125, in.hyp);
  }

  {
    double splaq = gf_avg_spatial_plaq(gf);
    double gplaq = gf_avg_plaq(gf);
    print0("spatial plaquette %.8e , plaquette %.8e \n", splaq, gplaq);
  }

  const long Dim = 3;
  const int  Nsrc = Dim;
  std::vector<colorFD > qlat_cf;qlat_cf.resize(Nsrc);
  for(int i=0;i<qlat_cf.size();i++){qlat_cf[i].init(geo);}

  //{
  //EigenV corr;fft_desc_basic fd(geo);
  //cf_simple_pion(qlat_cf, qlat_cf, corr, fd, 1, true, 4.0);
  //}

  quda_begin(mpi_layout);

  long V = geo.local_volume();
  qlat::vector<qlat::Complex > quda_gf;quda_gf.resize(V * 4 * 3*3);
  quda_convert_gauge(quda_gf, gf);

  quda_inverter qinv(geo, QUDA_PERIODIC_T);

  qinv.setup_link(quda_gf, 1);
  //qinv.setup_stagger(in.fermion_mass, 1e-10);
  qinv.setup_stagger();
  qinv.setup_eigen(in.fermion_mass, in.nvec, 1e-13);
  qinv.check_residue = 1;

  ////int setup_stagger = 1;
  ////===Invertion part
  {
  //// Vector construct START
  ////-----------------------------------------------------------------------------------
  //std::vector<quda::ColorSpinorField *> quda_in(Nsrc);
  //std::vector<quda::ColorSpinorField *> quda_out(Nsrc);
  //quda::ColorSpinorParam cs_param;
  /////constructWilsonTestSpinorParam(&cs_param, &inv_param, &gauge_param);

  for (int i = 0; i < Nsrc; i++) {
    // Populate the host spinor with random numbers.
    //quda_in[i]->Source(QUDA_RANDOM_SOURCE);

    //////set point src at zero
    qlat::Complex* res = (qlat::Complex*) (qinv.csrc->V());
    long Vh = V / 2;
    #pragma omp parallel for
    for (long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
      const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
      int quda_idx = eo * Vh + qlat_idx_4d / 2;

      for(int dc=0;dc<Dim;dc++){
      if(xg[0] == sp[0] and xg[1] == sp[1] and xg[2] == sp[2] and xg[3] == sp[3] and dc == i){
        res[quda_idx*Dim + dc] = qlat::Complex(1.0, 0.0);
      }
      else{
        res[quda_idx*Dim + dc] = qlat::Complex(0.0, 0.0);
      }
      }
    }
    //////set point src at zero

    qinv.do_inv(qinv.cres->V(), qinv.csrc->V(), in.fermion_mass, in.cg_err, in.niter);


    quda_cf_to_qlat_cf(qlat_cf[i], (qlat::Complex*) qinv.cres->V());
  }

  }

  EigenV corr;fft_desc_basic fd(geo);
  cf_simple_pion(qlat_cf, qlat_cf, corr, fd, 1, true, 4.0);

  //Propagator4d qlat_prop;
  //Fermion_to_prop4d(qlat_prop, qlat_ff);
  ////qlat_prop.init(geo);
  ////////qlat_prop *= (2*in.kappa);
  //char namew[500];
  //sprintf(namew,in.Pname.c_str(), in.icfg);
  ////////save_qlat_prop(namew, qlat_prop , false);
  //save_gwu_prop(namew, qlat_prop);


  //host_free(milc_fatlink);

  qinv.free_mem();

  quda_end();

  //qlat::vector<double > va;
  //va.resize(100);

  //qlat::vector_acc<double > vb;
  //vb.resize(100);
  fflush_MPI();
  qlat::Timer::display();
  qlat::end();
  return 0;
}

