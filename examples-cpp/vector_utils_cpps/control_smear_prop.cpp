#include <sys/sysinfo.h>
#include "utils_io_vec.h"
#include "general_funs.h"
#include "utils_smear_vecs.h"
#include "utils_check_fun.h"

int main(int argc, char* argv[])
{
  using namespace qlat;
  {

  inputpara in;begin_Lat(&argc, &argv, in);

  int icfg  = in.icfg;
  //int ionum = in.ionum;

  ////int vini  = 0;
  //int n_vec = in.nvec;

  omp_set_num_threads(omp_get_max_threads());
  qmessage("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());

  Coordinate total_site = Coordinate(in.nx, in.ny, in.nz, in.nt);
  Geometry geo;
  geo.init(total_site); 
  fflush_MPI();

  char rbc_conf[500],prop_name[500],namep[500];
  GaugeField gf;
  //GaugeFieldT<Complexq> gf;


  gf.init(geo);

  //Propagator4d propS;
  //Propagator4d prop_s0;
  //Propagator4d prop_s1;
  Propagator4dT<Ftype > propS;
  Propagator4dT<Ftype > prop_s0;
  Propagator4dT<Ftype > prop_s1;

  propS.init(geo);
  prop_s0.init(geo);prop_s1.init(geo);

  sprintf(rbc_conf,in.Ename.c_str(), icfg);
  /////load_gauge_field(gf_gwu,rbc_conf,true);
  load_gwu_link(rbc_conf, gf);

  if(in.nx == 24)twist_boundary_at_boundary(gf, -0.5, 3 );
  //if(nx == 24)twist_boundary_at_boundary(gf, EIGEN_PI/1.0, 3 );
  //qacc_for(isp, Long(geo1.local_volume_expanded()), {
  //  auto* res = gfF.get_elem(isp).p;
  //  auto* src =  gf.get_elem(isp).p;
  //  for(int m=0;m<3*3*4;m++)res[isp*3*3*4 + m] = src[isp*3*3*4 + m];
  //});
  //gfF = gf1;


  sprintf(prop_name,in.Pname.c_str(), icfg);
  load_gwu_prop(prop_name, propS);

  qmessage("%s \n",prop_name);

  if(in.paraI != "None"){
    std::vector<std::string > Li = stringtolist(in.paraI);
    qmessage("Li %s, size %d \n", in.paraI.c_str(),int(Li.size()) );
    fflush_MPI();

    for(int si=0;si<1;si++)
    {
      int nsmear   = stringtonum(   Li[si*2+0]);
      double width = stringtodouble(Li[si*2+1]);
      qmessage("sn%03dsk%6.4f \n", nsmear, width);

      smear_propagator_gwu_convension(propS, gf, width, nsmear);

      //GaugeField gf1;
      ////GaugeFieldT<Complexq> gf1;
      /////gf1.init(geo);
      //set_left_expanded_gauge_field(gf1, gf);

      //smear_propagator_gpu4(propS, gf1, width, nsmear);
      qlat::WilsonMatrixT<Ftype >& v0 =  propS.get_elem(0);
      qmessage("check %.3e %.3e \n", v0(0,0).real(), v0(0,0).imag() );


      sprintf(namep, "%s.sn%03dsk%6.4f", prop_name, nsmear, width);
      qmessage("%s \n",namep);
      load_gwu_prop(namep, prop_s1);
      fflush_MPI();

      diff_prop(propS,prop_s1, 1e-7);
      fflush_MPI();
    

    }
  }

  return end_Lat();
  }
}



