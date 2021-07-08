#include "io_gwu.h"
#include "check_fun.h"
//#include "meson_contra.h"

#include "../load-select-data/compute-check-prop.h"
#include "../load-select-data/compute-chvp.h"
#include "../load-select-data/compute-meson-chvp.h"
#include "../load-select-data/compute-meson-snk-src.h"
#include "../load-select-data/compute-meson-vv-meson.h"
#include "../load-select-data/compute-meson-vv.h"
#include "../load-select-data/compute-psel-fsel-distribution.h"
#include "../load-select-data/compute-three-point-func.h"
#include "../load-select-data/compute-two-point-func.h"
////#include "../load-select-data/compute-wall-src-prop-norm-ratio.h"

int main(int argc, char* argv[])
{
  using namespace qlat;
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 1));
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 4));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 2, 16));
  size_node_list.push_back(Coordinate(1, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 4, 16));
  size_node_list.push_back(Coordinate(2, 4, 4, 16));
  size_node_list.push_back(Coordinate(4, 4, 4, 16));
  begin(&argc, &argv, size_node_list);
  display_geometry_node();


  //inputpara in;
  ////in.load_para("input.txt");
  //in.load_para(argc,argv);

  int nx,ny,nz,nt;
  nx = 24;
  ny = 24;
  nz = 24;
  nt = 64;

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  //geo.init(get_total_site(job_tag), 1); 
  geo.init(total_site, 1);

  char filename[500],rbc_conf[500];
  sprintf(filename,"/global/homes/g/genwang/cscratch/24DH/confs/rbc_2464_24IDC_330_000330");
  sprintf(rbc_conf,"/global/homes/g/genwang/cscratch/24DH/confs/ckpoint_lat.330");
  char prop_n[500],prop_s[500];
  sprintf(prop_n,"/global/homes/g/genwang/cscratch/24DH/confs/f.000330.tem.No01.prop");
  sprintf(prop_s,"/global/homes/g/genwang/cscratch/24DH/confs/f.000330.tem.No01.src");
  char prop_nP[500],prop_sP[500];
  sprintf(prop_nP,"/global/homes/g/genwang/cscratch/24DH/confs/Perodic/f.000330.tem.No00.prop");
  sprintf(prop_sP,"/global/homes/g/genwang/cscratch/24DH/confs/Perodic/f.000330.tem.No00.src");


  char lat[500];
  int traj = 330;
  int countsave = 0;
  sprintf(lat,"24DH");

  int ionum = 8;
  io_gwu io_use(geo,ionum);

  qlat::FieldM<qlat::Complex,1> noi,grid;
  noi.init(geo);grid.init(geo);
  qlat::set_zero(noi);qlat::set_zero(grid);


  std::vector<qlat::FermionField4dT<qlat::Complex> > prop_qlat;
  prop_qlat.resize(12);for(int iv=0;iv<12;iv++){prop_qlat[iv].init(geo);qlat::set_zero(prop_qlat[iv]);}

  sprintf(filename,"data/dwf_gwu/%s/rbc.%s.%06d.N%06d.S",lat,lat,traj,countsave);
  load_gwu_noi(filename,noi ,io_use);
  sprintf(filename,"data/dwf_gwu/%s/rbc.%s.%06d.N%06d.G",lat,lat,traj,countsave);
  load_gwu_noi(filename,grid,io_use);
  sprintf(filename,"data/dwf_gwu/%s/rbc.%s.%06d.N%06d.prop",lat,lat,traj,countsave);

  ////load_gwu_prop(filename,prop_qlat,io_use);
  load_gwu_prop(prop_nP,prop_qlat,io_use);

  std::vector<qlat::FermionField4dT<qlat::Complex> > prop_inv;
  prop_inv.resize(12);for(int iv=0;iv<12;iv++){prop_inv[iv].init(geo);qlat::set_zero(prop_inv[iv]);}
  qlat::FieldM<qlat::Complex,1> noi_1;
  noi_1.init(geo);
  qlat::set_zero(noi_1);
  load_gwu_noi(prop_s,noi_1,io_use);
  load_gwu_prop(prop_n,prop_inv,io_use);

  //////////Convert prop from grid to qlat
  ////Propagator4d prop;
  ////Propagator4d prop_wm;prop_wm.init(geo);
  ////prop4d_to_Fermion(prop,prop_inv,0);
  ////convert_wm_from_mspincolor(prop_wm,prop);
  ////prop4d_to_Fermion(prop_wm,prop_inv,1);
  //////////Convert prop from grid to qlat


  print_src(noi);
  MPI_Barrier(qlat::get_comm());fflush(stdout);

  print_src(noi_1);
  MPI_Barrier(qlat::get_comm());fflush(stdout);

  int countp = 0;
  double diffp = 0;
  for (long index = 0; index < noi.geo().local_volume(); ++index) {
    Coordinate xl0 = noi.geo().coordinate_from_index(index);
    Coordinate xg0 = noi.geo().coordinate_g_from_l(xl0);

    qlat::Complex tem = grid.get_elems(index)[0];
    double sum = 0.0;
    sum += std::fabs(tem.real());
    sum += std::fabs(tem.imag());
    if(sum > 1e-6){
      {
      double p0 = 0.0;
      double p1 = 0.0;
      for(int d0=0;d0<12;d0++)
      {
        qlat::Complex* s0 = (qlat::Complex*)&(prop_qlat[d0].get_elem(index));
        qlat::Complex* s1 = (qlat::Complex*)&( prop_inv[d0].get_elem(index));
        for(int d1=0;d1<12;d1++)
        {
          p0 += (s0[d1]*std::conj(s0[d1])).real();
          p1 += (s1[d1]*std::conj(s1[d1])).real();
        }
      }
      double diff = 0.0;
      if(fabs(p0) > 1e-28){diff = std::fabs((p0-p1)/p0);}
      diffp += diff;
      if(diff > 1e-13 and countp < 24)
      {
        printf("x %3d, y %3d, z %3d, t %3d, value %.5e %.5e \n",xg0[0],xg0[1],xg0[2],xg0[3],p0,p1);
        countp += 1;
      }
      }

      //for(int d0=0;d0<12;d0++){
      //  double* s0 = (double*)&(prop_qlat[d0].get_elem(index));
      //  double* s1 = (double*)&( prop_inv[d0].get_elem(index));
      //  for(int d1=0;d1<24;d1++)
      //  {
      //    double p0 = s0[d1];
      //    double p1 = s1[d1];
      //    double diff = 0.0;
      //    if(fabs(p0) > 1e-28){diff = std::fabs((p0-p1)/p0);}
      //    diffp += diff;
      //    if(diff > 1e-13 and countp < 24)
      //    {
      //      printf("x %3d, y %3d, z %3d, t %3d, value %.5e %.5e \n",xg0[0],xg0[1],xg0[2],xg0[3],s0[d1],s1[d1]);
      //      countp += 1;
      //    }
      //  }
      //}
    }

  }
  sum_all_size(&diffp,1);
  MPI_Barrier(qlat::get_comm());fflush(stdout);
  print0("==prop diff %.5e \n",diffp);

  Coordinate xg;xg[0] = 0;xg[1] = 0;xg[2] = 0;
  xg[3] = 35;
  Coordinate xg1;xg1[0] = 0;xg1[1] = 0;xg1[2] = 0;
  xg1[3] = 0;

  ////int nt = 64;
  std::vector<double > write;write.resize(2*nt);
  std::vector<double > writ0;writ0.resize(2*nt);
  for(int ti=0;ti<write.size();ti++){write[ti]=0.0;}
  for(int ti=0;ti<writ0.size();ti++){writ0[ti]=0.0;}
  get_corr_pion(prop_qlat,xg,write);
  get_corr_pion(prop_inv,xg1,writ0);
  for(int ti=0;ti<nt;ti++)
  {
    print0("pion t %3d, %.9e %.9e \n",ti, write[ti*2+0], writ0[ti*2+0]);
  }


  //GaugeField gf;
  //GaugeField gf_gwu;
  //gf.init(geo);
  //gf_gwu.init(geo);
  //load_gauge_field(gf,rbc_conf,true);
  //load_gwu_link(filename,gf_gwu);

  //diff_gauge(gf,gf_gwu);


  end();
  return 0;
}


