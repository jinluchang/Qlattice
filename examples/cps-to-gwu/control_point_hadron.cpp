#include <sys/sysinfo.h>
#include "general_funs.h"
#include "io_vec.h"
#include "utils_construction.h"
#include "utils_eigensys.h"
#include "check_fun.h"
#include "utils_smear_vecs.h"
#include "utils_lms_funs.h"


int main(int argc, char* argv[])
{
  using namespace qlat;

  inputpara in;
  //int mode_dis = 0;
  int mode_dis = -1;
  begin_Lat(&argc, &argv, in, mode_dis);

  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;

  int icfg  = in.icfg;
  int ionum = in.ionum;

  int n_vec = in.nvec;
  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 
  fflush_MPI();

  io_vec io_use(geo,ionum);

  double length = (geo.local_volume()*pow(0.5,30))*12*sizeof(Complexq);
  size_t freeM = 0;size_t totalM = 0;
  #ifdef QLAT_USE_ACC
  cudaMemGetInfo(&freeM,&totalM);
  #endif
  double freeD = freeM*pow(0.5,30);double totalD = totalM*pow(0.5,30); 
  struct sysinfo s_info;
  sysinfo(&s_info);
  print0("Eign system vector size %.3e GB, Total %.3e GB; CPU free %.3e GB, total %.3e GB; GPU free %.3e GB, total %.3e GB. \n"
          ,length,n_vec*length, s_info.freeram*pow(0.5,30),s_info.totalram*pow(0.5,30),freeD,totalD);


  /////========load links
  double width = 0.0; int step = 0;
  get_smear_para(in.src_smear_para, width, step);
  GaugeField gf;GaugeField gfD;
  if(step > 0){
    gf.init(geo);
    char rbc_conf[500];
    sprintf(rbc_conf,in.Link_name.c_str(), icfg);
    load_gwu_link(rbc_conf, gf);
    set_left_expanded_gauge_field(gfD, gf);
  }
  /////========load links
  int nmass = in.nmass;
  std::vector<double> massL = in.masses;

  ////===load eigen
  fflush_MPI();
  fft_desc_basic fd(geo);
  eigen_ov ei(fd, n_vec, in.bini, in.nmass + 1);

  fflush_MPI();
  int mode_sm = 0;
  //ei.load_eigen(icfg, in.Ename,  io_use, 1, 0.2, in.Eerr);
  ei.load_eigen(icfg, in.Ename,  io_use);
  if(in.Ename_Sm != std::string("NONE"))
  {
    ei.load_eigen_Mvec(icfg, in.Ename_Sm, io_use, 1);
    mode_sm = 2;
  }
  //ei.random_eigen();ei.random_eigen(1);mode_sm = 2;

  ei.initialize_mass(massL, 12);
  ei.print_info();
  print0("Low eigen done. \n");
  ////===load eigen


  ////size_t Nvol = geo.local_volume();

  char key_T[1000], dimN[1000];
  sprintf(key_T, "%d   %d   %d  %d  %d %d", in.nsource, 3, 32, int(massL.size()), in.nt,2);
  sprintf(dimN , "src LHF operator masses nt complex");

  std::string ktem(key_T);
  std::string dtem(dimN);
  corr_dat res(ktem, dtem);
  res.print_info();
  
  /////===load noise and prop
  char names[450],namep[500];
  for(int si = 0; si < in.nsource; si++)
  {
    sprintf(names, in.srcN[si].c_str(),icfg);
    std::vector<qnoi > noi;noi.resize(1);
    load_gwu_noi(names, noi[0], io_use);

    std::vector<qprop > FpropV;FpropV.resize(nmass);
    for(int im=0;im<nmass;im++)
    {
      sprintf(names, in.propN[si].c_str(),icfg);
      sprintf(namep, "%s.m%8.6f", names, massL[im]);
      load_gwu_prop(namep, FpropV[im], io_use);
      //if(step > 0)
      //{
      //  smear_propagator_gwu_convension(propS[im], gfD, width, step);
      //}
    }
    /////===load noise and prop

    point_corr(noi, FpropV, massL, ei, fd, res, mode_sm, 1, in.lms);
  }

  sprintf(names, in.output.c_str(),icfg);
  res.print_info();
  res.write_dat(names);


  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

