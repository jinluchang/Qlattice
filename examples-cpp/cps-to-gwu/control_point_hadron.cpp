#include <sys/sysinfo.h>
#include "general_funs.h"
#include "utils_io_vec.h"
#include "utils_construction.h"
#include "utils_eigensys.h"
#include "check_fun.h"
#include "utils_smear_vecs.h"
#include "utils_lms_funs.h"

int main(int argc, char* argv[])
{
  using namespace qlat;

  inputpara in;begin_Lat(&argc, &argv, in);

  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;

  int ckpoint = 1;
  in.find_para(std::string("ckpoint"), ckpoint);

  int icfg  = in.icfg;
  int ionum = in.ionum;

  int n_vec = in.nvec;
  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 
  fflush_MPI();

  momentum_dat mdat(geo, in.mom_cut);

  {
  io_vec io_use(geo,ionum);

  print_mem_info();

  /////========load links
  double src_width = 0.0; int src_step = 0;
  double sink_width = 0.0; int sink_step = 0;
  get_smear_para(in.src_smear_para , src_width, src_step);
  get_smear_para(in.sink_smear_para, sink_width, sink_step);
  GaugeField gf;////GaugeField gfD;
  if(src_step != 0 or sink_step != 0){
    gf.init(geo);
    char rbc_conf[500];
    sprintf(rbc_conf,in.Link_name.c_str(), icfg);
    load_gwu_link(rbc_conf, gf);
    /////set_left_expanded_gauge_field(gfD, gf);
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
  char ename[500];
  //char enamew[550];
  sprintf(ename, in.Ename.c_str(),icfg);
  ei.load_eigen(std::string(ename),  io_use);

  //sprintf(enamew, "%s.w", ename);
  //ei.save_eigen_Mvec(std::string(enamew), io_use);

  if(src_step != 0){
    sprintf(ename, in.Ename_Sm.c_str(),icfg);
    ei.smear_eigen(std::string(ename), io_use, gf, src_width, src_step);
    mode_sm = 2;
  }

  //if(in.Ename_Sm != std::string("NONE"))
  //{
  //  sprintf(ename, in.Ename_Sm.c_str(),icfg);
  //  ei.load_eigen_Mvec(std::string(ename), io_use, 1);
  //  if(src_step != 0){mode_sm = 2;}
  //  ////mode_sm = 2;
  //}
  //ei.random_eigen();ei.random_eigen(1);mode_sm = 2;

  ei.initialize_mass(massL, 12);
  ei.print_info();
  print0("Low eigen done. \n");
  ////===load eigen


  ////size_t Nvol = geo.local_volume();

  char key_T[1000], dimN[1000];
  
  if(sink_step!=0){sprintf(key_T, "%d   %d    %d   %d  %d  %d %d", in.nsource, 2, 3, 32, int(massL.size()), in.nt,2);}
  else{sprintf(key_T, "%d   %d    %d   %d  %d  %d %d", in.nsource, 1, 3, 32, int(massL.size()), in.nt,2);}
  sprintf(dimN , "src  sm  LHF operator masses nt complex");

  std::string  INFO_Mass = mass_to_string(massL);

  std::string ktem(key_T);
  std::string dtem(dimN);
  corr_dat<Ftype > res(ktem, dtem);
  res.print_info();
  
  /////===load noise and prop
  char names[450],namep[500];
  std::vector<qprop > FpropV;FpropV.resize(nmass);
  Propagator4dT<Complexq > tmp;tmp.init(geo);
  for(int si = 0; si < in.nsource; si++)
  {
    sprintf(names, in.srcN[si].c_str(),icfg);
    //std::vector<qnoi > noi;noi.resize(1);
    qnoi noi;noi.init(geo);
    load_gwu_noi(names, noi, io_use);
    print0("%s \n", names);

    for(int im=0;im<nmass;im++)
    {
      sprintf(names, in.propN[si].c_str(),icfg);
      sprintf(namep, "%s.m%8.6f", names, massL[im]);

      load_gwu_prop(namep, tmp);
      prop4d_to_qprop(FpropV[im], tmp);
      //qprop_to_prop4d(tmp, FpropV[im]);
      //prop4d_to_qprop(FpropV[im], tmp);

      //load_gwu_prop(namep, FpropV[im], io_use);
      ////smear_propagator_gwu_convension(FpropV[im], gf, width, step);
    }
    /////===load noise and prop

    bool save_vecs_lms = false;
    if(in.output_vec != std::string("NONE")){
      sprintf(names, in.output_vec.c_str(), icfg, si);
      save_vecs_lms = true;
    }
    lms_para<Complexq > srcI;srcI.init();
    srcI.ionum = in.ionum;
    srcI.do_all_low = in.do_all_low;
    srcI.lms      = in.lms;
    srcI.combineT =  in.combineT;
    srcI.mode_eig_sm = mode_sm;
    srcI.SRC_PROP_WITH_LOW = in.SRC_PROP_WITH_LOW;
    srcI.INFOA.push_back(INFO_Mass);
    srcI.mom_cut = in.mom_cut;
    srcI.ckpoint = ckpoint;
    if(save_vecs_lms){
      sprintf(namep, "%s.pt", names);
      srcI.name_mom_vecs = std::string(namep);
      sprintf(namep, "%s.pt.zero.vec", names);
      srcI.name_zero_vecs = std::string(namep);
    }

    if(in.output == std::string("NONE")){
      srcI.save_zero_corr = 0;
    }

    point_corr(noi, FpropV, massL, ei, fd, res, srcI, mdat, 1);

    if(sink_step != 0){
      srcI.lms      = in.lms;
      srcI.combineT = in.combineT;
      srcI.mode_eig_sm = 3;

      for(int im=0;im<nmass;im++){
      smear_propagator_gwu_convension(FpropV[im], gf, sink_width, sink_step);
      //copy_noise_to_prop(FpropV[im], prop4d, 1);

      //smear_propagator_gwu_convension(prop4d, gf, width, step);
      //smear_propagator_gwu_convension(FpropV[im], gf, width, step);

      //copy_noise_to_prop(FpropV[im], prop4dS, 1);
      //diff_prop(prop4dS, prop4d);
      }

      if(save_vecs_lms){
        sprintf(namep, "%s.sm", names);
        srcI.name_mom_vecs  = std::string(namep); 
        sprintf(namep, "%s.sm.zero.vec", names);
        srcI.name_zero_vecs = std::string(namep);
      }

      point_corr(noi, FpropV, massL, ei, fd, res, srcI, mdat, 1);
    }
  }

  if(in.output != std::string("NONE")){
    sprintf(names, in.output.c_str(),icfg);
    res.INFOA.push_back(INFO_Mass);
    res.print_info();
    res.write_dat(names);
  }

  }

  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

