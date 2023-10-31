#include <sys/sysinfo.h>
#include "general_funs.h"
#include "utils_io_vec.h"
#include "utils_construction.h"
#include "utils_eigen_ov.sh"
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

  int prop_type = 1;
  in.find_para(std::string("prop_type"), prop_type);
  std::string out_vec_sink  = std::string("NONE");
  std::string out_vec_mom   = std::string("NONE");
  in.find_para(std::string("out_vec_sink"), out_vec_sink);
  in.find_para(std::string("out_vec_mom"), out_vec_mom);

  int prop_smear_inv_factor = 1; ////1, if use modified factors within inverstion ; 0 --> normal factors
  in.find_para(std::string("prop_smear_inv_factor"), prop_smear_inv_factor);
  int save_full_vec = 0;
  in.find_para(std::string("save_full_vec"), save_full_vec);

  int icfg  = in.icfg;

  int n_vec = in.nvec;
  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 
  fflush_MPI();

  {

  momentum_dat mdat(geo, in.mom_cut);
  print_mem_info("momentum dat");

  print_mem_info("io_vec");

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
    if(in.anti_peri == 1){twist_boundary_at_boundary(gf, -0.5, 3);}
    /////set_left_expanded_gauge_field(gfD, gf);
  }
  /////========load links
  int nmass = in.nmass;qassert(nmass > 0);
  std::vector<double> massL = in.masses;

  ////qprop tmpa;tmpa.init(geo);
  ////print0("vol %ld %ld \n", geo.local_volume(), long(qlat::get_data_size(tmpa)));

  ////===load eigen
  print_time();
  fflush_MPI();
  ////fft_desc_basic fd(geo);
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  eigen_ov ei(geo, n_vec, in.bini, in.nmass + 1);

  fflush_MPI();
  int mode_sm = 0;
  char ename[500];
  sprintf(ename, in.Ename.c_str(),icfg);
  ei.load_eigen(std::string(ename));


  if(src_step != 0){
    sprintf(ename, in.Ename_Sm.c_str(),icfg);
    ei.smear_eigen(std::string(ename), gf, src_width, src_step);
    mode_sm = 2;
  }
  print_time();

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

  corr_dat<Ftype > res(std::string(""));
  if(in.output != std::string("NONE")){
    std::string ktem(key_T);
    std::string dtem(dimN);
    res.create_dat(ktem, dtem);
    res.print_info();
  }
  
  /////===load noise and prop
  char names[450],namep[500], names_vec[450], names_mom[450];
  std::vector<qprop > FpropV;FpropV.resize(nmass);
  Propagator4dT<Ftype > tmp;tmp.init(geo);
  lms_para<Complexq > srcI;/////buffers and parameters for lms
  print_time();
  for(int si = 0; si < in.nsource; si++)
  {
    if(out_vec_sink != std::string("NONE") and ckpoint == 1){
      const size_t vol = size_t(fd.nx) * fd.ny * fd.nz * fd.nt;
      int flag_do_job = 0;
      sprintf(names, out_vec_sink.c_str(), icfg, si);
      sprintf(namep, "%s.pt.zero.vec", names);
      if(get_file_size_MPI(namep) > vol * 32 * nmass * 8){ 
        flag_do_job += 1;
      }   
      sprintf(namep, "%s.sm.zero.vec", names);
      if(get_file_size_MPI(namep) > vol * 32 * nmass * 8){ 
        flag_do_job += 1;
      } 
      if(flag_do_job == 2){
        print0("Pass %s \n", names);
        continue ;
      }   
    }

    sprintf(names, in.srcN[si].c_str(),icfg);
    //std::vector<qnoi > noi;noi.resize(1);
    if(get_file_size_MPI(names) == 0){
      print0("Src pass %s \n", names );
      continue;
    }

    qnoi noi;noi.init(geo);
    load_gwu_noi(names, noi);
    print0("%s \n", names);

    for(int im=0;im<nmass;im++)
    {
      sprintf(names, in.propN[si].c_str(),icfg);
      sprintf(namep, "%s.m%8.6f", names, massL[im]);

      if(prop_type == 0){load_gwu_prop( namep, tmp);}
      if(prop_type == 1){load_qlat_prop(namep, tmp);}

      if(src_width < 0 and prop_smear_inv_factor == 0){
        const Complexq sm_factor = std::pow( 1.0 * (-1.0*src_width), 3);
        prop4D_factor(tmp, sm_factor);
      }

      prop4d_to_qprop(FpropV[im], tmp);
      ////double sum = check_sum_prop(FpropV[im]);
      ////print0("===checksum %s %.8e \n", namep, sum);
    }
    /////===load noise and prop

    bool save_vecs_vec = false;
    bool save_vecs_mom = false;
    if(out_vec_sink != std::string("NONE")){
      sprintf(names_vec, out_vec_sink.c_str(), icfg, si);
      save_vecs_vec = true;
    }
    if(out_vec_mom != std::string("NONE")){
      sprintf(names_mom, out_vec_mom.c_str(), icfg, si);
      save_vecs_mom = true;
    }

    srcI.init();
    srcI.do_all_low = in.do_all_low;
    srcI.lms      = in.lms;
    srcI.combineT =  in.combineT;
    srcI.mode_eig_sm = mode_sm;
    srcI.SRC_PROP_WITH_LOW = in.SRC_PROP_WITH_LOW;
    srcI.INFOA.push_back(INFO_Mass);
    srcI.mom_cut = in.mom_cut;
    srcI.ckpoint = ckpoint;
    srcI.save_full_vec = save_full_vec;
    if(save_vecs_vec){
      sprintf(namep, "%s.pt.zero.vec", names_vec);
      srcI.name_zero_vecs = std::string(namep);
    }

    if(save_vecs_mom){
      sprintf(namep, "%s.pt", names_mom);
      srcI.name_mom_vecs = std::string(namep);
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

      if(save_vecs_vec){
        sprintf(namep, "%s.sm.zero.vec", names_vec);
        srcI.name_zero_vecs = std::string(namep);
      }

      if(save_vecs_mom){
        sprintf(namep, "%s.sm", names_mom);
        srcI.name_mom_vecs = std::string(namep);
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

  srcI.free_buf();
  ei.clear_GPU_mem(1);

  }

  return end_Lat();
}

