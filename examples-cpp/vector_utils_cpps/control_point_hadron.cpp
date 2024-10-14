#include <sys/sysinfo.h>
#include "general_funs.h"
#include "utils_io_vec.h"
#include "utils_construction.h"
#include "utils_eigen_ov.h"
#include "utils_check_fun.h"
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
  int mass_split = -1;
  in.find_para(std::string("prop_type"), prop_type);
  std::string out_vec_sink  = std::string("NONE");
  std::string out_vec_mom   = std::string("NONE");
  in.find_para(std::string("out_vec_sink"), out_vec_sink);
  in.find_para(std::string("out_vec_mom"), out_vec_mom);
  in.find_para(std::string("mass_split"), mass_split);

  int check_prop_norm = 0;
  in.find_para(std::string("check_prop_norm"), check_prop_norm);

  int prop_smear_inv_factor = 1; ////1, if use modified factors within inverstion ; 0 --> normal factors
  in.find_para(std::string("prop_smear_inv_factor"), prop_smear_inv_factor);
  int save_full_vec = 0;
  in.find_para(std::string("save_full_vec"), save_full_vec);

  int sleep_for_final = 0;
  in.find_para(std::string("sleep_for_final"), sleep_for_final);

  int do_point_sink   = 1;
  in.find_para(std::string("do_point_sink"), do_point_sink);

  std::string mom_shift_ = std::string("0 0 0 0");
  in.find_para(std::string("mom_shift"), mom_shift_);
  std::vector<Coordinate > mom_shiftL;
  {
    std::vector<std::string > Li = stringtolist(mom_shift_);
    Qassert(Li.size() % 4 == 0);
    const int Nmom = Li.size() / 4;
    mom_shiftL.resize(Nmom);
    for(int momi=0;momi<Nmom;momi++)
    {
      for(int i=0;i<4;i++){mom_shiftL[momi][i] = stringtonum(Li[momi*4 + i]);}
    }
  }
  //const Coordinate mom_shift = string_to_Coordinate(mom_shift_);

  std::string mom_smear_ = std::string("NONE");
  in.find_para(std::string("mom_smear"), mom_smear_);
  CoordinateD mom_smear;for(int i=0;i<4;i++){mom_smear[i] = 0;}
  if(mom_smear_ != std::string("NONE"))
  {
    std::vector<std::string > Li = stringtolist(mom_smear_);
    Qassert(Li.size() == 4);
    for(int i=0;i<4;i++){mom_smear[i] = stringtodouble(Li[i]);}
  }

  int icfg  = in.icfg;

  int n_vec = in.nvec;
  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site); 
  fflush_MPI();

  {


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
  int nmass_group = in.nmass;qassert(nmass_group > 0);
  if(mass_split != -1 and mass_split < nmass_group){
    nmass_group = mass_split;
  }
  std::vector<Long > mass_jobA = job_create(in.nmass, nmass_group);

  ////qprop tmpa;tmpa.init(geo);
  ////print0("vol %ld %ld \n", geo.local_volume(), Long(qlat::get_data_size(tmpa)));

  ////===load eigen
  print_time();
  fflush_MPI();
  ////fft_desc_basic fd(geo);
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  eigen_ov ei(geo, n_vec, in.bini, nmass_group + 1);

  if(src_step != sink_step)
  {
    do_point_sink = 0;
    //Qassert(do_point_sink == 0);// eigen system too large!
  }

  fflush_MPI();
  int mode_sm = 0;
  char ename[500];
  sprintf(ename, in.Ename.c_str(),icfg);

  /////ei.load_eigen(std::string(ename));
  {
    char enamev[600];
    ////sprintf(ename, "%s", ov_evecname);
    sprintf(enamev,"%s.eigvals", ename);
    print0("Vector File name: %s \n", ename );
    print0("Values File name: %s \n", enamev);
    //////Load eigen values
    const double kappa= 0.2;
    const double rho_tem = 4 - 1.0/(2*kappa);
    const double eigenerror = 1e-11;
    const int nini = 0;
    const int checknorm = 1;
    ei.load_eivals(std::string(enamev), rho_tem, eigenerror, nini);

    if(src_step == sink_step){
    ei.load_eigen_Mvec_smear(std::string(ename), gf, nini, checknorm,
      src_width , src_step , 0.0, 0);}

    if(src_step != sink_step){
    ei.load_eigen_Mvec_smear(std::string(ename), gf, nini, checknorm,
      src_width , src_step , sink_width, sink_step);}

    if(src_step != 0){
      mode_sm = 2;
    }
  }

  //if(src_step != 0)
  //{
  //  sprintf(ename, in.Ename_Sm.c_str(),icfg);
  //  ei.smear_eigen(std::string(ename), gf, src_width, src_step, mom_smear);
  //  mode_sm = 2;
  //}
  print_time();

  qnoi noi;noi.init(geo);
  Propagator4dT<Ftype > tmp;tmp.init(geo);

  for(LInt jobi=0;jobi < mass_jobA.size()/2; jobi++)
  {
    Long bji = mass_jobA[jobi*2 + 0]; Long bcut = mass_jobA[jobi*2+1];

    std::string outputG = in.output;
    std::string out_vec_sinkG = out_vec_sink;
    std::string out_vec_momG  = out_vec_mom ;

    if(nmass_group != in.nmass){
      if(outputG != std::string("NONE")){
        outputG = ssprintf("%s.mgroup%02d", outputG.c_str(), jobi);
      }

      if(out_vec_sinkG != std::string("NONE")){
        out_vec_sinkG = ssprintf("%s.mgroup%02d", out_vec_sinkG.c_str(), jobi);
      }

      if(out_vec_momG != std::string("NONE")){
        out_vec_momG = ssprintf("%s.mgroup%02d", out_vec_momG.c_str(), jobi);
      }
    }

    momentum_dat mdat(geo, in.mom_cut, mom_shiftL);
    print_mem_info("momentum dat");

    std::vector<double> massL;massL.resize(0);
    for(int mi=0;mi<bcut;mi++){
      massL.push_back(in.masses[mi + bji]);
    }
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
    if(outputG != std::string("NONE")){
      std::string ktem(key_T);
      std::string dtem(dimN);
      res.create_dat(ktem, dtem);
      res.print_info();
    }
    
    /////===load noise and prop
    char names[450],namep[500], names_vec[450], names_mom[450];
    std::vector<qprop > FpropV;FpropV.resize(0);FpropV.resize(bcut);
    lms_para<Complexq > srcI;/////buffers and parameters for lms
    print_time();
    for(int si = 0; si < in.nsource; si++)
    {
      if(jobi == mass_jobA.size()/2 - 1 and si == in.nsource -1 and sleep_for_final == 1)
      {
        print0("SLEEPPING!!!\n");
        sleep(30);// some weild writting error at final loop
      }
      fflush_MPI();

      if(out_vec_sinkG != std::string("NONE") and ckpoint == 1){
        const size_t vol = size_t(fd.nx) * fd.ny * fd.nz * fd.nt;
        int flag_do_job = 0;
        sprintf(names, out_vec_sinkG.c_str(), icfg, si);
        sprintf(namep, "%s.pt.zero.vec", names);
        if(get_file_size_MPI(namep) > vol * 32 * bcut * 8){ 
          flag_do_job += 1;
        }   
        sprintf(namep, "%s.sm.zero.vec", names);
        if(get_file_size_MPI(namep) > vol * 32 * bcut * 8){ 
          flag_do_job += 1;
        } 

        int flag_need = 0;
        if(do_point_sink == 1){flag_need = 2;}
        if(do_point_sink == 0){flag_need = 1;}
        Qassert(flag_need != 0);

        if(flag_do_job == flag_need){
          print0("Pass %s \n", names);
          continue ;
        }   
      }

      sprintf(names, in.srcN[si].c_str(),icfg);
      if(get_file_size_MPI(names) == 0){
        print0("Src pass %s \n", names );
        continue;
      }

      load_gwu_noi(names, noi);
      print0("%s \n", names);

      for(int im=0;im<bcut;im++)
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

        if(check_prop_norm == 1){
          print0("mgroup %2d, source %3d, mass %3d ", int(jobi), si, im);
          print_norm2(FpropV[im]);
        }
        ////double sum = check_sum_prop(FpropV[im]);
        ////print0("===checksum %s %.8e \n", namep, sum);
      }
      /////===load noise and prop

      bool save_vecs_vec = false;
      bool save_vecs_mom = false;
      if(out_vec_sinkG != std::string("NONE")){
        sprintf(names_vec, out_vec_sinkG.c_str(), icfg, si);
        save_vecs_vec = true;
      }
      if(out_vec_momG != std::string("NONE")){
        sprintf(names_mom, out_vec_momG.c_str(), icfg, si);
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
      srcI.check_prop_norm = check_prop_norm;
      if(save_vecs_vec){
        sprintf(namep, "%s.pt.zero.vec", names_vec);
        srcI.name_zero_vecs = std::string(namep);
      }

      if(save_vecs_mom){
        sprintf(namep, "%s.pt", names_mom);
        srcI.name_mom_vecs = std::string(namep);
      }

      if(outputG == std::string("NONE")){
        srcI.save_zero_corr = 0;
      }

      // point sink
      if(do_point_sink == 1){
        point_corr(noi, FpropV, massL, ei, fd, res, srcI, mdat, 1);
      }

      if(sink_step != 0 or do_point_sink == 0)
      {
        srcI.lms      = in.lms;
        srcI.combineT = in.combineT;
        srcI.mode_eig_sm = 3;

        //replace sink with the modified eigen
        if(src_step != sink_step){
          srcI.mode_eig_sm = 2;
        }

        for(int im=0;im<bcut;im++){
        smear_propagator_gwu_convension(FpropV[im], gf, sink_width, sink_step, mom_smear);
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

    if(outputG != std::string("NONE")){
      std::string Name = ssprintf(outputG.c_str(), icfg);
      res.INFOA.push_back(INFO_Mass);
      res.print_info();
      res.write_dat(Name);
    }

    srcI.free_buf();
  }

  if(sleep_for_final == 1){
    print0("SLEEPPING!!!\n");
    sleep(30);// some weild writting error at final loop
  }
  fflush_MPI();

  ei.clear_GPU_mem(1);

  }

  return end_Lat();
}

