#include <sys/sysinfo.h>
#include "general_funs.h"
#include "utils_io_vec.h"
#include "utils_construction.h"
#include "utils_eigensys.h"
#include "check_fun.h"
#include "utils_smear_vecs.h"

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

  int nmass = in.nmass;
  std::vector<double> massL = in.masses;

  fflush_MPI();
  fft_desc_basic fd(geo);
  eigen_ov ei(fd, n_vec, in.bini, in.nmass + 1);

  /////ei.ncut0 = in.ncut0;
  /////ei.ncut1 = in.ncut1;

  fflush_MPI();
  int mode_sm = 0;
  ei.load_eigen(icfg, in.Ename,  io_use);
  if(in.Ename_Sm != std::string("NONE"))
  {
    ei.load_eigen_Mvec(icfg, in.Ename_Sm, io_use, 1);
    mode_sm = 2;
  }

  ei.initialize_mass(massL, 12);
  ei.print_info();

  print0("Low eigen done. \n");

  size_t Nvol = geo.local_volume();
  
  /////EigenV props;
  EigenV prop1;
  EigenM s0;

  char names[450],namep[500];

  std::vector<Propagator4d > propS;propS.resize(1);for(unsigned int i=0;i<propS.size();i++)propS[i].init(geo);

  sprintf(names, in.Sname.c_str(),icfg);
  fflush_MPI();

  std::vector<qlat::FieldM<Complexq, 1> > noi;noi.resize(1);
  load_gwu_noiP(names, propS[0]);
  load_gwu_noi(names, noi[0], io_use);
  fflush_MPI();

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


  //if(step > 0)
  //{
  //  smear_propagator_gwu_convension(propS[0], gfD, width, step);
  //}

  copy_propE(propS, s0, fd);

  qlat::vector_gpu<Complexq > stmp, ptmp, bufP;
  stmp.resize(s0.size()*s0[0].size());
  ptmp.resize(massL.size()*s0.size()*s0[0].size());
  bufP.resize(massL.size()*s0.size()*s0[0].size());
  ptmp.set_zero();bufP.set_zero();

  std::vector<qlat::FieldM<Complexq , 12*12> > noi_prop;
  FieldM_src_to_FieldM_prop(noi, noi_prop, true);
  copy_eigen_src_to_FieldM(stmp, noi_prop, ei.b_size, fd, 1, true, false);

  //Complexq* psrc = stmp.data();
  //Complexq* pres = ptmp.data();

  //copy_eigen_prop_to_EigenM(psrc, s0, ei.b_size, 1, fd, 1);

  //size_t tmp_Lp = 2*massL.size()*ei.bfac * 12*ei.b_size;
  //zero_Ty(pres, tmp_Lp, 0, true);

  /////props.resize(nmass*12* 12*Nvol);zeroE(props);
  //qlat::set_zero(ei.ptmp);
  //ei.ptmp.set_zero();

  //prop_L_device(ei, psrc, pres, 12, massL, mode_sm);

  prop_L_device(ei, stmp.data(), ptmp.data(), 12, massL, mode_sm);

  ////print_sum((Complexq*)qlat::get_data(ei.ptmp).data(), qlat::get_data(ei.ptmp).data_size()/sizeof(Complexq),
  ////    "=====Final prop check");

  ////print_sum(pres, tmp_Lp, "=====Final prop check");

  std::vector< qlat::FieldM<qlat::Complex , 12*12> > FpropV;FpropV.resize(nmass);
  ///for(unsigned int iv=0;iv<FpropV.size();iv++){FpropV[iv].init(geo);}

  copy_eigen_src_to_FieldM(ptmp, FpropV, ei.b_size, fd, 0, 1, false);

  //qlat::Complex* tmp = (qlat::Complex*) qlat::get_data(FpropV[0]).data();
  //print0("size0 %ld, size1 %ld \n", 
  //  long(geo.local_volume()*12*12), long(qlat::get_data(FpropV[0]).data_size()/sizeof(qlat::Complex)));
  //print_sum(tmp, geo.local_volume(), "sum FpropV 0", 1);

  //print_sum(tmp, geo.local_volume(), "sum FpropV 2", 1);
  //print_sum((qlat::Complex*) qlat::get_data(FpropV[0]).data(), qlat::get_data(FpropV[0]).data_size()/sizeof(qlat::Complex), 
  //  "sum FpropV 2", 1);


  copy_eigen_src_to_FieldM(bufP, FpropV, ei.b_size, fd, 1, 1, false);

  //print_sum(tmp, geo.local_volume(), "sum FpropV 1", 1);
  //print_sum((qlat::Complex*) qlat::get_data(FpropV[0]).data(), qlat::get_data(FpropV[0]).data_size()/sizeof(qlat::Complex), 
  //  "sum FpropV 1", 1);

  copy_eigen_prop_to_EigenM(bufP.data(), s0, ei.b_size, nmass, fd, 0);

  //copy_eigen_prop_to_EigenM(ptmp.data(), s0, ei.b_size, nmass, fd, 0);

  std::vector<Propagator4d > propS_new;propS_new.resize(nmass);for(unsigned int i=0;i<propS_new.size();i++)propS_new[i].init(geo);
  copy_propE(propS_new, s0, fd, 1);

  //copy_eigen_prop_to_EigenM(ptmp.data(), s0, ei.b_size, nmass, fd, 0);
  //std::vector<Propagator4d > propS_new;propS_new.resize(nmass);for(unsigned int i=0;i<propS_new.size();i++)propS_new[i].init(geo);
  //copy_propE(propS_new, s0, fd, 1);

  //std::vector< qlat::FieldM<qlat::Complex , 12*12> > FpropV;FpropV.resize(nmass);
  //for(int im=0;im<nmass;im++){
  //  sprintf(names, in.Pname.c_str(),icfg);
  //  sprintf(namep, "%s.m%8.6f", names, massL[im]);
  //  load_gwu_prop(namep, FpropV[im], io_use);
  //}
  //copy_eigen_src_to_FieldM(bufP, FpropV, ei.b_size, fd, 1, 1, false);
  //std::vector<Propagator4d > bufS;bufS.resize(nmass);for(unsigned int i=0;i<propS_new.size();i++)bufS[i].init(geo);
  //copy_eigen_prop_to_EigenM(bufP.data(), s0, ei.b_size, nmass, fd, 0);
  //copy_propE(bufS, s0, fd, 1);

  //for(int im=0;im<nmass;im++){
  //print_sum((Complex*)qlat::get_data(propS_new[im]).data(), qlat::get_data(propS_new[im]).data_size()/sizeof(qlat::Complex),
  //    "=====Final prop");
  //}

  propS.resize(0);
  propS.resize(nmass);for(unsigned int i=0;i<propS.size();i++)propS[i].init(geo);
  //////for(unsigned int i=0;i<propS.size();i++){propS[i] = bufS[i];}

  print0("Prop size %5d \n", int(propS.size()));
  for(int im=0;im<nmass;im++)
  {
    sprintf(names, in.Pname.c_str(),icfg);
    sprintf(namep, "%s.m%8.6f"
        , names, massL[im]);
    load_gwu_prop(namep, propS[im]);
    //if(step > 0)
    //{
    //  smear_propagator_gwu_convension(propS[im], gfD, width, step);
    //}
  }

  for(int mi=0;mi<nmass;mi++){
  diff_prop(propS[mi], propS_new[mi]);}


  if(in.debuga != 0)
  {
    EigenV src;
    src.resize(12 * 12*Nvol);
    for(int i=0;i<in.debuga;i++)
    {
      random_EigenM(src);
      /////ei.prop_L(src, props, massL);
      //#ifdef QLAT_USE_ACC
      //cudaMemcpy(  ei.stmp, &src[0]  , ei.stmp_size*sizeof(Complexq),cudaMemcpyDeviceToDevice);
      //#else
      //memcpy(  ei.stmp, &src[0]  , ei.stmp_size*sizeof(Complexq));
      //#endif

      //prop_L_device(ei, ei.stmp.data(), ei.ptmp.data(), 12, massL);
    }
  }



  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

