#include <sys/sysinfo.h>
#include "general_funs.h"
#include "utils_io_vec.h"
#include "utils_construction.h"
#include "utils_eigen_ov.h"
#include "utils_check_fun.h"
#include "utils_smear_vecs.h"

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

  int icfg  = in.icfg;
  int ionum = in.ionum;

  int n_vec = in.nvec;
  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  const Geometry& geo = get_geo_cache(total_site);
  fflush_MPI();

  io_vec io_use(geo,ionum);

  double length = (geo.local_volume()*pow(0.5,30))*12*sizeof(Complexq);
  size_t freeM = 0;size_t totalM = 0;
  #ifdef QLAT_USE_ACC
  qlat_GPU_MemGetInfo(&freeM,&totalM);
  #endif
  double freeD = freeM*pow(0.5,30);double totalD = totalM*pow(0.5,30); 
  struct sysinfo s_info;
  sysinfo(&s_info);
  qmessage("Eign system vector size %.3e GB, Total %.3e GB; CPU free %.3e GB, total %.3e GB; GPU free %.3e GB, total %.3e GB. \n"
          ,length,n_vec*length, s_info.freeram*pow(0.5,30),s_info.totalram*pow(0.5,30),freeD,totalD);

  int nmass = in.nmass;
  std::vector<double> massL = in.masses;

  fflush_MPI();
  fft_desc_basic fd(geo);
  eigen_ov ei(geo, n_vec, in.bini, in.nmass + 1);

  /////ei.ncut0 = in.ncut0;
  /////ei.ncut1 = in.ncut1;

  fflush_MPI();
  int mode_sm = 0;
  char ename[500];
  sprintf(ename, in.Ename.c_str(),icfg);
  ei.load_eigen(ename);
  if(in.Ename_Sm != std::string("NONE"))
  {
    char ename[500];
    sprintf(ename, in.Ename_Sm.c_str(),icfg);
    ei.load_eigen_Mvec(ename, 1);
    mode_sm = 2;
  }

  ei.initialize_mass(massL, 12);
  ei.print_info();

  qmessage("Low eigen done. \n");

  size_t Nvol = geo.local_volume();
  
  /////EigenV props;
  ////EigenV prop1;
  ////EigenM s0;

  char names[450],namep[500];

  std::vector<Propagator4d > propS;propS.resize(1);for(unsigned int i=0;i<propS.size();i++)propS[i].init(geo);

  sprintf(names, in.Sname.c_str(),icfg);
  fflush_MPI();

  std::vector<qlat::FieldM<Complexq, 1> > noi;noi.resize(1);noi[0].init(geo);
  load_gwu_noiP(names, propS[0]);
  load_gwu_noi(names, noi[0]);
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

  const Long Size_prop = fd.get_prop_size();
  //copy_propE(propS, s0, fd);

  qlat::vector_gpu<Complexq > stmp, ptmp, bufP;
  stmp.resize(Size_prop);
  ptmp.resize(massL.size()*Size_prop);
  bufP.resize(massL.size()*Size_prop);
  ptmp.set_zero();bufP.set_zero();

  std::vector<qlat::FieldM<Complexq , 12*12> > noi_prop;
  FieldM_src_to_FieldM_prop(noi, noi_prop, true);
  copy_eigen_src_to_FieldM(stmp, noi_prop, ei.b_size, fd, 1, true, false);

  //qmessage("===src norm ");stmp.print_norm2();
  prop_L_device(ei, stmp.data(), ptmp.data(), 12, massL, mode_sm);
  //qmessage("===res norm0 ");ptmp.print_norm2();

  std::vector< qlat::FieldM<qlat::ComplexD , 12*12> > FpropV;FpropV.resize(nmass);

  //qmessage("===res norm ");ptmp.print_norm2();
  //copy_eigen_src_to_FieldM(ptmp, FpropV, ei.b_size, fd, 0, 1, false);
  copy_eigen_src_to_FieldM(ptmp, FpropV, ei.b_size, fd, 0, 1, false);

  //copy_eigen_src_to_FieldM(ptmp, FpropV, ei.b_size, fd, 1, 1, false);
  //qmessage("===res norm1 ");ptmp.print_norm2();
  //qmessage("===res norm ");ptmp.print_norm2();

  std::vector<Propagator4d > propS_new;propS_new.resize(nmass);
  for(unsigned int i=0;i<propS_new.size();i++){
    propS_new[i].init(geo);
    qprop_to_prop4d(propS_new[i], FpropV[i]);
  }

  propS.resize(0);
  propS.resize(nmass);for(unsigned int i=0;i<propS.size();i++)propS[i].init(geo);

  qmessage("Prop size %5d \n", int(propS.size()));
  for(int im=0;im<nmass;im++)
  {
    sprintf(names, in.Pname.c_str(),icfg);
    sprintf(namep, "%s.m%8.6f"
        , names, massL[im]);
    load_gwu_prop(namep, propS[im]);
  }

  for(int mi=0;mi<nmass;mi++){
    diff_prop(propS[mi], propS_new[mi]);
  }

  if(in.debuga != 0)
  {
    EigenV src;
    src.resize(12 * 12*Nvol);
    for(int i=0;i<in.debuga;i++)
    {
      random_EigenM(src);
      /////ei.prop_L(src, props, massL);
      //#ifdef QLAT_USE_ACC
      //qlat_GPU_Memcpy(  ei.stmp, &src[0]  , ei.stmp_size*sizeof(Complexq),qlat_GPU_MemcpyDeviceToDevice);
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

