#include <qlat/qlat.h>
#include <sys/sysinfo.h>
#include "io_gwu.h"
////#include "utils_low_rho.h"
#include "utils_eigensys.h"
#include "utils_construction.h"

inline Complexq inv_self(const Complexq& lam, double m, double rho,int one_minus_halfD=1)
{
  //Complexq tem = (one_minus_halfD>0)?(1-lam/2)/(rho*lam+m*(1-lam/2)):1.0/(rho*lam+m*(1-lam/2));
  std::complex<double > tem(lam.real(),lam.imag());
  std::complex<double > v0 = (one_minus_halfD>0)?(1.0-tem/2.0)/(rho*tem+m*(1.0-tem/2.0)):1.0/(rho*tem+m*(1.0-tem/2.0));
  Complexq res(v0.real(),v0.imag());
  return res;
}

int main(int argc, char* argv[])
{
  using namespace qlat;

  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1,  1));
  size_node_list.push_back(Coordinate(1, 1, 1,  2));
  size_node_list.push_back(Coordinate(1, 1, 1,  4));
  size_node_list.push_back(Coordinate(1, 1, 1,  6));
  size_node_list.push_back(Coordinate(1, 1, 1,  8));
  size_node_list.push_back(Coordinate(1, 1, 1, 12));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 1, 24));
  size_node_list.push_back(Coordinate(1, 1, 1, 32));
  size_node_list.push_back(Coordinate(1, 1, 1, 64));
  //size_node_list.push_back(Coordinate(1, 1, 1, 48));
  //size_node_list.push_back(Coordinate(1, 1, 1, 96));
  //size_node_list.push_back(Coordinate(1, 1, 1,128));
  size_node_list.push_back(Coordinate(4, 4, 8, 16));
  size_node_list.push_back(Coordinate(4, 8, 8, 16));
  //size_node_list.push_back(Coordinate(1, 2, 2, 16));
  //size_node_list.push_back(Coordinate(1, 1, 2, 16));
  //size_node_list.push_back(Coordinate(1, 2, 2, 16));
  //size_node_list.push_back(Coordinate(2, 2, 2, 16));
  //size_node_list.push_back(Coordinate(2, 2, 4, 16));
  //size_node_list.push_back(Coordinate(2, 4, 4, 16));
  //size_node_list.push_back(Coordinate(4, 4, 4, 16));

  //begin_thread(&argc, &argv, size_node_list);
  begin(&argc, &argv, size_node_list);
  set_GPU();
  ///set_GPU_threads();

  //fft_desc_basic fd();

  //Coordinate size_node = Coordinate(fd.mx, fd.my, fd.mz, fd.mt);
  //begin(fd.rank, size_node);
  //begin(MPI_COMM_WORLD, size_node);

  inputpara in;
  in.load_para("input.txt");

  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;

  int icfg  = in.icfg;
  int ionum = in.ionum;

  ////int vini  = 0;
  int n_vec = in.nvec;

  omp_set_num_threads(omp_get_max_threads());
  print0("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 
  fflush_MPI();

  char ename[500],enamev[500];

  sprintf(ename,in.Ename.c_str(), icfg);
  sprintf(enamev,"%s.eigvals",ename);

  io_gwu io_use(geo,ionum);

  double length = (geo.local_volume()*pow(0.5,30))*12*sizeof(Complexq);
  size_t freeM = 0;size_t totalM = 0;
  #ifdef QLAT_USE_ACC
  cudaMemGetInfo(&freeM,&totalM);
  #endif
  double freeD = freeM*pow(0.5,30);double totalD = totalM*pow(0.5,30); 
  struct sysinfo s_info;
  sysinfo(&s_info);
  print0("Eign system vector size %.3e GB, Total %.3e GB; CPU free %.3e GB, total %.3e GB; GPU free %.3e GB, total %.3e GB. \n"
          ,length,n_vec*length,s_info.totalram*pow(0.5,30), s_info.freeram*pow(0.5,30),freeD,totalD);

  fflush_MPI();
  fft_desc_basic fd(geo);
  eigen_ov ei(fd, n_vec, in.bini, in.nmass*12+12);
  ei.ncut0 = in.ncut0;
  ei.ncut1 = in.ncut1;
  ////ei.ncutgpu = in.ncut0;

  //if(in.paraI != "None"){
  //  std::vector<std::string > Li = stringtolist(in.paraI);
  //  print0("Li %s, size %d \n", in.paraI.c_str(),int(Li.size()) );
  //  fflush_MPI();
  //  int gpuf0    = stringtonum(Li[0]);
  //  ei.buffGPU = gpuf0;
  //  //////if(Li.size() > 1)ei.ncutgpu = stringtonum(Li[1]);
  //}


  fflush_MPI();
  ei.load_eigen(icfg, in.Ename,  io_use);
  print0("Low eigen done");

  /////size_t Nvol = ei.Mvec[0].size();
  size_t Nvol = geo.local_volume();
  double diff = 0.0;
  
  EigenV src;EigenV props;
  EigenV prop1;
  EigenM s0;
  int nmass = 7;
  std::vector<double> massL;massL.resize(nmass);
  massL[0] = 0.008090;
  massL[1] = 0.010200;
  massL[2] = 0.013500;
  massL[3] = 0.016000;
  massL[4] = 0.020300;
  massL[5] = 0.057600;
  massL[6] = 0.063000;

  if(in.nmass != 0){nmass = in.nmass;massL.resize(nmass);}

  //for(unsigned int i=0;i<massL.size();i++){massL[i] = 0.8;}


  char names[500],namep[500];
  /////namep[500];

  std::vector<Propagator4d > propS;propS.resize(1);for(unsigned int i=0;i<propS.size();i++)propS[i].init(geo);

  //sprintf(names, "/global/homes/g/genwang/cscratch/24IDc/Prop/Test/rbc.24D.%06d.N000004.S",icfg);
  sprintf(names, in.Pname.c_str(),icfg);
  fflush_MPI();

  load_gwu_noiP(names, propS[0]);
  fflush_MPI();

  copy_propE(propS, s0, fd);

  ei.initiallize_mass(massL, 12);

  #ifdef QLAT_USE_ACC
  cudaMemset(ei.ptmp, 0, ei.ptmp_size*sizeof(Complexq));
  #else
  memset(    ei.ptmp, 0, ei.ptmp_size*sizeof(Complexq));
  #endif


  ///resize_EigenM(src  , 12, 12*Nvol);
  src.resize(12 * 12*Nvol);
  int NTt  = fd.Nv[3];
  Complexq* useS = ei.stmp;
  Complexq* useP = ei.ptmp;
  int b_size = ei.b_size; int bfac = ei.bfac;
  qacc_for(c, long(Nvol),{
    int  ti = c/(Nvol/NTt);
    long vi = c%(Nvol/NTt);

    for(unsigned int d0=0;d0<12;d0++)
    for(unsigned int d1=0;d1<12;d1++)
    {
    int is = d0;int Ns = 12;
    int chi = d1/6;
    long xi = (d1%6)*Nvol + ti*(Nvol/NTt) + vi;
    long bi = xi/b_size;
    long bj = xi%b_size;

    //src[(chi*ei.bfac+bi)*Ns*ei.b_size  + is*ei.b_size + bj] = s0[(d0*12 + d1)*NTt+ti][vi];
    useS[(chi*bfac+bi)*Ns*b_size  + is*b_size + bj] = s0[(d0*12 + d1)*NTt+ti][vi];
    }
  });

  //qacc_for(c, long(12*12*Nvol),{
  //  int  d0 = c/(12*Nvol);
  //  int  d1 = (c%(12*Nvol))/Nvol;
  //  long vI = c%Nvol;
  //  int  ti = vI/(Nvol/NTt);
  //  long vi = vI%(Nvol/NTt);

  //  int is = d0;int Ns = 12;
  //  int chi = d1/6;
  //  long xi = (d1%6)*Nvol + ti*(Nvol/NTt) + vi;
  //  long bi = xi/b_size;
  //  long bj = xi%b_size;
  //  //src[(chi*ei.bfac+bi)*Ns*ei.b_size  + is*ei.b_size + bj] = s0[(d0*12 + d1)*NTt+ti][vi];
  //  useS[(chi*bfac+bi)*Ns*b_size  + is*b_size + bj] = s0[(d0*12 + d1)*NTt+ti][vi];
  //
  //});
  //for(unsigned int d0=0;d0<12;d0++)
  //for(unsigned int d1=0;d1<12;d1++)
  //for(int ti=0;ti<NTt;ti++)
  //for(int vi=0;vi<int(Nvol/NTt);vi++)
  //{
  //  //src[d0][d1*Nvol + ti*(Nvol/NTt) + vi] = s0[(d0*12 + d1)*NTt+ti][vi];
  //  int is = d0;int Ns = 12;
  //  int chi = d1/6;
  //  long xi = (d1%6)*Nvol + ti*(Nvol/NTt) + vi;
  //  long bi = xi/ei.b_size;
  //  long bj = xi%ei.b_size;
  //  //src[(chi*ei.bfac+bi)*Ns*ei.b_size  + is*ei.b_size + bj] = s0[(d0*12 + d1)*NTt+ti][vi];
  //  ei.stmp[(chi*ei.bfac+bi)*Ns*ei.b_size  + is*ei.b_size + bj] = s0[(d0*12 + d1)*NTt+ti][vi];
  //}

  //resize_EigenM(prop1, nmass*12, 12*Nvol);
  //resize_EigenM(props, nmass*12, 12*Nvol);
  prop1.resize(nmass*12* 12*Nvol);
  props.resize(nmass*12* 12*Nvol);
  zeroE(prop1);zeroE(props);



  prop_L_device(ei, ei.stmp, ei.ptmp, 12, massL);
  //ei.prop_L(src, props, massL);

  ////resize_EigenM(src  , 12, 12*Nvol);
  ////resize_EigenM(props, nmass*12, 12*Nvol);
  ////ei.prop_L(src, props, massL);

  /////for(unsigned int i=0;i<20;i++){ for(int j=0;j<12;j++)src[i][j] = 1.0;}

  ////{ei.setup_bfac(  10);
  ////TIMER("Low bfac  10");
  ////ei.copy_Mvec(ei.Mvec, ei.Mvec0);
  ////ei.prop_L(src, props, massL);}

  ////{ei.setup_bfac(  30);
  ////TIMER("Low bfac  30");
  ////ei.copy_Mvec(ei.Mvec, ei.Mvec0);
  ////ei.prop_L(src, props, massL);}

  ////{ei.setup_bfac( 10);
  ////TIMER("Low bfac 10");
  ////ei.copy_Mvec(ei.Mvec, ei.Mvec0);
  ////ei.prop_L(src, props, massL);}

  ////{ei.setup_bfac( 300);
  ////TIMER("Low bfac 300");
  ////ei.copy_Mvec(ei.Mvec, ei.Mvec0);
  ////ei.prop_L(src, props, massL);}

  ////{ei.setup_bfac( 500);
  ////TIMER("Low bfac 500");
  ////ei.copy_Mvec(ei.Mvec, ei.Mvec0);
  ////ei.prop_L(src, props, massL);}

  ////{ei.setup_bfac( 1000);
  ////TIMER("Low bfac 1000");
  ////ei.copy_Mvec(ei.Mvec, ei.Mvec0);
  ////ei.prop_L(src, props, massL);}


  //////ei.ncut0 = in.ncut0;

  ////ei.ncut0  = 30;
  ////ei.ncut1  = 300;
  ////for(int il=0;il<20;il++){
  ////TIMER("Low   30 300");
  ////resize_EigenM(props, nmass*12, 12*Nvol);
  ////ei.prop_L(src, props, massL);
  ////}

  //////ei.ncut0  = 20;
  //////ei.ncut1  = 800;
  ////for(int il=0;il<20;il++){
  ////TIMER("Low   20 800");
  ////resize_EigenM(props, nmass*12, 12*Nvol);
  ////ei.prop_L(src, props, massL);
  ////}


  ////for(int il=0;il<10;il++){
  ////{TIMER("Low   10");
  ////ei.ncut0  =   10;
  ////resize_EigenM(props, nmass*12, 12*Nvol);
  ////ei.prop_L(src, props, massL);}
  ////{TIMER("Low   50");
  ////ei.ncut0  =   50;
  ////resize_EigenM(props, nmass*12, 12*Nvol);
  ////ei.prop_L(src, props, massL);}

  ////{TIMER("Low  100");
  ////ei.ncut0  =  100;
  ////resize_EigenM(props, nmass*12, 12*Nvol);
  ////ei.prop_L(src, props, massL);}


  ////{TIMER("Low  200");
  ////ei.ncut0  =  200;
  ////resize_EigenM(props, nmass*12, 12*Nvol);
  ////ei.prop_L(src, props, massL);}

  ////{TIMER("Low  300");
  ////ei.ncut0  =  300;
  ////resize_EigenM(props, nmass*12, 12*Nvol);
  ////ei.prop_L(src, props, massL);}

  ////{TIMER("Low  400");
  ////ei.ncut0  =  400;
  ////resize_EigenM(props, nmass*12, 12*Nvol);
  ////ei.prop_L(src, props, massL);}

  ////{TIMER("Low  500");
  ////ei.ncut0  =  500;
  ////resize_EigenM(props, nmass*12, 12*Nvol);
  ////ei.prop_L(src, props, massL);}

  ////{TIMER("Low  600");
  ////ei.ncut0  =  600;
  ////resize_EigenM(props, nmass*12, 12*Nvol);
  ////ei.prop_L(src, props, massL);}

  ////{TIMER("Low  800");
  ////ei.ncut0  =  800;
  ////resize_EigenM(props, nmass*12, 12*Nvol);
  ////ei.prop_L(src, props, massL);}

  ////{TIMER("Low 1000");
  ////ei.ncut0  = 1000;
  ////resize_EigenM(props, nmass*12, 12*Nvol);
  ////ei.prop_L(src, props, massL);}
  ////}


  propS.resize(0);
  propS.resize(nmass);for(unsigned int i=0;i<propS.size();i++)propS[i].init(geo);
  print0("Prop size %5d \n", int(propS.size()));
  /////ga_matrices_milc ga_milc;
  for(int im=0;im<nmass;im++)
  {
    ////sprintf(names, "/global/homes/g/genwang/cscratch/24IDc/Prop/Test/rbc.24D.%06d.N000004.S.m%8.6f",icfg,massL[im]);
    //sprintf(names, "/global/homes/g/genwang/cscratch/24IDc/Prop/Test/rbc.24D.%06d.N000004.S.double.v%05d.m%8.6f"
    //    ,icfg, n_vec, massL[im]);
    //sprintf(names, "/global/homes/g/genwang/cscratch/24IDc/Prop/Test/rbc.24D.%06d.N000004.S.v%05d.m%8.6f"
    //    ,icfg, n_vec,massL[im]);

    sprintf(namep, "%s.double.v%05d.m%8.6f"
        ,names, n_vec, massL[im]);

    print0("%s \n",names);
    load_gwu_prop(namep, propS[im]);
    ////prop4d_src_gamma( propS[im], ga_milc.ga[0][5]);
    ////prop4d_sink_gamma(propS[im], ga_milc.ga[0][5]);
  }

  copy_propE(propS, s0, fd);

  for(int d0=0;d0<nmass*12;d0++)
  for(unsigned int d1=0;d1<12;d1++)
  for(int ti=0;ti<NTt;ti++)
  for(int vi=0;vi<int(Nvol/NTt);vi++)
  {
    //prop1[d0][d1*Nvol + ti*(Nvol/NTt) + vi] = s0[(d0*12 + d1)*NTt+ti][vi];

    int is = d0;int Ns = nmass*12;
    int chi = d1/6;
    long xi = (d1%6)*Nvol + ti*(Nvol/NTt) + vi;
    long bi = xi/ei.b_size;
    long bj = xi%ei.b_size;
    prop1[(chi*ei.bfac+bi)*Ns*ei.b_size  + is*ei.b_size + bj] = s0[(d0*12 + d1)*NTt+ti][vi];

  }

  diff = 0.0;
  qlat::vector<double > diffL;diffL.resize(Nvol);

  qacc_for(c, long(Nvol),{
    int ti = c/(Nvol/NTt);
    int vi = c%(Nvol/NTt);
    for(int d0=0;d0<nmass*12;d0++)
    for(unsigned int d1=0;d1<12;d1++)
    {
      Complexq tem = useP[d0*12*Nvol + d1*Nvol + ti*(Nvol/NTt) + vi] - prop1[d0*12*Nvol + d1*Nvol + ti*(Nvol/NTt) + vi];
      //Complexq tem = useP[d0*12*Nvol + d1*Nvol + ti*(Nvol/NTt) + vi];
      diffL[c] += (tem.real()*tem.real() + tem.imag()*tem.imag());
    }
  });
  for(int i=0;i<Nvol;i++){diff+=diffL[i];}

  //for(int d0=0;d0<nmass*12;d0++)
  //for(unsigned int d1=0;d1<12;d1++)
  //for(int ti=0;ti<NTt;ti++)
  //for(int vi=0;vi<int(Nvol/NTt);vi++)
  //{
  //  //Complexq tem = props[d0*12*Nvol + d1*Nvol + ti*(Nvol/NTt) + vi] - prop1[d0*12*Nvol + d1*Nvol + ti*(Nvol/NTt) + vi];
  //  Complexq tem = ei.ptmp[d0*12*Nvol + d1*Nvol + ti*(Nvol/NTt) + vi] - prop1[d0*12*Nvol + d1*Nvol + ti*(Nvol/NTt) + vi];
  //  diff += (tem.real()*tem.real() + tem.imag()*tem.imag());
  //}

  sum_all_size(&diff,1);
  print0("Prop diff %.5e \n",diff);

  if(in.debuga != 0)
  {
    //resize_EigenM(src  , 12, 12*Nvol);
    //resize_EigenM(props, nmass*12, 12*Nvol);
    //src.resize(12 * 12*Nvol);
    //props.resize(nmass*12* 12*Nvol);
    for(int i=0;i<in.debuga;i++)
    {
      ran_EigenM(src);
      /////ei.prop_L(src, props, massL);
      #ifdef QLAT_USE_ACC
      cudaMemcpy(  ei.stmp, &src[0]  , ei.stmp_size*sizeof(Complexq),cudaMemcpyDeviceToDevice);
      #else
      memcpy(  ei.stmp, &src[0]  , ei.stmp_size*sizeof(Complexq));
      #endif
      prop_L_device(ei, ei.stmp, ei.ptmp, 12, massL);
    }
  }



  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

