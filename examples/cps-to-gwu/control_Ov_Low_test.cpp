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
  size_node_list.push_back(Coordinate(1, 1, 3,  1));
  size_node_list.push_back(Coordinate(1, 1, 1,  4));
  size_node_list.push_back(Coordinate(1, 1, 3,  2));
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
  /////ei.ncutgpu = in.ncut0;

  int Ns = 1;
  int vfac = 8;
  if(in.paraI != "None"){
    std::vector<std::string > Li = stringtolist(in.paraI);
    print0("Li %s, size %d \n", in.paraI.c_str(),int(Li.size()) );
    fflush_MPI();
    //int gpuf0    = stringtonum(Li[0]);
    //ei.buffGPU = gpuf0;
    if(Li.size() > 0)Ns   = stringtonum(Li[0]);
    if(Li.size() > 1)vfac = stringtonum(Li[1]);
    ////if(Li.size() > 1)ei.ncutgpu = stringtonum(Li[1]);
  }

  eigen_ov ei(fd, n_vec, in.bini, (in.nmass*12 + 12)*Ns);

  ei.ncutgpu = vfac;
  ei.setup_bfac(in.bini, (in.nmass*12 + 12)*Ns);
  ei.ncut0 = in.ncut0;
  ei.ncut1 = in.ncut1;


  fflush_MPI();
  ei.random_eigen();
  print0("Low eigen done");

  /////size_t Nvol = ei.Mvec[0].size();
  size_t Nvol = geo.local_volume();
  double diff = 0.0;
  
  EigenV src;EigenV props;
  int nmass = 12;
  std::vector<double> massL;massL.resize(nmass);
  massL[ 0] = 0.008090;
  massL[ 1] = 0.010200;
  massL[ 2] = 0.013500;
  massL[ 3] = 0.016000;
  massL[ 4] = 0.020300;
  massL[ 5] = 0.057600;
  massL[ 6] = 0.063000;
  massL[ 7] = 0.073000;
  massL[ 8] = 0.083000;
  massL[ 9] = 0.093000;
  massL[10] = 0.163000;
  massL[11] = 0.173000;

  if(in.nmass != 0){nmass = in.nmass;massL.resize(in.nmass);}


  ei.initiallize_mass(massL, 12*Ns);
  #ifdef QLAT_USE_ACC
  cudaMemset(ei.ptmp, 0, ei.ptmp_size*sizeof(Complexq));
  #else
  memset(    ei.ptmp, 0, ei.ptmp_size*sizeof(Complexq));
  #endif


  if(in.debuga != 0)
  {
    src.resize(   Ns*12*12*Nvol);
    //props.resize(nmass*12*12*Nvol);zeroE(props);
    for(int i=0;i<in.debuga;i++)
    {
      ran_EigenM(src, 1);
      #ifdef QLAT_USE_ACC
      cudaMemcpy(  ei.stmp, &src[0]  , ei.stmp_size*sizeof(Complexq),cudaMemcpyDeviceToDevice);
      #else
      memcpy(  ei.stmp, &src[0]  , ei.stmp_size*sizeof(Complexq));
      #endif
      /////ei.prop_L(src, props, massL);
      prop_L_device(ei, ei.stmp, ei.ptmp, 12, massL);

    }
  }


  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

