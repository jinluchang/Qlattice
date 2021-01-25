#include <qlat/qlat.h>
#include "io_gwu.h"
#include "utils_low_rho.h"

int main(int argc, char* argv[])
{
  using namespace qlat;

  /////namespace qcd
  //init_machine_thread(argc, argv,false);
  //timer walltime;walltime.start("over all");

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

  //fft_desc_basic fd();

  //Coordinate size_node = Coordinate(fd.mx, fd.my, fd.mz, fd.mt);
  //begin(fd.rank, size_node);
  //begin(MPI_COMM_WORLD, size_node);

  int nx,ny,nz,nt;
  nx = 24;
  ny = 24;
  nz = 24;
  nt = 64;

  int icfg  = 1040;
  int ionum = 8;

  int vini  = 0;
  int n_vec = 30;

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 

  //char namew[500],namer[500],prop_tmp[500],name[500],name_tem[500];
  //char name0[500],name1[500],filename[500];

  char ename[500],enamev[500];


  sprintf(ename ,"/global/homes/g/genwang/cscratch/24IDc/Eigen/f.rbc_conf_2464_m0.00107_0.0850_%06d_hyp.eta_zero.half.overlap.eigensystem",icfg);
  sprintf(enamev,"/global/homes/g/genwang/cscratch/24IDc/Eigen/f.rbc_conf_2464_m0.00107_0.0850_%06d_hyp.eta_zero.half.overlap.eigensystem.eigvals",icfg);

  io_gwu io_use(geo,ionum);

  std::vector<qlat::FermionField4dT<qlat::ComplexF> > eigen;
  {TIMER("new io read");load_gwu_eigen(ename,eigen,io_use,vini,n_vec,true);}

  std::vector<double > values,errors;
  load_gwu_eigenvalues(values,errors,enamev);

  std::vector<Ftype > Mres;

  get_low_rho(eigen,values,Mres,geo);

  fflush_MPI();
  qlat::Timer::display();



  qlat::end();
  return 0;
}

