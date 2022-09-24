#include <sys/sysinfo.h>
#include "utils_io_vec.h"
#include "utils_low_rho.h"

#define Cfield qlat::FieldM<qlat::MvectorT<3,Complexq > ,1> 

int main(int argc, char* argv[])
{
  using namespace qlat;
  ////int n_node = init_mpi(&argc, &argv);
  //int n_node = init_mpi_thread(&argc, &argv, 1);

  //inputpara in; 
  //in.load_para(argc, argv);
  //int nx,ny,nz,nt;
  //nx = in.nx;
  //ny = in.ny;
  //nz = in.nz;
  //nt = in.nt;
  //Coordinate Lat(in.nx, in.ny, in.nz, in.nt);
  //Coordinate spreadT = guess_nodeL(n_node, Lat);

  //std::vector<Coordinate > size_node_list;
  //size_node_list.push_back(spreadT);
  //begin_comm(MPI_COMM_WORLD , spreadT);
  //set_GPU();

  inputpara in;int mode_dis = 0;
  begin_Lat(&argc, &argv, in, mode_dis);
  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;


  int icfg  = in.icfg;
  int ionum = in.ionum;

  int vini  = 0;
  int n_vec = in.nvec;

  ////Eigen::initParallel();
  omp_set_num_threads(omp_get_max_threads());
  print0("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 

  int GPUFM = 0;

  #ifdef QLAT_USE_ACC
  GPUFM = 1;
  {
    int Nt = geo.node_site[3];
    int Nmpi = qlat::get_num_node();
    if(nt/Nt != Nmpi){
      print0("Number of nodes not sparsed on time. \n");
      GPUFM = 2;
      ///qassert(false);
    }
  }
  #endif

  char ename[500],enamev[510];

  sprintf(ename,in.Ename.c_str(), icfg);
  sprintf(enamev,"%s.eigvals",ename);


  print_mem_info();

  double length = (geo.local_volume()*pow(0.5,30))*12*sizeof(Complexq);
  print0("Eign system vector size %.3e GB, total %.3e GB; \n", length,n_vec*length);

  ////std::vector<qlat::FermionField4dT<Complexq > > eigen;
  std::vector<qlat::FieldM<Complexq, 12> > eigen;
  ////eigen.resize(n_vec);

  //int threadio = false;
  //int threadio = true;
  {
    TIMER("new io read");
    int Nv = omp_get_max_threads();
    int useio = ionum;
    if(useio > qlat::get_num_node()){useio = qlat::get_num_node();}
    if(Nv*useio > 64){Nv = 64/useio;}
    io_vec io_use(geo,ionum, true, Nv);

    //load_gwu_eigen(ename,eigen,io_use,vini,n_vec,true, true);
    inputpara in_read_eigen;
    FILE* file_read  = open_eigensystem_file(ename, vini, n_vec + vini, true , io_use , in_read_eigen ); 
    load_eigensystem_vecs(file_read ,   eigen, io_use , in_read_eigen , vini, n_vec + vini);
    close_eigensystem_file(file_read , io_use , in_read_eigen );
  }
  

  std::vector<double > eval_self,errors;
  load_gwu_eigenvalues(eval_self,errors,enamev);

  qlat::vector_acc<Ftype > Mres;

  ////====Set up eigen_chi
  //if(ordermem == 1)
  {
  TIMER("Rotate eigen d,c");
  /////Switch d,c to outter side
  size_t Nvol = geo.local_volume();
  int Nt = geo.node_site[3];
  qlat::FermionField4dT<Complexq > tmp;tmp.init(geo);
  for(LInt iv=0;iv<eigen.size();iv++){
    Complexq* src = (Complexq* ) qlat::get_data(eigen[iv]).data();
    Complexq* buf = (Complexq* ) qlat::get_data(tmp).data();
    memcpy(buf,src, Nvol*12*sizeof(Complexq));
    if(GPUFM == 1 or GPUFM == 2)
    {
      #pragma omp parallel for
      for(LInt isp=0;isp<Nvol;isp++)
      {
        for(int d=0;d<4;d++)
        for(int c=0;c<3;c++){
          src[(d*3 + c)*Nvol + isp ] = buf[isp*12 + d*3 + c];
        }
      }
    }

    ////////t in the outside
    if(GPUFM == 0)
    {
      LInt Nsum = Nvol/Nt;
      #pragma omp parallel for
      for(LInt isp=0;isp<Nsum;isp++)
      for(int it=0;it<Nt;it++)
      {
        for(int d=0;d<4;d++)
        for(int c=0;c<3;c++){
          src[it*12*Nsum + (d*3 + c)*Nsum + isp ] = buf[(it*Nsum+isp)*12 + d*3 + c];
        }
      }
    }

  }
  }
  //====Set up eigen_chi

  std::vector<double > mass;mass = in.masses;
  /////====Set up values
  qlat::vector_acc<Complexq > values_sys,values;
  values_sys.resize(n_vec);values.resize(n_vec);
  for(int i=0;i<n_vec;i++){
    values_sys[i] = Complexq(eval_self[i*2+0],eval_self[i*2+1]);
  }

  double kappa = 0.2;double rho = 4 - 1.0/(2*kappa);int one_minus_halfD = 1;
  //double Eerr = 1e-11;
  int n_zero = 0;
  Complexq rhoc = Complexq(rho,0.0);
  for(int j=0; j<n_vec; ++j){
    values_sys[j] = values_sys[j]/rhoc;
    if(qnorm(values_sys[j]) < in.Eerr) n_zero += 1;
  }
  int Nmass = mass.size();
  values.resize(Nmass*n_vec);
  for(int mi=0;mi<Nmass;mi++)
  for(int iv=0;iv<n_vec;iv++){values[iv*Nmass + mi] = inv_self(values_sys[iv], mass[mi], rho,one_minus_halfD);}
  print0("===num zero %d \n", n_zero);
  /////====Set up values
  fflush_MPI();

  {TIMER("Kernal functions");get_low_rho(eigen,values, n_zero, Mres,geo, GPUFM);}

  fflush_MPI();
  {
    //TIMER("Save corr");
    //int Nmpi = qlat::get_num_node();
    //std::vector<double > write;write.resize(Mres.size());
    //for(int iv=0;iv<Mres.size();iv++){
    //  write[iv] = Mres[iv];
    //}
    //char filename[500];
    //////sprintf(filename,"res/rho_N%06d_%06d.Nmpi%02d.dat",in.nvec,icfg,Nmpi);
    //sprintf(filename, in.output.c_str(),icfg);
    //print0("%s",filename);
    //write_data(write,filename);


    char key_T[1000], dimN[1000];
    sprintf(key_T, "%d  %d   %d  %d  %d", in.nmass, 16, in.nt, in.nt, 1);
    sprintf(dimN , "masses operator t0 nt complex");

    std::string ktem(key_T);
    std::string dtem(dimN);
    corr_dat res(ktem, dtem);
    res.print_info();

    char names[500];
    res.write_corr((Ftype*) Mres.data(), Mres.size());
    sprintf(names, in.output.c_str(),icfg);
    res.write_dat(names);

  }

  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

