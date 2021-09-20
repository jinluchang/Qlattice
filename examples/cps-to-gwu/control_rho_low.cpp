#include <sys/sysinfo.h>
#include "io_vec.h"
#include "utils_low_rho.h"

#define Cfield qlat::FieldM<qlat::MvectorT<3,Complexq > ,1> 

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
  //int n_node = init_mpi(&argc, &argv);
  int n_node = init_mpi_thread(&argc, &argv, 1);

  inputpara in; 
  in.load_para(argc, argv);
  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;
  Coordinate Lat(in.nx, in.ny, in.nz, in.nt);
  Coordinate spreadT = guess_nodeL(n_node, Lat);

  std::vector<Coordinate > size_node_list;
  size_node_list.push_back(spreadT);
  begin_comm(MPI_COMM_WORLD , spreadT);
  set_GPU();

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

  char ename[500],enamev[500];
  char Ename[500];

  sprintf(ename,in.Ename.c_str(), icfg);
  sprintf(enamev,"%s.eigvals",ename);


  print_mem_info();

  double length = (geo.local_volume()*pow(0.5,30))*12*sizeof(Complexq);
  print0("Eign system vector size %.3e GB, total %.3e GB; \n", length,n_vec*length);

  std::vector<qlat::FermionField4dT<Complexq > > eigen;
  //int threadio = false;
  //int threadio = true;
  {
    TIMER("new io read");
    int Nv = omp_get_max_threads();
    int useio = ionum;
    if(useio > qlat::get_num_node()){useio = qlat::get_num_node();}
    if(Nv*useio > 64){Nv = 64/useio;}
    io_vec io_use(geo,ionum, true, Nv);
    load_gwu_eigen(ename,eigen,io_use,vini,n_vec,true, true);

    //if(threadio == false){
    //  io_vec io_use(geo,ionum);
    //  load_gwu_eigen(ename,eigen,io_use,vini,n_vec,true);
    //}

    //if(threadio == true){
    //  int Nvec = n_vec - vini;
    //  eigen.resize(0);
    //  eigen.resize(Nvec);
    //  for(int iv=0;iv<Nvec;iv++){eigen[iv].init(geo);}

    //  int Nv = omp_get_max_threads();
    //  ////int Nv = 1;
    //  int Group = (Nvec-1)/Nv + 1;
    //  print0("===Nvec %d, Nv %d, Group %d \n", Nvec, Nv, Group);
    //  #pragma omp parallel for
    //  for(int tid=0;tid<Nv;tid++)
    //  {
    //    io_vec io_use(geo, ionum);
    //    /////int tid = omp_get_thread_num();
    //    int currN = Group; if((tid+1)*Group > Nvec){currN = Nvec - tid*Group;}
    //    if(currN > 0){
    //      int iniN  = tid*Group;int endN = tid*Group + currN;
    //      std::vector<Ftype* > resp;resp.resize(currN);
    //      for(int iv=0;iv<currN;iv++){
    //        resp[iv]=(Ftype*)(qlat::get_data(eigen[iniN + iv]).data());
    //      }
    //      load_gwu_eigen(ename, resp, io_use, vini+iniN, vini+endN, true, true);
    //    }
    //  }
    //}

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
  for(int iv=0;iv<eigen.size();iv++){
    Complexq* src = (Complexq* ) &(eigen[iv].get_elem(0));
    Complexq* buf = (Complexq* ) &(tmp.get_elem(0));
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
  //std::vector<double > mass;mass.resize(16);
  //mass[0 ]=0.005200;
  //mass[1 ]=0.006000;
  //mass[2 ]=0.009000;
  //mass[3 ]=0.014000;
  //mass[4 ]=0.022000;
  //mass[5 ]=0.032000;
  //mass[6 ]=0.052000;
  //mass[7 ]=0.080000;
  //mass[8 ]=0.100000;
  //mass[9 ]=0.150000;
  //mass[10]=0.170000;
  //mass[11]=0.190000;
  //mass[12]=0.300000;
  //mass[13]=0.400000;
  //mass[14]=0.500000;
  //mass[15]=0.800000;

  /////====Set up values
  qlat::vector_acc<Complexq > values_sys,values;
  values_sys.resize(n_vec);values.resize(n_vec);
  for(int i=0;i<n_vec;i++){
    values_sys[i] = Complexq(eval_self[i*2+0],eval_self[i*2+1]);
  }

  double kappa = 0.2;double rho = 4 - 1.0/(2*kappa);int one_minus_halfD = 1;double Eerr = 1e-10;
  int n_zero = 0;
  Complexq rhoc = Complexq(rho,0.0);
  for(int j=0; j<n_vec; ++j){
    values_sys[j] = values_sys[j]/rhoc;
    if(qnorm(values_sys[j]) < Eerr) n_zero += 1;
  }
  int Nmass = mass.size();
  values.resize(Nmass*n_vec);
  for(int mi=0;mi<Nmass;mi++)
  for(int iv=0;iv<n_vec;iv++){values[iv*Nmass + mi] = inv_self(values_sys[iv], mass[mi], rho,one_minus_halfD);}
  /////====Set up values
  fflush_MPI();

  {TIMER("Kernal functions");get_low_rho(eigen,values, n_zero, Mres,geo, GPUFM);}

  fflush_MPI();
  {
    TIMER("Save corr");
    int Nmpi = qlat::get_num_node();
    std::vector<double > write;write.resize(Mres.size());
    for(int iv=0;iv<Mres.size();iv++){
      write[iv] = Mres[iv];
    }
    char filename[500];
    sprintf(filename,"res/rho_N%06d_%06d.Nmpi%02d.dat",in.nvec,icfg,Nmpi);
    print0("%s",filename);
    write_data(write,filename);
  }

  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

