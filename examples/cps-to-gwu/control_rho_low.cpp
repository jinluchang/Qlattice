#include <qlat/qlat.h>
#include <sys/sysinfo.h>
#include "io_gwu.h"
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


  //std::cout << "executable= " << argv[0] << std::endl;
  //for (int i=1; i<argc; i++) {
  //  std::string s(argv[i]); //put char array into a string
  //  std::cout << "arg["<<i<<"]="<<s<<std::endl;
  //  for (int j=0; j<6; j+=2) {
  //    std::string byteString = s.substr(j, 2);
  //    char byte = (char) strtol(byteString.c_str(), NULL, 16);
  //    std::cout << "byteString= "<<byteString << " as integer= "<<(int)byte<<std::endl;
  //  }
  //}

  //in.bSize = 16;
  //if(argc >= 1){
  //  std::string s(argv[1]);in.bSize = std::stoi(s);
  //  print0("===input bSize %d \n",in.bSize);
  //}

  /////namespace qcd
  //init_machine_thread(argc, argv,false);
  //timer walltime;walltime.start("over all");

  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 1));
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 4));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
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
  nx = 24;
  ny = 24;
  nz = 24;
  nt = 64;

  int icfg  = 1040;
  /////int icfg  = 1320;
  ////int icfg  = 1320;
  int ionum = 16;

  int vini  = 0;
  int n_vec = 30;
  //int n_vec = 200;
  //int n_vec = 350;
  //int n_vec = 400;
  //int n_vec = 1000;

  omp_set_num_threads(omp_get_max_threads());
  print0("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());
  /////int Nv = omp_get_num_threads();
  //#pragma omp parallel
  //for(unsigned long index=0;index<omp_get_max_threads();index++){
  //  print0("====index %8d, thread number %8d \n",index, omp_get_thread_num());
  //}   

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 
  {
    int Nt = geo.node_site[3];
    int noden = qlat::get_num_node();
    if(nt/Nt != noden){
      print0("Number of nodes not sparsed on time. \n");
      qassert(false);
    }
  }

  //////Bind cpu and gpu
  //{
  //  int rank, local_rank, local_size;
  //  MPI_Comm local_comm;
  //  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank,  MPI_INFO_NULL, &local_comm);
  //  MPI_Comm_size(local_comm, &local_size);
  //  MPI_Comm_rank(local_comm, &local_rank);
  //  cudaSetDevice(local_rank%local_size);
  //}


  //////===Test gammas
  ///std::vector<Complexq > Av;Av.resize(16);
  ///std::vector<Complexq > a;a.resize(4);
  ///std::vector<Complexq > b;b.resize(4);
  ///a[0] = Complexq(1.0,0.2);
  ///a[1] = Complexq(0.0,1.2);
  ///a[2] = Complexq(8.0,0.0);
  ///a[3] = Complexq(1.0,1.0);

  ///b[0] = Complexq(3.0,0.2);
  ///b[1] = Complexq(0.0,5.2);
  ///b[2] = Complexq(9.0,0.0);
  ///b[3] = Complexq(1.0,7.0);

  ///for(int i=0;i<16;i++){
  ///  Av[i] = a[i/4]*b[i%4];
  ///}

  ///std::vector<ga_M > gL;gL.resize(16);
  ///{int o=0;
  ///for(int i=0;i<6;i++){gL[o] = ga_cps.ga[0][i];o+=1;}
  ///for(int i=2;i<6;i++){gL[o] = ga_cps.ga[1][i];o+=1;}
  ///for(int i=3;i<6;i++){gL[o] = ga_cps.ga[2][i];o+=1;}
  ///for(int i=4;i<6;i++){gL[o] = ga_cps.ga[3][i];o+=1;}
  ///for(int i=5;i<6;i++){gL[o] = ga_cps.ga[4][i];o+=1;}}

  ///for(int i=0;i<16;i++){
  ///  Complexq res = reduce_gamma(Av,gL[i]);
  ///  print0("i %3d, %.3f %.3f \n",i,res.real(),res.imag());
  ///}


  /////char namew[500],namer[500],prop_tmp[500],name[500],name_tem[500];
  /////char name0[500],name1[500],filename[500];

  char ename[500],enamev[500];

  sprintf(ename ,"/global/homes/g/genwang/cscratch/24IDc/Eigen/f.rbc_conf_2464_m0.00107_0.0850_%06d_hyp.eta_zero.half.overlap.eigensystem",icfg);
  sprintf(enamev,"/global/homes/g/genwang/cscratch/24IDc/Eigen/f.rbc_conf_2464_m0.00107_0.0850_%06d_hyp.eta_zero.half.overlap.eigensystem.eigvals",icfg);

  io_gwu io_use(geo,ionum);

  double length = (geo.local_volume()*pow(0.5,30))*12*sizeof(Complexq);
  size_t freeM = 0;size_t totalM = 0;
  #ifdef QLAT_USE_ACC
  cudaMemGetInfo(&freeM,&totalM);
  #endif
  double freeD = freeM*pow(0.5,30);double totalD = totalM*pow(0.5,30); 
  struct sysinfo s_info;
  sysinfo(&s_info);
  print0("Eign system vector size %.3e GB, total %.3e GB; CPU free %.3e GB, total %.3e GB; GPU free %.3e GB, total %.3e GB. \n"
          ,n_vec*length,length,s_info.totalram*pow(0.5,30), s_info.freeram*pow(0.5,30),freeD,totalD);

  std::vector<qlat::FermionField4dT<Complexq > > eigen;
  {TIMER("new io read");load_gwu_eigen(ename,eigen,io_use,vini,n_vec,true);}

  std::vector<double > eval_self,errors;
  load_gwu_eigenvalues(eval_self,errors,enamev);

  //test_sum_vec(eigen,eval_self,geo);

  qlat::vector<Ftype > Mres;

  ////====Set up eigen_chi
  //if(ordermem == 1)
  {
  TIMER("Rotate eigen d,c");
  /////Switch d,c to outter side
  size_t Nvol = geo.local_volume();
  qlat::FermionField4dT<Complexq > tmp;tmp.init(geo);
  for(int iv=0;iv<eigen.size();iv++){
    //Complexq* src = (Complexq* ) qlat::get_data(eigen[iv]).data();
    //Complexq* buf = (Complexq* ) qlat::get_data(tmp).data();
    Complexq* src = (Complexq* ) &(eigen[iv].get_elem(0));
    Complexq* buf = (Complexq* ) &(tmp.get_elem(0));
    memcpy(buf,src, Nvol*12*sizeof(Complexq));
    ///std::vector<Complexq* > rp;rp.resize(4);
    ///for(int ir=0;ir<4;ir++){rp[ir] = (Complexq* ) qlat::get_data(eigen_chi[iv*4+ir]).data();}
    #pragma omp parallel for
    for(LInt isp=0;isp<Nvol;isp++){
      for(int d=0;d<4;d++)
      for(int c=0;c<3;c++){
        src[(d*3 + c)*Nvol + isp ] = buf[isp*12 + d*3 + c];
      }
    }
  }
  }
  //====Set up eigen_chi

  std::vector<double > mass;mass.resize(16);
  mass[0 ]=0.005200;
  mass[1 ]=0.006000;
  mass[2 ]=0.009000;
  mass[3 ]=0.014000;
  mass[4 ]=0.022000;
  mass[5 ]=0.032000;
  mass[6 ]=0.052000;
  mass[7 ]=0.080000;
  mass[8 ]=0.100000;
  mass[9 ]=0.150000;
  mass[10]=0.170000;
  mass[11]=0.190000;
  mass[12]=0.300000;
  mass[13]=0.400000;
  mass[14]=0.500000;
  mass[15]=0.800000;


  /////====Set up values
  qlat::vector<Complexq > values_sys,values;
  values_sys.resize(n_vec);values.resize(n_vec);
  for(int i=0;i<n_vec;i++){
    values_sys[i] = Complexq(eval_self[i*2+0],eval_self[i*2+1]);
  }

  double kappa = 0.2;double rho = 4 - 1.0/(2*kappa);int one_minus_halfD = 1;double Eerr = 1e-9;
  //double kappa = 0.2;double rho = 4 - 1.0/(2*kappa);int one_minus_halfD = 1;double Eerr = 1e-11;
  int n_zero = 0;
  Complexq rhoc = Complexq(rho,0.0);
  for(int j=0; j<n_vec; ++j){
    values_sys[j] = values_sys[j]/rhoc;
    if(qnorm(values_sys[j]) < Eerr) n_zero += 1;
  }
  int Nmass = mass.size();
  values.resize(Nmass*n_vec);
  for(int mi=0;mi<Nmass;mi++)
  //for(int iv=0;iv<n_vec;iv++)values[mi*n_vec+iv] = inv_self(values_sys[iv], mass[mi], rho,one_minus_halfD);
  for(int iv=0;iv<n_vec;iv++)values[iv*Nmass + mi] = inv_self(values_sys[iv], mass[mi], rho,one_minus_halfD);
  /////====Set up values


  {TIMER("Kernal functions");get_low_rho(eigen,values, n_zero,Mres,geo, in);}

  {
    int noden = qlat::get_num_node();
    std::vector<double > write;write.resize(Mres.size());
    for(int iv=0;iv<Mres.size();iv++){
      write[iv] = Mres[iv];
      ///write[iv*2+0] = Mres[iv].real();
      ///write[iv*2+1] = Mres[iv].imag();
    }
    char filename[500];
    sprintf(filename,"res/rho_test_%06d.dat",icfg);write_data(write,filename);
    sprintf(filename,"res/rho_test_%06d.Nmpi%02d.dat",icfg,noden);write_data(write,filename);
  }


  //int nmass = mass.size();
  //for(int op=0;op<16;op++)
  //for(int mi=0;mi<nmass;mi++)
  //for(int t=0;t<nt;t++)
  //{
  //  int offM = op*nmass*nt*nt + (mi*nt+0)*nt+t;
  //  print0("op %3d, mi %3d,t %5d, %13.8e  %13.8e \n",op,mi,t,Mres[offM].real(),Mres[offM].imag());
  //}


  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

