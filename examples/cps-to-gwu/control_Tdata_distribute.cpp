#include <sys/sysinfo.h>
#include "io_vec.h"
#include "general_funs.h"
////#include "utils_smear_src.h"
#include "check_fun.h"
#include "utils_Vec_redistribute.h"
////#include "utils_low_rho.h"

int main(int argc, char* argv[])
{
  using namespace qlat;
  int mode_dis = 1;

  inputpara in; 
  if(mode_dis == 0)
  { 
  int n_node = init_mpi(&argc, &argv);
  in.load_para("input.txt");
  Coordinate Lat(in.nx, in.ny, in.nz, in.nt);
  Coordinate spreadT = guess_nodeL(n_node, Lat);
  begin_comm(MPI_COMM_WORLD , spreadT);
  } 

  if(mode_dis == 1)
  {
  std::vector<Coordinate> size_node_list;
  add_nodeL(size_node_list);
  begin(&argc, &argv, size_node_list);
  in.load_para("input.txt");
  }

  set_GPU();

  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;

  ////int vini  = 0;

  omp_set_num_threads(omp_get_max_threads());
  print0("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 
  fflush_MPI();

  fft_desc_basic fd(geo);
  //Vec_redistribute vec_large(fd, 1);
  bool test_original = false;
  bool test_gather3D = 1;
  bool test_gather4D = 1;

  if(test_gather4D){
  Rotate_vecs rot(fd, 1);

  std::vector<int > secT;secT.resize(fd.Nmpi);
  for(int i=0;i<secT.size();i++){secT[i] = fd.Nt;}

  int Nt = fd.Nt;

  int nvol = in.nx*in.ny*in.nz;
  qlat::vector_acc<qlat::Complex > sendbuf;////recvbuf;
  qlat::vector_acc<qlat::Complex > databuf;
  /////sendbuf.set_acc(true);databuf.set_acc(true);
  int civ = 1;int Nvec = 3;
  int biva = Nvec * (fd.vol*fd.nt)/(fd.Nvol);
  int b0 = Nvec; 
  int c0 = 2*civ;

  sendbuf.resize(biva*fd.Nvol*civ);
  ///recvbuf.resize(b0*civ*nvol*fd.Nt);
  databuf.resize(biva*fd.Nvol*civ);
  for(long bi=0;bi<biva;bi++)
  for(long long vi=0;vi<fd.Nvol;vi++)
  for(long ci=0;ci<civ;ci++)
  {
    long long off = (bi*fd.Nvol+vi)*civ+ci;
    std::vector<int > p = fd.coordinate_g_from_index(vi, fd.rank);
    databuf[off] = ((p[3]*900+p[2])*900+p[1])*900+p[0] + std::cos(ci) + qlat::Complex(0.0,std::cos(bi));
  }

  memcpy(&sendbuf[0], &databuf[0], 2*sizeof(double)*sendbuf.size());
  //
  //qlat::Complex* tem_buf;
  //gpuMalloc(tem_buf, sendbuf.size(), qlat::Complex);


  ////cpy_data_thread(&tem_buf[0], &sendbuf[0], sendbuf.size(), 0);
  ////qacc_for(i, sendbuf.size(), {sendbuf[i]=sendbuf[i];});

  //rot.set_mem<double >(b0, civ, 1);
  rot.set_mem<qlat::Complex >(b0, civ, 1);

  rot.reorder(true, (qlat::Complex*) &sendbuf[0]);

  //cpy_data_thread(tem_buf, &sendbuf[0], sendbuf.size(),0);
  //rot.reorder(true, (double*) &tem_buf[0]);
  //cpy_data_thread(&sendbuf[0], tem_buf, sendbuf.size(),0);

  double diff = 0.0;
  for(long ni=0;ni<b0;ni++)
  for(long long vi=0;vi<fd.nt*nvol;vi++)
  for(long ci=0;ci<civ;ci++)
  {
    ////int t = Nti + fd.init;
    long long off = (ni*fd.nt*nvol+vi)*civ+ci;
    std::vector<int > p = fd.coordinate_g_from_g_index(vi);
    /////Becarefule about the order
    int bi = fd.get_mi_curr(4)*b0 + ni;

    qlat::Complex tem = sendbuf[off] - (((p[3]*900+p[2])*900+p[1])*900+p[0] + std::cos(ci) + qlat::Complex(0.0,std::cos(bi)));
    diff += qnorm(tem);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 0 %.3e .\n",diff);

  rot.reorder(false, (qlat::Complex*) &sendbuf[0]);
  //cpy_data_thread(tem_buf, &sendbuf[0], sendbuf.size(),0);
  //rot.reorder(false, (double*) &tem_buf[0]);
  //cpy_data_thread(&sendbuf[0], tem_buf, sendbuf.size(),0);


  diff = 0.0;
  for(size_t i = 0;i< sendbuf.size();i++)
  {
    diff += qnorm(sendbuf[i] - databuf[i]);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 1 %.3e .\n",diff);

  qlat::Complex* sendT;
  gpuMalloc(sendT, sendbuf.size() , qlat::Complex);
  for(int i=0;i<100;i++){TIMER("=====Rotate 4D");rot.reorder(false, (double*) &sendT[0]);}
  gpuFree(sendT);

  //for(int ni=0;ni<in.nvec;ni++){
  //rot.reorder(false, (qlat::Complex*) &sendbuf[0]);}
  //gpuFree(tem_buf);

  }


  if(test_gather3D){
  Rotate_vecs rot(fd, 0);

  /////Need to understand
  std::vector<int > secT;secT.resize(fd.Nmpi);
  for(int i=0;i<secT.size();i++){secT[i] = fd.Nt;}

  int Nt = fd.Nt;

  int nvol = in.nx*in.ny*in.nz;
  qlat::vector_acc<qlat::Complex > sendbuf;////recvbuf;
  qlat::vector_acc<qlat::Complex > databuf;
  ////sendbuf.set_acc(true);databuf.set_acc(true);
  int civ = 6;int Nvec = 6;
  int biva = Nvec * (fd.vol*fd.Nt)/(fd.Nvol);
  int b0 = Nvec; 
  int c0 = 2*civ;

  sendbuf.resize(b0*civ*nvol*fd.Nt);
  ///recvbuf.resize(b0*civ*nvol*fd.Nt);
  databuf.resize(b0*civ*nvol*fd.Nt);
  for(long bi=0;bi<biva;bi++)
  for(long long vi=0;vi<fd.Nvol;vi++)
  for(long ci=0;ci<civ;ci++)
  {
    long long off = (bi*fd.Nvol+vi)*civ+ci;
    std::vector<int > p = fd.coordinate_g_from_index(vi, fd.rank);
    databuf[off] = ((p[3]*900+p[2])*900+p[1])*900+p[0] + std::cos(ci) + qlat::Complex(0.0,std::cos(bi));
  }

  memcpy(&sendbuf[0], &databuf[0], 2*sizeof(double)*sendbuf.size());
  //qlat::Complex* tem_buf;

  //gpuMalloc(tem_buf, sendbuf.size(), qlat::Complex);

  rot.set_mem<double >(b0, c0, 0);
  rot.reorder(true, (double*) &sendbuf[0]);
    
  //cpy_data_thread(tem_buf, &sendbuf[0], sendbuf.size(),0);
  //rot.reorder(true, (double*) &tem_buf[0]);
  //cpy_data_thread(&sendbuf[0], tem_buf, sendbuf.size(),0);

  double diff = 0.0;
  for(long ni=0;ni<b0;ni++)
  for(int Nti=0;Nti<fd.Nt;Nti++)
  for(long long vi=0;vi<nvol;vi++)
  for(long ci=0;ci<civ;ci++)
  {
    int t = Nti + fd.init;
    long long off = ((ni*fd.Nt+Nti)*nvol+vi)*civ+ci;
    std::vector<int > p = fd.coordinate_g_from_g_index(t*nvol + vi);
    /////Becarefule about the order
    int bi = fd.get_mi_curr()*b0 + ni;

    qlat::Complex tem = sendbuf[off] - (((p[3]*900+p[2])*900+p[1])*900+p[0] + std::cos(ci) + qlat::Complex(0.0,std::cos(bi)));
    diff += qnorm(tem);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 0 %.3e .\n",diff);

  rot.reorder(false, (double*) &sendbuf[0]);

  //cpy_data_thread(tem_buf, &sendbuf[0], sendbuf.size(),0);
  //rot.reorder(false, (double*) &tem_buf[0]);
  //cpy_data_thread(&sendbuf[0], tem_buf, sendbuf.size(),0);

  diff = 0.0;
  for(size_t i = 0;i< sendbuf.size();i++)
  {
    diff += qnorm(sendbuf[i] - databuf[i]);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 1 %.3e .\n",diff);

  qlat::Complex* sendT;
  gpuMalloc(sendT, sendbuf.size() , qlat::Complex);
  for(int i=0;i<100;i++){TIMER("=====Rotate 3D");rot.reorder(false, (double*) &sendT[0]);}
  gpuFree(sendT);

  //for(int ni=0;ni<in.nvec;ni++){
  //rot.reorder(false, (double*) &tem_buf[0]);}
  //gpuFree(tem_buf);

  }


  if(test_original){
  Vec_redistribute vec_large(fd);


  /////Need to understand
  std::vector<int > secT;secT.resize(fd.Nmpi);
  for(int i=0;i<secT.size();i++){secT[i] = fd.Nt;}


  int Nt = fd.Nt;

  int nvol = in.nx*in.ny*in.nz;
  qlat::vector_acc<qlat::Complex > sendbuf,recvbuf;
  qlat::vector_acc<qlat::Complex > databuf;
  ////sendbuf.set_acc(true);recvbuf.set_acc(true);databuf.set_acc(true);
  int biva = 2*(fd.vol*fd.Nt)/(fd.Nvol);int civ = 6;
  int Nvec = biva*fd.Nvol/(fd.vol*fd.Nt);
  int b0 = Nvec; 
  int c0 = 2*civ;

  sendbuf.resize(b0*civ*nvol*fd.Nt);
  recvbuf.resize(b0*civ*nvol*fd.Nt);
  databuf.resize(b0*civ*nvol*fd.Nt);
  for(long bi=0;bi<biva;bi++)
  for(long long vi=0;vi<fd.Nvol;vi++)
  for(long ci=0;ci<civ;ci++)
  {
    long long off = (bi*fd.Nvol+vi)*civ+ci;
    //////Need a function in fd to get the coordinates?
    //int ti = fd.Pos0[fd.rank][3] +  vi/(fd.Nz*fd.Ny*fd.Nx);
    //int zi = fd.Pos0[fd.rank][2] + (vi%(fd.Nz*fd.Ny*fd.Nx))/(fd.Ny*fd.Nx);
    //int yi = fd.Pos0[fd.rank][1] + (vi%(fd.Ny*fd.Nx))/fd.Nx;
    //int xi = fd.Pos0[fd.rank][0] + (vi%fd.Nx);
    std::vector<int > p = fd.coordinate_g_from_index(vi, fd.rank);
    databuf[off] = ((p[3]*900+p[2])*900+p[1])*900+p[0] + std::cos(ci) + qlat::Complex(0.0,std::cos(bi));

  }

  ///for(int i=0;i<biva*12*12/civ;i++){memcpy(&sendbuf[i*Nt*N0*N1*N2*civ + 0],&srcE_civ[i](0),sizeof(Ftype)*2*Nt*N0*N1*N2*civ);}
  memcpy(&sendbuf[0], &databuf[0], 2*sizeof(double)*sendbuf.size());

  /////Data will be modified for sendbuf and recvbuf, results on sendbuf
  {
  TIMER("Single Vec redistribute");
  vec_large.reorder((double*) &sendbuf[0],(double*) &recvbuf[0],b0,c0 ,  0);
  } 

  double diff = 0.0;
  for(long ni=0;ni<b0;ni++)
  for(int Nti=0;Nti<fd.Nt;Nti++)
  for(long long vi=0;vi<nvol;vi++)
  for(long ci=0;ci<civ;ci++)
  {
    int t = Nti + fd.init;
    long long off = ((ni*fd.Nt+Nti)*nvol+vi)*civ+ci;
    ////////Need a function in fd to get the coordinates?
    //int ti = fd.Pos0[fd.rank][3] + Nti;
    ///////Order by ti 
    //int zi = (vi%(fd.nz*fd.ny*fd.nx))/(fd.ny*fd.nx);
    //int yi = (vi%(fd.ny*fd.nx))/fd.nx;
    //int xi = (vi%fd.nx);
    std::vector<int > p = fd.coordinate_g_from_g_index(t*nvol + vi);
    ///////Need a function for this offset from coordinate

    //int bi = ni*fd.mz*fd.my*fd.mx + fd.get_mi_curr();
    /////Becarefule about the order
    int bi = fd.get_mi_curr()*b0 + ni;

    /////databuf[off] = ((p[3]*900+p[2])*900+p[1])*900+p[0] + std::cos(ci) + qlat::Complex(0.0,std::cos(bi));
    qlat::Complex tem = sendbuf[off] - (((p[3]*900+p[2])*900+p[1])*900+p[0] + std::cos(ci) + qlat::Complex(0.0,std::cos(bi)));
    diff += qnorm(tem);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 0 %.3e .\n",diff);


  vec_large.reorder((double*) &sendbuf[0],(double*) &recvbuf[0],b0,c0 ,100);

  diff = 0.0;
  for(size_t i = 0;i< sendbuf.size();i++)
  {
    diff += qnorm(sendbuf[i] - databuf[i]);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 1 %.3e .\n",diff);
  }


  //qlat::Complex* sendT;
  //qlat::Complex* sbufT;
  ////gpuErrchk( cudaMalloc(&sendT    , b0*civ*nvol*fd.Nt * sizeof(qlat::Complex)));
  ////gpuErrchk( cudaMalloc(&sbufT    , b0*civ*nvol*fd.Nt * sizeof(qlat::Complex)));

  //gpuMalloc(sendT, b0*civ*nvol*fd.Nt , qlat::Complex);
  //gpuMalloc(sbufT, b0*civ*nvol*fd.Nt , qlat::Complex);

  //for(int i=0;i<20;i++){
  //TIMER("Single Vec redistribute GPU");
  //vec_large.reorder((double*) &sendT[0],(double*) &sbufT[0],b0,c0 ,  0);
  ////vec_large.reorder((double*) &sendbuf[0],(double*) &recvbuf[0],b0,c0 ,  0);
  //} 

  //gpuFree(sendT);
  //gpuFree(sbufT);


  //Propagator4d prop;prop.init(geo);
  //random_point_src(prop, 0);

  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

