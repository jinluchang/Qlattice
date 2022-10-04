#include <sys/sysinfo.h>
#include "utils_io_vec.h"
#include "general_funs.h"
////#include "utils_smear_src.h"
#include "check_fun.h"
#include "utils_Vec_redistribute.h"
////#include "utils_low_rho.h"

int main(int argc, char* argv[])
{
  using namespace qlat;
  {

  inputpara in;begin_Lat(&argc, &argv, in);
  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;

  set_GPU();

  ////int vini  = 0;

  omp_set_num_threads(omp_get_max_threads());
  print0("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 
  fflush_MPI();

  fft_desc_basic fd(geo);
  //Vec_redistribute vec_large(fd, 1);
  bool test_simple   = false;
  bool test_gather3D = false;
  bool test_gather4D = false;
  bool test_original = false;
  std::vector<std::string > Li = stringtolist(in.paraI);
  if(Li.size() >= 4){
    if(stringtonum(Li[0]) == 1){test_simple = true;}
    if(stringtonum(Li[1]) == 1){test_gather3D = true;}
    if(stringtonum(Li[2]) == 1){test_gather4D = true;}
    if(stringtonum(Li[3]) == 1){test_original = true;}
  }else{print0("======Need paraI to set test sets. \n");}

  if(test_simple){

  Vec_redistribute vec_large(fd, 1);

  qlat::vector_gpu<qlat::Complex > prop;
  qlat::vector_gpu<qlat::Complex > prop_buf;
  int NVmpi = fd.mv[0]*fd.mv[1]*fd.mv[2];
  long Nvol = geo.local_volume();
  prop.resize(NVmpi * Nvol * 12);
  prop_buf.resize(prop.size());
  random_Ty(prop.data(), prop.size(), 1, in.seed);

  vec_large.reorder(prop.data(), prop_buf.data(), 1, 12 ,   0);
  }

  if(test_gather4D){
  Rotate_vecs rot(fd, 1);

  ////int Nt = fd.Nt;

  int nvol = in.nx*in.ny*in.nz;
  qlat::vector_acc<qlat::Complex > sendbuf;////recvbuf;
  qlat::vector_gpu<qlat::Complex > sendbufG;////recvbuf;
  qlat::vector_acc<qlat::Complex > databuf;
  /////sendbuf.set_acc(true);databuf.set_acc(true);
  int civ = 1;int Nvec = 3;
  int biva = Nvec * (fd.vol*fd.nt)/(fd.Nvol);
  int b0 = Nvec; 
  ////int c0 = 2*civ;

  sendbuf.resize(biva*fd.Nvol*civ);
  ///recvbuf.resize(b0*civ*nvol*fd.Nt);
  databuf.resize(biva*fd.Nvol*civ);
  for(long bi=0;bi<biva;bi++)
  for(unsigned long long vi=0;vi<fd.Nvol;vi++)
  for(long ci=0;ci<civ;ci++)
  {
    long long off = (bi*fd.Nvol+vi)*civ+ci;
    const Coordinate& p = fd.coordinate_g_from_index(vi, fd.rank);
    databuf[off] = ((p[3]*900+p[2])*900+p[1])*900+p[0] + std::cos(ci) + qlat::Complex(0.0,std::cos(bi));
  }

  memcpy(&sendbuf[0], &databuf[0], 2*sizeof(double)*sendbuf.size());
  rot.set_mem<qlat::Complex >(b0, civ, 1);

  ////rot.reorder(true, (qlat::Complex*) &sendbuf[0]);
  sendbufG.copy_from(sendbuf, 1, 1);
  rot.reorder(true, (qlat::Complex*) &sendbufG[0]);
  sendbufG.copy_to(sendbuf, 1);

  double diff = 0.0;
  for(long ni=0;ni<b0;ni++)
  for(long long vi=0;vi<fd.nt*nvol;vi++)
  for(long ci=0;ci<civ;ci++)
  {
    ////int t = Nti + fd.init;
    long long off = (ni*fd.nt*nvol+vi)*civ+ci;
    const Coordinate& p = fd.coordinate_g_from_g_index(vi);
    /////Be careful about the order
    int bi = fd.get_mi_curr(4)*b0 + ni;

    qlat::Complex tem = sendbuf[off] - (((p[3]*900+p[2])*900+p[1])*900+p[0] + std::cos(ci) + qlat::Complex(0.0,std::cos(bi)));
    diff += qnorm(tem);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 0 %.3e .\n",diff);

  //rot.reorder(false, (qlat::Complex*) &sendbuf[0]);

  sendbufG.copy_from(sendbuf, 1, 1);
  rot.reorder(false, (qlat::Complex*) &sendbufG[0]);
  sendbufG.copy_to(sendbuf, 1);


  diff = 0.0;
  for(long i = 0;i< sendbuf.size();i++)
  {
    diff += qnorm(sendbuf[i] - databuf[i]);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 1 %.3e .\n",diff);

  qlat::vector_gpu<qlat::Complex > sendT;sendT.resize(sendbuf.size());
  for(int i=0;i<100;i++){TIMER("=====Rotate 4D");rot.reorder(false, (double*) sendT.data());}

  }

  if(test_gather3D){
  Rotate_vecs rot(fd, 1);

  /////Need to understand
  ////int Nt = fd.Nt;

  int nvol = in.nx*in.ny*in.nz;
  qlat::vector_acc<qlat::Complex > sendbuf;////recvbuf;
  qlat::vector_gpu<qlat::Complex > sendbufG;////recvbuf;
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
  for(unsigned long long vi=0;vi<fd.Nvol;vi++)
  for(long ci=0;ci<civ;ci++)
  {
    long long off = (bi*fd.Nvol+vi)*civ+ci;
    const Coordinate& p = fd.coordinate_g_from_index(vi, fd.rank);
    databuf[off] = ((p[3]*900+p[2])*900+p[1])*900+p[0] + std::cos(ci) + qlat::Complex(0.0,std::cos(bi));
  }

  memcpy(&sendbuf[0], &databuf[0], 2*sizeof(double)*sendbuf.size());

  rot.set_mem<double >(b0, c0, 0);

  sendbufG.copy_from(sendbuf, 1, 1);
  rot.reorder(true, (double*) sendbufG.data());
  sendbufG.copy_to(sendbuf, 1);

  //rot.reorder(true, (double*) &sendbuf[0]);
  //rot.reorder(true, (double*) sendbuf.v.p);
    
  double diff = 0.0;
  for(long ni=0;ni<b0;ni++)
  for(int Nti=0;Nti<fd.Nt;Nti++)
  for(long long vi=0;vi<nvol;vi++)
  for(long ci=0;ci<civ;ci++)
  {
    int t = Nti + fd.init;
    long long off = ((ni*fd.Nt+Nti)*nvol+vi)*civ+ci;
    const Coordinate& p = fd.coordinate_g_from_g_index(t*nvol + vi);
    /////Be careful about the order
    int bi = fd.get_mi_curr()*b0 + ni;

    qlat::Complex tem = sendbuf[off] - (((p[3]*900+p[2])*900+p[1])*900+p[0] + std::cos(ci) + qlat::Complex(0.0,std::cos(bi)));
    diff += qnorm(tem);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 2 %.3e .\n",diff);

  sendbufG.copy_from(sendbuf, 1, 1);
  rot.reorder(false, (double*) sendbufG.data());
  sendbufG.copy_to(sendbuf, 1);

  //rot.reorder(false, (double*) sendbuf.v.p);
  //double* tmp = (double*) sendbuf.v.p;
  //rot.reorder(false, tmp);

  diff = 0.0;
  for(long i = 0;i< sendbuf.size();i++)
  {
    diff += qnorm(sendbuf[i] - databuf[i]);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 3 %.3e .\n",diff);

  qlat::vector_gpu<qlat::Complex > sendT;sendT.resize(sendbuf.size());
  for(int i=0;i<100;i++){TIMER("=====Rotate 3D");rot.reorder(false, (double*) sendT.data());}

  }


  if(test_original){
  Vec_redistribute vec_large(fd);

  ////int Nt = fd.Nt;

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
  for(unsigned long long vi=0;vi<fd.Nvol;vi++)
  for(long ci=0;ci<civ;ci++)
  {
    long long off = (bi*fd.Nvol+vi)*civ+ci;
    //////Need a function in fd to get the coordinates?
    //int ti = fd.Pos0[fd.rank][3] +  vi/(fd.Nz*fd.Ny*fd.Nx);
    //int zi = fd.Pos0[fd.rank][2] + (vi%(fd.Nz*fd.Ny*fd.Nx))/(fd.Ny*fd.Nx);
    //int yi = fd.Pos0[fd.rank][1] + (vi%(fd.Ny*fd.Nx))/fd.Nx;
    //int xi = fd.Pos0[fd.rank][0] + (vi%fd.Nx);
    const Coordinate& p = fd.coordinate_g_from_index(vi, fd.rank);
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
    const Coordinate& p = fd.coordinate_g_from_g_index(t*nvol + vi);
    ///////Need a function for this offset from coordinate

    //int bi = ni*fd.mz*fd.my*fd.mx + fd.get_mi_curr();
    /////Be careful about the order
    int bi = fd.get_mi_curr()*b0 + ni;

    /////databuf[off] = ((p[3]*900+p[2])*900+p[1])*900+p[0] + std::cos(ci) + qlat::Complex(0.0,std::cos(bi));
    qlat::Complex tem = sendbuf[off] - (((p[3]*900+p[2])*900+p[1])*900+p[0] + std::cos(ci) + qlat::Complex(0.0,std::cos(bi)));
    diff += qnorm(tem);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 4 %.3e .\n",diff);


  vec_large.reorder((double*) &sendbuf[0],(double*) &recvbuf[0],b0,c0 ,100);

  diff = 0.0;
  for(long i = 0;i< sendbuf.size();i++)
  {
    diff += qnorm(sendbuf[i] - databuf[i]);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 5 %.3e .\n",diff);
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
  }

  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

