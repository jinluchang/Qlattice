#include <sys/sysinfo.h>
#include "utils_io_vec.h"
#include "general_funs.h"
#include "check_fun.h"
#include "utils_FFT_GPU.h"

#define TyD qlat::Complex
#define TyF qlat::Complex 

int main(int argc, char* argv[])
{
  using namespace qlat;

  ///int mode_dis = 0;
  //inputpara in;int mode_dis = 0;
  //begin_Lat(&argc, &argv, in, mode_dis);
  inputpara in;begin_Lat(&argc, &argv, in);

  //if(mode_dis == 0)
  //{
  //int n_node = init_mpi(&argc, &argv);
  //in.load_para("input.txt");
  //Coordinate Lat(in.nx, in.ny, in.nz, in.nt);
  //Coordinate spreadT = guess_nodeL(n_node, Lat);
  /////3D begin
  //////begin_comm(MPI_COMM_WORLD , spreadT);

  /////4D begin
  //int id_node, n;
  //MPI_Comm_size(MPI_COMM_WORLD, &n);
  //MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
  //int t =  id_node/(spreadT[0]*spreadT[1]*spreadT[2]);
  //int z = (id_node/(spreadT[0]*spreadT[1]))%(spreadT[2]);
  //int y = (id_node/(spreadT[0]))%(spreadT[1]);
  //int x = (id_node%(spreadT[0]));
  /////int new_id = ((z*spreadT[1] + y)*spreadT[0] + x)*spreadT[3] + t;
  ////int new_id = ((x*spreadT[1] + y)*spreadT[2] + z)*spreadT[3] + t;
  ////begin(new_id, spreadT);
  //begin(id_node, spreadT);
  //}

  //if(mode_dis == 1)
  //{
  //std::vector<Coordinate> size_node_list;
  //add_nodeL(size_node_list);
  //begin(&argc, &argv, size_node_list);
  //in.load_para("input.txt");
  //}
  //begin_Lat(&argc, &argv, "input.txt", in, mode_dis);


  fflush_MPI();

  set_GPU();

  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;

  omp_set_num_threads(omp_get_max_threads());
  print0("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());
  fflush_MPI();
  print_mem_info();

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 
  //char nameP0[500]
  //char nameP1[500];

  bool check_fft_with_qlat   = 1;

  bool GPU = true;
  bool ft4D = true;
  bool checkdiff = true;
  const int Nvec =    10;
  int bfac = in.bfac;
  int mode_FFT_MPI = in.mode_FFT_MPI;
  if(in.debuga == 3){ft4D = false;}
  if(sizeof(TyF) == sizeof(ComplexF)){checkdiff = false;}
  //if(in.nvec >= 2){checkdiff = false;}
  //if(checkdiff){bfac = 1;}
  if(bfac != 1){checkdiff = false;}

  /////==========test correct
  if(check_fft_with_qlat){

  std::vector<int > nv,Nv,mv;
  geo_to_nv(geo, nv,Nv,mv);
  fft_desc_basic fd(geo);
  std::vector<int > dimN;

  std::vector<qlat::FieldM<TyF, Nvec> > src;src.resize(in.nvec);
  for(int iv=0;iv<in.nvec;iv++)src[iv].init(geo);
  TyF* P0 = (TyF*) (qlat::get_data(src[0]).data());


  qlat::RngState rs(qlat::get_id_node() + 134 );
  double ini = qlat::u_rand_gen(rs);
  for(int iv=0;iv<in.nvec;iv++){
  P0 = (TyF*) (qlat::get_data(src[iv]).data());
  qacc_for(isp, geo.local_volume(),{
  for(int di=0;di<Nvec;di++){
      P0[isp*Nvec + di ] = TyF(std::cos((ini+isp+di)*0.5) , (5.0/(isp+1))*(ini+di)*0.1);
    }
  });}

  fft_schedule fft3D(fd, GPU);
  dimN.resize(3);dimN[0] = nv[2];dimN[1] = nv[1];dimN[2] = nv[0];
  if(!ft4D){
  fft3D.set_mem<TyF >(src.size(), Nvec, dimN, mode_FFT_MPI,  Nvec);
  //fft3D.set_mem<TyF >(src.size(), Nvec, dimN, -3,  Nvec);
  fft3D.print_info();}


  fft_schedule fft4D(fd, GPU);
  dimN.resize(4);dimN[0] = nv[3];dimN[1] = nv[2];dimN[2] = nv[1];dimN[3] = nv[0];
  if( ft4D){
  fft4D.set_mem<TyF >(src.size(), Nvec, dimN, mode_FFT_MPI,   Nvec);
  fft4D.print_info();}

  qlat::FieldM<TyD, Nvec> srcF;srcF.init(geo);
  qlat::FieldM<TyF, Nvec> srcT;srcT.init(geo);
  srcT = src[0];

  P0 = (TyF*) (qlat::get_data(src[0]).data());
  TyD* Pa = (TyD*) (qlat::get_data(srcF  ).data());
  for(long isp=0;isp<geo.local_volume(); isp++)
  for(int di=0;di<Nvec;di++){
    Pa[isp*Nvec + di] = P0[isp*Nvec + di];
  }

  if( ft4D){
  if(checkdiff){TIMER("qlat fft");qlat::fft_complex_field(srcF, false);}
  {TIMER("new fft ");fft_fieldM(fft4D, src, false);}}

  if(!ft4D){
  if(checkdiff){
    TIMER("qlat fft");
    qlat::fft_complex_field_spatial(srcF, false);
    //if(int iv=0;iv<in.nvec;iv++){qlat::fft_complex_field_spatial(srcF, false);}
  }
  {
    TIMER("new fft ");
    for(int iv=0;iv<bfac;iv++){
      fft_fieldM(fft3D, src, false);
    }
    //fft_fieldM(fft3D, src, false);

    //fft_fieldM(src, false, ft4D);

    //FFTGPUPlanKey fkey = get_fft_gpu_plan_key(src, ft4D);
    //////fft_fieldM(*((fft_schedule*) get_fft_gpu_plan(fkey).fftP), src, false);
    ////get_fft_gpu_plan(fkey).fftP->print_info();
    ////fft_gpu_copy ft = make_fft_gpu_plan(fkey);
    ////ft.set();
    ////ft.fftP->print_info();
    //fft_fieldM(*(get_fft_gpu_plan(fkey).fftP), src, false);
    //fft_fieldM(*(get_fft_gpu_plan(fkey).fftP), src, true );
    //fft_fieldM(*(get_fft_gpu_plan(fkey).fftP), src, false);


    //for(int iv=0;iv<in.nvec;iv++)
    //{
    //P0 = (TyF*) (qlat::get_data(src[iv]).data());
    //qacc_for(isp, geo.local_volume(),{
    //for(int di=0;di<Nvec;di++){
    //    P0[isp*Nvec + di ] /= (nx*ny*nz);
    //  }
    //});}


  }}

  double fftdiff = 0.0;
  P0 = (TyF*) (qlat::get_data(src[0]).data());
  Pa = (TyD* ) (qlat::get_data(srcF  ).data());
  for(long isp=0;isp<geo.local_volume(); isp++)
  for(int di=0;di<Nvec;di++){
    TyD tem = Pa[isp*Nvec + di];
    fftdiff += qnorm(TyF(tem.real(),tem.imag()) - P0[isp*Nvec + di]);
  }
  sum_all_size(&fftdiff , 1);
  print0("===fft total diff %.3e \n", fftdiff);

  std::vector<TyD >  dat0;
  std::vector<TyD >  dat1;
  std::vector<TyD >  dat2;
  std::vector<int > mom(4);
  mom[0] = 0;mom[1] = 0;mom[2] = 3;mom[3] = 32;
  get_mom_apply(srcT , mom, dat0, false, ft4D);
  get_mom_fft(src[0] , mom, dat1, ft4D);
  get_mom_fft(srcF   , mom, dat2, ft4D);

  if(checkdiff){
  double diff[2];
  diff[0] = 0.0;diff[1] = 0.0;
  for(long di=0;di<long(dat0.size());di++){
    diff[0] += qlat::qnorm(dat0[di] - dat1[di]);
    diff[1] += qlat::qnorm(dat0[di] - dat2[di]);
    if(qlat::qnorm(dat0[di] - dat1[di]) > 1e-5 or qlat::qnorm(dat0[di] - dat2[di]) > 1e-5)
    print0("diff di %5ld, expect %.3e %.3e, new %.3e %.3e, qlat %.3e %.3e \n",
      di, dat0[di].real(), dat0[di].imag(),
          dat1[di].real(), dat1[di].imag(),
          dat2[di].real(), dat2[di].imag() );
  }
  print0("===diff 0 %.3e , 1 %.3e \n", diff[0], diff[1]);
  }

  //for(unsigned int i=0;i<nvec;i++){
  //  TIMER("Test runs");
  //  if( ft4D)fft_fieldM(fft4D, src, false);
  //  if(!ft4D)fft_fieldM(fft3D, src, false);
  //}

  }
  /////==========test correct

  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

