#include <sys/sysinfo.h>
#include "utils_io_vec.h"
#include "general_funs.h"
#include "utils_smear_vecs.h"
#include "check_fun.h"

int main(int argc, char* argv[])
{
  using namespace qlat;
  inputpara in;begin_Lat(&argc, &argv, in);

  set_GPU();

  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;

  omp_set_num_threads(omp_get_max_threads());
  print0("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 
  fflush_MPI();

  ////if(nx == 24)twist_boundary_at_boundary(gf, EIGEN_PI/1.0, 3 );

  fft_desc_basic fd(geo);
  Vec_redistribute vec_large(fd);
  long Nvol = geo.local_volume();

  if(in.paraI != "None"){
    std::vector<std::string > Li = stringtolist(in.paraI);
    print0("Li %s, size %d \n", in.paraI.c_str(),int(Li.size()) );
    fflush_MPI();
    qassert(Li.size()%3 == 0);

    for(int si=0;si<Li.size()/3;si++)
    {
      int step     = stringtonum(   Li[si*3+0]);
      double width = stringtodouble(Li[si*3+1]);
      int  Nprop   = stringtodouble(Li[si*3+2]);
      print0("sn%03dsk%6.4f \n", step, width);

      unsigned int NVmpi = fd.mz*fd.my*fd.mx;
      int groupP = (12*Nprop + NVmpi-1)/NVmpi;
      int repeat = 1;
      if(groupP > 12){repeat = (groupP+12-1)/12;groupP=12;}
      print0("====Vec redistribute setup, repeat %d, NVmpi %d, groupP %d \n", repeat, NVmpi, groupP);

      ////EigenV propT;EigenV propT_buf;
      ////propT.resize(repeat*NVmpi*Nvol*groupP*12);
      ////propT_buf.resize(repeat*NVmpi*Nvol*groupP*12);
      ////EigenV gfET;gfET.resize(    NVmpi*6*Nvol*9);

      EigenV propT_tmp; propT_tmp.resize(repeat*NVmpi*Nvol*groupP*12);
      EigenV gfET_tmp;  gfET_tmp.resize(    NVmpi*6*Nvol*9);

      ////EigenV propT;EigenV propT_buf;
      Complexq* propT=NULL;
      Complexq* propT_buf=NULL;
      Complexq* gfET=NULL;

      //////size_t Np = repeat*NVmpi*Nvol*groupP*12;
      //gpuMalloc(((void**)&propT    ),propT_tmp.size()*sizeof(Complexq));
      //gpuMalloc(((void**)&propT_buf),propT_tmp.size()*sizeof(Complexq));
      //gpuMalloc(((void**)&gfET     ),gfET_tmp.size()*sizeof(Complexq));

      gpuMalloc(propT    ,propT_tmp.size(), Complexq);
      gpuMalloc(propT_buf,propT_tmp.size(), Complexq);
      gpuMalloc(gfET     ,gfET_tmp.size() , Complexq);


      //gpuErrchk(cudaMalloc(&propT, propT_tmp.size()*sizeof(Complexq)));
      //gpuErrchk(cudaMalloc(&propT_buf, propT_tmp.size()*sizeof(Complexq)));
      //gpuErrchk(cudaMalloc(&gfET, gfET_tmp.size()*sizeof(Complexq)));

      cudaMemcpy(&propT[0], &propT_tmp[0], propT_tmp.size()*sizeof(Complexq),cudaMemcpyDeviceToDevice);
      cudaMemcpy(&gfET[0], &gfET_tmp[0], gfET_tmp.size()*sizeof(Complexq),cudaMemcpyDeviceToDevice);

      ////random_numbers(propT, Np);

      smear_fun smf;
      smf.init_distribute(geo);

      long long Tfloat = 0;
      double mem       = 0.0;

      {long long Lat = geo.local_volume();
      int nsrc = groupP * NVmpi * repeat; 
      long long vGb = Lat *nsrc*4;
      int Fcount = 3*(3*6 + 2*2); 
      int direction   = 6;
      Tfloat = step*direction*vGb*Fcount;
      mem = (Lat*nsrc*12 + Lat*4*9)*8.0;}
      ////timer.flops += Tfloat;
      print0("Memory size %.3e GB, %.3e Gflop \n",
        mem/(1024.0*1024*1024), Tfloat/(1024.0*1024*1024));

      {
        TIMER_FLOPS("==compute time");

        {TIMER("Vec prop");vec_large.reorder(&propT[0],&propT_buf[0], repeat, groupP*12 ,   0);}
        for(int i =0; i< repeat; i++ )
        {
          long off = i*NVmpi*Nvol*groupP*12;
          int bfac = groupP; int d0 = 4;
          smear_propagator4(&propT[off], &gfET[0], width, step, &propT_buf[0], smf, bfac, d0);
        }
        {TIMER("Vec prop");vec_large.reorder(&propT[0],&propT_buf[0], repeat, groupP*12 , 100);}
        timer.flops += Tfloat;
      }

      gpuFree(propT);propT = NULL;
      gpuFree(propT_buf);propT_buf = NULL;
      gpuFree(gfET);gfET = NULL;


    }
  }

  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

