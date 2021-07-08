#include "io_gwu.h"
//#include "utils_low_rho.h"
#include "cach_reduce.h"
#include "reduce_V.h"

////#define Complexq qlat::ComplexF
#define Complexq qlat::Complex
#define Cfield qlat::FieldM<qlat::MvectorT<3,Complexq > ,1> 

inline Complexq inv_self(const Complexq& lam, double m, double rho,int one_minus_halfD=1)
{
  Complexq tem = (one_minus_halfD>0)?(1-lam/2)/(rho*lam+m*(1-lam/2)):1.0/(rho*lam+m*(1-lam/2));
  return tem;
}

void printDevProp(cudaDeviceProp devProp)
{
  printf("%s\n", devProp.name);
  printf("Major revision number:         %d\n", devProp.major);
  printf("Minor revision number:         %d\n", devProp.minor);
  printf("Total global memory:           %u", devProp.totalGlobalMem);
  printf(" bytes\n");
  printf("Number of multiprocessors:     %d\n", devProp.multiProcessorCount);
  printf("Max Number of threads per multiprocessors:    %d\n", devProp.maxThreadsPerMultiProcessor);
  printf("Total amount of shared memory per block: %u\n",devProp.sharedMemPerBlock);
  printf("Total registers per block:     %d\n", devProp.regsPerBlock);
  printf("Warp size:                     %d\n", devProp.warpSize);
  printf("Maximum memory pitch:          %u\n", devProp.memPitch);
  printf("Total amount of constant memory:         %u\n",   devProp.totalConstMem);
  return;
}


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

  int n_vec = 1;

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 

  //char namew[500],namer[500],prop_tmp[500],name[500],name_tem[500];
  //char name0[500],name1[500],filename[500];
  int nDevices;
  cudaGetDeviceCount(&nDevices);
  if(qlat::get_id_node() == 0)
  for (int i = 0; i < nDevices; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    printf("Device Number: %d\n", i);
    printf("  Device name: %s\n", prop.name);
    printf("  Memory Clock Rate (KHz): %d\n",
           prop.memoryClockRate);
    printf("  Memory Bus Width (bits): %d\n",
           prop.memoryBusWidth);
    printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
           2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
    printDevProp(prop);
  }

  int mode = 0;

  if(mode == 0){
    qlat::FieldM<Complexq ,1> fd;fd.init(geo);

    int Nt = geo.node_site[3];
    long Nsum = geo.local_volume()/geo.node_site[3];
    Complexq* p = &fd.get_elem(0);

    qlat::vector<Complexq > res0;res0.resize(Nt);set_zero(res0);
    qlat::vector<Complexq > res1;res1.resize(Nt);set_zero(res1);
    qlat::vector<Complexq > res2;res2.resize(Nt);set_zero(res2);

    ////for(size_t isp=0;isp<size_t(geo.local_volume());isp++){fint.get_elems(isp)[0] = 1.0*std::cos(isp);}
    ////get_elems return pointer V
    qacc_for(index, geo.local_volume(), {
      ///////Make a search of time 
      fd.get_elem(index) = 1.0*cos(index);
    });

  {
    {TIMER("Reduce CPU 0");
      for(int it=0;it<Nt;it++){
      p = &fd.get_elem(it*Nsum);
      Complexq& Ntem = res0[it];
      for(long isp=0;isp<Nsum;isp++)Ntem += p[isp];
    }}


    set_zero(res0);
    for(int iv=0;iv<n_vec;iv++)
    {TIMER("Reduce CPU 1");
      for(int it=0;it<Nt;it++){
      p = &fd.get_elem(it*Nsum);
      Complexq& Ntem = res0[it];
      for(long isp=0;isp<Nsum;isp++)Ntem += p[isp];
    }}

    //set_zero(res1);
    //{TIMER("Reduce GPU 0");
    //  p = &fd.get_elem(0);
    //  for(int it=0;it<Nt;it++){
    //  reduce_gpu((Complexq*) &p[it*Nsum],res1[it],Nsum,1,1,16, true);
    //}}

    set_zero(res1);
    {{TIMER("Reduce GPU 1");
      reduce_gpu2d_6(&fd.get_elem(0),&res1[0],Nsum,Nt);
    }}

    set_zero(res1);
    for(int iv=0;iv<n_vec;iv++)
    {
      //qacc_for(index, geo.local_volume(), {
      //  ///////Make a search of time 
      //  fd.get_elem(index) = iv*cos(index);
      //});

      {TIMER("Reduce GPU 2");
      reduce_gpu2d_6(&fd.get_elem(0),&res1[0],Nsum,Nt);
      }
    }

    //set_zero(res0);
    //for(int iv=0;iv<n_vec;iv++)
    //{
    //  {TIMER("Reduce GPU 2");
    //  reduce_gpu2d_6((double*)&fd.get_elem(0),(double* )&res0[0],Nsum,Nt);
    //  }
    //}

    //set_zero(res1);
    //for(int iv=0;iv<n_vec;iv++)
    //{{TIMER("Reduce GPU 3");
    //  reduce_gpu((double*)&fd.get_elem(0),(double*)&res1[0],Nsum,Nt);
    //}}

    ////////Cuda cach reduce
    {
    {TIMER("Reduce GPU 4");
    set_zero(res2);
    //reduce_gpu2d_6_cach((thrust::complex<float>*)&fd.get_elem(0),(thrust::complex<float>*)&res2[0] , Nsum, Nt);
    //reduce_gpu2d_6_cach((float*)&fd.get_elem(0),(float*)&res2[0] , Nsum, Nt, 2);
    //reduce_gpu2d_6_cach((double*)&fd.get_elem(0),(double*)&res2[0] , Nsum, Nt, 1);
    reduce_gpu2d_6_cach((double*)&fd.get_elem(0),(double*)&res2[0] , Nsum, Nt, 2);
    ////////Cuda cach reduce
    }
    }

    //for(int it=0;it<Nt;it++){
    //  print0("it%2d, %.6e %.6e, %.6e %.6e, \n",it,
    //      res0[it].real(),res0[it].imag(),res2[it].real(),res2[it].imag());
    //}

    double diff = 0.0;
    double diff2= 0.0;
    double sum0=0.0;double sum1=0.0;double sum2=0.0;
    for(int it=0;it<Nt;it++){
      diff += (res0[it].real() - res1[it].real())*(res0[it].real() - res1[it].real());
      diff += (res0[it].imag() - res1[it].imag())*(res0[it].imag() - res1[it].imag());

      diff2+= (res0[it].real() - res2[it].real())*(res0[it].real() - res2[it].real());
      diff2+= (res0[it].imag() - res2[it].imag())*(res0[it].imag() - res2[it].imag());

      sum0 += (res0[it].real())*(res0[it].real()) + (res0[it].imag())*(res0[it].imag());
      sum1 += (res1[it].real())*(res1[it].real()) + (res1[it].imag())*(res1[it].imag());
      sum2 += (res2[it].real())*(res2[it].real()) + (res2[it].imag())*(res1[it].imag());
    }
    print0("cpu %.6e, gpu %.6e %.6e, diff %.6e %.6e \n",sum0,sum1,sum2,diff,diff2);
  }


  //{
  //  qlat::vector<Complexq > Vres;Vres.resize(15);
  //  set_zero(Vres);
  //  
  //  //////set_zero(res1);

  //  reduce_buf<Complexq > Vbp(1,Nt*Nsum,16);
  //  Complexq *a = &fd.get_elem(0);
  //  {
  //    reduce_gpu2d(a,(Complexq*)&Vres[1],Vbp);
  //  }
  //  reduce_gpu(a,Vres[1],geo.local_volume(),128, 8,16);

  //  set_zero(Vres);
  //  for(int iv=0;iv<n_vec;iv++){
  //    {TIMER("Reduce CPU V ");reduce_cpu(a,Vres[0],geo.local_volume());}
  //  }

  //  //{TIMER("Reduce GPU V 0 ");
  //  //  Vbp.src[0] = a;
  //  //  reduce_gpu2d((Complexq*)&Vres1,Vbp);
  //  //}

  //  for(int iv=0;iv<n_vec;iv++){
  //    TIMER("Reduce GPU V 7");
  //    reduce_gpu(a,Vres[3],geo.local_volume(),128, 8,16);
  //  }

  //  for(int iv=0;iv<n_vec;iv++){
  //    TIMER("Reduce GPU V6 warm up");
  //    reduce_gpu2d_6(a,(Complexq*)&Vres[14],geo.local_volume(),1,16);
  //  }

  //  for(int iv=0;iv<n_vec;iv++){
  //    TIMER("Reduce GPU V6 0");
  //    reduce_gpu2d_6(a,(Complexq*)&Vres[4],geo.local_volume(),1,16);
  //  }

  //  for(int iv=0;iv<n_vec;iv++){
  //    //{TIMER("Reduce Global  V1");reduce_gpu(a,Vres1,geo.local_volume(), 1, 1,16, true);}
  //    //Vres1 = 0;
  //    {TIMER("Reduce GPU V 1");
  //      reduce_gpu2d(a,(Complexq*)&Vres[1],Vbp);
  //    }
  //  }

  //  for(int iv=0;iv<n_vec;iv++){
  //    TIMER("Reduce GPU V 2");
  //    reduce_gpu2d(a,(Complexq*)&Vres[2],geo.local_volume(),1,16);
  //  }


  //  ////Vres0 = 0;
  //  ////{TIMER("Reduce GPU V 8");}
  //  ////Vres0 = 0;
  //  ////{TIMER("Reduce GPU V 9");reduce_gpu(a,Vres0,geo.local_volume(),32, 32,32);}
  //  ////Vres0 = 0;
  //  ////{TIMER("Reduce GPU V 10");reduce_gpu(a,Vres0,geo.local_volume());}

  //  for(int i=0;i<10;i++){
  //    print0("i%2d, %.6e %.6e, \n",i,Vres[i].real(),Vres[i].imag());
  //  }
  //  ////print0("cpu %.6e %.6e, gpu %.6e %.6e .\n",Vres0.real(),Vres0.imag(),Vres1.real(),Vres1.imag());

  //}

  }

  //if(mode == 1)
  //{
  //  qlat::FieldM<double ,1> fint;fint.init(geo);
  //  for(size_t isp=0;isp<size_t(geo.local_volume());isp++){fint.get_elems(isp)[0]=1.0*std::cos(isp);}

  //  ///touchv(fint);
  //  //int *a = fint.field.v.p;
  //  double *a = (double*) &(fint.get_elems(0)[0]);
  //  double res1 = 0;
  //  double res0 = 0;
  //  double res2 = 0;
  //  double res3 = 0;

  //  for(int iv=0;iv<n_vec;iv++)
  //  {
  //    res1 = 0;
  //    {TIMER("Reduce CPU  ");reduce_cpu(a,res1,geo.local_volume());}
  //    res0 = 0;
  //    {TIMER("Reduce GPU 0");reduce_gpu(a,res0,geo.local_volume());}
  //    res0 = 0;
  //    {TIMER("Reduce GPU 1");reduce_gpu(a,res0,geo.local_volume(),128, 8, 4);}
  //    res0 = 0;
  //    {TIMER("Reduce GPU 2");reduce_gpu(a,res0,geo.local_volume(),256, 16,2);}
  //    res0 = 0;
  //    {TIMER("Reduce GPU 3");reduce_gpu(a,res0,geo.local_volume(),128, 16,2);}
  //    res0 = 0;
  //    {TIMER("Reduce GPU 4");reduce_gpu(a,res0,geo.local_volume(),64,  16,2);}
  //    res0 = 0;
  //    {TIMER("Reduce GPU 5");reduce_gpu(a,res0,geo.local_volume(),16, 16,2);}
  //    res0 = 0;
  //    {TIMER("Reduce GPU 6");reduce_gpu(a,res0,geo.local_volume(), 8, 8, 2);}
  //    res0 = 0;
  //    {TIMER("Reduce GPU 7");reduce_gpu(a,res0,geo.local_volume(),128, 8,16);}
  //    res0 = 0;
  //    {TIMER("Reduce GPU 8");reduce_gpu(a,res0,geo.local_volume(),128, 8,32);}
  //    res0 = 0;
  //    {TIMER("Reduce GPU 9");reduce_gpu(a,res0,geo.local_volume(),32, 32,32);}
  //    res0 = 0;
  //    {TIMER("Reduce GPU 10");reduce_gpu(a,res0,geo.local_volume());}

  //    res3 = 0;
  //    {TIMER("Reduce Global 1");reduce_gpu(a,res3,geo.local_volume(), 1, 1,16, true);}


  //    //res2 = 0;
  //    //{TIMER("Reduce Global 0");reduce_global(a,res2,geo.local_volume());}

  //  }

  //  double length = n_vec*sizeof(double)*(geo.local_volume()/(1024*1024*1024.0));
  //  print0("Data size %.6e GB \n",length);
  //  print0("gpu %.6e, cpu %.6e, global %.6e %.6e, geoL %8d \n",res0,res1,res2,res3,geo.local_volume());
  //}


  fflush_MPI();
  qlat::Timer::display();



  qlat::end();
  return 0;
}

